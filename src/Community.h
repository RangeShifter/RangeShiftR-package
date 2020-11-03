/*----------------------------------------------------------------------------
 *	
 *	Copyright (C) 2020 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell 
 *	
 *	This file is part of RangeShifter.
 *	
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *	
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *	
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *	
 --------------------------------------------------------------------------*/
 
 
/*------------------------------------------------------------------------------

RangeShifter v2.0 Community

Implements the Community class

There is ONLY ONE instance of a Community in an individual replicate simulation.
It holds a SubCommunity for each Patch in the Landscape (including the matrix),
and is thus the highest-level entity accessed for most processing concerned with
simulated populations.

Optionally, the Community maintains a record of the occupancy of suitable cells
or patches during the course of simulation of multiple replicates.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Peer G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 7 January 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef CommunityH
#define CommunityH

#include <vector>
#include <algorithm>
using namespace std;

#include "Version.h"
#include "SubCommunity.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "Species.h"

//---------------------------------------------------------------------------
struct commStats {
int ninds,nnonjuvs,suitable,occupied;
int minX,maxX,minY,maxY;
};

class Community {

public:
	Community(Landscape*);
	~Community(void);
	SubCommunity* addSubComm(Patch*,int);
	// functions to manage populations occurring in the community
	void initialise(
		Species*,	// pointer to Species
		int				// year (relevent only for seedType == 2)
	);
	void addManuallySelected(void);
	void resetPopns(void);
	void localExtinction(int);
	void patchChanges(void);
	void reproduction(
		int				// year
	);
	void emigration(void);
	void dispersal(
		short,	// landscape change index
		short		// season / year
	);

	void survival(
		short,	// part:		0 = determine survival & development,
						//		 			1 = apply survival changes to the population
		short,	// option0:	0 = stage 0 (juveniles) only         )
						//					1 = all stages                       ) used by part 0 only
						//					2 = stage 1 and above (all non-juvs) )
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
	);
	void ageIncrement(void);
	int totalInds(void);
	Population* findPop( // Find the population of a given species in a given patch
		Species*, // pointer to Species
		Patch*		// pointer to Patch
	);
	commStats getStats(void);
	void createOccupancy(
		int,	// no. of rows = (no. of years / interval) + 1
		int		// no. of replicates
	);
	void updateOccupancy(
		int,	// row = (no. of years / interval)
		int		// replicate
	);
	void deleteOccupancy(
		int		// no. of rows (as above)
	);

	bool outRangeHeaders( // Open range file and write header record
		Species*,	// pointer to Species
		int				// Landscape number (-999 to close the file)
	);
	void outRange( // Write record to range file
		Species*, // pointer to Species
		int,			// replicate
		int,			// year
		int				// generation
	);
	bool outPopHeaders( // Open population file and write header record
		Species*, // pointer to Species
		int       // option: -999 to close the file
	);
	void outPop( // Write records to population file
		int,	// replicate
		int,	// year
		int		// generation
	);

	void outInds( // Write records to individuals file
		int,	// replicate
		int,	// year
		int,	// generation
		int		// Landscape number (>= 0 to open the file, -999 to close the file
					//									 -1 to write data records)
	);
	void outGenetics( // Write records to genetics file
		int,	// replicate
		int,	// year
		int,	// generation
		int		// Landscape number (>= 0 to open the file, -999 to close the file
					//									 -1 to write data records)
	);
	// Open occupancy file, write header record and set up occupancy array
	bool outOccupancyHeaders(
		int		// option: -999 to close the file
	);
	void outOccupancy(void);
	void outOccSuit(
		bool	// TRUE if occupancy graph is to be viewed on screen
	);
	void viewOccSuit( // Update the occupancy graph on the screen
										// NULL for the batch version
		int,		// year
		double,	// mean occupancy
		double	// standard error of occupancy
	);
	bool outTraitsHeaders( // Open traits file and write header record
		Species*,	// pointer to Species
		int				// Landscape number (-999 to close the file)
	);
	bool outTraitsRowsHeaders( // Open trait rows file and write header record
		Species*, // pointer to Species
		int       // Landscape number (-999 to close the file)
	);
	void outTraits( // Write records to traits file
		traitCanvas,// pointers to canvases for drawing variable traits
//		emigCanvas,	// pointers to canvases for drawing emigration traits
//		trfrCanvas, // pointers to canvases for drawing emigration traits
								// see SubCommunity.h
								// in the batch version, these are replaced by integers set to zero
		Species*,		// pointer to Species
		int,				// replicate
		int,				// year
		int					// generation
	);
	void writeTraitsRows( // Write records to trait rows file
		Species*,	// pointer to Species
		int,			// replicate
		int,			// year
		int,			// generation
		int,			// row number (Y cell co-ordinate)
		traitsums	// structure holding sums of trait genes for dispersal (see Population.h)
	);
	void draw(	// Draw the Community on the landscape map and optionally save the map
							// NULL for the batch version
		int,	// replicate
		int,	// year
		int,	// generation
		int		// Landscape number
	);
    Rcpp::IntegerMatrix addYearToPopList(int,int);

private:
	Landscape *pLandscape;
	int indIx;				// index used to apply initial individuals
	float **occSuit;	// occupancy of suitable cells / patches
	std::vector <SubCommunity*> subComms;

};

extern paramSim *paramsSim;
extern paramInit *paramsInit;


//---------------------------------------------------------------------------
#endif
