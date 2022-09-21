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
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 25 June 2021 by Anne-Kathleen Malchow

------------------------------------------------------------------------------*/

#ifndef CommunityH
#define CommunityH

#if VCL
//#include <System.Classes.hpp>
#include <VCLTee.Chart.hpp>
#endif

#include <vector>
#include <algorithm>
using namespace std;

#include "Version.h"
#include "SubCommunity.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "Species.h"
#if RS_ABC
#include "ABC.h"
#endif
#if PEDIGREE
#include "Group.h"
#include "Pedigree.h"
#endif
#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN 

//---------------------------------------------------------------------------
struct commStats {
int ninds,nnonjuvs,suitable,occupied;
#if GOBYMODEL
int nsocial,nasocial;
#endif
int minX,maxX,minY,maxY;
};

class Community {

public:
	Community(Landscape*);
	~Community(void);
	SubCommunity* addSubComm(Patch*,int);
#if RS_CONTAIN
	void setHabIndex(
		Species*,	// pointer to Species
		short			// landscape change index
	);
#endif // RS_CONTAIN 
	// functions to manage populations occurring in the community
#if PEDIGREE
	void initialise(
		Species*,		// pointer to Species
		Pedigree*,	// pointer to Pedigree
		int					// year (relevent only for seedType == 2)
	);
#else
	void initialise(
		Species*,	// pointer to Species
		int				// year (relevent only for seedType == 2)
	);
#endif
	void addManuallySelected(void);
	void resetPopns(void);
	void localExtinction(int);
	void patchChanges(void);
#if SEASONAL
	void reproduction(
		int,			// year
		short			// season
	);
#else
#if GROUPDISP
	void reproduction(
		Species*,	// pointer to Species
		int				// year
	);
#else
#if BUTTERFLYDISP
	void reproduction(
		Species*,	// pointer to Species
		int,			// year
		short			// option: 0 = default (all reproduction before dispersal),
							// 1 = mating only (before dispersal),
							// 2 = parturition only (after dispersal)
	);
	void fledge(void);
#else
	void reproduction(
		int				// year
	);
#endif // BUTTERFLYDISP 
#endif // GROUPDISP 
#endif // SEASONAL
#if RS_DISEASE
	void emigration(
		Species*,	// pointer to Species
		short			// season
	);
#else
#if SEASONAL
	void emigration(
		short		// season
	);
#else
	void emigration(void);
#endif // SEASONAL
#endif // RS_DISEASE  
#if RS_ABC
	void dispersal(
		short,	// landscape change index
		bool		// TRUE if there are any connectivity observations
	);
#else
#if PEDIGREE
	void dispersal(
		Pedigree*,	// pointer to Pedigree
		int,        // replicate
		int,        // year
		int,        // generation
		short				// landscape change index
	);
#else
#if SEASONAL || RS_RCPP
	void dispersal(
		short,	// landscape change index
		short		// season / year
	);
#else
	void dispersal(
		short		// landscape change index
	);
#endif // SEASONAL || RS_RCPP
#endif // PEDIGREE
#endif // RS_ABC

#if SPATIALMORT
	void survival(
		short,	// part:		0 = determine survival & development,
						//		 			1 = apply survival changes to the population
		short,	// spatial mortality period (0 or 1)
		short		// option:	0 = stage 0 (juveniles) only         )
						//					1 = all stages                       ) used by part 0 only
						//					2 = stage 1 and above (all non-juvs) )
	);
#else
#if SEASONAL
	void survival(
		short,	// season
		short,	// part:		0 = determine survival & development,
						//		 			1 = apply survival changes to the population
		short,	// option0:	0 = stage 0 (juveniles) only         )
						//					1 = all stages                       ) used by part 0 only
						//					2 = stage 1 and above (all non-juvs) )
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
	);
#else
#if PEDIGREE
	void survival(
		Pedigree*,	// pointer to Pedigree
		short,			// part:		0 = determine survival & development,
								//		 			1 = apply survival changes to the population
		short,			// option0:	0 = stage 0 (juveniles) only         )
								//					1 = all stages                       ) used by part 0 only
								//					2 = stage 1 and above (all non-juvs) )
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
	);
#else
	void survival(
		short,	// part:		0 = determine survival & development,
						//		 			1 = apply survival changes to the population
		short,	// option0:	0 = stage 0 (juveniles) only         )
						//					1 = all stages                       ) used by part 0 only
						//					2 = stage 1 and above (all non-juvs) )
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
	);
#endif // PEDIGREE
#endif // SEASONAL 
#endif // SPATIALMORT
#if RS_CONTAIN
	int findCullTargets(Cull*,int,int);
	void cullAllTargets(Cull*);
	void cullRandomTargets(Cull*,int);
	void cullTargets(Cull*);
//	void cullClosestTargets(Cull*,int);
	void resetCullTargets(void);
	void resetCull(void); 
	void updateDamage(Species*,Cull*);
#endif // RS_CONTAIN 
	void ageIncrement(void);
	int totalInds(void);
	Population* findPop( // Find the population of a given species in a given patch
		Species*, // pointer to Species
		Patch*		// pointer to Patch
	);
#if SEASONAL
	commStats getStats(
		short		// season
	);
#else
	commStats getStats(void);
#endif // SEASONAL 
	void createOccupancy(
		int,	// no. of rows = (no. of years / interval) + 1
		int		// no. of replicates
	);
#if SEASONAL
	void updateOccupancy(
		int,	// row = (no. of years / interval)
		int,	// replicate
		short	// season
	);
#else
	void updateOccupancy(
		int,	// row = (no. of years / interval)
		int		// replicate
	);
#endif // SEASONAL 
	void deleteOccupancy(
		int		// no. of rows (as above)
	);

	bool outRangeHeaders( // Open range file and write header record
		Species*,	// pointer to Species
		int				// Landscape number (-999 to close the file)
	);
#if RS_ABC
	void outRange( // Write record to range file
		Species*, 	// pointer to Species
		int,				// replicate
		int,				// year
		int,				// generation
		ABCmaster*,	// pointer to ABC master object
		bool				// TRUE if ABC observations in current year
	);
#else
	void outRange( // Write record to range file
		Species*, // pointer to Species
		int,			// replicate
		int,			// year
		int				// generation
	);
#endif
	bool outPopHeaders( // Open population file and write header record
		Species*, // pointer to Species
		int       // option: -999 to close the file
	);
#if RS_ABC
	void outPop( // Write records to population file
		Species*, 	// pointer to Species
		int,				// replicate
		int,				// year
		int,				// generation
		ABCmaster*,	// pointer to ABC master object
		bool,				// TRUE if ABC observations in current year
		bool				// TRUE if normal population output in current year
	);
#else
	void outPop( // Write records to population file
		int,	// replicate
		int,	// year
		int		// generation
	);
#endif

#if RS_CONTAIN
	bool outCullHeaders( // Open cull file and write header record
		Species*, // pointer to Species
		int       // option: -999 to close the file
	);
	void outCull( // Write records to cull file
		int,	// replicate
		int,	// year
		int		// generation
	);
#endif // RS_CONTAIN 

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
#if RS_ABC
	void outABCpreds( // Write predictions for ABC
//		Species*, 	// pointer to Species
		int,				// replicate
		int,				// year
//		int,				// generation
		ABCmaster*,	// pointer to ABC master object
		float				// landscape resolution
	);
#endif
#if PEDIGREE
	bool outGroupHeaders( // Open groups file and write header record
		int       // option: -999 to close the file
	);
#endif
#if RS_RCPP && !R_CMD
    Rcpp::IntegerMatrix addYearToPopList(int,int);
#endif

private:
	Landscape *pLandscape;
	int indIx;				// index used to apply initial individuals
	float **occSuit;	// occupancy of suitable cells / patches
	std::vector <SubCommunity*> subComms;

};

extern paramSim *paramsSim;
extern paramInit *paramsInit;

#if VCL
extern bool stopRun;
#endif


//---------------------------------------------------------------------------
#endif
