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

RangeShifter v2.0 SubCommunity

Implements the SubCommunity class

There is ONE instance of a SubCommunity for each Patch in the Landscape
(including the matrix). The SubCommunity holds a number of Populations, one for
each Species represented in the simulation.
CURRENTLY the number of Populations withn a SubCommunity is LIMITED TO ONE.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Peer G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 17 June 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef SubCommunityH
#define SubCommunityH

#include <vector>
#include <algorithm>
using namespace std;

#include "Version.h"
#include "Parameters.h"
#include "Landscape.h"
#include "Population.h"

//---------------------------------------------------------------------------

struct traitCanvas { // canvases for drawing variable traits
	int *pcanvas[NTRAITS]; // dummy variables for batch version
};

class SubCommunity {

public:
	SubCommunity(Patch*,int);
	~SubCommunity(void);
	int getNum(void);
	Patch* getPatch(void);
	locn getLocn(void);

	// functions to manage populations occurring in the SubCommunity
	popStats getPopStats(void);
	void setInitial(bool);
	void initialise(Landscape*,Species*);
	void initialInd(Landscape*,Species*,Patch*,Cell*,int);
	Population* newPopn( // Create a new population, and return its address
		Landscape*,	// pointer to Landscape
		Species*,		// pointer to Species
		Patch*,			// pointer to Patch
		int					// no. of Individuals
	);
	void resetPopns(void);
	void resetPossSettlers(void);
	void localExtinction( // Extirpate all populations
		int	// option: 	0 - random local extinction probability
				//					1 - local extinction probability gradient
	);
	void patchChange(void);
	void reproduction(
		int,		// Landscape resolution
		float,	// epsilon - global stochasticity value
		short,	// raster type (see Landscape)
		bool		// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void emigration(void);
	// Remove emigrants from their natal patch and add to patch 0 (matrix)
	void initiateDispersal(
		SubCommunity*	// pointer to matrix SubCommunity
	);
// Add an individual into the local population of its species in the patch
	void recruit(
		Individual*,	// pointer to Individual
		Species*			// pointer to Species
	);
	int transfer( // Transfer through matrix - run for matrix SubCommunity only
		Landscape*,	// pointer to Landscape
		short,			// landscape change index
		short				// season / year
	);
	// Remove emigrants from patch 0 (matrix) and transfer to SubCommunity in which
	// their destination co-ordinates fall (executed for the matrix patch only)
	void completeDispersal(
		Landscape*,	// pointer to Landscape
		bool				// TRUE to increment connectivity totals
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
	// Find the population of a given species in a given patch
	Population* findPop(Species*,Patch*);
	void createOccupancy(
		int	// no. of rows = (no. of years / interval) + 1
	);
	void updateOccupancy(
		int	// row = (no. of years / interval)
	);
	int getOccupancy(
		int	// row = (no. of years / interval)
	);
	void deleteOccupancy(void);

	bool outPopHeaders( // Open population file and write header record
		Landscape*,	// pointer to Landscape
		Species*,		// pointer to Species
		int					// option: -999 to close the file
	);
	void outPop( // Write records to population file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int					// generation
	);

	void outInds( // Write records to individuals file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int,				// generation
		int					// Landscape number (>= 0 to open the file, -999 to close the file
								//									 -1 to write data records)
	);
	void outGenetics( // Write records to genetics file
		int,				// replicate
		int,				// year
		int,				// generation
		int					// Landscape number (>= 0 to open the file, -999 to close the file
								//									 -1 to write data records)
	);
	bool outTraitsHeaders( // Open traits file and write header record
		Landscape*,	// pointer to Landscape
		Species*,		// pointer to Species
		int					// Landscape number (-999 to close the file)
	);
	traitsums outTraits( // Write records to traits file and return aggregated sums
		traitCanvas,	// pointers to canvases for drawing variable traits
									// in the batch version, these are replaced by integers set to zero
		Landscape*, 	// pointer to Landscape
		int,					// replicate
		int,					// year
		int,					// generation
		bool					// true if called to summarise data at community level
	);
	int stagePop( // Population size of a specified stage
		int	// stage
	);

private:
	intptr subCommNum;	// SubCommunity number
		// 0 is reserved for the SubCommunity in the inter-patch matrix
//	intptr *occupancy;	// pointer to occupancy array
	Patch *pPatch;
	int *occupancy;	// pointer to occupancy array
	std::vector <Population*> popns;
	bool initial; 	// WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES ...

};

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;

//---------------------------------------------------------------------------
#endif
