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
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 25 June 2021 by Greta Bocedi

------------------------------------------------------------------------------*/

#ifndef SubCommunityH
#define SubCommunityH

#if VCL
#include <VCLTee.Chart.hpp>
#endif

#include <vector>
#include <algorithm>
using namespace std;

#include "Version.h"
#include "Parameters.h"
#include "Landscape.h"
#include "Population.h"
#if PEDIGREE
#include "Pedigree.h"
#endif
#if ROBFITT
#include "PatchSelection.h"
#endif
#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN 

//---------------------------------------------------------------------------

struct traitCanvas { // canvases for drawing variable traits
#if VCL
	TCanvas *pcanvas[NTRAITS];
#else
	int *pcanvas[NTRAITS]; // dummy variables for batch version
#endif
};

class SubCommunity {

public:
	SubCommunity(Patch*,int);
	~SubCommunity(void);
	intptr getNum(void);
	Patch* getPatch(void);
	locn getLocn(void);
#if RS_CONTAIN
	void setHabIndex(
		Species*,	// pointer to Species
		short,		// landscape raster type
		short			// landscape change index
	);
#endif // RS_CONTAIN 

	// functions to manage populations occurring in the SubCommunity
	popStats getPopStats(void);
	void setInitial(bool);
#if PEDIGREE
	void initialise(Landscape*,Species*,Pedigree*);
	void initialInd(Landscape*,Species*,Pedigree*,Patch*,Cell*,int);
	Population* newPopn( // Create a new population, and return its address
		Landscape*,	// pointer to Landscape
		Species*,		// pointer to Species
		Pedigree*,	// pointer to Pedigree
		Patch*,			// pointer to Patch
		int					// no. of Individuals
	);
#else
	void initialise(Landscape*,Species*);
	void initialInd(Landscape*,Species*,Patch*,Cell*,int);
	Population* newPopn( // Create a new population, and return its address
		Landscape*,	// pointer to Landscape
		Species*,		// pointer to Species
		Patch*,			// pointer to Patch
		int					// no. of Individuals
	);
#endif
	void resetPopns(void);
	void resetPossSettlers(void);
	void localExtinction( // Extirpate all populations
		int	// option: 	0 - random local extinction probability
				//					1 - local extinction probability gradient
	);
	void patchChange(void);
#if SEASONAL
	void reproduction(
		int,		// Landscape resolution
		float,	// epsilon - global stochasticity value
		short,	// season
		short,	// raster type (see Landscape)
		bool		// TRUE for a patch-based model, FALSE for a cell-based model
	);
#else
#if GROUPDISP
	Individual* getFather(int,int);
	void reproduction(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const int,		// minimum breeding stage
		const std::vector <Individual*> *, // pointer to list of global 'fathers'
		const int,		// length of vector
		const int,		// Landscape resolution
		const locn,		// minimum breeding range
		const locn,		// maximum breeding range
		const float,	// epsilon - global stochasticity value
		const short,	// raster type (see Landscape)
		const bool		// TRUE for a patch-based model, FALSE for a cell-based model
	);
#else
#if BUTTERFLYDISP
	void reproduction(
		int,		// Landscape resolution
		float,	// epsilon - global stochasticity value
		short,	// dispersal timing: 0 = during reprodn, 1 = after reprodn
		short,	// option: 0 = default (all reproduction before dispersal),
						// 1 = mating only (before dispersal),
						// 2 = parturition only (after dispersal)
		short,	// raster type (see Landscape)
		bool		// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void fledge(void);
#else
	void reproduction(
		int,		// Landscape resolution
		float,	// epsilon - global stochasticity value 
		short,	// raster type (see Landscape)
		bool		// TRUE for a patch-based model, FALSE for a cell-based model
	);
#endif // BUTTERFLYDISP
#endif // GROUPDISP
#endif // SEASONAL
#if RS_DISEASE
	void emigration(
		Species*,		// pointer to Species
		short				// season
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
	// Remove emigrants from their natal patch and add to patch 0 (matrix)
	void initiateDispersal(
		SubCommunity*	// pointer to matrix SubCommunity
	);
// Add an individual into the local population of its species in the patch
#if GROUPDISP
	void recruit(
		Patch*,				// pointer to Patch
		Individual*,	// pointer to Individual
		Species*,			// pointer to Species
		bool					// new dispersal group
	);
#if PEDIGREE
	void outGroups(Pedigree*,int,int,int,bool);
#endif // PEDIGREE 

#endif // GROUPDISP 
	void recruit(
		Individual*,	// pointer to Individual
		Species*			// pointer to Species
	);
#if SEASONAL || RS_RCPP
	int transfer( // Transfer through matrix - run for matrix SubCommunity only
		Landscape*,	// pointer to Landscape
		short,			// landscape change index
		short				// season / year
	);
#else
	int transfer( // Transfer through matrix - run for matrix SubCommunity only
		Landscape*,	// pointer to Landscape
		short				// landscape change index
	);
#endif // SEASONAL || RS_RCPP
	// Remove emigrants from patch 0 (matrix) and transfer to SubCommunity in which
	// their destination co-ordinates fall (executed for the matrix patch only)
#if PEDIGREE
	void completeDispersal(
		Landscape*,	// pointer to Landscape
		Pedigree*,	// pointer to Pedigree
		bool				// TRUE to increment connectivity totals
	);
#else
	void completeDispersal(
		Landscape*,	// pointer to Landscape
		bool				// TRUE to increment connectivity totals
	);
#endif
#if GROUPDISP
	// Delete dispersal groups once dispersal has finished
	// This function is executed for the matrix patch only
	void deleteGroups(void);
#endif
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
#if SPATIALMORT
	void survival(
		short,	// part:		0 = determine survival & development,
						//		 			1 = apply survival changes to the population
		short,	// spatial mortality period (0 or 1)
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
		short 			// option1:	0 - development only (when survival is annual)
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
#endif // SPATIALMORT
#endif // SEASONAL
#if RS_CONTAIN
	short findCullTarget(Cull*,int,int,int);
	bool isCullTarget(void);
	int initialYear(void);
	double damageIndex(void);
	void resetCullTarget(void);
	void resetCull(void); 
//	void cullPatch(Cull*,int,float);
	void cullPatch(Cull*,int,int);
	void updateDamage(Landscape*,Species*,Cull*);
	int getCullCount(void);
	double prevDamage(void);
#endif // RS_CONTAIN 
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
#if RS_ABC
	void outPop( // Write records to population file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int,				// generation
		ABCmaster*,	// pointer to ABC master object
		bool,				// TRUE if ABC observations in current year
		bool				// TRUE if normal population output in current year
	);
#else
	void outPop( // Write records to population file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int					// generation
	);
#endif

#if RS_CONTAIN
	bool outCullHeaders( // Open cull file and write header record
		Landscape*,	// pointer to Landscape
		Species*,		// pointer to Species
		int					// option: -999 to close the file
	);
	void outCull( // Write records to cull file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int					// generation
	);
#endif // RS_CONTAIN 

	void outInds( // Write records to individuals file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int,				// generation
		int					// Landscape number (>= 0 to open the file, -999 to close the file
								//									 -1 to write data records)
	);
#if GROUPDISP || ROBFITT
	void outGenetics( // Write records to genetics file
		int,				// replicate
		int,				// year
		int,				// generation
		int,				// Landscape number (>= 0 to open the file, -999 to close the file
								//									 -1 to write data records)
		bool				// patch-based landscape
	);
#else
	void outGenetics( // Write records to genetics file
		int,				// replicate
		int,				// year
		int,				// generation
		int					// Landscape number (>= 0 to open the file, -999 to close the file
								//									 -1 to write data records)
	);
#endif
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
#if RS_ABC
dispstats getDispStats(float);
#endif

#if !LINUX_CLUSTER && VCL
	void draw( // Draw the SubCommunity on the landscape map - NULL for the batch version
		TCanvas*,		// pointer to canvas
		Landscape*	// pointer to Landscape
	);
#endif

private:
	intptr subCommNum;	// SubCommunity number
		// 0 is reserved for the SubCommunity in the inter-patch matrix
//	intptr *occupancy;	// pointer to occupancy array
	Patch *pPatch;
	int *occupancy;	// pointer to occupancy array
	std::vector <Population*> popns;
	bool initial; 	// WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES ...
#if RS_CONTAIN
	int habIndex;		// current habitat index for the patch in which the sub-community occurs
		// NOTE: the habitat is based on a randomly chosen cell within the patch and
		// therefore habitat-dependent demography should be applied only for a landscape
		// in which patches are homogeneous
#endif // RS_CONTAIN 
#if RS_CONTAIN
	bool cullTarget;
	int firstYear;	// first year qualified as a cull target
	int cullCount;	// estimated count for purpose of cull
#endif // RS_CONTAIN 

};

#if RS_CONTAIN
extern DamageParams *pDamageParams;
#endif // RS_CONTAIN 

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;

//---------------------------------------------------------------------------
#endif
