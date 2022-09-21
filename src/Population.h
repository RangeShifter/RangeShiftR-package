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

RangeShifter v2.0 Population

Implements the Population class

There is ONE instance of a Population for each Species within each SubCommunity
(including the matrix). The Population holds a list of all the Individuals in
the Population.

The matrix Population(s) hold(s) Individuals which are currently in the process
of transfer through the matrix.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe?er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species? responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 25 June 2021 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef PopulationH
#define PopulationH

#include <vector>
#include <algorithm>
//#include <Math.hpp>
//#include <math.h>
//#include <stdlib.h>
using namespace std;

#include "Parameters.h"
#if GROUPDISP
#include "Group.h"
#endif
#include "Individual.h"
#include "Species.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#if RS_ABC
#include "ABC.h"
#endif
#if PEDIGREE
#include "Pedigree.h"
#endif
#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN

//---------------------------------------------------------------------------

struct popStats {
	Species *pSpecies; Patch *pPatch; int spNum,nInds,nNonJuvs,nAdults; bool breeding;
#if GOBYMODEL
	int nSocial,nAsocial;
#endif
#if SOCIALMODEL
	int nSocial,nAsocial;
#endif
};
#if GROUPDISP
struct dispgroup {
	Group *pGroup; Cell *pCell; int groupsize; bool yes;
};
#endif
#if RS_ABC
struct dispstats {
	int nPhilo,nDisp,nSuccess; double sumDist,sumDist2;
};
#endif
struct disperser {
	Individual *pInd; Cell *pCell; bool yes;
};
struct traitsums { // sums of trait genes for dispersal
	int ninds[NSEXES];				// no. of individuals
	double sumD0[NSEXES];			// sum of maximum emigration probability
	double ssqD0[NSEXES];			// sum of squares of maximum emigration probability
	double sumAlpha[NSEXES];	// sum of slope of emigration dens-dep reaction norm
	double ssqAlpha[NSEXES];	// sum of squares of slope of emigration den-dep reaction norm
	double sumBeta[NSEXES]; 	// sum of inflection point of emigration reaction norm
	double ssqBeta[NSEXES]; 	// sum of squares of inflection point of emigration reaction norm
	double sumDist1[NSEXES]; 	// sum of kernel I mean
	double ssqDist1[NSEXES]; 	// sum of squares of kernel I mean
	double sumDist2[NSEXES]; 	// sum of kernel II mean
	double ssqDist2[NSEXES]; 	// sum of squares of kernel II mean
	double sumProp1[NSEXES]; 	// sum of propn using kernel I
	double ssqProp1[NSEXES]; 	// sum of squares of propn using kernel I
	double sumDP[NSEXES]; 		// sum of SMS directional persistence
	double ssqDP[NSEXES]; 		// sum of squares of SMS directional persistence
	double sumGB[NSEXES]; 		// sum of SMS goal bias
	double ssqGB[NSEXES]; 		// sum of squares of SMS goal bias
	double sumAlphaDB[NSEXES];	// sum of SMS dispersal bias decay rate
	double ssqAlphaDB[NSEXES]; 	// sum of squares of SMS dispersal bias decay rate
	double sumBetaDB[NSEXES];		// sum of SMS dispersal bias decay infl. pt.
	double ssqBetaDB[NSEXES]; 	// sum of squares of SMS dispersal bias decay infl. pt.
	double sumStepL[NSEXES]; 	// sum of CRW step length
	double ssqStepL[NSEXES]; 	// sum of squares of CRW step length
	double sumRho[NSEXES]; 		// sum of CRW correlation coefficient
	double ssqRho[NSEXES]; 		// sum of squares of CRW correlation coefficient
	double sumS0[NSEXES];			// sum of maximum settlement probability
	double ssqS0[NSEXES];			// sum of squares of maximum settlement probability
	double sumAlphaS[NSEXES];	// sum of slope of settlement den-dep reaction norm
	double ssqAlphaS[NSEXES];	// sum of squares of slope of settlement den-dep reaction norm
	double sumBetaS[NSEXES]; 	// sum of inflection point of settlement reaction norm
	double ssqBetaS[NSEXES]; 	// sum of squares of inflection point of settlement reaction norm
};

class Population {

public:
	Population(void); // default constructor
#if PEDIGREE
	Population( // constructor for a Population of a specified size
		Species*,		// pointer to Species
		Pedigree*,	// pointer to Pedigree
		Patch*,			// pointer to Patch
		int,				// no. of Individuals
		int					// Landscape resolution
	);
#else
	Population( // constructor for a Population of a specified size
		Species*,	// pointer to Species
		Patch*,		// pointer to Patch
		int,			// no. of Individuals
		int				// Landscape resolution
	);
#endif
	~Population(void);
	traitsums getTraits(Species*);
#if RS_CONTAIN
	popStats getStats(
		short		// habitat index
	);
#else
#if SPATIALDEMOG
	popStats getStats(
			std::vector <float>
	);
#else
	popStats getStats(void);
#endif // SPATIALDEMOG
#endif // RS_CONTAIN
	Species* getSpecies(void);
#if GROUPDISP
	int getNGroups(void);
#endif
	int getNInds(void);
#if VIRTUALECOLOGIST
	Individual* getInd(int);
#endif
	int totalPop(void);
	int stagePop( // return no. of Individuals in a specified stage
		int	// stage
	);
#if RS_ABC
	int sexPop( // return no. of Individuals of a specified sex
		short		// sex
	);
	int stageSexPop( // return no. of Individuals in a specified stage/sex
		short,	// stage
		short		// sex
	);
#endif // RS_ABC
	void extirpate(void); // Remove all individuals
#if RS_CONTAIN
	void cull(Cull*,double); // Remove individuals according to cull rate
	void resetCull(void);
#endif // RS_CONTAIN
#if RS_CONTAIN
#if SEASONAL
	void reproduction(
		const int,		// habitat index
		const int,		// season
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
#else
	void reproduction(
		const int,		// habitat index
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
#endif // SEASONAL
#else
#if SEASONAL
	void reproduction(
		const int,		// season
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
#else
#if GROUPDISP
	Individual* getFather(int,int);
	void reproduction(
		const std::vector <Individual*> *, // pointer to list of global 'fathers'
		const int,		// length of global 'fathers' vector
		const std::vector <Individual*> *, // pointer to list of neighbourhood 'fathers'
		const int,		// length of neighbourhood 'fathers' vector
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
#else
#if BUTTERFLYDISP
	void reproduction(
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int,		// Landscape resolution
		const short		// option: 0 = default (all reproduction before dispersal),
									// 1 = mating only (before dispersal),
									// 2 = parturition only (after dispersal)
	);
#else
#if SPATIALDEMOG
	void reproduction(
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int,		// Landscape resolution
		std::vector <float>    // local demographic scaling
	);
#else
	void reproduction(
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
#endif // SPATIALDEMOG
#endif // BUTTERFLYDISP
#endif // GROUPDISP
#endif // SEASONAL
#endif // RS_CONTAIN
	// Following reproduction of ALL species, add juveniles to the population
	void fledge(void);
#if SEASONAL
	void emigration( // Determine which individuals will disperse
		float,  // local carrying capacity
		short	  // season
	);
#else
	void emigration( // Determine which individuals will disperse
		float   // local carrying capacity
	);
#endif
	void allEmigrate(void); // All individuals emigrate after patch destruction
	// If an individual has been identified as an emigrant, remove it from the Population
	disperser extractDisperser(
		int		// index no. to the Individual in the inds vector
	);
#if GROUPDISP
	// For a group identified as being in the matrix population:
	// if it is a settler, return its new location and remove it from the current population
	// otherwise, leave it in the matrix population for possible reporting before deletion
	dispgroup extractGroupSettler(
		int   // index no. to the Group in the group vector
	);
//	void setGroupStatus(int,short);
	void deleteGroup(int);
#endif
	// For an individual identified as being in the matrix population:
	// if it is a settler, return its new location and remove it from the current population
	// otherwise, leave it in the matrix population for possible reporting before deletion
	disperser extractSettler(
		int   // index no. to the Individual in the inds vector
	);
#if GROUPDISP
	void recruit( // Add a specified individual to the population
		Patch*,				// pointer to Patch
		Individual*,	// pointer to Individual
		bool,					// TRUE if movement model
		short,				// 1 = SMS, 2 = CRW
		bool					// new dispersal group
	);
#endif
	void recruit( // Add a specified individual to the population
		Individual*	// pointer to Individual
	);
#if GROUPDISP
	int grouptransfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short				// landscape change index
	);
	void deleteGroups(void); // Delete dispersal groups once dispersal has finished
#if PEDIGREE
	void outGroups(Pedigree*,int,int,int,bool);
#endif
#endif // GROUPDISP
#if SEASONAL
	int transfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short,			// landscape change index
		short				// season
	);
	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool matePresent(
		Cell*,	// pointer to the Cell which the potential settler has reached
		short,	// sex of the required mate (0 = female, 1 = male)
		short		// season
	);
#else
#if RS_RCPP
	int transfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short,				// landscape change index
		short				// year
	);
	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool matePresent(
		Cell*,	// pointer to the Cell which the potential settler has reached
		short		// sex of the required mate (0 = female, 1 = male)
	);
#else
	int transfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short				// landscape change index
	);
	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool matePresent(
		Cell*,	// pointer to the Cell which the potential settler has reached
		short		// sex of the required mate (0 = female, 1 = male)
	);
#endif // RS_RCPP
#endif // SEASONAL
	// Determine survival and development and record in individual's status code
	// Changes are NOT applied to the Population at this stage
#if RS_CONTAIN
#if SEASONAL
	void survival0(
		float,	// local carrying capacity
		short,	// habitat index
		short,	// season
		short,	// option0:	0 - stage 0 (juveniles) only
						//	  			1 - all stages
						//					2 - stage 1 and above (all non-juveniles)
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
						//	  	 		2 - survival only (when survival is annual)
	);
#else
	void survival0(
		float,	// local carrying capacity
		short,	// habitat index
		short,	// option0:	0 - stage 0 (juveniles) only
						//	  			1 - all stages
						//					2 - stage 1 and above (all non-juveniles)
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
						//	  	 		2 - survival only (when survival is annual)
	);
#endif // SEASONAL
#else
#if SEASONAL
	void survival0(
		float,	// local carrying capacity
		short,	// season
		short,	// option0:	0 - stage 0 (juveniles) only
						//	  			1 - all stages
						//					2 - stage 1 and above (all non-juveniles)
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
						//	  	 		2 - survival only (when survival is annual)
	);
#else
#if SPATIALMORT
	void survival0(
		float,	// local carrying capacity
		float,	// spatially varying mortality
		short,	// option0:	0 - stage 0 (juveniles) only
						//	  			1 - all stages
						//					2 - stage 1 and above (all non-juveniles)
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
						//	  	 		2 - survival only (when survival is annual)
	);
#else
#if PEDIGREE
	void survival0(
		Pedigree*,	// pointyer to Pedigree
		float,			// local carrying capacity
		short,			// option0:	0 - stage 0 (juveniles) only
								//	  			1 - all stages
								//					2 - stage 1 and above (all non-juveniles)
		short 			// option1:	0 - development only (when survival is annual)
								//	  	 		1 - development and survival
								//	  	 		2 - survival only (when survival is annual)
	);
#else
#if SPATIALDEMOG
	void survival0(
			float,	// local carrying capacity
			short,	// option0:	0 - stage 0 (juveniles) only
							//	  			1 - all stages
							//					2 - stage 1 and above (all non-juveniles)
			short, 	// option1:	0 - development only (when survival is annual)
							//	  	 		1 - development and survival
							//	  	 		2 - survival only (when survival is annual)
			std::vector <float> // local demographic scaling
		);
#else
	void survival0(
		float,	// local carrying capacity
		short,	// option0:	0 - stage 0 (juveniles) only
						//	  			1 - all stages
						//					2 - stage 1 and above (all non-juveniles)
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
						//	  	 		2 - survival only (when survival is annual)
	);
#endif // SPATIALDEMOG
#endif // PEDIGREE
#endif // SPATIALMORT
#endif // SEASONAL
#endif // RS_CONTAIN
#if PEDIGREE
	void survival1(				// Apply survival changes to the population
		Pedigree*
	);
#else
	void survival1(void); // Apply survival changes to the population
#endif
#if SEASONAL && PARTMIGRN
	void extremeEvent(float);
#endif // SEASONAL && PARTMIGRN
	void ageIncrement(void);
	bool outPopHeaders( // Open population file and write header record
		int,	// Landscape number (-999 to close the file)
		bool	// TRUE for a patch-based model, FALSE for a cell-based model
	);
#if RS_ABC
	void outPopulation( // Write record to population file
		int,				// replicate
		int,				// year
		int,				// generation
		float,			// epsilon - global stochasticity value
		bool,				// TRUE for a patch-based model, FALSE for a cell-based model
		bool,				// TRUE to write environmental data
		bool,				// TRUE if there is a gradient in carrying capacity
		bool,				// TRUE if ABC observations in current year
		ABCmaster*	// pointer to ABC master object
	);
#else
	void outPopulation( // Write record to population file
		int,		// replicate
		int,		// year
		int,		// generation
		float,	// epsilon - global stochasticity value
		bool,		// TRUE for a patch-based model, FALSE for a cell-based model
		bool,		// TRUE to write environmental data
		bool		// TRUE if there is a gradient in carrying capacity
	);
#endif // RS_ABC

#if RS_CONTAIN

	bool outCullHeaders( // Open cull file and write header record
		Landscape*,	// pointer to Landscape
		int,				// Landscape number (-999 to close the file)
		bool				// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void outCullData( // Write record to cull file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int,				// generation
		bool				// TRUE for a patch-based model, FALSE for a cell-based model
	);

#endif // RS_CONTAIN

	void outIndsHeaders( // Open individuals file and write header record
		int,	// replicate
		int,	// Landscape number (-999 to close the file)
		bool	// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void outIndividual( // Write records to individuals file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int,				// generation
		int					// Patch number
	);
#if RS_ABC
dispstats getDispStats(float);
#endif
#if GROUPDISP || ROBFITT
	void outGenetics( // Write records to genetics file
		const int,		// replicate
		const int,		// year
		const int,		// landscape number
		const bool		// patch-based landscape
	);
	void groupclean(void); // Remove zero pointers to dead or dispersed individuals
#else
	void outGenetics( // Write records to genetics file
		const int,		// replicate
		const int,		// year
		const int	 		// landscape number
	);
#endif
	void clean(void); // Remove zero pointers to dead or dispersed individuals

private:
	short nStages;
	short nSexes;
	Species *pSpecies;	// pointer to the species
	Patch *pPatch;			// pointer to the patch
	int nInds[NSTAGES][NSEXES];		// no. of individuals in each stage/sex

	std::vector <Individual*> inds; // all individuals in population except ...
	std::vector <Individual*> juvs; // ... juveniles until reproduction of ALL species
																	// has been completed
#if GROUPDISP
	std::vector <Group*> groups;		// dispersal groups (matrix only)
	Group *pGroup;			// pointer to the current group
#endif

#if RS_CONTAIN
	int nCulled;				// no. of individuals culled
	bool selectedForCull;
#endif // RS_CONTAIN

};

//---------------------------------------------------------------------------

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;
extern paramSim *paramsSim;
extern RSrandom *pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

//---------------------------------------------------------------------------
#endif

