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

RangeShifter v2.0 Individual

Implements the Individual class

Various optional attributes (genes for traits, movement parameters, etc.) are
allocated dynamically and accessed by pointers if required.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 28 July 2021 by Greta Bocedi

------------------------------------------------------------------------------*/

#ifndef IndividualH
#define IndividualH


#include <queue>
#include <algorithm>
using namespace std;

//#include "mathlib.h"
#include "Parameters.h"
#include "Species.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "Genome.h"

#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN 

#define NODATACOST 100000 // cost to use in place of nodata value for SMS
#define ABSNODATACOST 100 // cost to use in place of nodata value for SMS
													// when boundaries are absorbing

//---------------------------------------------------------------------------

struct indStats {
	short stage; short sex; short age; short status; short fallow;
#if SEASONAL
#if PARTMIGRN
	short migrnstatus;
#endif // PARTMIGRN 
#endif
	bool isDeveloping;
#if GOBYMODEL
	bool asocial;
#endif
#if SOCIALMODEL
	bool asocial;
#endif
};
struct pathData { // to hold path data common to SMS and CRW models
	int year, total, out; // nos. of steps
#if SEASONAL
	int season;
#endif
	Patch* pSettPatch;		// pointer to most recent patch tested for settlement
	short settleStatus; 	// whether ind may settle in current patch
												// 0 = not set, 1 = debarred through density dependence rule
												// 2 = OK to settle subject to finding a mate
//	bool leftNatalPatch;	// individual has moved out of its natal patch
#if RS_RCPP
	short pathoutput;
#endif
};
struct pathSteps { // nos. of steps for movement model
	int year, total, out;
#if SEASONAL
	int season;
#endif
};
struct settlePatch {
	Patch* pSettPatch; short settleStatus;
};
struct crwParams { // to hold data for CRW movement model
	float prevdrn;	// direction of previous step (UNITS)
	float xc,yc;		// continuous cell co-ordinates
	float stepL;		// phenotypic step length (m)
	float rho;			// phenotypic step correlation coefficient
};
struct array3x3d { double cell[3][3]; };
struct movedata { float dist; float cost; };
struct smsdata {
	locn prev;			// location of previous cell
	locn goal;			// location of goal
	float dp;				// directional persistence
	float gb;				// goal bias
	float alphaDB;	// dispersal bias decay rate
	int betaDB;			// dispersal bias decay inflection point (no. of steps)
#if PARTMIGRN
	short goalType;
#endif  // PARTMIGRN 
};

#if SEASONAL
struct patchlist { Patch *pPatch; short season; bool breeding; bool fixed; };
#endif

class Individual {

public:
	static int indCounter; // used to create ID, held by class, not members of class
#if GROUPDISP
	Individual(); // Default constructor
#endif
#if RS_CONTAIN
	Individual( // Individual constructor
		Cell*,	// pointer to Cell
		Patch*,	// pointer to patch
		short,	// stage
		short,	// age
		short,	// reproduction interval (no. of years/seasons between breeding attempts)
		short,	// mother's stage
		float,	// probability that sex is male
		bool,		// TRUE for a movement model, FALSE for kernel-based transfer
		short		// movement type: 1 = SMS, 2 = CRW
	);
#else
#if PARTMIGRN
	Individual( // Individual constructor
		Species*,	// pointer to Species
		Cell*,		// pointer to Cell
		Patch*,		// pointer to patch
		short,		// stage
		short,		// age
		short,		// reproduction interval (no. of years/seasons between breeding attempts)
		float,		// probability that sex is male
		bool,			// TRUE for a movement model, FALSE for kernel-based transfer
		short			// movement type: 1 = SMS, 2 = CRW
	);
#else
	Individual( // Individual constructor
		Cell*,	// pointer to Cell
		Patch*,	// pointer to patch
		short,	// stage
		short,	// age
		short,	// reproduction interval (no. of years/seasons between breeding attempts)
		float,	// probability that sex is male
		bool,		// TRUE for a movement model, FALSE for kernel-based transfer
		short		// movement type: 1 = SMS, 2 = CRW
	);
#endif // PARTMIGRN 
#endif // RS_CONTAIN 
	~Individual(void);
#if BUTTERFLYDISP
void setMate(Individual*);
//void setMated(short,Individual*);
//int getNJuvs(void);
Individual* getMate(void);
#endif
	void setGenes( // Set genes for individual variation from species initialisation parameters
		Species*,			// pointer to Species
		int						// Landscape resolution
	);
	void setGenes( // Inherit genome from parents
		Species*,			// pointer to Species
		Individual*,	// pointer to mother
		Individual*,	// pointer to father (must be 0 for an asexual Species)
		int						// Landscape resolution
	);
#if VIRTUALECOLOGIST
	Genome* getGenome(void);
#endif
	void setEmigTraits( // Set phenotypic emigration traits
		Species*,	// pointer to Species
		short,		// location of emigration genes on genome
		short,		// number of emigration genes
		bool			// TRUE if emigration is sex-dependent
	);
	emigTraits getEmigTraits(void); // Get phenotypic emigration traits

	void setKernTraits( // Set phenotypic transfer by kernel traits
		Species*,	// pointer to Species
		short,		// location of kernel genes on genome
		short,		// number of kernel genes
		int,			// Landscape resolution
		bool			// TRUE if transfer is sex-dependent
	);
	trfrKernTraits getKernTraits(void); // Get phenotypic transfer by kernel traits

	void setSMSTraits( // Set phenotypic transfer by SMS traits
		Species*,	// pointer to Species
		short,		// location of SMS genes on genome
		short,		// number of SMS genes
		bool			// TRUE if transfer is sex-dependent
	);
	trfrSMSTraits getSMSTraits(void); // Get phenotypic transfer by SMS traits
	void setCRWTraits( // Set phenotypic transfer by CRW traits
		Species*,	// pointer to Species
		short,		// location of CRW genes on genome
		short,		// number of CRW genes
		bool			// TRUE if transfer is sex-dependent
	);
	trfrCRWTraits getCRWTraits(void); // Get phenotypic transfer by CRW traits

	void setSettTraits( // Set phenotypic settlement traits
		Species*,	// pointer to Species
		short,		// location of settlement genes on genome
		short,		// number of settlement genes
		bool			// TRUE if settlement is sex-dependent
	);
	settleTraits getSettTraits(void); // Get phenotypic settlement traits

	// Identify whether an individual is a potentially breeding female -
	// if so, return her stage, otherwise return 0
	int breedingFem(void);
	int getId(void);
#if GROUPDISP
	int getParentId(short);
#if PEDIGREE
	Individual* getParent(short);
#endif // PEDIGREE
	void setGroupId(int);
	int getGroupId(void);
#endif // GROUPDISP 
#if PEDIGREE
	void setMatPosn(unsigned int);	// Set position in relatedness matrix
	unsigned int getMatPosn(void);	// Get position in relatedness matrix
#endif
	int getSex(void);
	int getStatus(void);
#if SEASONAL
#if PARTMIGRN
	int getMigrnStatus(void);
#endif // PARTMIGRN 
#endif
	indStats getStats(void);
#if GOBYMODEL
	bool isAsocial(void);
#endif
#if SOCIALMODEL
	bool isAsocial(void);
#endif
	Cell* getLocn( // Return location (as pointer to Cell)
		const short	// option: 0 = get natal locn, 1 = get current locn
	); //
	Patch* getNatalPatch(void);
	void setYearSteps(int);
	pathSteps getSteps(void);
	settlePatch getSettPatch(void);
	void setSettPatch(const settlePatch);
#if SEASONAL
	void resetPathSeason(void);
#endif
	void setStatus(short);
#if SEASONAL
#if PARTMIGRN
	void setMigrnStatus(short);
#endif // PARTMIGRN 
	void setPrevPatch(Patch*);
#if PARTMIGRN
//	void setNpatches( // set max. no. of patches to be held in memory
//		const short		// no. of patches
//	);
	void addPatch( // add patch to memory
		patchlist				// patch data
	);
	patchlist getPatch( // get specified patch from memory
		const int		// patch list no.
	);
	void setGoal( // set goal type and location
		const locn,			// goal location
		const short,		// goal type
		const bool			// breeding season
	);
#endif // PARTMIGRN 
#endif // SEASONAL
	void developing(void);
	void develop(void);
	void ageIncrement( // Age by one year
		short	// maximum age - if exceeded, the Individual dies
	);
	void incFallow(void); // Inrement no. of reproductive seasons since last reproduction
	void resetFallow(void);
	void moveto( // Move to a specified neighbouring cell
		Cell*	// pointer to the new cell
	);
#if GROUPDISP
	void moveTo( // Move to any specified cell
		Cell*	// pointer to the new cell
	);
#endif
	// Move to a new cell by sampling a dispersal distance from a single or double
	// negative exponential kernel
	// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
#if SEASONAL
	int moveKernel(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// reproduction type (see Species)
		const short,	// season
		const bool    // absorbing boundaries?
	);
#else
	int moveKernel(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// reproduction type (see Species)
		const bool    // absorbing boundaries?
	);
#endif // SEASONAL 
	// Make a single movement step according to a mechanistic movement model
	// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
#if SEASONAL
	int moveStep(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// landscape change index
		const short,	// season
		const bool    // absorbing boundaries?
	);
#else
	int moveStep(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// landscape change index
		const bool    // absorbing boundaries?
	);
#endif // SEASONAL 
	void drawMove(	// Visualise paths resulting from movement simulation model
									// NULL for the batch version
		const float,	// initial x co-ordinate
		const float,	// initial y co-ordinate
		const float,	// final x co-ordinate
		const float		// final y co-ordinate
	);
	movedata smsMove( // Move to a neighbouring cell according to the SMS algorithm
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// landscape change index
		const bool,		// TRUE if still in (or returned to) natal patch
		const bool,   // individual variability?
		const bool    // absorbing boundaries?
	);
	array3x3d getSimDir( // Weight neighbouring cells on basis of current movement direction
		const int,	// current x co-ordinate
		const int,	// current y co-ordinate
		const float	// directional persistence value
	);
	array3x3d getGoalBias( // Weight neighbouring cells on basis of goal bias
		const int,	// current x co-ordinate
		const int,	// current y co-ordinate
		const int,	// goal type: 0 = none, 1 = towards goal (NOT IMPLEMENTED), 2 = dispersal bias
		const float	// GOAL BIAS VALUE
	);
	array3x3d calcWeightings( // Calculate weightings for neighbouring cells
		const double,	// base for power-law (directional persistence or goal bias value)
		const double	// direction in which lowest (unit) weighting is to be applied
	);
	array3x3f getHabMatrix( // Weight neighbouring cells on basis of (habitat) costs
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const int,		// current x co-ordinate
		const int,		// current y co-ordinate
		const short,	// perceptual range (cells)
		const short,	// perceptual range evaluation method (see Species)
		const short,	// landscape change index
		const bool    // absorbing boundaries?
	);
#if GROUPDISP || ROBFITT
	void outGenetics( // Write records to genetics file
		const int,		 	// replicate
		const int,		 	// year
		const int,		 	// species number
		const int,	 	 	// landscape number
		const bool,  		// patch-based landscape
		const bool	 		// output as cross table?
	);
#else
	void outGenetics( // Write records to genetics file
		const int,		 	// replicate
		const int,		 	// year
		const int,		 	// species number
		const int,		 	// landscape number
		const bool	 		// output as cross table?
	);
#endif
#if RS_RCPP
	void outMovePath( // Write records to movement paths file
		const int		 	// year
	);
#endif

#if GROUPDISP

protected:
	short stage;
	short sex;
	Cell *pPrevCell;						// pointer to previous Cell
	Cell *pCurrCell;						// pointer to current Cell
	Patch *pNatalPatch;					// pointer to natal Patch
	short status;	// 0 = initial status in natal patch / philopatric recruit
								// 1 = disperser
								// 2 = disperser awaiting settlement in possible suitable patch
								// 3 = waiting between dispersal events
								// 4 = completed settlement
								// 5 = completed settlement in a suitable neighbouring cell
								// 6 = died during transfer by failing to find a suitable patch
								//     (includes exceeding maximum number of steps or crossing
								//			absorbing boundary)
								// 7 = died during transfer by constant, step-dependent,
								//     habitat-dependent or distance-dependent mortality
								// 8 = failed to survive annual (demographic) mortality
								// 9 = exceeded maximum age
	pathData *path; 						// pointer to path data for movement model
	crwParams *crw;     				// pointer to CRW traits and data
	smsdata *smsData;						// pointer to variables required for SMS
	emigTraits *emigtraits;			// pointer to emigration traits
	trfrKernTraits *kerntraits;	// pointers to transfer by kernel traits
	settleTraits *setttraits;		// pointer to settlement traits
	Genome *pGenome;

private:
	int indId;
#if PEDIGREE
	Individual *pParent[2];
#endif
	int parentId[2];
	int groupId;
#if PEDIGREE
	unsigned int matPosn;	// position in relatedness matrix
#endif
	short age;
	short fallow; // reproductive seasons since last reproduction
	bool isDeveloping;
#if GOBYMODEL
	bool asocial;
#endif
#if SOCIALMODEL
	bool asocial;
#endif
	std::queue <locn> memory;		// memory of last N squares visited for SMS

#else

private:
	int indId;
	short stage;
	short sex;
	short age;
	short status;	// 0 = initial status in natal patch / philopatric recruit
								// 1 = disperser
								// 2 = disperser awaiting settlement in possible suitable patch
								// 3 = waiting between dispersal events
								// 4 = completed settlement
								// 5 = completed settlement in a suitable neighbouring cell
								// 6 = died during transfer by failing to find a suitable patch
								//     (includes exceeding maximum number of steps or crossing
								//			absorbing boundary)
								// 7 = died during transfer by constant, step-dependent,
								//     habitat-dependent or distance-dependent mortality
								// 8 = failed to survive annual (demographic) mortality
								// 9 = exceeded maximum age
#if RS_CONTAIN
	short motherstage;	// mother's stage for purpose of implementing WALD kernel
#endif // RS_CONTAIN 
#if SEASONAL
#if PARTMIGRN
	short migrnstatus;	// 0 = not yet determined
											// 1 = philopatric resident, breeding and non-breeding fixed
											// 2 = philopatric migrant, non-breeding site fixed
											// 3 = philopatric migrant, non-breeding site not fixed
											// 4 = dispersed resident, breeding and non-breeding fixed
											// 5 = dispersed migrant, breeding fixed, non-breeding fixed
											// 6 = dispersed migrant, breeding fixed, non-breeding not fixed
#endif // PARTMIGRN 
	short npatches; 		// no. of patches held in memory
#endif // SEASONAL 
	short fallow; // reproductive seasons since last reproduction
#if BUTTERFLYDISP
//	short nJuvs;
	Individual *pMate;
#endif
	bool isDeveloping;
#if GOBYMODEL
	bool asocial;
#endif
#if SOCIALMODEL
	bool asocial;
#endif
	Cell *pPrevCell;						// pointer to previous Cell
	Cell *pCurrCell;						// pointer to current Cell
	Patch *pNatalPatch;					// pointer to natal Patch
#if SEASONAL
	Patch *pPrevPatch;						// pointer to previous Patch
#endif
	emigTraits *emigtraits;			// pointer to emigration traits
	trfrKernTraits *kerntraits;	// pointers to transfer by kernel traits
	pathData *path; 						// pointer to path data for movement model
	crwParams *crw;     				// pointer to CRW traits and data
	smsdata *smsData;						// pointer to variables required for SMS
	settleTraits *setttraits;		// pointer to settlement traits
	std::queue <locn> memory;		// memory of last N squares visited for SMS
#if SEASONAL
	std::vector <patchlist> patches;		// memory of patches used
#endif

	Genome *pGenome;

#endif // GROUPDISP

};


//---------------------------------------------------------------------------

double cauchy(double location, double scale) ;
double wrpcauchy (double location, double rho = exp(double(-1)));

extern RSrandom *pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

#if RS_RCPP
extern ofstream outMovePaths;
#endif

//---------------------------------------------------------------------------
#endif
