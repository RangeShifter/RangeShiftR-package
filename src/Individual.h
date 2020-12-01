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

Last updated: 27 November 2020 by Anne-Kathleen Malchow, Potsdam University

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

#define NODATACOST 100000 // cost to use in place of nodata value for SMS
#define ABSNODATACOST 100 // cost to use in place of nodata value for SMS
													// when boundaries are absorbing

//---------------------------------------------------------------------------

struct indStats {
	short stage; short sex; short age; short status; short fallow;
	bool isDeveloping;
};
struct pathData { // to hold path data common to SMS and CRW models
	int year, total, out; // nos. of steps
	Patch* pSettPatch;		// pointer to most recent patch tested for settlement
	short settleStatus; 	// whether ind may settle in current patch
												// 0 = not set, 1 = debarred through density dependence rule
												// 2 = OK to settle subject to finding a mate
//	bool leftNatalPatch;	// individual has moved out of its natal patch
	short pathoutput;
};
struct pathSteps { // nos. of steps for movement model
	int year, total, out;
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
struct array3x3d { float cell[3][3]; };
struct movedata { float dist; float cost; };
struct smsdata {
	locn prev;			// location of previous cell
	locn goal;			// location of goal
	float dp;				// directional persistence
	float gb;				// goal bias
	float alphaDB;	// dispersal bias decay rate
	int betaDB;			// dispersal bias decay inflection point (no. of steps)
};

class Individual {

public:
	static int indCounter; // used to create ID, held by class, not members of class
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
	~Individual(void);
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
	int getSex(void);
	int getStatus(void);
	indStats getStats(void);
	Cell* getLocn( // Return location (as pointer to Cell)
		const short	// option: 0 = get natal locn, 1 = get current locn
	); //
	Patch* getNatalPatch(void);
	void setYearSteps(int);
	pathSteps getSteps(void);
	settlePatch getSettPatch(void);
	void setSettPatch(const settlePatch);
	void setStatus(short);
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
	// Move to a new cell by sampling a dispersal distance from a single or double
	// negative exponential kernel
	// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
	int moveKernel(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// reproduction type (see Species)
		const bool    // absorbing boundaries?
	);
	// Make a single movement step according to a mechanistic movement model
	// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
	int moveStep(
		Landscape*,		// pointer to Landscape
		Species*,			// pointer to Species
		const short,	// landscape change index
		const bool    // absorbing boundaries?
	);
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
		const float,	// base for power-law (directional persistence or goal bias value)
		const float	// direction in which lowest (unit) weighting is to be applied
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
	void outGenetics( // Write records to genetics file
		const int,		 	// replicate
		const int,		 	// year
		const int,		 	// species number
		const int,		 	// landscape number
		const bool	 		// output as cross table?
	);
	void outMovePath( // Write records to movement paths file
		const int		 	// year
	);

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
	short fallow; // reproductive seasons since last reproduction
	bool isDeveloping;
	Cell *pPrevCell;						// pointer to previous Cell
	Cell *pCurrCell;						// pointer to current Cell
	Patch *pNatalPatch;					// pointer to natal Patch
	emigTraits *emigtraits;			// pointer to emigration traits
	trfrKernTraits *kerntraits;	// pointers to transfer by kernel traits
	pathData *path; 						// pointer to path data for movement model
	crwParams *crw;     				// pointer to CRW traits and data
	smsdata *smsData;						// pointer to variables required for SMS
	settleTraits *setttraits;		// pointer to settlement traits
	std::queue <locn> memory;		// memory of last N squares visited for SMS

	Genome *pGenome;

};


//---------------------------------------------------------------------------

double cauchy(double location, double scale) ;
double wrpcauchy (double location, double rho = exp(double(-1)));

extern RSrandom *pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

extern ofstream outMovePaths;

//---------------------------------------------------------------------------
#endif
