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

RangeShifter v2.0 Parameters

Implements the following classes:

paramGrad  - Environmental gradient parameters
paramInit  - Initialisation (seeding) parameters
paramSim   - Simulation parameters
paramStoch - Environmental stochasticity parameters

Also declares some structures and functions used throughout the program.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 26 November 2020 by Greta Bocedi

------------------------------------------------------------------------------*/

#ifndef ParametersH
#define ParametersH

//#if LINUX_CLUSTER
//#include <string.h>
//#else
#include <string>
//#endif
#include <fstream>
#include <iostream>
//#include <io.h>
#include <iomanip>
#include <stdlib.h>
#include <vector>
using namespace std;

#include "Version.h"
#include "RSrandom.h"

#define NSTAGES 10		// maximum number of stages permitted
#define NSEXES 2			// maximum number of sexes permitted
#define PARAMDEBUG 0
#define NTRAITS 18		// maximum number of variable traits which can be displayed
											// in GUI (VCL version)
#define NSD 3.0				// no. of s.d. to use to control range for displaying traits

typedef intptr_t intptr;

    #ifndef R_EXT_CONSTANTS_H_  // the R headers define PI as a macro, so that the 'else' line results in an error
        #define M_2PI 6.283185307179586
        const double PI = 3.141592653589793238462643383279502884197169399375;
    #endif

const double SQRT2 = std::sqrt(double(2.0)); // more efficient than calculating every time

//---------------------------------------------------------------------------

// Common declarations

struct locn { int x; int y; };
struct rgb { // colour scheme for drawing maps
	int r,g,b;
};

const string Int2Str(const int);
const string Int2Str(const int, unsigned int);
const string Float2Str(const float);
const string Double2Str(const double);
const rgb draw_wheel(int);

//---------------------------------------------------------------------------

// Environmental gradient parameters

// SHOULD THIS BE PART OF LANDSCAPE OBJECT OR A SEPARATE OBJECT?????????????

struct envGradParams {
	bool gradient; bool shifting;
	int gradType; float grad_inc; float opt_y; float factor; float extProbOpt;
	float shift_rate; int shift_begin; int shift_stop;
};

class paramGrad {

public:
	paramGrad(void);
	~paramGrad(void);
	void setGradient(
		int,		// gradient type
		float,	// gradient steepness
		float,	// optimum row (Y dimension)
		float,	// local scaling factor
		float		// local extinction probability at optimum
	);
	void setShifting(
		float,	// shifting rate (rows/year)
		int,		// first year of shifting
		int			// last year of shifting
	);
	void noGradient(void);
	void noShifting(void);
	envGradParams getGradient(void);
	void incrOptY(void);
	void resetOptY(void);

private:
	bool gradient;		// there a gradient
	bool shifting;		// the gradient is shifting
	int gradType;			// 0 = none, 1  = carrying capacity,
										// 2 = growth rate, 3 = local extinction probability
	float grad_inc;		// gradient steepness
	float opt_y;			// optimum row (Y dimension)
	float opt_y0;			// optimum row at year 0 (internal use only)
	float factor;			// local scaling factor
	float extProbOpt;	// local extinction probability at optimum (if gradient = 4, otherwise 0)
	float shift_rate;	// rows/year
	int shift_begin;	// first year of shifting
	int shift_stop;		// last year of shifting
};

//---------------------------------------------------------------------------

// Environmental stochasticity parameters

// SHOULD THIS BE PART OF LANDSCAPE OBJECT OR A SEPARATE OBJECT?????????????

struct envStochParams {
	bool stoch; bool local; bool inK; bool localExt;
	double ac; double std; 
	double locExtProb;
};

class paramStoch {

public:
	paramStoch(void);
	~paramStoch(void);
	void setStoch(envStochParams);
	bool envStoch(void);
	envStochParams getStoch(void);

private:
	bool stoch;				// stochasticity applied
	bool local;				// applied locally (if not, application is global)
	bool inK;					// in carrying capacity (if not, in growth rate)
	bool localExt;		// local extinction applied
	double ac;					// temporal autocorrelation coefficient		
	double std;				// amplitude of fluctuations: sampled from N(0,std)
	double locExtProb;	// local extinction probability
};

//---------------------------------------------------------------------------

// Initialisation (seeding) parameters

struct initParams {
	short seedType; short freeType; short spDistType; short initDens;
	short initAge; int initFrzYr; bool restrictRange;
	int restrictRows; int restrictFreq; int finalFrzYr;
	int indsCell; float indsHa;
	int minSeedX; int maxSeedX; int minSeedY; int maxSeedY;
	int nSeedPatches; int nSpDistPatches;
	string indsFile;
};

struct initInd {
	int year,patchID,x,y; short species,sex,age,stage;
};

class paramInit {

public:
	paramInit(void);
	~paramInit(void);
	void setInit(initParams);
	initParams getInit(void);
	void setProp(
		short,	// stage
		float		// initial proportion
	);
	float getProp(
		short		// stage
	);
	void addInitInd(initInd);
	initInd getInitInd(int);
	void resetInitInds(void);
	int numInitInds(void);

private:
	short seedType;		 	// initialisation type: 0 = free, 1 = from species distn,
											// 2 = initial individuals, 3 = from file
	short freeType;		 	// free initialisation type:
											// 0 = random (given no.)
											// 1 = all suitable cells/patches
											// 2 = manually selected cells/patches
	short spDistType;	 	// species distribution initialisation type:
											// 0 = all suitable cells/patches,
											// 1 = some randomly chosen suitable cells/patches,
											// 2 = all cells/patches within selected sp. dist. cells
	short initDens;		 	// initialisation density:
											// 0 = at carrying capacity
											// 1 = at half carrying capacity
											// 2 = specified no. per cell or density
	short initAge;		 	// initial age distribution within each stage:
											// 0 = lowest possible age
											// 1 = randomised
											// 2 = quasi-equilibrium
	int initFrzYr;		 	// year until which initial range is frozen
	bool restrictRange;	// restrict range to northern front
	int restrictRows;		// no. of rows to retain behind front
	int restrictFreq;		// frequency of range restriction
	int finalFrzYr;		 	// year after which range is frozen
	int indsCell;			 	// initial individuals / cell (cell-based model)
	float indsHa;			 	// initial density (patch-based model)
	int minSeedX;			 	// )
	int maxSeedX;			 	// ) min. and max. of area to initialise (cell numbers)
	int minSeedY;			 	// ) only applied if seedType is 0
	int maxSeedY;			 	// )
	int nSeedPatches;	 	// no. of cells/patches to initialise
	int nSpDistPatches;	// no. of species distribution cells to initialise
	string indsFile;		// no. of species distribution cells to initialise
	float initProp[NSTAGES];	// initial stage proportions (structured population only)

	vector <initInd> initinds;	// individuals to be initialised

};

//---------------------------------------------------------------------------

// Simulation parameters

struct simParams {
	int batchNum;
	int simulation; int reps; int years;
//	int outStartRange;
//	int outStartOcc;
	int outStartPop; int outStartInd; int outStartGenetic;
	int outStartTraitCell; int outStartTraitRow; int outStartConn;
	int outIntRange; int outIntOcc; int outIntPop; int outIntInd; int outIntGenetic;
	int outIntTraitCell; int outIntTraitRow; int outIntConn;
	int mapInt; int traitInt;
	bool batchMode; bool absorbing;
	bool outRange; bool outOccup; bool outPop; bool outInds;
	bool outGenetics; short outGenType; bool outGenXtab;
	bool outTraitsCells; bool outTraitsRows; bool outConnect;
	bool saveMaps;
	bool drawLoaded; bool saveTraitMaps;
	bool saveVisits;
	int outStartPaths; int outIntPaths;
	bool outPaths;	bool ReturnPopRaster; bool CreatePopFile;
};

struct simView {
	bool viewLand; bool viewPatch; bool viewGrad; bool viewCosts;
	bool viewPop; bool viewTraits; bool viewPaths; bool viewGraph;
	int slowFactor;
};

class paramSim {

public:
	paramSim(void);
	~paramSim(void);
	void setSim(simParams);
	simParams getSim(void);
	int getSimNum(void);
	void setViews(simView);
	simView getViews(void);
	void setDir(string);
	string getDir(int);
	bool getReturnPopRaster(void);
	bool getCreatePopFile(void);

private:
	int batchNum;						// batch number
	int simulation;					// simulation no.
	int reps;								// no. of replicates
	int years;							// no. of years
//	int outStartRange;			// output start year for range file
//	int outStartOcc;				// output start year for occupancy file
	int outStartPop;				// output start year for population file
	int outStartInd;				// output start year for individuals file
	int outStartGenetic; 		// output start year for genetics file
	int outStartTraitCell;	// output start year for traits by cell file
	int outStartTraitRow;		// output start year for traits by row file
	int outStartConn;				// output start year for connectivity matrix
	int outIntRange;				// output interval for range file
	int outIntOcc;					// output interval for occupancy file
	int outIntPop;					// output interval for population file
	int outIntInd;					// output interval for individuals file
	int outIntGenetic;			// output interval for genetics file
	int outIntTraitCell;		// output interval for traits by cell file
	int outIntTraitRow;			// output interval for traits by row file
	int outIntConn;					// output interval for connectivity matrix
	int mapInt;							// output interval for maps
	int traitInt;						// output interval for evolving traits maps
	int slowFactor;					// to reduce speed of movement paths on screen
	bool batchMode;					//
	bool absorbing; 				// landscape boundary and no-data regions are
													// absorbing boundaries
	bool outRange;					// produce output range file?
	bool outOccup;					// produce output occupancy file?
	bool outPop;						// produce output population file?
	bool outInds;						// produce output individuals file?
	bool outGenetics;				// produce output genetics file?
	short outGenType;				// produce output genetics for: 0 = juveniles only
													// 1 = all individuals, 2 = adults (i.e. final stage) only
	bool outGenXtab;				// produce output genetics as a cross table?
	bool outTraitsCells;		// produce output summary traits by cell file?
	bool outTraitsRows;			// produce output summary traits by row (y) file?
	bool outConnect;				// produce output connectivity file?
	bool saveMaps;					// save landscape/population maps?
	bool saveVisits;        // save dispersal visits heat maps?
	int outStartPaths;
	int outIntPaths;
	bool outPaths;
	bool ReturnPopRaster;
	bool CreatePopFile;
	bool drawLoaded;				// draw initial distribution on landscape/population maps?
	bool saveTraitMaps;			// save summary traits maps?
	bool viewLand;					// view landscape map on screen?
	bool viewPatch;					// view map of landscape patches on screen?
	bool viewGrad;					// view gradient map on screen?
	bool viewCosts;					// view costs map on screen?
	bool viewPop;						// view population density on landscape map on screen?
	bool viewTraits;				// view summary traits map(s) on screen?
	bool viewPaths;					// view individual movement paths on screen?
	bool viewGraph;					// view population/occupancy graph on screen?
	string dir;							// full name of working directory

};

#if RSDEBUG
extern ofstream DEBUGLOG;
void DebugGUI(string);
#endif

//---------------------------------------------------------------------------
#endif
