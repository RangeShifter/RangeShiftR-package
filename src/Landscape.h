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

RangeShifter v2.0 Landscape

Implements the following classes:

InitDist  - Initial species distribution

Landscape - Landscape grid

The Landscape is a rectangular array of Cells which are grouped together in
Patches. As far as the user is aware, the Landscape is either patch-based or
cell-based (having no Patches), but internally Patches are always present (they
each comprise only one Cell for a cell-based model). The grain of the Landscape
may be any positive integer, and is nominally in metres.

The Landscape is either input from one or more text files in ArcGIS raster export
format, or it is generated artificially as a random or fractal binary array (in
which case, it must be cell-based). An input 'real' Landscape may hold within each
Cell either discrete habitat classes, or percent cover of habitat classes, or a
continuous quality index (1 to 100%).

The Landscape may be dynamic, in which case the user specifies a set of changes
to be applied at certain years during each simulation. The changes may be to
habitat only, patches only (if a patch-based model) or habitats and patches.
Although the changes must be supplied as entire habitat / patch files (which
must match the original Landscape in terms of cell size and extent), internally
they are recorded as lists of dynamic habitat and patch changes.

The initial species distribution is a rectangular array if distribution cells
(DistCell) covering the same spatial extent at the Landscape. Its grain may be
either that of the Landscape or an integer multiple of it.

The Landscape can record a list (in the vector initcells) of Cells or Patches
to be intialised, which are specified by the user in FormSeeding. This option is
available in the GUI version only.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 28 July 2021 by Greta Bocedi
------------------------------------------------------------------------------*/

#ifndef LandscapeH
#define LandscapeH

#include "Version.h"

#include <algorithm>
#include <fstream>
#include <vector>

/*
#include <locale>
#include <codecvt>
#include <Rcpp.h>
#include <armadillo>
*/
#if RS_RCPP
#include <locale>
#if !RSWIN64
#include <codecvt>
#endif //!RSWIN64
#if !R_CMD
#include <Rcpp.h>
#if SPATIALDEMOG
#include <armadillo>
#endif //SPATIALDEMOG
#endif //!R_CMD
#endif //RS_RCPP

using namespace std;

#include "Parameters.h"
#include "Patch.h"
#include "Cell.h"
#include "Species.h"
#include "FractalGenerator.h"

#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN 

//---------------------------------------------------------------------------


struct landOrigin {
	double minEast; double minNorth;
};


// Initial species distribution

class InitDist{
public:
	InitDist(Species*);
	~InitDist();
	int readDistribution(
#if RS_THREADSAFE
		Rcpp::NumericMatrix,
		landOrigin,
		int
#else
		string // name of species distribution file
#endif
	);
	void setDistribution(
		int	// no. of distribution cells to be initialised (0 for all cells)
	);
	void setDistCell( // Set a specified cell (by position in cells vector)
		int,	// index no. of DistCell in cells vector
		bool  // value to be set
	);
	void setDistCell( // Set a specified cell (by co-ordinates)
		locn,	// structure holding x (column) and y (row) co-ordinates
		bool
	);
	bool inInitialDist( // Specified location is within the initial distribution?
		locn  // structure holding x (column) and y (row) co-ordinates
	);
	int cellCount(void);
	locn getCell( // Return the co-ordinates of a specified initial distribution cell
		int  // index no. of DistCell in cells vector
	);
	locn getSelectedCell( // Return the co-ordinates of a specified initial distribution
												// cell if it has been selected
												// otherwise return negative co-ordinates
		int  // index no. of DistCell in cells vector
	);
	locn getDimensions(void);
	void resetDistribution(void);

private:
	Species *pSpecies;		// pointer to species
	int resol;						// species distribution cell size (m)
	int maxX, maxY;				// dimensions
	double minEast;				// ) real world min co-ordinates
	double minNorth;			// ) read from raster file

	// list of cells in the initial distribution
	// cells MUST be loaded in the sequence ascending x within descending y
	std::vector <DistCell*> cells;

};


//---------------------------------------------------------------------------

struct landParams {       
	bool patchModel; bool spDist; bool generated;
	bool dynamic;
#if RS_CONTAIN
	bool dmgLoaded;
#endif // RS_CONTAIN 
#if SPATIALDEMOG
	bool spatialdemog;
#endif // SPATIALDEMOG 
	int landNum; int resol; int spResol; int nHab; int nHabMax;
	int dimX,dimY,minX,minY,maxX,maxY;
	short rasterType;
};
struct landData {
	int resol; int dimX,dimY,minX,minY,maxX,maxY;
};
struct genLandParams {
	bool fractal; bool continuous;
	float minPct,maxPct; float propSuit; float hurst; int maxCells;
};
struct landPix {
	int pix; float gpix;
};
struct rasterHdr {
	bool ok;
	int errors,ncols,nrows,cellsize;
	double xllcorner,yllcorner;
};
struct rasterdata {
	bool ok;
	int errors,ncols,nrows,cellsize;
	double xllcorner,yllcorner;
#if RS_RCPP
	bool utf;
#endif
};
struct patchData {
	Patch *pPatch; int patchNum,nCells; int x,y;
};
struct landChange {
	int chgnum,chgyear; string habfile,pchfile,costfile;
};
struct patchChange {
	int chgnum, x, y, oldpatch, newpatch;
};
struct costChange {
	int chgnum,x,y,oldcost,newcost;
};

#if SEASONAL
//#if PARTMIGRN
struct extEvent { // extreme event
	int year,season,patchID,x,y; float probMort;
};
//#endif // PARTMIGRN 
#endif // SEASONAL

class Landscape{
public:
	Landscape();
	~Landscape();
	void resetLand(void);

	// functions to set and return parameter values

	void setLandParams(
		landParams,	// structure holding landscape parameters
		bool				// batch mode
	);
	landParams getLandParams(void);
	landData getLandData(void);
	void setGenLandParams(genLandParams);
	genLandParams getGenLandParams(void);
	void setLandLimits(
		int,	// minimum available X
		int,	// minimum available Y
		int,	// maximum available X
		int		// maximum available Y
	);
	void resetLandLimits(void);
	void setLandPix(landPix);

	landPix getLandPix(void);
	void setOrigin(landOrigin);
	landOrigin getOrigin(void);

	// functions to handle habitat codes

	bool habitatsIndexed(void);
	void listHabCodes(void);
	void addHabCode(int);
	int findHabCode(int);
	int getHabCode(int);
	void clearHabitats(void);
	void addColour(rgb);
	void changeColour(int,rgb);
	rgb getColour(int);
	int colourCount(void);

	// functions to handle patches and cells

	void setCellArray(void);
	void addPatchNum(int);
	void generatePatches(void); 		// create an artificial landscape
#if SEASONAL
	void allocatePatches( // create patches for a cell-based landscape
		Species*,		// pointer to Species
		short				// no. of seasons
	);	
	Patch* newPatch(
		int,		// patch sequential no. (id no. is set to equal sequential no.)
		short		// no. of seasons
	);
	Patch* newPatch(
		int,	  // patch sequential no.
		int,		// patch id no.
		short		// no. of seasons
	);
#else
	void allocatePatches(Species*);	// create patches for a cell-based landscape
	Patch* newPatch(
		int		// patch sequential no. (id no. is set to equal sequential no.)
	);
	Patch* newPatch(
		int,  // patch sequential no.
		int		// patch id no.
	);
#endif // SEASONAL
	void resetPatches(void);
	void addNewCellToLand(
		int,    // x co-ordinate
		int,    // y co-ordinate
		float   // habitat quality value
	);
	void addNewCellToLand(
		int,    // x co-ordinate
		int,    // y co-ordinate
		int     // habitat class no.
	);
	void addCellToPatch(
		Cell*,	// pointer to Cell
		Patch*	// pointer to Patch
	);
	void addCellToPatch(
		Cell*,  // pointer to Cell
		Patch*, // pointer to Patch
		float		// habitat quality value
	);
	void addCellToPatch(
		Cell*,	// pointer to Cell
		Patch*, // pointer to Patch
		int     // habitat class no.
	);
	void addNewCellToPatch(
		Patch*, // pointer to Patch
		int,    // x co-ordinate
		int,    // y co-ordinate
		int     // habitat class no.
	);
	void addNewCellToPatch(
		Patch*, // pointer to Patch
		int,    // x co-ordinate
		int,    // y co-ordinate
		float   // habitat quality value
	);
#if RS_CONTAIN
	void setDamage(
//		Patch*, // pointer to Patch
		int,    // x co-ordinate
		int,    // y co-ordinate
		intptr,	// pointer (cast as integer) to the Patch in which cell lies (if any)
		int			// damage index 
	);
	void updateDamage(
		int,    // x co-ordinate
		int,    // y co-ordinate
		intptr	// pointer (cast as integer) to the Patch in which cell lies (if any)
	);
	void updateDamageIndices(void);
	void setAlpha(double);
	double getAlpha(void);
	void resetDamageLocns(void);
//	void updateDamageLocns(Species*);
	double totalDamage(
		bool		// transfer by SMS?
	);
	void viewDamage(	// Update the damage graph on the screen
										// NULL for the batch version
		int,		// year
		double,	// mean damage
		double,	// standard error of damage
		bool		// show standard error?
	);
	void createTotDamage(int nrows,int reps);
	void updateTotDamage(unsigned short row,unsigned short rep,float damage); 
	void deleteTotDamage(int nrows);
	void outTotDamage(bool view);
//	void resetPrevDamage(void);
#endif // RS_CONTAIN 
	patchData getPatchData(
		int		// index no. of Patch in patches vector
	);
	bool existsPatch(
		int		// Patch id no.
	);
	Patch* findPatch(
		int   // Patch id no.
	);
	int checkTotalCover(void);
	void resetPatchPopns(void);
	void updateCarryingCapacity(
		Species*,	// pointer to Species
		int,			// year
		short			// landscape change index (always zero if not dynamic)
	);
	Cell* findCell(
		int,		// x co-ordinate
		int			// y co-ordinate
	);
	int patchCount(void);
	void updateHabitatIndices(void);
#if SEASONAL
	void setEnvGradient(
		Species*,	// pointer to Species
		short,		// no. of seasons
		bool      // TRUE for initial instance that gradient is set
	);
#else
	void setEnvGradient(
		Species*, // pointer to Species
		bool      // TRUE for initial instance that gradient is set
	);
#endif // SEASONAL 
	void setGlobalStoch(
		int		// no. of years
	);
#if BUTTERFLYDISP
	void readGlobalStoch(
		int,		// no. of years
		string	// filename
	);
#endif // BUTTERFLYDISP 
	float getGlobalStoch(
		int		// year
	);
	void updateLocalStoch(void);
	void resetCosts(void);
	void resetEffCosts(void);

	// functions to handle dynamic changes

	void setDynamicLand(bool);
	void addLandChange(
		landChange	// structure holding landscape change data
	);
	int numLandChanges(void);
	landChange getLandChange(
		short	// change number
	);
	void deleteLandChanges(void);
#if RS_THREADSAFE
	int readLandChange(
		int,		// change number
		Rcpp::NumericMatrix,// habitat raster
		Rcpp::NumericMatrix,// patch raster
		Rcpp::NumericMatrix	// cost raster
#if SPATIALDEMOG
		,Rcpp::NumericVector// array of demographic scaling layers
#endif
	);
#else
#if RS_RCPP && !R_CMD
	int readLandChange(
	    int,		// change file number
		bool,		// change SMS costs?
		wifstream&, // habitat file stream
		wifstream&, // patch file stream
		wifstream&, // cost file stream
		int,		// habnodata
		int,		// pchnodata
		int			// costnodata
	);
#else
	int readLandChange(
		int,	// change file number
		bool	// change SMS costs?
	);
#endif
#endif // RS_THREADSAFE
	void createPatchChgMatrix(void);
	void recordPatchChanges(int);
	void deletePatchChgMatrix(void);
	int numPatchChanges(void);
	patchChange getPatchChange(
		int	// patch change number
	);
	void createCostsChgMatrix(void);
	void recordCostChanges(int);
	void deleteCostsChgMatrix(void);
	int numCostChanges(void);
	costChange getCostChange(
		int	// cost change number
	);
#if SPATIALDEMOG
	void updateDemoScalings(short);
#endif // SPATIALDEMOG 

#if SPATIALMORT
	// functions to handle spatial mortality

	int readMortalityFiles(
		string,	// mortality file name period 0
		string	// mortality file name period 1
	);
#endif // SPATIALMORT 

	// functions to handle species distributions

	int newDistribution(
		Species*,	// pointer to Species
#if RS_THREADSAFE
		Rcpp::NumericMatrix,
		int
#else
		string		// name of initial distribution file
#endif
	);
	void setDistribution(
		Species*, // pointer to Species
		int				// no. of distribution squares to initialise
	);
	bool inInitialDist( // Specified cell matches one of the distn cells to be initialised?
		Species*, // pointer to Species
		locn			// structure holding co-ordinates of Cell
	);
	void deleteDistribution(
		Species*  // pointer to Species
	);
	int distnCount(void);	// Return no. of initial distributions in the Landscape
	int distCellCount(    // Return no. of distribution cells in an initial distribution
		int // index no. of InitDist in distns vector
	);
	locn getDistnCell( // Get co-ordinates of a specified cell in a specified initial distn
		int,	// index no. of InitDist in distns vector
		int		// index no. of DistCell in cells vector
	);
	locn getSelectedDistnCell(	// Get co-ordinates of a specified cell in a specified initial distn
															// Returns negative co-ordinates if the cell is not selected
		int,  // index no. of InitDist in distns vector
		int   // index no. of DistCell in cells vector
	);
	locn getDistnDimensions(	// Get the dimensions of a specified initial distribution
		int   // index no. of InitDist in distns vector
	);
	void setDistnCell( // Set a cell in a specified init distn (by posn in cells vector)
		int,  // index no. of InitDist in distns vector
		int,  // index no. of DistCell in cells vector
		bool	// value to be set
	);
	void setDistnCell( // Set a cell in a specified init distn (by given co-ordinates)
		int,  // index no. of InitDist in distns vector
		locn, // structure holding co-ordinates of DistCell
		bool  // value to be set
	);
	void resetDistribution(
		Species*	// pointer to Species
	);

	// functions to handle initialisation cells

	int initCellCount(void);
	void addInitCell( // Create a new DistCell and add to the initcells vector
		int,	// x co-ordinate
		int   // y co-ordinate
	);
	locn getInitCell(
		int   // index no. of DistCell in initcells vector
	);
	void clearInitCells(void);

	// functions to handle connectivity matrix

	void createConnectMatrix(void);
	void resetConnectMatrix(void);
	void incrConnectMatrix(
		int,	// sequential no. of origin Patch
		int   // sequential no. of settlement Patch
	);
	void deleteConnectMatrix(void);
	bool outConnectHeaders( // Write connectivity file headers
		int		// option - set to -999 to close the connectivity file
	);
#if RS_RCPP
	void outPathsHeaders(int, int);
#endif
#if SEASONAL
	void outConnect(
		int,	// replicate no.
		int,  // year
		short // season
	);
#else
	void outConnect(
		int,	// replicate no.
		int   // year
	);
#endif // SEASONAL 

	// functions to handle input and output

#if RS_THREADSAFE
	int readLandscape(
		int, 				// no. of seasonss
		Rcpp::NumericMatrix,// habitat raster
		Rcpp::NumericMatrix,// patch raster
		Rcpp::NumericMatrix	// cost raster
#if SPATIALDEMOG
		,Rcpp::NumericVector // array of demographic scaling layers
#endif
	);
#else //RS_THREADSAFE
#if RS_CONTAIN
#if SEASONAL
	int readLandscape(
		int,		// no. of seasonss
		int,		// fileNum == 0 for (first) habitat file and optional patch file
						// fileNum > 0  for subsequent habitat files under the %cover option
		string,	// habitat file name
		string,	// patch file name
		string,	// cost file name (may be NULL)
		string	// damage file name (may be NULL)
	);
#else
	int readLandscape(
		int,		// fileNum == 0 for (first) habitat file and optional patch file
						// fileNum > 0  for subsequent habitat files under the %cover option
		string,	// habitat file name
		string,	// patch file name
		string,	// cost file name (may be NULL)
		string	// damage file name (may be NULL)
	);
#endif // SEASONAL 
	bool outSummDmgHeaders( // Open summary damage file and write header record
		int				// Landscape number (-999 to close the file)
	);
	void outSummDmg( // Write record to summary damage file
		int,			// replicate
		int,			// year
		bool,			// transfer by SMS?
		bool			// view damage on screen
	);
	bool outDamageHeaders( // Open damage file and write header record
		int				// Landscape number (-999 to close the file)
	);
	void outDamage( // Write record to damage file
		int,			// replicate
		int,			// year
		bool			// transfer by SMS?
	);
#else
#if SEASONAL
	int readLandscape(
		int,		// no. of seasonss
		int,		// fileNum == 0 for (first) habitat file and optional patch file
						// fileNum > 0  for subsequent habitat files under the %cover option
		string,	// habitat file name
		string,	// patch file name
		string	// cost file name (may be NULL)
	);
#else
	int readLandscape(
		int,		// fileNum == 0 for (first) habitat file and optional patch file
						// fileNum > 0  for subsequent habitat files under the %cover option
		string,	// habitat file name
		string,	// patch file name
		string	// cost file name (may be NULL)
	);
#endif // SEASONAL
#endif // RS_CONTAIN
#endif // RS_THREADSAFE
	void listPatches(void);
	int readCosts(
		string	// costs file name
	);
	// the following four functions are implemented for the GUI version only
	// in the batch version, they are defined, but empty
	void setLandMap(void);
	void drawLandscape(
		int,	// replicate no.
		int,	// landscape index number (always 0 if landscape is not dynamic)
		int		// landscape no.
	);
	void drawGradient(void); // Draw environmental gradient map
	void drawGlobalStoch(	// Draw environmental stochasticity time-series
		int		// no. of years
	);

	void resetVisits(void);
#if VCL
	void saveVisits(int,int); // save SMS path visits map to .bmp file
#endif
	void outVisits(int,int);	// save SMS path visits map to raster text file

#if RS_ABC
// Returns connectivity (no. of successful dispersers) for given start and end patches
int outABCconnect(int,int);
#endif

#if SEASONAL
//#if PARTMIGRN
	// extreme events
	void addExtEvent(extEvent);
	extEvent getExtEvent(int);
	void resetExtEvents(void);
	int numExtEvents(void);
//#endif // PARTMIGRN 
#endif // SEASONAL

private:
	bool generated;				// artificially generated?
	bool patchModel;			//
	bool spDist;					// initial species distribution loaded
	bool fractal;					//
	bool continuous;			//
	bool dynamic;					// landscape changes during simulation
	bool habIndexed;			// habitat codes have been converted to index numbers
#if RS_CONTAIN
	bool dmgLoaded;				// economic / environmental damage values have been input
#endif // RS_CONTAIN 
#if SPATIALDEMOG
	bool spatialdemog;			// are there spatially varying demographic rates?
#endif // SPATIALDEMOG 
	short rasterType;			// 0 = habitat codes 1 = % cover 2 = quality 9 = artificial landscape
	int landNum;					// landscape number
	int resol;						// cell size (m)
	int spResol;					// species distribution cell size (m)
	int nHab;							// no. of habitats
	int nHabMax;					// max. no. of habitats (used for batch input only)
	int dimX,dimY;				// dimensions
	int minX,minY;				// minimum available X and Y co-ordinates
	int maxX,maxY;				// maximum available X and Y co-ordinates
	float minPct,maxPct;  // min and max percentage of habitat in a cell
	float propSuit;				// proportion of suitable cells
	float hurst;					// Hurst exponent
	int maxCells;					// max. cells per patch (artificial landscapes)
	int pix;							// image display ratio
	float gpix;						// image display ratio for gradient map
	double minEast;				// ) real world min co-ordinates
	double minNorth;			// ) read from habitat raster
#if RS_CONTAIN
	double alpha;					// economic / environmental damage distance decay coefficient
#endif // RS_CONTAIN 

	// list of cells in the landscape
	// cells MUST be loaded in the sequence ascending x within descending y
	Cell ***cells;

	// list of patches in the landscape - can be in any sequence
	std::vector <Patch*> patches;

	// list of patch numbers in the landscape
	std::vector <int> patchnums;

	// list of habitat codes
	std::vector <int> habCodes;

	// list of colours for habitat codes
	std::vector <rgb> colours;

	// list of dynamic landscape changes
	std::vector <landChange> landchanges;
	std::vector <patchChange> patchchanges;
	std::vector <costChange> costschanges;

	// list of initial individual species distributions
	std::vector <InitDist*> distns;

	// list of cells to be initialised for ALL species
	std::vector <DistCell*> initcells;

	// patch connectivity matrix
	// indexed by [start patch seq num][end patch seq num]
	int **connectMatrix;

	// global environmental stochasticity (epsilon)
	float *epsGlobal;	// pointer to time-series	

	// patch and costs change matrices (temporary - used when reading dynamic landscape)
	// indexed by [descending y][x][period]
	// where there are three periods, 0=original 1=previous 2=current
	int ***patchChgMatrix;
	int ***costsChgMatrix;

#if SEASONAL
//#if PARTMIGRN
	// extreme events
	std::vector <extEvent> extevents;
//#endif // PARTMIGRN 
#endif // SEASONAL 

#if RS_CONTAIN
	std::vector <DamageLocn*> dmglocns;
	float **totDamage;	// total damage record to view on screen
#endif // RS_CONTAIN 
	
};

// NOTE: the following function is not a behaviour of Landscape, as it is run by the
// batch routine to check raster files before any Landscape has been initiated
rasterdata CheckRasterFile(string);

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;
extern paramSim *paramsSim;
extern RSrandom *pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
extern void DebugGUI(string);
#endif

#if VCL
extern void MemoLine(UnicodeString);
#else
extern void MemoLine(string);
#endif

#if RS_RCPP
extern rasterdata landraster,patchraster,spdistraster,costsraster;
extern void EOFerrorR(string);
extern void StreamErrorR(string);
#endif

//---------------------------------------------------------------------------
#endif
