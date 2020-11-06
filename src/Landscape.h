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

Last updated: 24 July 2020 by Anne-Kathleen Malchow, Potsdam University

------------------------------------------------------------------------------*/

#ifndef LandscapeH
#define LandscapeH

//#include <stdlib.h>
//#include <math.h>
#include <algorithm>
#include <fstream>
//#include <iostream.h>
//#include <stdio.h>
#include <vector>

using namespace std;

#include "Parameters.h"
#include "Patch.h"
#include "Cell.h"
#include "Species.h"
#include "FractalGenerator.h"
#include <locale>
#if !RSWIN64
#include <codecvt>
#endif
#include <Rcpp.h>

//---------------------------------------------------------------------------

// Initial species distribution

class InitDist{
public:
	InitDist(Species*);
	~InitDist();
	int readDistribution(
		string // name of species distribution file
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
struct landOrigin {
	double minEast; double minNorth;
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
	bool utf;
};
struct patchData {
	Patch *pPatch; int patchNum,nCells; int x,y;
};
struct landChange {
	int chgnum,chgyear; string habfile,pchfile,costfile;
};
struct patchChange {
	int chgnum,x,y,oldpatch,newpatch;
};
struct costChange {
	int chgnum,x,y,oldcost,newcost;
};

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
	void allocatePatches(Species*);	// create patches for a cell-based landscape
	Patch* newPatch(
		int		// patch sequential no. (id no. is set to equal sequential no.)
	);
	Patch* newPatch(
		int,  // patch sequential no.
		int		// patch id no.
	);
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
	void setEnvGradient(
		Species*, // pointer to Species
		bool      // TRUE for initial instance that gradient is set
	);
	void setGlobalStoch(
		int		// no. of years
	);
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

	// functions to handle species distributions

	int newDistribution(
		Species*,	// pointer to Species
		string		// name of initial distribution file
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
	void outPathsHeaders(int, int);
	void outConnect(
		int,	// replicate no.
		int   // year
	);

	// functions to handle input and output

	int readLandscape(
		int,		// fileNum == 0 for (first) habitat file and optional patch file
						// fileNum > 0  for subsequent habitat files under the %cover option
		string,	// habitat file name
		string,	// patch file name
		string	// cost file name (may be NULL)
	);
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
	void outVisits(int,int);	// save SMS path visits map to raster text file

private:
	bool generated;				// artificially generated?
	bool patchModel;			//
	bool spDist;					// initial species distribution loaded
	bool fractal;					//
	bool continuous;			//
	bool dynamic;					// landscape changes during simulation
	bool habIndexed;			// habitat codes have been converted to index numbers
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

extern void MemoLine(string);

extern rasterdata landraster,patchraster,spdistraster,costsraster;
extern void EOFerrorR(string);
extern void StreamErrorR(string);

//---------------------------------------------------------------------------
#endif
