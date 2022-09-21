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

RangeShifter v2.0 Model

Implements three functions which run the model and produce output common to both
GUI and batch version.

RunModel() handles looping through replicates, years and generations

Further functions are declared here, but defined differently in main function of
GUI and batch versions.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 28 July 2021 by Greta Bocedi
------------------------------------------------------------------------------*/

#ifndef ModelH
#define ModelH

#include <sys/types.h>
#include <sys/stat.h>
#if RS_RCPP
#include <Rcpp.h>
#endif // RS_RCPP 

#include "Version.h"
#include "Parameters.h"
#include "Landscape.h"
#include "Community.h"
#include "SubCommunity.h"
#include "Species.h"

#if !RS_EMBARCADERO && !LINUX_CLUSTER && !RS_RCPP
#include <filesystem>
using namespace std::filesystem;
#endif

#if VIRTUALECOLOGIST
#include "VirtualEcologist.h"
#endif
#if RS_ABC
#include "ABC.h"
#endif
#if PEDIGREE
#include "Pedigree.h"
#endif
#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN
#if RS_THREADSAFE
#include "Rinterface.h"
#endif // RS_THREADSAFE 

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

#if RS_ABC
int RunModel(
	Landscape*,	// pointer to Landscape
	int,				// sequential simulation number (always 0 for VCL version)
	ABCmaster*	// pointer to ABC master object
);
#else
#if RS_RCPP && !R_CMD
#if RS_THREADSAFE
Rcpp::List RunModel(
	Landscape*,	// pointer to Landscape
	int,		// sequential simulation number (always 0 for VCL version)
	Rcpp::S4	// parameter master to read initial conditions in each replicate simulation
);
#else
Rcpp::List RunModel(
	Landscape*,	// pointer to Landscape
	int					// sequential simulation number (always 0 for VCL version)
);
#endif // RS_THREADSAFE
#else
int RunModel(
	Landscape*,	// pointer to Landscape
	int					// sequential simulation number (always 0 for VCL version)
);
#endif // RS_RCPP && !R_CMD
#endif // RS_ABC
bool CheckDirectory(void);
void PreReproductionOutput(
	Landscape*,	// pointer to Landscape
	Community*, // pointer to Community
	int,				// replicate
	int,				// year
	int					// generation
);
#if RS_ABC
void RangePopOutput(
	Community*, // pointer to Community
	int,				// replicate
	int,				// year
	int,				// generation
	ABCmaster*,	// pointer to ABC master object
	bool				// TRUE if ABC observations in current year
);
#else
void RangePopOutput(
	Community*, // pointer to Community
	int,				// replicate
	int,				// year
	int					// generation
);
#endif
void Outputs_Visuals_B(
	int,	// replicate
	int,	// year
	int,	// generation
	int		// Landscape number
);
void RefreshVisualCost(void);
traitCanvas SetupTraitCanvas(void);
void SetupVisualOutput(void);
void ResetVisualOutput(void);
void DrawPopnGraph(
	Community*,	// pointer to Community
	int					// year
);
#if RS_CONTAIN
void ManagementCull(Landscape*,int,int);
#endif // RS_CONTAIN 
void OutParameters(
	Landscape*	// pointer to Landscape
);

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern Species *pSpecies;
extern paramSim *paramsSim;
extern paramInit *paramsInit;
extern Community *pComm;
#if VIRTUALECOLOGIST
extern VirtualEcologist *pVirt;
#endif
#if RS_CONTAIN
extern Cull *pCull;
extern DamageParams *pDamageParams;	
#endif // RS_CONTAIN 

#if VCL
extern bool batchMode;
#else
const bool batchMode = true;
#endif
extern string landFile;
extern vector <string> hfnames;
extern string habmapname;		// see FormLand.cpp (VCL) OR Main.cpp (batch)
extern string patchmapname;	// see FormLand.cpp (VCL) OR Main.cpp (batch)
extern string distnmapname;	// see FormLand.cpp (VCL) OR Main.cpp (batch)
extern string costmapname;	// see FormMove.cpp (VCL) OR Main.cpp (batch)
extern string genfilename;	// see FormGenetics.cpp (VCL) OR Main.cpp (batch)
#if RS_CONTAIN
extern string dmgmapname;		// see FormLand.cpp (VCL) OR Main.cpp (batch)
#endif // RS_CONTAIN 
#if SPATIALMORT
extern string mortmapname[2];	// see FormLand.cpp (VCL) OR Main.cpp (batch)
#endif // SPATIALMORT 
#if TEMPMORT
extern string mortfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // TEMPMORT  
#if VIRTUALECOLOGIST
extern string locfilename;		// see FormVirtEcol.cpp (VCL) OR Main.cpp (batch)
extern string patchfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // VIRTUALECOLOGIST 
#if BUTTERFLYDISP
extern string envstochfilename;
#endif // BUTTERFLYDISP 
extern RSrandom *pRandom;

// these functions to have different version for GUI and batch applications ...
#if BATCH
extern void MemoLine(string);
#endif
#if VCL
extern void MemoLine(UnicodeString);
#endif
void GUIsetLandScale(
	int,	// landscape image height (pixels)
	int		// landscape image width  (pixels)
);

#if RS_RCPP
extern std::uint32_t RS_random_seed;
extern string name_landscape, name_patch, name_costfile, name_sp_dist;
#endif
//---------------------------------------------------------------------------
#endif
