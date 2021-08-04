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

#include "Version.h"
#include "Parameters.h"
#include "Landscape.h"
#include "Community.h"
#include "SubCommunity.h"
#include "Species.h"

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

Rcpp::List RunModel(
	Landscape*,	// pointer to Landscape
	int					// sequential simulation number (always 0 for VCL version)
);
bool CheckDirectory(void);
void PreReproductionOutput(
	Landscape*,	// pointer to Landscape
	Community*, // pointer to Community
	int,				// replicate
	int,				// year
	int					// generation
);
void RangePopOutput(
	Community*, // pointer to Community
	int,				// replicate
	int,				// year
	int					// generation
);
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
void OutParameters(
	Landscape*	// pointer to Landscape
);

extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern Species *pSpecies;
extern paramSim *paramsSim;
extern paramInit *paramsInit;
extern Community *pComm;

const bool batchMode = true;
extern string landFile;
extern vector <string> hfnames;
extern string habmapname;		// see FormLand.cpp (VCL) OR Main.cpp (batch)
extern string patchmapname;	// see FormLand.cpp (VCL) OR Main.cpp (batch)
extern string distnmapname;	// see FormLand.cpp (VCL) OR Main.cpp (batch)
extern string costmapname;	// see FormMove.cpp (VCL) OR Main.cpp (batch)
extern string genfilename;	// see FormGenetics.cpp (VCL) OR Main.cpp (batch)
extern RSrandom *pRandom;

// these functions to have different version for GUI and batch applications ...
extern void MemoLine(string);
void GUIsetLandScale(
	int,	// landscape image height (pixels)
	int		// landscape image width  (pixels)
);

extern std::uint32_t RS_random_seed;
extern string name_landscape, name_patch, name_costfile, name_sp_dist;
//---------------------------------------------------------------------------
#endif
