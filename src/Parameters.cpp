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
 
 
//---------------------------------------------------------------------------

#include "Parameters.h"
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

// Environmental gradient parameters

paramGrad::paramGrad(void) {
gradient = false; gradType = 0; grad_inc = 0.05f;
opt_y0 = opt_y = factor = extProbOpt = 0.0;
shifting = false; shift_rate = 0.5; shift_begin = 0; shift_stop = 100;
}

paramGrad::~paramGrad() { }

void paramGrad::setGradient(int gtype, float inc, float y, float f, float p)
{
if (gtype > 0 && gtype < 4)
{ // valid gradient type
	gradient = true; gradType = gtype;
	if (inc >= 0.0 && inc <= 1.0) grad_inc = inc;
	if (y >= 0.0) opt_y0 = opt_y = y;
	if (f >= 0.0) factor = f;
	if (p > 0.0 && p < 1.0) extProbOpt = p;
}
else {
	gradient = false; gradType = 0;
}
}

void paramGrad::setShifting(float r, int begin, int end)
{
shifting = true;
if (r > 0.0) shift_rate = r;
if (begin >= 0) shift_begin = begin;
if (end > 0) shift_stop = end;
}

void paramGrad::noGradient(void) { gradient = false; gradType = 0; }

void paramGrad::noShifting(void) { shifting = false; }

envGradParams paramGrad::getGradient(void) {
envGradParams g;
g.gradient = gradient; g.gradType = gradType; g.grad_inc = grad_inc;
g.opt_y = opt_y; g.factor = factor; g.extProbOpt = extProbOpt;
g.shifting = shifting; g.shift_rate = shift_rate;
g.shift_begin = shift_begin; g.shift_stop = shift_stop;
return g;
}

void paramGrad::incrOptY(void)
{ if (gradient && shifting) opt_y += shift_rate; }

void paramGrad::resetOptY(void) { opt_y = opt_y0; }

//---------------------------------------------------------------------------

// Environmental stochasticity parameters

paramStoch::paramStoch(void) {
stoch = false; local = false; inK = false; localExt = false;
ac = 0.0; std = 0.25;
locExtProb = 0.1;
}

paramStoch::~paramStoch(void) {}


void paramStoch::setStoch(envStochParams e)
{
stoch = e.stoch; local = e.local; inK = e.inK; localExt = e.localExt;
if (e.ac >= 0.0 && e.ac < 1.0) ac = e.ac;
if (e.std > 0.0 && e.std <= 1.0) std = e.std;
locExtProb = e.locExtProb;
}

bool paramStoch::envStoch(void) { return stoch; }

envStochParams paramStoch::getStoch(void)
{
envStochParams e;
e.stoch = stoch; e.local = local; e.inK = inK; e.localExt = localExt;
e.ac = ac; e.std = std;
e.locExtProb = locExtProb;
return e;
}


//---------------------------------------------------------------------------

// Initialisation (seeding) parameters

paramInit::paramInit(void) {
seedType = freeType = spDistType = initDens = 0;
initAge = initFrzYr = 0;
restrictRange = false;
restrictRows = 100;
restrictFreq = 10;
finalFrzYr = 99999999;
indsCell = 1; indsHa = 0.0;
minSeedX = 0; maxSeedX = 99999999; minSeedY = 0; maxSeedY = 99999999;
nSeedPatches = 1; nSpDistPatches = 1;
indsFile = "NULL";
for (int i = 0; i < NSTAGES; i++) {
	initProp[i] = 0.0;
}
}

paramInit::~paramInit(void) {
initinds.clear();
}

void paramInit::setInit(initParams i) {
if (i.seedType >= 0 && i.seedType <= 3) seedType = i.seedType;
if (i.freeType >= 0 && i.freeType <= 2) freeType = i.freeType;
if (i.spDistType >= 0 && i.spDistType <= 2) spDistType = i.spDistType;
initDens = i.initDens;
initAge = i.initAge;
if (i.initFrzYr >= 0) initFrzYr = i.initFrzYr;
restrictRange = i.restrictRange;
if (i.restrictRows > 0) restrictRows = i.restrictRows;
if (i.restrictFreq > 0) restrictFreq = i.restrictFreq;
if (i.finalFrzYr > 0) finalFrzYr = i.finalFrzYr;
if (i.indsCell >= 1) indsCell = i.indsCell;
if (i.indsHa > 0.0) indsHa = i.indsHa;
if (i.minSeedX >= 0) minSeedX = i.minSeedX;
if (i.maxSeedX >= 0) maxSeedX = i.maxSeedX;
if (i.minSeedY >= 0) minSeedY = i.minSeedY;
if (i.maxSeedY >= 0) maxSeedY = i.maxSeedY;
if (i.nSeedPatches >= 1) nSeedPatches = i.nSeedPatches;
if (i.nSpDistPatches >= 1) nSpDistPatches = i.nSpDistPatches;
indsFile = i.indsFile;
}

initParams paramInit::getInit(void) {
initParams i;
i.seedType = seedType; i.freeType = freeType; i.spDistType = spDistType;
i.initDens = initDens; i.initAge = initAge;
i.initFrzYr = initFrzYr;
i.restrictRange = restrictRange;
i.restrictRows = restrictRows; i.restrictFreq = restrictFreq;
i.finalFrzYr = finalFrzYr;
i.indsCell = indsCell; i.indsHa = indsHa;
i.minSeedX = minSeedX; i.minSeedY = minSeedY;
i.maxSeedX = maxSeedX; i.maxSeedY = maxSeedY;
i.nSeedPatches = nSeedPatches; i.nSpDistPatches = nSpDistPatches;
i.indsFile = indsFile;
return i;
}

void paramInit::setProp(short stg,float p) {
if (stg >= 0 && stg < NSTAGES && p >= 0.0 && p <= 1.0) initProp[stg] = p;
}

float paramInit::getProp(short stg) {
float p = 0.0;
if (stg >= 0 && stg < NSTAGES) p = initProp[stg];
return p;
}

void paramInit::addInitInd(initInd iind) {
#if RSDEBUG
//DebugGUI(("paramInit::addInitInd(): iind.patchID=" + Int2Str(iind.patchID)
//	+ " iind.x=" + Int2Str(iind.x)
//	+ " iind.y=" + Int2Str(iind.y)
//	).c_str());
#endif
initinds.push_back(iind);
}

initInd paramInit::getInitInd(int ix) {
initInd iind;
if (ix >= 0 && ix < (int)initinds.size()) {
	iind = initinds[ix];
}
else {
	iind.year = iind.patchID = iind.x = iind.y = iind.sex = iind.age = iind.stage = 0;
	iind.species = -1;
}
#if RSDEBUG
//DEBUGLOG << "paramInit::getInitInd(): ix=" << ix << " size()=" << initinds.size()
//	<< " iind.patchID=" << iind.patchID
//	<< " iind.x=" << iind.x
//	<< " iind.y=" << iind.y
//	<< endl;
#endif
return iind;
}

void paramInit::resetInitInds(void) { initinds.clear(); }

int paramInit::numInitInds(void) { return (int)initinds.size(); }


//---------------------------------------------------------------------------

// Simulation parameters

paramSim::paramSim(void) {
	simulation = 0;
	reps = years = 1;
	outIntRange = 1;
//	outStartRange = outStartOcc = outStartPop = outStartInd = 0;
	outStartPop = outStartInd = outStartGenetic = 0;
	outStartTraitCell = outStartTraitRow = outStartConn = 0;
	outIntOcc = outIntPop = outIntInd = outIntGenetic = 10;
	outIntTraitCell = outIntTraitRow = outIntConn = 10;
	mapInt = traitInt = 10;
	slowFactor = 1;
	batchMode = absorbing = false;
	outRange = outOccup = outPop = outInds = false;
	outGenetics = outGenXtab = false; outGenType = 0;
	outTraitsCells = outTraitsRows = outConnect = false;
	saveMaps = false; saveTraitMaps = false;
	saveVisits = false;
	outStartPaths = 0; outIntPaths = 0;
	outPaths = false; ReturnPopRaster = false; CreatePopFile = true;
	drawLoaded = false;
	viewLand = false; viewPatch = false; viewGrad = false; viewCosts = false;
	viewPop = false; viewTraits = false; viewPaths = false; viewGraph = false;
	dir = ' ';
}

paramSim::~paramSim(void) { }

void paramSim::setSim(simParams s) {
if (s.batchNum >= 0) batchNum = s.batchNum;
if (s.simulation >= 0) simulation = s.simulation;
if (s.reps >= 1) reps	= s.reps;
if (s.years >= 1) years	= s.years;
if (s.mapInt >= 1) mapInt = s.mapInt;
if (s.traitInt >= 1) traitInt = s.traitInt;
batchMode = s.batchMode; absorbing = s.absorbing;
outRange	= s.outRange; outOccup	= s.outOccup;
outPop	= s.outPop; outInds	= s.outInds;
outGenetics = s.outGenetics;
if (s.outGenType >= 0 && s.outGenType <= 2) {
	outGenType = s.outGenType;
}
outGenXtab = s.outGenXtab;
outTraitsCells = s.outTraitsCells; outTraitsRows = s.outTraitsRows;
outConnect = s.outConnect;
//if (s.outStartRange >= 0) outStartRange =	s.outStartRange;
//if (s.outStartOcc >= 0) outStartOcc =	s.outStartOcc;
if (s.outStartPop >= 0) outStartPop =	s.outStartPop;
if (s.outStartInd >= 0) outStartInd =	s.outStartInd;
if (s.outStartGenetic >= 0) outStartGenetic =	s.outStartGenetic;
if (s.outStartTraitCell >= 0) outStartTraitCell =	s.outStartTraitCell;
if (s.outStartTraitRow >= 0) outStartTraitRow =	s.outStartTraitRow;
if (s.outStartConn >= 0) outStartConn =	s.outStartConn;
if (s.outIntRange >= 1) outIntRange = s.outIntRange;
if (s.outIntOcc >= 1) outIntOcc = s.outIntOcc;
if (s.outIntPop >= 1) outIntPop = s.outIntPop;
if (s.outIntInd >= 1) outIntInd = s.outIntInd;
if (s.outIntGenetic >= 1) outIntGenetic = s.outIntGenetic;
if (s.outIntTraitCell >= 1) outIntTraitCell = s.outIntTraitCell;
if (s.outIntTraitRow >= 1) outIntTraitRow = s.outIntTraitRow;
if (s.outIntConn >= 1) outIntConn = s.outIntConn;
saveMaps = s.saveMaps; saveTraitMaps = s.saveTraitMaps;
saveVisits = s.saveVisits;
outStartPaths = s.outStartPaths;
outIntPaths = s.outIntPaths;
outPaths = s.outPaths;
ReturnPopRaster = s.ReturnPopRaster;
CreatePopFile = s.CreatePopFile;
drawLoaded = s.drawLoaded;
}

simParams paramSim::getSim(void) {
simParams s;
s.batchNum = batchNum;
s.simulation = simulation; s.reps = reps; s.years = years;
s.outRange = outRange; s.outOccup = outOccup; s.outPop = outPop; s.outInds = outInds;
s.outGenetics = outGenetics; s.outGenType = outGenType; s.outGenXtab = outGenXtab;
s.outTraitsCells = outTraitsCells; s.outTraitsRows = outTraitsRows; s.outConnect = outConnect;
//s.outStartRange =	outStartRange;
//s.outStartOcc =	outStartOcc;
s.outStartPop =	outStartPop; s.outStartInd = outStartInd; s.outStartGenetic = outStartGenetic;
s.outStartTraitCell =	outStartTraitCell; s.outStartTraitRow = outStartTraitRow;
s.outStartConn = outStartConn;
s.outIntRange = outIntRange;
s.outIntOcc = outIntOcc; s.outIntPop = outIntPop;
s.outIntInd = outIntInd; s.outIntGenetic = outIntGenetic;
s.outIntTraitCell = outIntTraitCell;
s.outIntTraitRow = outIntTraitRow;
s.outIntConn = outIntConn;
s.batchMode = batchMode;
s.absorbing = absorbing;
s.saveMaps = saveMaps; s.saveTraitMaps = saveTraitMaps;
s.saveVisits = saveVisits;
s.mapInt = mapInt; s.traitInt = traitInt;
s.outStartPaths = outStartPaths;
s.outIntPaths = outIntPaths;
s.outPaths = outPaths;
s.ReturnPopRaster = ReturnPopRaster;
s.CreatePopFile = CreatePopFile;
s.drawLoaded = drawLoaded;
return s;
}

int paramSim::getSimNum(void) { return simulation; }

void paramSim::setViews(simView v) {
viewLand = v.viewLand; viewPatch = v.viewPatch;
viewGrad = v.viewGrad; viewCosts = v.viewCosts;
viewPop = v.viewPop; viewTraits = v.viewTraits;
viewPaths = v.viewPaths; viewGraph = v.viewGraph;
if (v.slowFactor > 0) slowFactor = v.slowFactor;
}

simView paramSim::getViews(void) {
simView v;
v.viewLand = viewLand; v.viewPatch = viewPatch;
v.viewGrad = viewGrad; v.viewCosts = viewCosts;
v.viewPop = viewPop; v.viewTraits = viewTraits;
v.viewPaths = viewPaths; v.viewGraph = viewGraph;
v.slowFactor = slowFactor;
return v;
}

void paramSim::setDir(string s) {
dir = s;
}

// return directory name depending on option specified
string paramSim::getDir(int option) {
string s;
switch (option) {
case 0: // working directory
	s = dir;
	break;
case 1: // Inputs folder
	s = dir + "Inputs/";
	break;
case 2: // Outputs folder
	s = dir + "Outputs/";
	break;
case 3: // Maps folder
	s = dir + "Output_Maps/";
	break;
default:
	s = "ERROR_ERROR_ERROR";
}
return s;
}

bool paramSim::getReturnPopRaster(void) { return ReturnPopRaster; }
bool paramSim::getCreatePopFile(void) { return CreatePopFile; }

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

