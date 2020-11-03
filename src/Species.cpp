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

#pragma hdrstop

#include "Species.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

Species::Species(void)
{
// initialise demographic parameters
repType = 0; nStages = 2;
stageStruct = false;
propMales = 0.5; harem = 1.0; bc	= 1.0; lambda	= 1.5; probRep = 1.0;
repSeasons = 1; 
repInterval = 0; maxAge	= 1000; survival	= 1;
fecDens = false; fecStageDens	 = false;
devDens	 = false; devStageDens	 = false;
survDens	 = false; survStageDens	 = false;
disperseOnLoss = false;
for (int i = 0; i < NSTAGES; i++) {
	for (int j = 0; j < NSEXES; j++) {
		fec[i][j] = 0.0; dev[i][j] = 0.0; surv[i][j] = 0.0; 
		minAge[i][j] = 0;
	}
}
devCoeff = survCoeff = 1.0;
ddwtFec = ddwtDev = ddwtSurv = 0; ddwtFecDim = ddwtDevDim = ddwtSurvDim = 0;
habK = 0; habDimK = 0;
minRK = 1.0; maxRK = 2.0;

// initialise genome attributes
nChromosomes = nTraits = 0;
emigTrait[0] = 0; emigTrait[1] = 0;
movtTrait[0] = 0; movtTrait[1] = 0;
settTrait[0] = 0; settTrait[1] = 0;
//genomeCanRecombine = false;
diploid = true;
neutralMarkers = false;
pleiotropic = false;
trait1Chromosome = true;
probMutn = 0.0001;
probCrossover = 0.0001;
alleleSD = mutationSD = 0.1;
nNLoci = 0;
nLoci = NULL;
traitdata = NULL;
traitnames = NULL;
nTraitNames = 0;

// initialise emigration parameters
densDepEmig = false; stgDepEmig = false; sexDepEmig = false; indVarEmig = false;
emigStage = 0;
for (int i = 0; i < NSTAGES; i++) {
	for (int j = 0; j < NSEXES; j++) {
		 d0[i][j] = 0.0; alphaEmig[i][j] = 0.0; betaEmig[i][j] = 1.0;
	}
}
for (int j = 0; j < NSEXES; j++) {
	d0Mean[0][j] = 0.0; alphaMean[0][j] = 0.0; betaMean[0][j] = 1.0;
	d0SD[0][j] = 0.0; alphaSD[0][j] = 0.0; betaSD[0][j] = 1.0;
}
d0Scale = alphaScale = betaScale = 0.0;

// initialise transfer parameters
moveModel = false; stgDepTrfr = false; sexDepTrfr = false; distMort = false; 
indVarTrfr = false;
twinKern = false;
habMort = false; 
costMap = false;
moveType = 1;
for (int i = 0; i < NSTAGES; i++) {
	for (int j = 0; j < NSEXES; j++) {
		meanDist1[i][j] = 100.0; meanDist2[i][j] = 1000.0; probKern1[i][j] = 0.99;
	}
}
for (int j = 0; j < NSEXES; j++) {
	dist1Mean[0][j] = 100.0; dist1SD[0][j] = 10.0;
	dist2Mean[0][j] = 1000.0; dist2SD[0][j] = 100.0;
	PKern1Mean[0][j] = 0.9; PKern1SD[0][j] = 0.01;
	stepLgthMean[0][j] = 10.0; stepLgthSD[0][j] = 1.0;
	rhoMean[0][j] = 0.9; rhoSD[0][j] = 0.01;
	dpMean[0][j] = 1.0; dpSD[0][j] = 0.1;
	gbMean[0][j] = 1.0; gbSD[0][j] = 0.1;
	alphaDBMean[0][j] = 1.0;  alphaDBSD[0][j] = 0.1;
	betaDBMean[0][j]  = 10.0; betaDBSD[0][j]  = 1.0;
}
pr = 1; prMethod = 1; memSize = 1; goalType = 0;
dp = 1.0; gb = 1.0; alphaDB = 1.0; betaDB = 100000;
stepMort = 0.0; stepLength = 10.0; rho = 0.9;
habStepMort = 0; habCost = 0;
//costMapFile = "NULL";
fixedMort = 0.0; mortAlpha = 0.0; mortBeta = 1.0;
dist1Scale = dist2Scale = PKern1Scale = stepLScale = rhoScale = 0.0;
dpScale = 0.1; gbScale = 0.1; alphaDBScale = 0.1; betaDBScale  = 1.0;
habDimTrfr = 0;
straigtenPath = false;
fullKernel = false;

// initialise settlement parameters
stgDepSett = false; sexDepSett = false; indVarSett = false;
minSteps = 0; maxSteps = 99999999;
for (int i = 0; i < NSTAGES; i++) {
	for (int j = 0; j < NSEXES; j++) {
		densDepSett[i][j] = false; wait[i][j] = false; go2nbrLocn[i][j] = false; findMate[i][j] = false;
		maxStepsYr[i][j] = 99999999;
		s0[i][j] = 1.0; alphaS[i][j] = 0.0; betaS[i][j] = 1.0;
	}
}

// initialise attributes
spNum = 0;

}

Species::~Species() {
// demographic parameters
if (habK != NULL) deleteHabK();
if (ddwtFec  != 0) deleteDDwtFec();
if (ddwtDev  != 0) deleteDDwtDev();
if (ddwtSurv != 0) deleteDDwtSurv();
// transfer parameters
if (habCost != 0 || habStepMort != 0) deleteHabCostMort();
if (nLoci != NULL) deleteLoci();
if (traitdata != NULL) deleteTraitData();
if (traitnames != NULL) deleteTraitNames();
}

short Species::getSpNum(void) { return spNum; }

//---------------------------------------------------------------------------

// Demographic functions

void Species::setDemogr(const demogrParams d) {
if (d.repType >= 0 && d.repType <= 2) repType = d.repType;
if (d.repSeasons >= 1) repSeasons = d.repSeasons;
stageStruct = d.stageStruct;
if (d.propMales > 0.0 && d.propMales < 1.0) propMales = d.propMales;
if (d.harem > 0.0) harem = d.harem;
if (d.bc > 0.0) bc = d.bc;
if (d.lambda > 0.0) lambda = d.lambda;
}

demogrParams Species::getDemogr(void) {
demogrParams d;
d.repType = repType;
d.repSeasons = repSeasons;
d.stageStruct = stageStruct;
d.propMales = propMales;
d.harem = harem;
d.bc = bc;
d.lambda = lambda;
return d;
}

short Species::getRepType(void) { return repType; }

bool Species::stageStructured(void) { return stageStruct; }

void Species::createHabK(short nhab) {
if (nhab >= 0) {
	habDimK = nhab;
	if (habK != 0) deleteHabK();
	habK = new float[nhab];
	for (int i = 0; i < nhab; i++) habK[i] = 0.0;
}
}

void Species::setHabK(short hx,float k) {
if (hx >= 0 && hx < habDimK) {
	if (k >= 0.0) habK[hx] = k;
}
}

float Species::getHabK(short hx) {
float k = 0.0;
if (hx >= 0 && hx < habDimK) k = habK[hx];
return k;
}

float Species::getMaxK(void) {
float k = 0.0;
for (int i = 0; i < habDimK; i++) {
	if (habK[i] > k) k = habK[i];
}
return k;
}

void Species::deleteHabK(void) {
if (habK != 0) {
	delete[] habK; habK = 0;
}
}

void Species::setStage(const stageParams s) {
if (s.nStages > 1) nStages = s.nStages;
if (s.repInterval >= 0) repInterval = s.repInterval;
if (s.maxAge >= 1) maxAge = s.maxAge;
if (s.survival >= 0 && s.survival <= 2) survival = s.survival;
if (s.probRep > 0.0 && s.probRep <= 1.0) probRep = s.probRep;
fecDens = s.fecDens;		fecStageDens = s.fecStageDens;
devDens = s.devDens;		devStageDens = s.devStageDens;
survDens = s.survDens;	survStageDens = s.survStageDens;
disperseOnLoss = s.disperseOnLoss;
}

stageParams Species::getStage(void) {
stageParams s;
s.nStages = nStages; s.repInterval = repInterval; s.maxAge = maxAge;
s.survival = survival; s.probRep = probRep;
s.fecDens = fecDens; s.fecStageDens = fecStageDens;
s.devDens = devDens; s.devStageDens = devStageDens;
s.survDens = survDens; s.survStageDens = survStageDens;
s.disperseOnLoss = disperseOnLoss;
return s;
}

void Species::setFec(short stg,short sex,float f) {
// NB fecundity for stage 0 must always be zero
if (stg > 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES && f >= 0)
	fec[stg][sex] = f;
}

float Species::getFec(short stg,short sex) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
	return fec[stg][sex];
else return 0.0;
}

float Species::getMaxFec(void) {
float maxfec = 0.0;
if (stageStruct) {
	for (int stg = 1; stg < NSTAGES; stg++) {
		if (fec[stg][0] > maxfec) maxfec = fec[stg][0];
	}
}
else maxfec = lambda;
return maxfec;
}

void Species::setDev(short stg,short sex,float d) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES && d >= 0)
	dev[stg][sex] = d;
}

float Species::getDev(short stg,short sex) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
	return dev[stg][sex];
else return 0.0;
}

void Species::setSurv(short stg,short sex,float s) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES && s >= 0)
	surv[stg][sex] = s;
}

float Species::getSurv(short stg,short sex) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
	return surv[stg][sex];
else return 0.0;
}

void Species::setMinAge(short stg,short sex,int age) {
// NB min age for stages 0 & 1 must always be zero
if (stg > 1 && stg < NSTAGES && sex >= 0 && sex < NSEXES && age >= 0)
	minAge[stg][sex] = age;
}

short Species::getMinAge(short stg,short sex) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
	return minAge[stg][sex];
else return 0;
}

void Species::setDensDep(float d, float s) {
if (d > 0.0) devCoeff = d;
if (s > 0.0) survCoeff = s;
}

densDepParams Species::getDensDep(void) {
densDepParams d;
d.devCoeff = devCoeff; d.survCoeff = survCoeff;
return d;
}

void Species::createDDwtFec(short mSize) {
if (mSize >= 0 && mSize < (NSTAGES * NSEXES)) {
	if (ddwtFec != 0) deleteDDwtFec();
	ddwtFecDim = mSize;
	ddwtFec = new float *[mSize];
	for (int i = 0; i < mSize; i++) {
		ddwtFec[i] = new float[mSize];
		for (int j = 0; j < mSize; j++) ddwtFec[i][j] = 1.0;
	}
}
}

void Species::setDDwtFec(short row,short col,float f) {
if (row >= 0 && row < ddwtFecDim && col >= 0 && col < ddwtFecDim)
	ddwtFec[row][col] = f;
}

float Species::getDDwtFec(short row,short col) {
if (row >= 0 && row < ddwtFecDim && col >= 0 && col < ddwtFecDim)
	return ddwtFec[row][col];
else return 0.0;
}

void Species::deleteDDwtFec(void) {
if (ddwtFec != 0) {
	for (int i = 0; i < ddwtFecDim; i++) if (ddwtFec[i] != 0) {
		delete[] ddwtFec[i];
	}
	delete[] ddwtFec; ddwtFec = 0;
}
}

void Species::createDDwtDev(short mSize) {
if (mSize >= 0 && mSize < (NSTAGES * NSEXES)) {
	if (ddwtDev != 0) deleteDDwtDev();
	ddwtDevDim = mSize;
	ddwtDev = new float *[mSize];
	for (int i = 0; i < mSize; i++) {
		ddwtDev[i] = new float[mSize];
		for (int j = 0; j < mSize; j++) ddwtDev[i][j] = 1.0;
	}
}
}

void Species::setDDwtDev(short row,short col,float f) {
if (row >= 0 && row < ddwtDevDim && col >= 0 && col < ddwtDevDim)
	ddwtDev[row][col] = f;
}

float Species::getDDwtDev(short row,short col) {
if (row >= 0 && row < ddwtDevDim && col >= 0 && col < ddwtDevDim)
	return ddwtDev[row][col];
else return 0.0;
}

void Species::deleteDDwtDev(void) {
if (ddwtDev != 0) {
	for (int i = 0; i < ddwtDevDim; i++) if (ddwtDev[i] != 0) {
		delete[] ddwtDev[i];
	}
	delete[] ddwtDev; ddwtDev = 0;
}
}

void Species::createDDwtSurv(short mSize) {
if (mSize >= 0 && mSize < (NSTAGES * NSEXES)) {
	if (ddwtSurv != 0) deleteDDwtSurv();
	ddwtSurvDim = mSize;
	ddwtSurv = new float *[mSize];
	for (int i = 0; i < mSize; i++) {
		ddwtSurv[i] = new float[mSize];
		for (int j = 0; j < mSize; j++) ddwtSurv[i][j] = 1.0;
	}
}
}

void Species::setDDwtSurv(short row,short col,float f) {
if (row >= 0 && row < ddwtSurvDim && col >= 0 && col < ddwtSurvDim)
	ddwtSurv[row][col] = f;
}

float Species::getDDwtSurv(short row,short col) {
if (row >= 0 && row < ddwtSurvDim && col >= 0 && col < ddwtSurvDim)
	return ddwtSurv[row][col];
else return 0.0;
}

void Species::deleteDDwtSurv(void) {
if (ddwtSurv != 0) {
	for (int i = 0; i < ddwtSurvDim; i++) if (ddwtSurv[i] != 0) {
		delete[] ddwtSurv[i];
	}
	delete[] ddwtSurv; ddwtSurv = 0;
}
}

// Functions to handle min/max R or K (under environmental stochasticity)
void Species::setMinMax(float min,float max) {
if (min >= 0.0 && max > min) {
	minRK = min; maxRK = max;
}
}
float Species::getMinMax(short opt) {
if (opt == 0) return minRK;
else return maxRK;
}

//---------------------------------------------------------------------------

// Genome functions

void Species::setGenomeData(genomeData d) {
diploid = d.diploid;
neutralMarkers = d.neutralMarkers;
trait1Chromosome = d.trait1Chromosome;
if (trait1Chromosome) {
	if (d.nLoci > 0) nLoci[0] = d.nLoci;
}
if (d.probMutn >= 0.0 && d.probMutn <= 1.0) probMutn = d.probMutn;
if (d.probCrossover >= 0.0 && d.probCrossover <= 1.0) probCrossover = d.probCrossover;
if (d.alleleSD > 0.0) alleleSD = d.alleleSD;
if (d.mutationSD > 0.0) mutationSD = d.mutationSD;
}

genomeData Species::getGenomeData(void) {
genomeData d;
d.diploid = diploid;
d.neutralMarkers = neutralMarkers;
d.pleiotropic = pleiotropic;
d.trait1Chromosome = trait1Chromosome;
if (nLoci != NULL) d.nLoci = nLoci[0];
else d.nLoci = 0;
d.probMutn	= probMutn;
d.probCrossover	= probCrossover;
d.alleleSD = alleleSD;
d.mutationSD = mutationSD;
return d;
}

bool Species::isDiploid(void) { return diploid; }

// Chromosome functions

void Species::setNChromosomes(int c) {
if (nLoci != NULL) deleteLoci();
if (c > 0) {
	nChromosomes = nNLoci = c;
	nLoci = new short [c];
	for (int i = 0; i < nNLoci; i++) nLoci[i] = 0;
}
else nChromosomes = nNLoci = 0;
}

int Species::getNChromosomes(void) { return nChromosomes; }

void Species::setNLoci(const short chr,const short nloc) {
if (chr >= 0 && chr < nNLoci) {
	if (nloc > 0) nLoci[chr] = nloc;
	else nLoci[chr] = 0;
}
}

int Species::getNLoci(const short chr) {
if (chr >= 0 && chr < nChromosomes) return nLoci[chr];
else return 0;
}

void Species::deleteLoci(void) {
if (nLoci != NULL) { delete[] nLoci; nLoci = NULL; }
}

// Trait functions

// Set 1:1 mapping of trait to chromosome
void Species::set1ChromPerTrait(const int nloc) {
nChromosomes = nTraits;
if (nLoci != NULL) deleteLoci();
nLoci = new short [1];
if (nloc > 0) nLoci[0] = nloc;
else nLoci[0] = 1;
}

bool Species::has1ChromPerTrait(void) { return trait1Chromosome; }

// Set trait attributes for the species
void Species::setTraits(void) {

emigTrait[0] = 0; emigTrait[1] = 0;
movtTrait[0] = 0; movtTrait[1] = 0;
settTrait[0] = 0; settTrait[1] = 0;
nTraits = 0;
#if RSDEBUG
DebugGUI("Species::setTraits(): 0000 nChromosomes=" + Int2Str(nChromosomes)
	+ " nTraits=" + Int2Str(nTraits)
	+ " indVarEmig=" + Int2Str((int)indVarEmig)
	+ " indVarTrfr=" + Int2Str((int)indVarTrfr)
	+ " indVarSett=" + Int2Str((int)indVarSett)
	);
#endif

if (indVarEmig) {
	if (sexDepEmig) {
		if (densDepEmig) nTraits += 6; else nTraits += 2;
	}
	else {
		if (densDepEmig) nTraits += 3; else nTraits += 1;
	}
	emigTrait[0] = 0; emigTrait[1] = nTraits;
}
#if RSDEBUG
//DebugGUI("Species::setTraits(): 1111 nTraits=" + Int2Str(nTraits));
#endif

int movttraits = 0;
if (indVarTrfr) {
	if (moveModel) {
		if (moveType == 1) { // SMS
			movttraits = 1;
			if (goalType == 2) movttraits += 3;
		}
		if (moveType == 2) movttraits = 2;
	}
	else {
		if (sexDepTrfr) {
			if (twinKern) movttraits = 6; else movttraits = 2;
		}
		else {
			if (twinKern) movttraits = 3; else movttraits = 1;
		}
	}
	movtTrait[0] = nTraits; movtTrait[1] = movttraits;
	nTraits += movttraits;
}
#if RSDEBUG
//DebugGUI("Species::setTraits(): 2222 nTraits=" + Int2Str(nTraits));
#endif

int setttraits = 0;
if (indVarSett) {
	if (sexDepSett) setttraits = 6; else setttraits = 3;
	settTrait[0] = nTraits; settTrait[1] = setttraits;
	nTraits += setttraits;
}

setTraitNames();

//if (trait1Chromosome) {
//	nChromosomes = nTraits;
//}
#if RSDEBUG
DebugGUI("Species::setTraits(): 9999 nChromosomes=" + Int2Str(nChromosomes)
	+ " nTraits=" + Int2Str(nTraits));
#endif

}

void Species::setTraitNames(void) {
#if RSDEBUG
//DebugGUI("Species::setTraitNames(): nTraits=" + Int2Str(nTraits)
//	+ " nTraitNames=" + Int2Str(nTraitNames)
//	+ " traitnames=" + Int2Str((int)traitnames)
//	);
//if (traitnames != NULL) {
//	DebugGUI("Species::setTraitNames(): traitnames[0]=" + traitnames[0]
//		);
//	if (nTraits > 1) {
//		DebugGUI("Species::setTraitNames(): traitnames[1]=" + traitnames[1]
//			);
//	}
//}
#endif
deleteTraitNames();
nTraitNames = nTraits;
traitnames = new string [nTraitNames];
int trait = 0;
if (indVarEmig) {
	if (sexDepEmig) {
		if (densDepEmig) {
			traitnames[trait++] = "d0_F";
			traitnames[trait++] = "d0_M";
			traitnames[trait++] = "alpha_F";
			traitnames[trait++] = "alpha_M";
			traitnames[trait++] = "beta_F";
			traitnames[trait++] = "beta_M";
		}
		else {
			traitnames[trait++] = "d0_F";
			traitnames[trait++] = "d0_M";
		}
	}
	else {
		traitnames[trait++] = "d0";
		if (densDepEmig) {
			traitnames[trait++] = "alpha";
			traitnames[trait++] = "beta";
		}
	}
}

if (indVarTrfr) {
	if (moveModel) {
		if (moveType == 1) { // SMS
			traitnames[trait++] = "DP";
			if (goalType == 2) {
				traitnames[trait++] = "GB";
				traitnames[trait++] = "alphaDB";
				traitnames[trait++] = "betaDB";
			}
		}
		if (moveType == 2) { // CRW
			traitnames[trait++] = "stepL";
			traitnames[trait++] = "rho";
		}
	}
	else {
		if (sexDepTrfr) {
			if (twinKern) 
			{
				traitnames[trait++] = "meanDistI_F";
				traitnames[trait++] = "meanDistI_M";
				traitnames[trait++] = "meanDistII_F";
				traitnames[trait++] = "meanDistII_M";
				traitnames[trait++] = "probKernI_F";
				traitnames[trait++] = "probKernI_M";
			}
			else {
				traitnames[trait++] = "meanDistI_F";
				traitnames[trait++] = "meanDistI_M";
			}
		}
		else {
			traitnames[trait++] = "meanDistI";
			if (twinKern) 
			{
				traitnames[trait++] = "meanDistII";
				traitnames[trait++] = "probKernI";
			}
		}
	}
}

if (indVarSett) {
	if (sexDepSett) {
		traitnames[trait++] = "s0_F";
		traitnames[trait++] = "s0_M";
		traitnames[trait++] = "alphaS_F";
		traitnames[trait++] = "alphaS_M";
		traitnames[trait++] = "betaS_F";
		traitnames[trait++] = "betaS_M";
	}
	else {
		traitnames[trait++] = "s0";
		traitnames[trait++] = "alphaS";
		traitnames[trait++] = "betaS";
	}
}
}

void Species::deleteTraitNames(void) {
if (traitnames != NULL) {
#if RSDEBUG
//DebugGUI("Species::deleteTraitNames(): traitnames=" + Int2Str((int)traitnames)
//	);
#endif
	delete[] traitnames;
	traitnames = NULL;
}
}

string Species::getTraitName(const int trait) {
string name = "not used";
if (traitnames != NULL) {
	if (trait >= 0 && trait < nTraits) {
		name = traitnames[trait];
	}
}
return name;
}

int Species::getNTraits(void) { return nTraits; }

void Species::setTraitData(const int ntraits) {
#if RSDEBUG
//DebugGUI(("Species::setTraitData(): traitdata=" + Int2Str((int)traitdata)
//	+ " ntraits=" + Int2Str(ntraits)
//	).c_str());
#endif
deleteTraitData();
traitdata = new traitData;
if (ntraits > 0) {
	traitdata->nTraitMaps = ntraits;
	traitdata->traitmaps = new traitMap *[ntraits];
	for (int i = 0; i < ntraits; i++) {
		traitdata->traitmaps[i] = new traitMap;
	}
}
else { // neutral markers only
	traitdata->nTraitMaps = 0;
}
traitdata->neutralloci = new traitMap;
traitdata->neutralloci->nAlleles = 0;
#if RSDEBUG
//DebugGUI(("Species::setTraitData(): traitdata=" + Int2Str((int)traitdata)
//	+ " nTraitMaps=" + Int2Str(traitdata->nTraitMaps)
//	).c_str());
#endif
}

void Species::deleteTraitData(void) {
if (traitdata != NULL) {
#if RSDEBUG
//DebugGUI(("Species::deleteTraitData(): traitdata=" + Int2Str((int)traitdata)
//	+ " nTraitMaps=" + Int2Str(traitdata->nTraitMaps)
//	).c_str());
#endif
	for (int i = 0; i < traitdata->nTraitMaps; i++) {
		if (traitdata->traitmaps[i]->traitalleles != 0) {
			for (int j = 0; j < traitdata->traitmaps[i]->nAlleles; j++) {
				delete traitdata->traitmaps[i]->traitalleles[j];
			}
		}
		delete[] traitdata->traitmaps[i];
	}
	deleteNeutralLoci();
	delete traitdata;
	traitdata = NULL;
}
}

int Species::getNTraitMaps(void) {
if (traitdata == NULL) return 0;
else return traitdata->nTraitMaps;
}

void Species::setTraitMap(const short trait,const short nalleles) {
traitdata->traitmaps[trait]->nAlleles = nalleles;
traitdata->traitmaps[trait]->traitalleles = new traitAllele *[nalleles];
for (int i = 0; i < nalleles; i++) {
	traitdata->traitmaps[trait]->traitalleles[i] = new traitAllele;
}
}

int Species::getNTraitAlleles(const int trait) {
int nalleles = 0;
if (traitdata != NULL) {
	if (trait >= 0 && trait < traitdata->nTraitMaps) {
		nalleles = traitdata->traitmaps[trait]->nAlleles;
	}
}
return nalleles;
}

void Species::setTraitAllele(const short trait,const short allele,
	const short chr,const short loc)
{
traitdata->traitmaps[trait]->traitalleles[allele] = new traitAllele;
if (chr >= 0 && loc >= 0) {
	traitdata->traitmaps[trait]->traitalleles[allele]->chromo = chr;
	traitdata->traitmaps[trait]->traitalleles[allele]->locus = loc;
}
else {
	traitdata->traitmaps[trait]->traitalleles[allele]->chromo = 0;
	traitdata->traitmaps[trait]->traitalleles[allele]->locus = 0;
}
}

traitAllele Species::getTraitAllele(const short trait,const short allele) {
traitAllele a; a.chromo = 0; a.locus = 0;
if (traitdata != NULL) {
	if (trait >= 0 && trait < traitdata->nTraitMaps) {
		if (allele >= 0 && allele < traitdata->traitmaps[trait]->nAlleles) {
			a = *traitdata->traitmaps[trait]->traitalleles[allele];
		}
	}
}
return a;
}

// Neutral loci functions

// Identify neutral loci and determine whether there is pleiotropy
void Species::setNeutralLoci(bool neutralMarkersOnly) {
bool neutral;
int nneutral = 0;
// find minimum of no. of defined / applied traits
int ntraits;
if (traitdata == 0 ) ntraits = 0;
else ntraits = traitdata->nTraitMaps;
if (ntraits > nTraits) ntraits = nTraits;
#if RSDEBUG
//DebugGUI("Species::setNeutralLoci(): neutralMarkersOnly=" + Int2Str((int)neutralMarkersOnly)
//	+ " nNLoci=" + Int2Str(nNLoci)
//	+ " nTraits=" + Int2Str(nTraits) + " ntraits=" + Int2Str(ntraits)
//);
#endif

// determine no. of neutral loci
deleteNeutralLoci();
for (int i = 0; i < nNLoci; i++) { // each chromosome
	for (int j = 0; j < nLoci[i]; j++) { // each locus
		neutral = true;
		for (int t = 0; t < ntraits; t++) { // each trait
			for (int a = 0; a < traitdata->traitmaps[t]->nAlleles; a++) {
#if RSDEBUG
//DebugGUI("Species::setNeutralLoci(): i=" + Int2Str(i)
//	+ " j=" + Int2Str(j) + " t=" + Int2Str(t) + " a=" + Int2Str(a)
//	+ " chromo=" + Int2Str(traitdata->traitmaps[t]->traitalleles[a]->chromo)
//	+ " locus=" + Int2Str(traitdata->traitmaps[t]->traitalleles[a]->locus)
//);
#endif
				if (i == traitdata->traitmaps[t]->traitalleles[a]->chromo
				&&  j == traitdata->traitmaps[t]->traitalleles[a]->locus) {
#if RSDEBUG
//DebugGUI("Species::setNeutralLoci(): FALSE");
#endif
					neutral = false; // as locus contributes to a trait
					a = 999999;
				}
			}
			if (!neutral) t = 999999;
		}
		if (neutral) nneutral++;
	}
}

traitdata->neutralloci = new traitMap;
traitdata->neutralloci->nAlleles = nneutral;
if (nneutral < 1) return;

// record neutral markers
traitdata->neutralloci->traitalleles = new traitAllele *[nneutral];
nneutral = 0;
for (int i = 0; i < nNLoci; i++) { // each chromosome
	for (int j = 0; j < nLoci[i]; j++) { // each locus
		neutral = true;
		for (int t = 0; t < ntraits; t++) { // each trait
			for (int a = 0; a < traitdata->traitmaps[t]->nAlleles; a++) {
				if (i == traitdata->traitmaps[t]->traitalleles[a]->chromo
				&&  j == traitdata->traitmaps[t]->traitalleles[a]->locus) {
					neutral = false; // as locus contributes to a trait
					a = 999999;
				}
			}
			if (!neutral) t = 999999;
		}
		if (neutral) {
			traitdata->neutralloci->traitalleles[nneutral] = new traitAllele;
			traitdata->neutralloci->traitalleles[nneutral]->chromo = i;
			traitdata->neutralloci->traitalleles[nneutral]->locus = j;
			nneutral++;
		}
	}
}

pleiotropic = false;
if (neutralMarkersOnly) return; // pleiotropy cannot apply

// determine whether there is pleiotropy
int chr,loc;
int nloci = 0; // maximum no. of loci on any one chromosome
for (int i = 0; i < nNLoci; i++) {
	if (nloci < nLoci[i]) nloci = nLoci[i];
}
int ***locfreq;
locfreq = new int **[nNLoci];
for (int i = 0; i < nNLoci; i++) {
	locfreq[i] = new int *[nloci];
	for (int j = 0; j < nloci; j++) {
		locfreq[i][j] = new int[ntraits];
		for (int t = 0; t < ntraits; t++) locfreq[i][j][t] = 0;
	}
}
for (int t = 0; t < ntraits; t++) { // each trait
	for (int a = 0; a < traitdata->traitmaps[t]->nAlleles; a++) {
		chr = traitdata->traitmaps[t]->traitalleles[a]->chromo;
		loc = traitdata->traitmaps[t]->traitalleles[a]->locus;
		locfreq[chr][loc][t]++;
	}
}
#if RSDEBUG
//for (int i = 0; i < nNLoci; i++) {
//	for (int j = 0; j < nloci; j++)
//		for (int t = 0; t < ntraits; t++)
//			DebugGUI("locfreq[" + Int2Str(i) + "][" + Int2Str(j) + "][" + Int2Str(t)
//				+ "]=" + Int2Str(locfreq[i][j][t]));
//}
#endif
for (int i = 0; i < nNLoci; i++) {
	for (int j = 0; j < nloci; j++) {
		// remove multiple contributions of a locus to a particular trait
		// (i.e. prevent recording of pseudo-pleiotropy)
		for (int t = 0; t < ntraits; t++) {
			if (locfreq[i][j][t] > 0) locfreq[i][j][t] = 1;
		}
		// sum at level of chromosome/locus
		for (int t = 1; t < ntraits; t++) {
			locfreq[i][j][0] += locfreq[i][j][t];
		}
		if (locfreq[i][j][0] > 1) pleiotropic = true;
	}
}
for (int i = 0; i < nNLoci; i++) {
	for (int j = 0; j < nloci; j++) {
		delete[] locfreq[i][j];
	}
	delete[] locfreq[i];
}
delete[] locfreq;

}

void Species::deleteNeutralLoci(void) {
if (traitdata->neutralloci != NULL) {
	for (int i = 0; i < traitdata->neutralloci->nAlleles; i++) {
		delete traitdata->neutralloci->traitalleles[i];
	}
	delete[] traitdata->neutralloci;
}
traitdata->neutralloci = NULL;
}

int Species::getNNeutralLoci(void) {
int nn = 0;
if (traitdata != NULL) {
	if (traitdata->neutralloci != NULL) {
		nn = traitdata->neutralloci->nAlleles;
	}
}
return nn;
}

traitAllele Species::getNeutralAllele(const short allele) {
traitAllele a; a.chromo = 0; a.locus = 0;
if (traitdata != NULL) {
	if (allele >= 0 && allele < traitdata->neutralloci->nAlleles) {
		a = *traitdata->neutralloci->traitalleles[allele];
	}
}
return a;
}

//---------------------------------------------------------------------------

// Emigration functions

void Species::setEmig(const emigRules e) {
#if RSDEBUG
//DebugGUI("Species::setEmig(): e.indVar=" + Int2Str((int)e.indVar));
#endif
densDepEmig = e.densDep; stgDepEmig = e.stgDep; sexDepEmig = e.sexDep;
indVarEmig = e.indVar;
if (e.emigStage >= 0) emigStage = e.emigStage;
//setGenome();
}

emigRules Species::getEmig(void) {
emigRules e;
e.densDep = densDepEmig; e.stgDep = stgDepEmig; e.sexDep = sexDepEmig;
e.indVar = indVarEmig; e.emigStage = emigStage;
e.emigTrait[0] = emigTrait[0]; e.emigTrait[1] = emigTrait[1];
return e;
}

void Species::setEmigTraits(const short stg,const short sex,const emigTraits e) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	if (e.d0 >= 0.0 && e.d0 <= 1.0) d0[stg][sex] = e.d0;
	alphaEmig[stg][sex] = e.alpha; betaEmig[stg][sex] = e.beta;
}
}

emigTraits Species::getEmigTraits(short stg,short sex) {
emigTraits e;
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	e.d0 = d0[stg][sex]; e.alpha = alphaEmig[stg][sex]; e.beta = betaEmig[stg][sex];
}
else {
	e.d0 = e.alpha = e.beta = 0.0;
}
return e;
}

float Species::getEmigD0(short stg,short sex) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	return d0[stg][sex];
}
else {
	return 0.0;
}
}

void Species::setEmigParams(const short stg,const short sex,const emigParams e) {
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < NSEXES) // implemented for stage 0 only
{
	if (e.d0Mean >= 0.0 && e.d0Mean < 1.0) d0Mean[stg][sex] = e.d0Mean;
	if (e.d0SD > 0.0 && e.d0SD < 1.0) d0SD[stg][sex] = e.d0SD;
//	if (e.d0MutnSize > 0.0 && e.d0MutnSize < 1.0) d0MutnSize = e.d0MutnSize;
	alphaMean[stg][sex] = e.alphaMean;
	if (e.alphaSD > 0.0) alphaSD[stg][sex] = e.alphaSD;
//	if (e.alphaMutnSize > 0.0) alphaMutnSize = e.alphaMutnSize;
	betaMean[stg][sex] = e.betaMean;
	if (e.betaSD > 0.0) betaSD[stg][sex] = e.betaSD;
//	if (e.betaMutnSize > 0.0) betaMutnSize = e.betaMutnSize;
}
}

emigParams Species::getEmigParams(short stg,short sex) {
emigParams e;
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < NSEXES) // implemented for stage 0 only
{
	e.d0Mean = d0Mean[stg][sex]; e.d0SD = d0SD[stg][sex];
	e.d0Scale = d0Scale;
	e.alphaMean = alphaMean[stg][sex]; e.alphaSD = alphaSD[stg][sex];
	e.alphaScale = alphaScale;
	e.betaMean = betaMean[stg][sex]; e.betaSD = betaSD[stg][sex];
	e.betaScale = betaScale;
}
else {
	e.d0Mean = e.alphaMean = e.betaMean = e.d0SD = e.alphaSD = e.betaSD = 0.0;
	e.d0Scale = d0Scale;
	e.alphaScale = alphaScale;
	e.betaScale = betaScale;
}
return e;
}

void Species::setEmigScales(const emigScales s) {
if (s.d0Scale >= 0.0 && s.d0Scale < 1.0 ) d0Scale = s.d0Scale;
if (s.alphaScale >= 0.0) alphaScale = s.alphaScale;
if (s.betaScale >= 0.0) betaScale = s.betaScale;
}

emigScales Species::getEmigScales(void) {
emigScales s;
s.d0Scale = d0Scale; s.alphaScale = alphaScale; s.betaScale = betaScale;
return s;
}

//---------------------------------------------------------------------------

// Transfer functions

void Species::setTrfr(const trfrRules t) {
#if RSDEBUG
//DebugGUI("Species::setTrfr(): t.indVar=" + Int2Str((int)t.indVar));
#endif
moveModel = t.moveModel; stgDepTrfr = t.stgDep; sexDepTrfr = t.sexDep;
distMort = t.distMort; indVarTrfr = t.indVar;
twinKern = t.twinKern; 
habMort = t.habMort;
moveType = t.moveType; costMap = t.costMap;
//setGenome();
}

trfrRules Species::getTrfr(void) {
trfrRules t;
t.moveModel = moveModel; t.stgDep = stgDepTrfr; t.sexDep = sexDepTrfr;
t.distMort = distMort; t.indVar = indVarTrfr;
t.twinKern = twinKern; 
t.habMort = habMort;
t.moveType = moveType; t.costMap = costMap;
t.movtTrait[0] = movtTrait[0]; t.movtTrait[1] = movtTrait[1];
return t;
}

void Species::setFullKernel(bool k) {
fullKernel = k;
}

bool Species::useFullKernel(void) { return fullKernel; }

void Species::setKernTraits(const short stg,const short sex,
	const trfrKernTraits k,const float resol)
{
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	if (k.meanDist1 > 0.0 && k.meanDist1 >= resol) meanDist1[stg][sex] = k.meanDist1;
	if (k.meanDist2 >= resol) meanDist2[stg][sex] = k.meanDist2;
	if (k.probKern1 > 0.0 && k.probKern1 < 1.0) probKern1[stg][sex] = k.probKern1;
}
}

trfrKernTraits Species::getKernTraits(short stg,short sex) {
trfrKernTraits k;
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	k.meanDist1 = meanDist1[stg][sex];
	k.meanDist2 = meanDist2[stg][sex];
	k.probKern1 = probKern1[stg][sex];
}
else {
	k.meanDist1 = 0.0; k.meanDist2 = 0.0; k.probKern1 = 1.0;
}
return k;
}

void Species::setMortParams(const trfrMortParams m) {
if (m.fixedMort >= 0.0 && m.fixedMort < 1.0) fixedMort = m.fixedMort;
mortAlpha = m.mortAlpha;
mortBeta = m.mortBeta;
}

trfrMortParams Species::getMortParams(void) {
trfrMortParams m;
m.fixedMort = fixedMort; m.mortAlpha = mortAlpha; m.mortBeta = mortBeta;
return m;
}

void Species::setMovtTraits(const trfrMovtTraits m) {
if (m.pr >= 1) pr = m.pr;
if (m.prMethod >= 1 && m.prMethod <= 3) prMethod = m.prMethod;
if (m.memSize >= 1 && m.memSize <= 14) memSize = m.memSize;
if (m.goalType >= 0 && m.goalType <= 2) goalType = m.goalType;
if (m.dp >= 1.0) dp = m.dp;
if (m.gb >= 1.0) gb = m.gb;
if (m.alphaDB > 0.0) alphaDB = m.alphaDB;
if (m.betaDB > 0) betaDB = m.betaDB;
if (m.stepMort >= 0.0 && m.stepMort < 1.0) stepMort = m.stepMort;
if (m.stepLength > 0.0) stepLength = m.stepLength;
if (m.rho > 0.0 && m.rho < 1.0) rho = m.rho;
straigtenPath = m.straigtenPath;
}

trfrMovtTraits Species::getMovtTraits(void) {
trfrMovtTraits m;
m.pr = pr; m.prMethod = prMethod; m.memSize = memSize; m.goalType = goalType;
m.dp = dp; m.gb = gb; m.alphaDB = alphaDB;  m.betaDB = betaDB;
m.stepMort = stepMort; m.stepLength = stepLength; m.rho = rho;
return m;
}

trfrCRWTraits Species::getCRWTraits(void) {
trfrCRWTraits m;
m.stepMort = stepMort; m.stepLength = stepLength; m.rho = rho;
m.straigtenPath = straigtenPath;
return m;
}

trfrSMSTraits Species::getSMSTraits(void) {
trfrSMSTraits m;
m.pr = pr; m.prMethod = prMethod; m.memSize = memSize; m.goalType = goalType;
m.dp = dp; m.gb = gb; m.alphaDB = alphaDB;  m.betaDB = betaDB; m.stepMort = stepMort;
m.straigtenPath = straigtenPath;
return m;
}

void Species::setKernParams(const short stg,const short sex,
	const trfrKernParams k,const float resol)
{
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < NSEXES) // implemented for stage 0 only
{
	if (k.dist1Mean > 0.0 && k.dist1Mean >= resol && k.dist1SD > 0.0) {
		dist1Mean[stg][sex] = k.dist1Mean;  dist1SD[stg][sex] = k.dist1SD;
	}
	if (k.dist2Mean > 0.0 && k.dist2Mean >= resol && k.dist2SD > 0.0) {
		dist2Mean[stg][sex] = k.dist2Mean;  dist2SD[stg][sex] = k.dist2SD;
	}
	if (k.PKern1Mean > 0.0 && k.PKern1Mean < 1.0 && k.PKern1SD > 0.0 && k.PKern1SD < 1.0 ) {
		PKern1Mean[stg][sex] = k.PKern1Mean; PKern1SD[stg][sex] = k.PKern1SD;
	}
}
}

trfrKernParams Species::getKernParams(short stg,short sex) {
trfrKernParams k;
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < NSEXES) // implemented for stage 0 only
{
	k.dist1Mean  = dist1Mean[stg][sex];  k.dist1SD  = dist1SD[stg][sex];
	k.dist2Mean  = dist2Mean[stg][sex];  k.dist2SD  = dist2SD[stg][sex];
	k.PKern1Mean = PKern1Mean[stg][sex]; k.PKern1SD = PKern1SD[stg][sex];
	k.dist1Scale = dist1Scale; k.dist2Scale = dist2Scale; k.PKern1Scale = PKern1Scale;
}
else {
	k.dist1Mean = 100000.0; k.dist1SD = 0.001; k.dist1Scale = 1.0;
	k.dist2Mean = 100000.0; k.dist2SD = 0.001; k.dist2Scale = 1.0;
	k.PKern1Mean = 0.5;  k.PKern1SD = 0.1; k.PKern1Scale = 0.1;
}
return k;
}

void Species::setSMSParams(const short stg,const short sex,const trfrSMSParams s) {
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < 1) // implemented for stage 0 & sex 0 only
{
	if (s.dpMean >= 1.0 && s.dpSD > 0.0) {
		dpMean[stg][sex] = s.dpMean; dpSD[stg][sex] = s.dpSD;
	}
	if (s.gbMean >= 1.0 && s.gbSD > 0.0) {
		gbMean[stg][sex] = s.gbMean; gbSD[stg][sex] = s.gbSD;
	}
	if (s.alphaDBMean > 0.0 && s.alphaDBSD > 0.0) {
		alphaDBMean[stg][sex] = s.alphaDBMean; alphaDBSD[stg][sex] = s.alphaDBSD;
	}
	if (s.betaDBMean >= 1.0 && s.betaDBSD > 0.0) {
		betaDBMean[stg][sex] = s.betaDBMean; betaDBSD[stg][sex] = s.betaDBSD;
	}
}
}

trfrSMSParams Species::getSMSParams(short stg,short sex) {
trfrSMSParams s;
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < 1) // implemented for stage 0 & sex 0 only
{
	s.dpMean = dpMean[stg][sex]; s.dpSD = dpSD[stg][sex];
	s.gbMean = gbMean[stg][sex]; s.gbSD = gbSD[stg][sex];
	s.alphaDBMean = alphaDBMean[stg][sex]; s.alphaDBSD = alphaDBSD[stg][sex];
	s.betaDBMean  = betaDBMean[stg][sex];  s.betaDBSD  = betaDBSD[stg][sex];
	s.dpScale = dpScale; s.gbScale = gbScale;
	s.alphaDBScale = alphaDBScale; s.betaDBScale = betaDBScale;
}
else {
	s.dpMean = 1.0; s.dpSD = 0.1; s.dpScale = 0.1;
	s.gbMean = 1.0; s.gbSD = 0.1; s.gbScale = 0.1;
	s.alphaDBMean = 1.0;  s.alphaDBSD = 0.1; s.alphaDBScale = 0.1;
	s.betaDBMean  = 10.0; s.betaDBSD  = 1.0; s.betaDBScale  = 1.0;
}
return s;
}

void Species::setCRWParams(const short stg,const short sex,const trfrCRWParams m) {
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < 1) // implemented for stage 0 & sex 0 only
{
	if (m.stepLgthMean > 0.0 && m.stepLgthSD > 0.0) {
		stepLgthMean[stg][sex] = m.stepLgthMean; stepLgthSD[stg][sex] = m.stepLgthSD;
	}
	if (m.rhoMean > 0.0 && m.rhoMean < 1.0 && m.rhoSD > 0.0 && m.rhoSD < 1.0 ) {
		rhoMean[stg][sex] = m.rhoMean; rhoSD[stg][sex] = m.rhoSD;
	}
}
}

trfrCRWParams Species::getCRWParams(short stg,short sex) {
trfrCRWParams m;
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < 1) // implemented for stage 0 & sex 0 only
{
	m.stepLgthMean = stepLgthMean[stg][sex]; m.stepLgthSD = stepLgthSD[stg][sex];
	m.rhoMean = rhoMean[stg][sex]; m.rhoSD = rhoSD[stg][sex];
	m.stepLScale = stepLScale; m.rhoScale = rhoScale;
}
else {
	m.stepLgthMean = 1.0; m.stepLgthSD = 0.1; m.stepLScale = 0.1;
	m.rhoMean = 0.5; m.rhoSD = 0.1; m.rhoScale = 0.1;
}
return m;
}

void Species::setTrfrScales(const trfrScales s) {
if (s.dist1Scale >= 0.0) dist1Scale = s.dist1Scale;
if (s.dist2Scale >= 0.0) dist2Scale = s.dist2Scale;
if (s.PKern1Scale > 0.0 && s.PKern1Scale < 1.0) PKern1Scale = s.PKern1Scale;
if (s.dpScale > 0.0) dpScale = s.dpScale;
if (s.gbScale > 0.0) gbScale = s.gbScale;
if (s.alphaDBScale > 0.0) alphaDBScale = s.alphaDBScale;
if (s.betaDBScale > 0.0)  betaDBScale  = s.betaDBScale;
if (s.stepLScale > 0.0) stepLScale = s.stepLScale;
if (s.rhoScale > 0.0 && s.rhoScale < 1.0) rhoScale = s.rhoScale;
}

trfrScales Species::getTrfrScales(void) {
trfrScales s;
s.dist1Scale = dist1Scale; s.dist2Scale = dist2Scale; s.PKern1Scale = PKern1Scale;
s.dpScale = dpScale; s.gbScale = gbScale;
s.alphaDBScale = alphaDBScale; s.betaDBScale = betaDBScale;
s.stepLScale = stepLScale; s.rhoScale = rhoScale;
return s;
}

short Species::getMovtHabDim() { return habDimTrfr; }

void Species::createHabCostMort(short nhab) {
if (nhab >= 0) {
	habDimTrfr = nhab;
	if (habCost != 0 || habStepMort != 0) deleteHabCostMort();
	habCost = new int[nhab];
	habStepMort = new double[nhab];
	for (int i = 0; i < nhab; i++) {
		habCost[i] = 1; habStepMort[i] = 0.0;
	}
}
}

void Species::setHabCost(short hab,int cost) {
if (hab >= 0 && hab < habDimTrfr) {
	if (cost >= 1) habCost[hab] = cost;
}
}

void Species::setHabMort(short hab,double mort) {
if (hab >= 0 && hab < habDimTrfr) {
	if (mort >= 0.0 && mort < 1.0) habStepMort[hab] = mort;
}
}

int Species::getHabCost(short hab) {
int cost = 0;
if (hab >= 0 && hab < habDimTrfr) cost = habCost[hab];
return cost;
}

double Species::getHabMort(short hab) {
double pmort = 0.0;
if (hab >= 0 && hab < habDimTrfr) pmort = habStepMort[hab];
return pmort;
}

void Species::deleteHabCostMort(void) {
if (habCost != 0) {
	delete[] habCost; habCost = 0;
}
if (habStepMort != 0) {
	delete[] habStepMort; habStepMort = 0;
}
}

//---------------------------------------------------------------------------

// Settlement functions

void Species::setSettle(const settleType s) {
stgDepSett = s.stgDep; sexDepSett = s.sexDep; indVarSett = s.indVar;
}

settleType Species::getSettle(void) {
settleType s;
s.stgDep = stgDepSett; s.sexDep = sexDepSett; s.indVar = indVarSett;
s.settTrait[0] = settTrait[0]; s.settTrait[1] = settTrait[1];
return s;
}

void Species::setSettRules(const short stg,const short sex,const settleRules s) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	densDepSett[stg][sex] = s.densDep; wait[stg][sex] = s.wait;
	go2nbrLocn[stg][sex] = s.go2nbrLocn; findMate[stg][sex] = s.findMate;
}
}

settleRules Species::getSettRules(short stg,short sex) {
settleRules s;
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	s.densDep = densDepSett[stg][sex]; s.wait = wait[stg][sex];
	s.go2nbrLocn = go2nbrLocn[stg][sex]; s.findMate = findMate[stg][sex];
}
return s;
}

void Species::setSteps(const short stg,const short sex,const settleSteps s) {
if (stg == 0 && sex == 0) {
	if (s.minSteps >= 0) minSteps	= s.minSteps;
	else minSteps = 0;
	if (s.maxSteps >= 1) maxSteps = s.maxSteps;
	else maxSteps = 99999999;
}
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	if (s.maxStepsYr >= 1) maxStepsYr[stg][sex] = s.maxStepsYr;
	else maxStepsYr[stg][sex] = 99999999;
}
}

settleSteps Species::getSteps(short stg,short sex) {
settleSteps s;
s.minSteps = minSteps;
s.maxSteps = maxSteps;
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) s.maxStepsYr = maxStepsYr[stg][sex];
else s.maxStepsYr = 99999999;
return s;
}

void Species::setSettTraits(const short stg,const short sex,const settleTraits dd) {
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	if (dd.s0 > 0.0 && dd.s0 <= 1.0 ) s0[stg][sex] = dd.s0;
	alphaS[stg][sex] = dd.alpha; betaS[stg][sex] = dd.beta;
}
}

settleTraits Species::getSettTraits(short stg,short sex) {
settleTraits dd;
if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES) {
	dd.s0 = s0[stg][sex]; dd.alpha = alphaS[stg][sex]; dd.beta = betaS[stg][sex];
}
else { dd.s0 = 1.0; dd.alpha = dd.beta = 0.0; }
return dd;
}

void Species::setSettParams(const short stg,const short sex,const settParams s) {
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < NSEXES) // implemented for stage 0 only
{
	if (s.s0Mean >= 0.0 && s.s0Mean < 1.0) s0Mean[stg][sex] = s.s0Mean;
	if (s.s0SD > 0.0 && s.s0SD < 1.0) s0SD[stg][sex] = s.s0SD;
	alphaSMean[stg][sex] = s.alphaSMean;
	if (s.alphaSSD > 0.0) alphaSSD[stg][sex] = s.alphaSSD;
	betaSMean[stg][sex] = s.betaSMean;
	if (s.betaSSD > 0.0) betaSSD[stg][sex] = s.betaSSD;
  if (sex == 0) {
		if (s.s0Scale > 0.0 && s.s0Scale < 1.0) s0Scale = s.s0Scale;
		if (s.alphaSScale > 0.0) alphaSScale = s.alphaSScale;
		if (s.betaSScale > 0.0) betaSScale = s.betaSScale;    
  }
}
}

settParams Species::getSettParams(short stg,short sex) {
settParams s;
//if (stg >= 0 && stg < NSTAGES && sex >= 0 && sex < NSEXES)
if (stg >= 0 && stg < 1 && sex >= 0 && sex < NSEXES) // implemented for stage 0 only
{
	s.s0Mean = s0Mean[stg][sex]; s.s0SD = s0SD[stg][sex];
	s.alphaSMean = alphaSMean[stg][sex]; s.alphaSSD = alphaSSD[stg][sex];
	s.betaSMean = betaSMean[stg][sex]; s.betaSSD = betaSSD[stg][sex];
}
else {
	s.s0Mean = s.alphaSMean = s.betaSMean = s.s0SD = s.alphaSSD = s.betaSSD = 0.0;
}
s.s0Scale = s0Scale;
s.alphaSScale = alphaSScale;
s.betaSScale = betaSScale;
return s;
}

void Species::setSettScales(const settScales s) {
if (s.s0Scale >= 0.0 && s.s0Scale < 1.0 ) s0Scale = s.s0Scale;
if (s.alphaSScale >= 0.0) alphaSScale = s.alphaSScale;
if (s.betaSScale >= 0.0) betaSScale = s.betaSScale;
}

settScales Species::getSettScales(void) {
settScales s;
s.s0Scale = s0Scale; s.alphaSScale = alphaSScale; s.betaSScale = betaSScale;
return s;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

