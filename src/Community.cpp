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

#include "Community.h"

//---------------------------------------------------------------------------


ofstream outrange;
ofstream outoccup,outsuit;
ofstream outtraitsrows;

//---------------------------------------------------------------------------

Community::Community(Landscape *pLand) {
pLandscape = pLand;
indIx = 0;
}

Community::~Community(void) {
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	delete subComms[i];
}
subComms.clear();
}

SubCommunity* Community::addSubComm(Patch *pPch,int num) {
int nsubcomms = (int)subComms.size();
subComms.push_back(new SubCommunity(pPch,num));
return subComms[nsubcomms];
}

void Community::initialise(Species *pSpecies,int year) 
{

int nsubcomms,npatches,ndistcells,spratio,patchnum,rr = 0;
locn distloc;
patchData pch;
patchLimits limits;
intptr ppatch,subcomm;
std::vector <intptr> subcomms;
std::vector <bool> selected;
SubCommunity *pSubComm;
Patch *pPatch;
Cell *pCell;
landParams ppLand = pLandscape->getLandParams();
initParams init = paramsInit->getInit();

nsubcomms = (int)subComms.size();

spratio = ppLand.spResol / ppLand.resol;

#if RSDEBUG
DEBUGLOG << endl << "Community::initialise(): this=" << this
	<< " seedType=" << init.seedType << " freeType=" << init.freeType
	<< " minSeedX=" << init.minSeedX << " minSeedY=" << init.minSeedY
	<< " maxSeedX=" << init.maxSeedX << " maxSeedY=" << init.maxSeedY
	<< " indsFile=" << init.indsFile
	<< " nsubcomms=" << nsubcomms << " spratio=" << spratio
	<< endl;
#endif

switch (init.seedType) {

case 0:	// free initialisation

	switch (init.freeType) {

	case 0:	// random
		// determine no. of patches / cells within the specified initialisation limits
		// and record their corresponding sub-communities in a list
		// parallel list records which have been selected
		npatches = pLandscape->patchCount();
		limits.xMin = init.minSeedX; limits.xMax = init.maxSeedX;
		limits.yMin = init.minSeedY; limits.yMax = init.maxSeedY;
		for (int i = 0; i < npatches; i++) {
			pch = pLandscape->getPatchData(i);
			if (pch.pPatch->withinLimits(limits)) {
				if (ppLand.patchModel) {
					if (pch.pPatch->getPatchNum() != 0) {
						subcomms.push_back(pch.pPatch->getSubComm());
						selected.push_back(false);
					}
				}
				else { // cell-based model - is cell(patch) suitable
					if (pch.pPatch->getK() > 0.0)
					{
						subcomms.push_back(pch.pPatch->getSubComm());
						selected.push_back(false);
					}
				}
			}
		}
		// select specified no. of patches/cells at random
		npatches = (int)subcomms.size();
		if (init.nSeedPatches > npatches/2) { // use backwards selection method
			for (int i = 0; i < npatches; i++) selected[i] = true;
			for (int i = 0; i < (npatches-init.nSeedPatches); i++) {
				do {
					rr = pRandom->IRandom(0,npatches-1);
				} while (!selected[rr]);
				selected[rr] = false;
			}
		}
		else { // use forwards selection method
			for (int i = 0; i < init.nSeedPatches; i++) {
				do {
					rr = pRandom->IRandom(0,npatches-1);
				} while (selected[rr]);
				selected[rr] = true;
			}
		}
		// selected sub-communities for initialisation
		for (int i = 0; i < nsubcomms; i++) { // all sub-communities
			subComms[i]->setInitial(false);
		}
		for (int i = 0; i < npatches; i++) {
			if (selected[i]) {
				pSubComm = (SubCommunity*)subcomms[i];
				pSubComm->setInitial(true);
			}
		}
		break;

	case 1:	// all suitable patches/cells
		npatches = pLandscape->patchCount();
		limits.xMin = init.minSeedX; limits.xMax = init.maxSeedX;
		limits.yMin = init.minSeedY; limits.yMax = init.maxSeedY;
		for (int i = 0; i < npatches; i++) {
			pch = pLandscape->getPatchData(i);
			if (pch.pPatch->withinLimits(limits)) {
				patchnum = pch.pPatch->getPatchNum();
				if (patchnum != 0) {
					if (pch.pPatch->getK() > 0.0) 
					{ // patch is suitable
						subcomm = pch.pPatch->getSubComm();
						if (subcomm == 0) {
							// create a sub-community in the patch
							pSubComm = addSubComm(pch.pPatch,patchnum);
						}
						else {
							pSubComm = (SubCommunity*)subcomm;
						}
						pSubComm->setInitial(true);
					}
				}
			}
		}

		break;

	case 2:	// manually selected patches/cells
		break;

	} // end of switch (init.freeType)
	nsubcomms = (int)subComms.size();
	for (int i = 0; i < nsubcomms; i++) { // all sub-communities
		subComms[i]->initialise(pLandscape,pSpecies);
	}                                   
	break;

case 1:	// from species distribution
	if (ppLand.spDist)
	{
		// deselect all existing sub-communities
		for (int i = 0; i < nsubcomms; i++) {
			subComms[i]->setInitial(false);
		}
		// initialise from loaded species distribution
		switch (init.spDistType) {
		case 0: // all presence cells
			pLandscape->setDistribution(pSpecies,0); // activate all patches
			break;
		case 1: // some randomly selected presence cells
			pLandscape->setDistribution(pSpecies,init.nSpDistPatches); // activate random patches
			break;
		case 2: // manually selected presence cells
			// cells have already been identified - no further action here
			break;
		}

		// THE FOLLOWING WILL HAVE TO BE CHANGED FOR MULTIPLE SPECIES...
		ndistcells = pLandscape->distCellCount(0);
		for (int i = 0; i < ndistcells; i++) {
			distloc = pLandscape->getSelectedDistnCell(0,i);
			if (distloc.x >= 0) { // distribution cell is selected
				// process each landscape cell within the distribution cell
				for (int x = 0; x < spratio; x++) {
					for (int y = 0; y < spratio; y++) {
						pCell = pLandscape->findCell(distloc.x*spratio+x,distloc.y*spratio+y);
						if (pCell != 0) { // not a no-data cell
							ppatch = pCell->getPatch();
							if (ppatch != 0) {
								pPatch = (Patch*)ppatch;
								if (pPatch->getSeqNum() != 0) { // not the matrix patch
									subcomm = pPatch->getSubComm();
									if (subcomm != 0) {
										pSubComm = (SubCommunity*)subcomm;
										pSubComm->setInitial(true);
									}									
								}								
							}
						}
					}
				}
			}
		}
		
		nsubcomms = (int)subComms.size();
		for (int i = 0; i < nsubcomms; i++) { // all sub-communities
			subComms[i]->initialise(pLandscape,pSpecies);
		}
	}
	else {
		// WHAT HAPPENS IF INITIAL DISTRIBUTION IS NOT LOADED ??? ....
		// should not occur - take no action - no initialisation will occur
	}
	break;

case 2:	// initial individuals in specified patches/cells
	if (year < 0) {
		// initialise matrix sub-community only
		subComms[0]->initialise(pLandscape,pSpecies);
		indIx = 0; // reset index for initial individuals
	}
	else { // add any initial individuals for the current year
		initInd iind; iind.year = 0;
		int ninds = paramsInit->numInitInds();
		while (indIx < ninds && iind.year <= year) {
			iind = paramsInit->getInitInd(indIx);
			while (iind.year == year) {
#if RSDEBUG
//DEBUGLOG << "Community::initialise(): year=" << year
//	<< " indIx=" << indIx << " iind.year=" << iind.year
//	<< " iind.patchID=" << iind.patchID << " iind.x=" << iind.x << " iind.y=" << iind.y
//	<< " iind.sex=" << iind.sex << " iind.age=" << iind.age << " iind.stage=" << iind.stage
//	<< endl;
#endif
				if (ppLand.patchModel) {
					if (pLandscape->existsPatch(iind.patchID)) {
						pPatch = pLandscape->findPatch(iind.patchID);
						if (pPatch->getK() > 0.0) 
						{ // patch is suitable
							subcomm = pPatch->getSubComm();
							if (subcomm == 0) {
								// create a sub-community in the patch
								pSubComm = addSubComm(pPatch,iind.patchID);
							}
							else {
								pSubComm = (SubCommunity*)subcomm;
							}
							pSubComm->initialInd(pLandscape,pSpecies,pPatch,pPatch->getRandomCell(),indIx);
						}
					}
				}
				else { // cell-based model
					pCell = pLandscape->findCell(iind.x,iind.y);
					if (pCell != 0) {
						intptr ppatch = pCell->getPatch();
						if (ppatch != 0) {
							pPatch = (Patch*)ppatch;
							if (pPatch->getK() > 0.0) 
							{ // patch is suitable
								subcomm = pPatch->getSubComm();
								if (subcomm == 0) {
									// create a sub-community in the patch
									pSubComm = addSubComm(pPatch,iind.patchID);
								}
								else {
									pSubComm = (SubCommunity*)subcomm;
								}
								pSubComm->initialInd(pLandscape,pSpecies,pPatch,pCell,indIx);
							}
						}
					}
				}
				indIx++;
				if (indIx < ninds) {
					iind = paramsInit->getInitInd(indIx);
				}
				else {
          iind.year = 99999999;
				}
			}
		}
	}
	break;

case 3:	// from file
	// this condition cannot occur here, as init.seedType will have been changed to 0 or 1
	// when the initialisation file was read
	break;

} // end of switch (init.seedType)

#if RSDEBUG
DEBUGLOG << "Community::initialise(): this=" << this
	<< " nsubcomms=" << nsubcomms
	<< endl;
#endif

}

// Add manually selected patches/cells to the selected set for initialisation
void Community::addManuallySelected(void) {
int npatches;
intptr subcomm,patch;
locn initloc;
Cell *pCell;
Patch *pPatch;
SubCommunity *pSubComm;

landParams ppLand = pLandscape->getLandParams();

npatches = pLandscape->initCellCount(); // no. of patches/cells specified
#if RSDEBUG
DEBUGLOG << "Community::addManuallySelected(): this = " << this
	<< " npatches = " << npatches << endl;
#endif
// identify sub-communities to be initialised
if (ppLand.patchModel) {
	for (int i = 0; i < npatches; i++) {
		initloc = pLandscape->getInitCell(i); // patch number held in x-coord of list
		pPatch = pLandscape->findPatch(initloc.x);
		if (pPatch != 0) {
			subcomm = pPatch->getSubComm();
			if (subcomm != 0) {
				pSubComm = (SubCommunity*)subcomm;
				pSubComm->setInitial(true);
			}
		}
	}
}
else { // cell-based model
	for (int i = 0; i < npatches; i++) {
		initloc = pLandscape->getInitCell(i);
		if (initloc.x >= 0 && initloc.x < ppLand.dimX
		&&  initloc.y >= 0 && initloc.y < ppLand.dimY) {
			pCell = pLandscape->findCell(initloc.x,initloc.y);
			if (pCell != 0) { // not no-data cell
				patch = pCell->getPatch();
#if RSDEBUG
DEBUGLOG << "Community::initialise(): i = " << i
	<< " x = " << initloc.x << " y = " << initloc.y
	<< " pCell = " << pCell << " patch = " << patch
	<< endl;
#endif
				if (patch != 0) {
					pPatch = (Patch*)patch;
					subcomm = pPatch->getSubComm();
#if RSDEBUG
DEBUGLOG << "Community::initialise(): i = " << i
	<< " pPatch = " << pPatch << " subcomm = " << subcomm
	<< endl;
#endif
					if (subcomm != 0) {
						pSubComm = (SubCommunity*)subcomm;
						pSubComm->setInitial(true);
#if RSDEBUG
DEBUGLOG << "Community::initialise(): i = " << i
	<< " pSubComm = " << pSubComm
	<< endl;
#endif
					}
				}
			}
		}
	}
}
}

void Community::resetPopns(void) {
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	subComms[i]->resetPopns();
}
// reset the individual ids to start from zero
Individual::indCounter = 0;
}

void Community::localExtinction(int option) {
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	if (subComms[i]->getNum() > 0) { // except in matrix
		subComms[i]->localExtinction(option);
	}
}
}

void Community::patchChanges(void) {
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	if (subComms[i]->getNum() > 0) { // except in matrix
		subComms[i]->patchChange();
	}
}
}

void Community::reproduction(int yr)
{
float eps = 0.0; // epsilon for environmental stochasticity
landParams land = pLandscape->getLandParams();
envStochParams env = paramsStoch->getStoch();
int nsubcomms = (int)subComms.size();
#if RSDEBUG
DEBUGLOG << "Community::reproduction(): this=" << this
	<< " nsubcomms=" << nsubcomms << endl;
#endif

for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	if (env.stoch) {
		if (!env.local) { // global stochasticty
			eps = pLandscape->getGlobalStoch(yr);
		}
	}
	subComms[i]->reproduction(land.resol,eps,land.rasterType,land.patchModel);
}
#if RSDEBUG
DEBUGLOG << "Community::reproduction(): finished" << endl;
#endif
}

void Community::emigration(void) 
{
int nsubcomms = (int)subComms.size();
#if RSDEBUG
DEBUGLOG << "Community::emigration(): this=" << this
	<< " nsubcomms=" << nsubcomms << endl;
#endif
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	subComms[i]->emigration();
}
#if RSDEBUG
DEBUGLOG << "Community::emigration(): finished" << endl;
#endif
}

void Community::dispersal(short landIx,short nextseason)
{
#if RSDEBUG
int t0,t1,t2;
t0 = time(0);
#endif

simParams sim = paramsSim->getSim();

int nsubcomms = (int)subComms.size();
// initiate dispersal - all emigrants leave their natal community and join matrix community
SubCommunity *matrix = subComms[0]; // matrix community is always the first
for (int i = 0; i < nsubcomms; i++) { // all populations
	subComms[i]->initiateDispersal(matrix);
}
#if RSDEBUG
t1 = time(0);
DEBUGLOG << "Community::dispersal(): this=" << this
	<< " nsubcomms=" << nsubcomms << " initiation time=" << t1-t0  << endl;
#endif

// dispersal is undertaken by all individuals now in the matrix patch
// (even if not physically in the matrix)
int ndispersers = 0;
do {
	for (int i = 0; i < nsubcomms; i++) { // all populations
		subComms[i]->resetPossSettlers();
	}
#if RSDEBUG
//DEBUGLOG << "Community::dispersal() 1111: ndispersers=" << ndispersers << endl;
#endif
	ndispersers = matrix->transfer(pLandscape,landIx,nextseason);
#if RSDEBUG
//DEBUGLOG << "Community::dispersal() 2222: ndispersers=" << ndispersers << endl;
#endif
	matrix->completeDispersal(pLandscape,sim.outConnect);
#if RSDEBUG
//DEBUGLOG << "Community::dispersal() 3333: ndispersers=" << ndispersers << endl;
#endif
} while (ndispersers > 0);

#if RSDEBUG
DEBUGLOG << "Community::dispersal(): matrix=" << matrix << endl;
t2 = time(0);
DEBUGLOG << "Community::dispersal(): transfer time=" << t2-t1 << endl;
#endif

#if RSDEBUG
//int t3 = time(0);
//DEBUGLOG << "Community::dispersal(): completion time = " << t3-t2 << endl;
#endif

}

void Community::survival(short part,short option0,short option1) 
{
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all communities (including in matrix)
	subComms[i]->survival(part,option0,option1);
}
}

void Community::ageIncrement(void) {
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all communities (including in matrix)
	subComms[i]->ageIncrement();
}
}

// Calculate total no. of individuals of all species
int Community::totalInds(void) {
popStats p;
int total = 0;
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all communities (including in matrix)
	p = subComms[i]->getPopStats();
	total += p.nInds;
}
return total;
}

// Find the population of a given species in a given patch
Population* Community::findPop(Species *pSp,Patch *pPch) {
Population *pPop = 0;
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all communities (including in matrix)
	pPop = subComms[i]->findPop(pSp,pPch);
	if (pPop != 0) break;
}
return pPop;
}

//---------------------------------------------------------------------------
void Community::createOccupancy(int nrows,int reps) {
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) {
	subComms[i]->createOccupancy(nrows);
}
// Initialise array for occupancy of suitable cells/patches
occSuit = new float *[nrows];
for (int i = 0; i < nrows; i++)
{
	occSuit[i] = new float[reps];
	for (int ii = 0; ii < reps; ii++) occSuit[i][ii] = 0.0;
}
}

void Community::updateOccupancy(int row,int rep) 
{
#if RSDEBUG
DEBUGLOG << "Community::updateOccupancy(): row=" << row << endl;
#endif
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) {
	subComms[i]->updateOccupancy(row);
}

commStats s = getStats();
occSuit[row][rep] = (float)s.occupied / (float)s.suitable;

}

void Community::deleteOccupancy(int nrows) {
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) {
	subComms[i]->deleteOccupancy();
}

for(int i = 0; i < nrows; i++)
	delete[] occSuit[i];
delete[] occSuit;

}

//---------------------------------------------------------------------------
// Count no. of sub-communities (suitable patches) and those occupied (non-zero populations)
// Determine range margins
commStats Community::getStats(void)
{
commStats s;
landParams ppLand = pLandscape->getLandParams();
s.ninds = s.nnonjuvs = s.suitable = s.occupied = 0;
s.minX = ppLand.maxX; s.minY = ppLand.maxY; s.maxX = s.maxY = 0;
float localK;
popStats patchPop;
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	patchPop = subComms[i]->getPopStats();
	s.ninds += patchPop.nInds;
	s.nnonjuvs += patchPop.nNonJuvs;
#if RSDEBUG
//DEBUGLOG << "Community::getStats(): i = " << i
//	<< " pSpecies = " << patchPop.pSpecies << " pPatch = " << patchPop.pPatch
//	<< " nInds = " << patchPop.nInds << endl;
#endif
	if (patchPop.pPatch != 0) { // not the matrix patch
#if RSDEBUG
//DEBUGLOG << "Community::getStats(): i = " << i
//	<< " patchNum = " << patchPop.pPatch->getPatchNum() << endl;
#endif
		if (patchPop.pPatch->getPatchNum() != 0) { // not matrix patch
			localK = patchPop.pPatch->getK();
#if RSDEBUG
//DEBUGLOG << "Community::getStats(): i= " << i
//	<< " pSpecies= " << patchPop.pSpecies << " pPatch= " << patchPop.pPatch
//	<< " patchNum= " << patchPop.pPatch->getPatchNum() << " localK= " << localK
//	<< " nInds= " << patchPop.nInds << " breeding= " << (int)patchPop.breeding
//	<< endl;
#endif
			if (localK > 0.0) s.suitable++;
			if (patchPop.nInds > 0 && patchPop.breeding) {
				s.occupied++;
				patchLimits pchlim = patchPop.pPatch->getLimits();
				if (pchlim.xMin < s.minX) s.minX = pchlim.xMin;
				if (pchlim.xMax > s.maxX) s.maxX = pchlim.xMax;
				if (pchlim.yMin < s.minY) s.minY = pchlim.yMin;
				if (pchlim.yMax > s.maxY) s.maxY = pchlim.yMax;
			}
		}
	}
}
return s;
}

//---------------------------------------------------------------------------

// Functions to control production of output files

// Open population file and write header record
bool Community::outPopHeaders(Species *pSpecies,int option) {
return subComms[0]->outPopHeaders(pLandscape,pSpecies,option);
}

// Write records to population file
void Community::outPop(int rep,int yr,int gen)
{
// generate output for each sub-community (patch) in the community
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	subComms[i]->outPop(pLandscape,rep,yr,gen);
}

}


// Write records to individuals file
void Community::outInds(int rep, int yr, int gen,int landNr) {

if (landNr >= 0) { // open the file
	subComms[0]->outInds(pLandscape,rep,yr,gen,landNr);
	return;
}
if (landNr == -999) { // close the file
	subComms[0]->outInds(pLandscape,rep,yr,gen,-999);
	return;
}
// generate output for each sub-community (patch) in the community
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	subComms[i]->outInds(pLandscape,rep,yr,gen,landNr);
}
}

// Write records to genetics file
void Community::outGenetics(int rep, int yr, int gen,int landNr) {
//landParams ppLand = pLandscape->getLandParams();
if (landNr >= 0) { // open the file
	subComms[0]->outGenetics(rep,yr,gen,landNr);
	return;
}
if (landNr == -999) { // close the file
	subComms[0]->outGenetics(rep,yr,gen,landNr);
	return;
}
// generate output for each sub-community (patch) in the community
int nsubcomms = (int)subComms.size();
for (int i = 0; i < nsubcomms; i++) { // all sub-communities
	subComms[i]->outGenetics(rep,yr,gen,landNr);
}
}

// Open range file and write header record
bool Community::outRangeHeaders(Species *pSpecies,int landNr)
{

if (landNr == -999) { // close the file
	if (outrange.is_open()) outrange.close();
	outrange.clear();
	return true;
}

string name;
landParams ppLand = pLandscape->getLandParams();
envStochParams env = paramsStoch->getStoch();
simParams sim = paramsSim->getSim();

// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
// ATTRIBUTES OF *ALL* SPECIES AS DETECTED AT MODEL LEVEL
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();

#if RSDEBUG
DEBUGLOG << "Community::outRangeHeaders(): simulation=" << sim.simulation
	<< " sim.batchMode=" << sim.batchMode
	<< " landNr=" << landNr << endl;
#endif

if (sim.batchMode) {
	name = paramsSim->getDir(2)
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation) + "_Land"
		+ Int2Str(landNr)
		+ "_Range.txt";
}
else {
	name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation) + "_Range.txt";
}
outrange.open(name.c_str());
outrange << "Rep\tYear\tRepSeason";
if (env.stoch && !env.local) outrange << "\tEpsilon";

outrange << "\tNInds";
if (dem.stageStruct) {
	for (int i = 1; i < sstruct.nStages; i++) outrange << "\tNInd_stage" << i;
	outrange << "\tNJuvs";
}
if (ppLand.patchModel) outrange << "\tNOccupPatches";
else outrange << "\tNOccupCells";
outrange << "\tOccup/Suit\tmin_X\tmax_X\tmin_Y\tmax_Y";

if (emig.indVar) {
	if (emig.sexDep) {
		if (emig.densDep) {
			outrange << "\tF_meanD0\tF_stdD0\tM_meanD0\tM_stdD0";
			outrange << "\tF_meanAlpha\tF_stdAlpha\tM_meanAlpha\tM_stdAlpha";
			outrange << "\tF_meanBeta\tF_stdBeta\tM_meanBeta\tM_stdBeta";
		}
		else {
			outrange << "\tF_meanEP\tF_stdEP\tM_meanEP\tM_stdEP";
		}
	}
	else {
		if (emig.densDep) {
			outrange << "\tmeanD0\tstdD0\tmeanAlpha\tstdAlpha";
			outrange << "\tmeanBeta\tstdBeta";
		}
		else {
			outrange << "\tmeanEP\tstdEP";
		}
	}
}
if (trfr.indVar) {
	if (trfr.moveModel) {
		if (trfr.moveType == 1) {
			outrange << "\tmeanDP\tstdDP\tmeanGB\tstdGB";
				outrange << "\tmeanAlphaDB\tstdAlphaDB\tmeanBetaDB\tstdBetaDB";
		}
		if (trfr.moveType == 2) {
			outrange << "\tmeanStepLength\tstdStepLength\tmeanRho\tstdRho";
		}
	}
	else {
		if (trfr.sexDep) {
			outrange << "\tF_mean_distI\tF_std_distI\tM_mean_distI\tM_std_distI";
			if (trfr.twinKern)
				outrange << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
		}
		else {
			outrange << "\tmean_distI\tstd_distI";
			if (trfr.twinKern)
				outrange << "\tmean_distII\tstd_distII\tmeanPfirstKernel\tstdPfirstKernel";
		}
	}
}
if (sett.indVar) {
	if (sett.sexDep) {
		outrange << "\tF_meanS0\tF_stdS0\tM_meanS0\tM_stdS0";
		outrange << "\tF_meanAlphaS\tF_stdAlphaS\tM_meanAlphaS\tM_stdAlphaS";
		outrange << "\tF_meanBetaS\tF_stdBetaS\tM_meanBetaS\tM_stdBetaS";

	}
	else {
		outrange << "\tmeanS0\tstdS0";
		outrange << "\tmeanAlphaS\tstdAlphaS";
		outrange << "\tmeanBetaS\tstdBetaS";
	}
}
outrange << endl;

#if RSDEBUG
DEBUGLOG << "Community::outRangeHeaders(): finished" << endl;
#endif

return outrange.is_open();
}

// Write record to range file
void Community::outRange(Species *pSpecies,int rep,int yr,int gen)
{
#if RSDEBUG
DEBUGLOG << "Community::outRange(): rep=" << rep
	<< " yr=" << yr << " gen=" << gen << endl;
#endif

landParams ppLand = pLandscape->getLandParams();
envStochParams env = paramsStoch->getStoch();

// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
// ATTRIBUTES OF *ALL* SPECIES AS DETECTED AT MODEL LEVEL
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();

outrange << rep << "\t" << yr << "\t" << gen;
if (env.stoch && !env.local) // write global environmental stochasticity
	outrange << "\t" << pLandscape->getGlobalStoch(yr);

commStats s = getStats();

if (dem.stageStruct) {
	outrange << "\t" << s.nnonjuvs;
	int stagepop;
	int nsubcomms = (int)subComms.size();
	// all non-juvenile stages
	for (int stg = 1; stg < sstruct.nStages; stg++) {
		stagepop = 0;
		for (int i = 0; i < nsubcomms; i++) { // all sub-communities
			stagepop += subComms[i]->stagePop(stg);
		}
		outrange <<"\t" << stagepop;
	}
	// juveniles born in current reproductive season
	stagepop = 0;
	for (int i = 0; i < nsubcomms; i++) { // all sub-communities
		stagepop += subComms[i]->stagePop(0);
	}
	outrange <<"\t" << stagepop;
}
else { // non-structured species
	outrange << "\t" << s.ninds;
}

float occsuit = 0.0;
if (s.suitable > 0) occsuit = (float)s.occupied / (float)s.suitable;
outrange << "\t" << s.occupied << "\t" << occsuit;
// RANGE MINIMA AND MAXIMA NEED TO BECOME A PROPERTY OF THE SPECIES
if (s.ninds > 0) {
	landOrigin origin = pLandscape->getOrigin();
	outrange << "\t" << (float)s.minX * (float)ppLand.resol + origin.minEast
		<< "\t" << (float)(s.maxX+1) * (float)ppLand.resol + origin.minEast
		<< "\t" << (float)s.minY * (float)ppLand.resol + origin.minNorth
		<< "\t" << (float)(s.maxY+1) * (float)ppLand.resol + origin.minNorth;
}
else
	outrange <<"\t0\t0\t0\t0";

if (emig.indVar || trfr.indVar || sett.indVar) { // output trait means
	traitsums ts;   
	traitsums scts; // sub-community traits
	traitCanvas tcanv; 
	int ngenes,popsize;

	tcanv.pcanvas[0] = NULL; 

	for (int i = 0; i < NSEXES; i++) {
		ts.ninds[i] = 0;
		ts.sumD0[i] = ts.ssqD0[i] = 0.0;
		ts.sumAlpha[i] = ts.ssqAlpha[i] = 0.0; ts.sumBeta[i] = ts.ssqBeta[i] = 0.0;
		ts.sumDist1[i] = ts.ssqDist1[i] = 0.0; ts.sumDist2[i] = ts.ssqDist2[i] = 0.0;
		ts.sumProp1[i] = ts.ssqProp1[i] = 0.0;
		ts.sumDP[i] = ts.ssqDP[i] = 0.0;
		ts.sumGB[i] = ts.ssqGB[i] = 0.0;
		ts.sumAlphaDB[i] = ts.ssqAlphaDB[i] = 0.0;
		ts.sumBetaDB[i]  = ts.ssqBetaDB[i]  = 0.0;
		ts.sumStepL[i] = ts.ssqStepL[i] = 0.0; ts.sumRho[i] = ts.ssqRho[i] = 0.0;
		ts.sumS0[i] = ts.ssqS0[i] = 0.0;
		ts.sumAlphaS[i] = ts.ssqAlphaS[i] = 0.0; ts.sumBetaS[i] = ts.ssqBetaS[i] = 0.0;
	}

	int nsubcomms = (int)subComms.size();
	for (int i = 0; i < nsubcomms; i++) { // all sub-communities (incl. matrix)
		scts = subComms[i]->outTraits(tcanv,pLandscape,rep,yr,gen,true); 
		for (int j = 0; j < NSEXES; j++) {
			ts.ninds[j]     += scts.ninds[j];
			ts.sumD0[j]     += scts.sumD0[j];     ts.ssqD0[j]     += scts.ssqD0[j];
			ts.sumAlpha[j]  += scts.sumAlpha[j];  ts.ssqAlpha[j]  += scts.ssqAlpha[j];
			ts.sumBeta[j]   += scts.sumBeta[j];   ts.ssqBeta[j]   += scts.ssqBeta[j];
			ts.sumDist1[j]  += scts.sumDist1[j];  ts.ssqDist1[j]  += scts.ssqDist1[j];
			ts.sumDist2[j]  += scts.sumDist2[j];  ts.ssqDist2[j]  += scts.ssqDist2[j];
			ts.sumProp1[j]  += scts.sumProp1[j];  ts.ssqProp1[j]  += scts.ssqProp1[j];
			ts.sumDP[j]     += scts.sumDP[j];     ts.ssqDP[j]     += scts.ssqDP[j];
			ts.sumGB[j]     += scts.sumGB[j];     ts.ssqGB[j]     += scts.ssqGB[j];
			ts.sumAlphaDB[j] += scts.sumAlphaDB[j]; ts.ssqAlphaDB[j] += scts.ssqAlphaDB[j];
			ts.sumBetaDB[j]  += scts.sumBetaDB[j];  ts.ssqBetaDB[j]  += scts.ssqBetaDB[j];
			ts.sumStepL[j]  += scts.sumStepL[j];  ts.ssqStepL[j]  += scts.ssqStepL[j];
			ts.sumRho[j]    += scts.sumRho[j];    ts.ssqRho[j]    += scts.ssqRho[j];
			ts.sumS0[j]     += scts.sumS0[j];     ts.ssqS0[j]     += scts.ssqS0[j];
			ts.sumAlphaS[j] += scts.sumAlphaS[j]; ts.ssqAlphaS[j] += scts.ssqAlphaS[j];
			ts.sumBetaS[j]  += scts.sumBetaS[j];  ts.ssqBetaS[j]  += scts.ssqBetaS[j];
#if RSDEBUG
//DEBUGLOG << "Community::outRange(): i=" << i << " j=" << j
//	<< " scts.ninds[j]=" << scts.ninds[j]
//	<< " scts.sumD0[j]=" << scts.sumD0[j]
//	<< " scts.ssqD0[j]=" << scts.ssqD0[j]
//	<< endl;
//DEBUGLOG << "Community::outRange(): i=" << i << " j=" << j
//	<< " ts.ninds[j]=" << ts.ninds[j]
//	<< " ts.sumD0[j]=" << ts.sumD0[j] << " ts.ssqD0[j]=" << ts.ssqD0[j]
//	<< endl;
#endif
		}
	}

	if (emig.indVar) {
		if (emig.sexDep) { // must be a sexual species
			ngenes = 2;
		}
		else {
			if (dem.repType == 0) { // asexual reproduction
				ngenes = 1;
			}
			else { // sexual reproduction
				ngenes = 1;
			}
		}
		double mnD0[2],mnAlpha[2],mnBeta[2],sdD0[2],sdAlpha[2],sdBeta[2];
		for (int g = 0; g < ngenes; g++) {
			mnD0[g] = mnAlpha[g] = mnBeta[g] = sdD0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
			// individuals may have been counted by sex if there was
			// sex dependency in another dispersal phase
			if (ngenes == 2) popsize = ts.ninds[g];
			else popsize = ts.ninds[0] + ts.ninds[1];
			if (popsize > 0) {
				mnD0[g] = ts.sumD0[g] / (double)popsize;
				mnAlpha[g] = ts.sumAlpha[g] / (double)popsize;
				mnBeta[g] = ts.sumBeta[g] / (double)popsize;
				if (popsize > 1) {
					sdD0[g] = ts.ssqD0[g]/(double)popsize	- mnD0[g]*mnD0[g];
					if (sdD0[g] > 0.0) sdD0[g] = sqrt(sdD0[g]); else sdD0[g] = 0.0;
					sdAlpha[g] = ts.ssqAlpha[g]/(double)popsize	- mnAlpha[g]*mnAlpha[g];
					if (sdAlpha[g] > 0.0) sdAlpha[g] = sqrt(sdAlpha[g]); else sdAlpha[g] = 0.0;
					sdBeta[g] = ts.ssqBeta[g]/(double)popsize	- mnBeta[g]*mnBeta[g];
					if (sdBeta[g] > 0.0) sdBeta[g] = sqrt(sdBeta[g]); else sdBeta[g] = 0.0;
				}
				else {
					sdD0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
				}
			}
#if RSDEBUG
//DEBUGLOG << "Community::outRange(): ngenes=" << ngenes << " g=" << g
//	<< " ts.ninds[g]=" << ts.ninds[g]
//	<< " ts.sumD0[g]=" << ts.sumD0[g]
//	<< " ts.ssqD0[g]=" << ts.ssqD0[g]
//	<< endl;
//DEBUGLOG << "Community::outRange(): popsize=" << popsize
//	<< " mnD0[g]" << mnD0[g] << " sdD0[g]" << sdD0[g]
//	<< endl;
#endif
		}
		if (emig.sexDep) {
			outrange << "\t" << mnD0[0] << "\t" << sdD0[0];
			outrange << "\t" << mnD0[1] << "\t" << sdD0[1];
			if (emig.densDep) {
				outrange << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
				outrange << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
				outrange << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
				outrange << "\t" << mnBeta[1]  << "\t" << sdBeta[1];
			}
		}
		else { // sex-independent
			outrange << "\t" << mnD0[0] << "\t" << sdD0[0];
			if (emig.densDep) {
				outrange << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
				outrange << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
			}
		}
	}

	if (trfr.indVar) {
		if (trfr.moveModel) {
			// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
			ngenes = 1;
		}
		else {
			if (trfr.sexDep) { // must be a sexual species
				ngenes = 2;
			}
			else {
				ngenes = 1;
			}
		}
		double mnDist1[2],mnDist2[2],mnProp1[2],mnStepL[2],mnRho[2];
		double sdDist1[2],sdDist2[2],sdProp1[2],sdStepL[2],sdRho[2];
		double mnDP[2],mnGB[2],mnAlphaDB[2],mnBetaDB[2];
		double sdDP[2],sdGB[2],sdAlphaDB[2],sdBetaDB[2];
		for (int g = 0; g < ngenes; g++) {
			mnDist1[g] = mnDist2[g] = mnProp1[g] = mnStepL[g] = mnRho[g] = 0.0;
			sdDist1[g] = sdDist2[g] = sdProp1[g] = sdStepL[g] = sdRho[g] = 0.0;
			mnDP[g] = mnGB[g] = mnAlphaDB[g] = mnBetaDB[g] = 0.0;
			sdDP[g] = sdGB[g] = sdAlphaDB[g] = sdBetaDB[g] = 0.0;
			// individuals may have been counted by sex if there was
			// sex dependency in another dispersal phase
			if (ngenes == 2) popsize = ts.ninds[g];
			else popsize = ts.ninds[0] + ts.ninds[1];
			if (popsize > 0) {
				mnDist1[g] = ts.sumDist1[g] / (double)popsize;
				mnDist2[g] = ts.sumDist2[g] / (double)popsize;
				mnProp1[g] = ts.sumProp1[g] / (double)popsize;
				mnStepL[g] = ts.sumStepL[g] / (double)popsize;
				mnRho[g] =   ts.sumRho[g]   / (double)popsize;
				mnDP[g] = ts.sumDP[g] / (double)popsize;
				mnGB[g] = ts.sumGB[g] / (double)popsize;
				mnAlphaDB[g] = ts.sumAlphaDB[g] / (double)popsize;
				mnBetaDB[g]  = ts.sumBetaDB[g]  / (double)popsize;
				if (popsize > 1) {
					sdDist1[g] = ts.ssqDist1[g]/(double)popsize	- mnDist1[g]*mnDist1[g];
					if (sdDist1[g] > 0.0) sdDist1[g] = sqrt(sdDist1[g]); else sdDist1[g] = 0.0;
					sdDist2[g] = ts.ssqDist2[g]/(double)popsize	- mnDist2[g]*mnDist2[g];
					if (sdDist2[g] > 0.0) sdDist2[g] = sqrt(sdDist2[g]); else sdDist2[g] = 0.0;
					sdProp1[g] = ts.ssqProp1[g]/(double)popsize	- mnProp1[g]*mnProp1[g];
					if (sdProp1[g] > 0.0) sdProp1[g] = sqrt(sdProp1[g]); else sdProp1[g] = 0.0;
					sdStepL[g] = ts.ssqStepL[g]/(double)popsize	- mnStepL[g]*mnStepL[g];
					if (sdStepL[g] > 0.0) sdStepL[g] = sqrt(sdStepL[g]); else sdStepL[g] = 0.0;
					sdRho[g] = ts.ssqRho[g]/(double)popsize	- mnRho[g]*mnRho[g];
					if (sdRho[g] > 0.0) sdRho[g] = sqrt(sdRho[g]); else sdRho[g] = 0.0;
					sdDP[g] = ts.ssqDP[g]/(double)popsize	- mnDP[g]*mnDP[g];
					if (sdDP[g] > 0.0) sdDP[g] = sqrt(sdDP[g]); else sdDP[g] = 0.0;
					sdGB[g] = ts.ssqGB[g]/(double)popsize	- mnGB[g]*mnGB[g];
					if (sdGB[g] > 0.0) sdGB[g] = sqrt(sdGB[g]); else sdGB[g] = 0.0;
					sdAlphaDB[g] = ts.ssqAlphaDB[g]/(double)popsize	- mnAlphaDB[g]*mnAlphaDB[g];
					if (sdAlphaDB[g] > 0.0) sdAlphaDB[g] = sqrt(sdAlphaDB[g]); else sdAlphaDB[g] = 0.0;
					sdBetaDB[g]  = ts.ssqBetaDB[g]/(double)popsize	- mnBetaDB[g]*mnBetaDB[g];
					if (sdBetaDB[g] > 0.0) sdBetaDB[g] = sqrt(sdBetaDB[g]); else sdBetaDB[g] = 0.0;
				}
			}
#if RSDEBUG
//DEBUGLOG << "Community::outRange(): ngenes=" << ngenes << " g=" << g
//	<< " ts.ninds[g]=" << ts.ninds[g]
//	<< " ts.sumDP[g]=" << ts.sumDP[g] << " ts.ssqDP[g]=" << ts.ssqDP[g]
//	<< " ts.sumGB[g]=" << ts.sumGB[g] << " ts.ssqGB[g]=" << ts.ssqGB[g]
//	<< endl;
//DEBUGLOG << "Community::outRange(): popsize=" << popsize
//	<< " mnDP[g]" << mnDP[g] << " sdDP[g]" << sdDP[g]
//	<< " mnGB[g]" << mnGB[g] << " sdGB[g]" << sdGB[g]
//	<< endl;
#endif
		}
		if (trfr.moveModel) {
			if (trfr.moveType == 1) {
				outrange << "\t" << mnDP[0] << "\t" << sdDP[0];
				outrange << "\t" << mnGB[0] << "\t" << sdGB[0];
				outrange << "\t" << mnAlphaDB[0] << "\t" << sdAlphaDB[0];
				outrange << "\t" << mnBetaDB[0] << "\t" << sdBetaDB[0];
			}
			if (trfr.moveType == 2) {
				outrange << "\t" << mnStepL[0] << "\t" << sdStepL[0];
				outrange << "\t" << mnRho[0] << "\t" << sdRho[0];
			}
		}
		else {
			if (trfr.sexDep) {
				outrange << "\t" << mnDist1[0] << "\t" << sdDist1[0];
				outrange << "\t" << mnDist1[1] << "\t" << sdDist1[1];
				if (trfr.twinKern) 
				{
					outrange << "\t" << mnDist2[0] << "\t" << sdDist2[0];
					outrange << "\t" << mnDist2[1] << "\t" << sdDist2[1];
					outrange << "\t" << mnProp1[0] << "\t" << sdProp1[0];
					outrange << "\t" << mnProp1[1] << "\t" << sdProp1[1];
				}
			}
			else { // sex-independent
				outrange << "\t" << mnDist1[0] << "\t" << sdDist1[0];
				if (trfr.twinKern) 
				{
					outrange << "\t" << mnDist2[0] << "\t" << sdDist2[0];
					outrange << "\t" << mnProp1[0] << "\t" << sdProp1[0];
				}
			}
		}
	}

	if (sett.indVar) {
		if (sett.sexDep) { // must be a sexual species
			ngenes = 2;
		}
		else {
			if (dem.repType == 0) { // asexual reproduction
				ngenes = 1;
			}
			else { // sexual reproduction
				ngenes = 1;
			}
		}
		// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
		double mnS0[2],mnAlpha[2],mnBeta[2],sdS0[2],sdAlpha[2],sdBeta[2];
		for (int g = 0; g < ngenes; g++) {
			mnS0[g] = mnAlpha[g] = mnBeta[g] = sdS0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
			// individuals may have been counted by sex if there was
			// sex dependency in another dispersal phase
			if (ngenes == 2) popsize = ts.ninds[g];
			else popsize = ts.ninds[0] + ts.ninds[1];
			if (popsize > 0) {
				mnS0[g] = ts.sumS0[g] / (double)popsize;
				mnAlpha[g] = ts.sumAlphaS[g] / (double)popsize;
				mnBeta[g] = ts.sumBetaS[g] / (double)popsize;
				if (popsize > 1) {
					sdS0[g] = ts.ssqS0[g]/(double)popsize	- mnS0[g]*mnS0[g];
					if (sdS0[g] > 0.0) sdS0[g] = sqrt(sdS0[g]); else sdS0[g] = 0.0;
					sdAlpha[g] = ts.ssqAlphaS[g]/(double)popsize	- mnAlpha[g]*mnAlpha[g];
					if (sdAlpha[g] > 0.0) sdAlpha[g] = sqrt(sdAlpha[g]); else sdAlpha[g] = 0.0;
					sdBeta[g] = ts.ssqBetaS[g]/(double)popsize	- mnBeta[g]*mnBeta[g];
					if (sdBeta[g] > 0.0) sdBeta[g] = sqrt(sdBeta[g]); else sdBeta[g] = 0.0;
				}
				else {
					sdS0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
				}
			}
		}
		if (sett.sexDep) {
			outrange << "\t" << mnS0[0] << "\t" << sdS0[0];
			outrange << "\t" << mnS0[1] << "\t" << sdS0[1];
			outrange << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
			outrange << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
			outrange << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
			outrange << "\t" << mnBeta[1]  << "\t" << sdBeta[1];
		}
		else {
			outrange << "\t" << mnS0[0] << "\t" << sdS0[0];
			outrange << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
			outrange << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
		}
	}

}

outrange << endl;
}

// Open occupancy file, write header record and set up occupancy array
bool Community::outOccupancyHeaders(int option)
{
if (option == -999) { // close the files
	if (outsuit.is_open()) outsuit.close();
	if (outoccup.is_open()) outoccup.close();
	outsuit.clear(); outoccup.clear();
	return true;
}

string name, nameI;
simParams sim = paramsSim->getSim();
//demogrParams dem = pSpecies->getDemogr();
landParams ppLand = pLandscape->getLandParams();
int outrows = (sim.years/sim.outIntOcc) + 1;

name = paramsSim->getDir(2);
if (sim.batchMode) {
	name += "Batch" + Int2Str(sim.batchNum) + "_";
	name += "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(ppLand.landNum);
}
else
	name += "Sim" + Int2Str(sim.simulation);
name += "_Occupancy_Stats.txt";
outsuit.open(name.c_str());
outsuit << "Year\tMean_OccupSuit\tStd_error" << endl;

name = paramsSim->getDir(2);
if (sim.batchMode) {
	name += "Batch" + Int2Str(sim.batchNum) + "_";
	name += "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(ppLand.landNum);
}
else
	name += "Sim" + Int2Str(sim.simulation);
name += "_Occupancy.txt";
outoccup.open(name.c_str());
if (ppLand.patchModel) {
	outoccup << "PatchID";
}
else {
	outoccup << "X\tY";
}
for (int i = 0; i < outrows; i++)
	outoccup << "\t" << "Year_" << i*sim.outIntOcc;
outoccup << endl;

// Initialise cells/patches occupancy array
createOccupancy(outrows,sim.reps);

return outsuit.is_open() && outoccup.is_open();
}

void Community::outOccupancy(void) {
landParams ppLand = pLandscape->getLandParams();
simParams sim = paramsSim->getSim();
//streamsize prec = outoccup.precision();
locn loc;

int nsubcomms = (int)subComms.size();
for (int i = 1; i < nsubcomms; i++) { // all except matrix sub-community
	if (ppLand.patchModel) {
		outoccup << subComms[i]->getPatch()->getPatchNum();
	}
	else {
		loc = subComms[i]->getLocn();
		outoccup << loc.x << "\t" << loc.y;
	}
	for (int row = 0; row <= (sim.years/sim.outIntOcc); row++)
	{
		outoccup << "\t" << (double)subComms[i]->getOccupancy(row)/(double)sim.reps;
	}
	outoccup << endl;
}
}

void Community::outOccSuit(bool view) {
double sum,ss,mean,sd,se;
simParams sim = paramsSim->getSim();
//streamsize prec = outsuit.precision();

#if RSDEBUG
//DEBUGLOG << "Community::outOccSuit(): sim.reps=" << sim.reps
//	<< " sim.years=" << sim.years << " sim.outInt=" << sim.outInt << endl;
#endif
for (int i = 0; i < (sim.years/sim.outIntOcc)+1; i++) {
	sum = ss = 0.0;
	for (int rep = 0; rep < sim.reps; rep++) {
		sum += occSuit[i][rep];
		ss  += occSuit[i][rep] * occSuit[i][rep];
#if RSDEBUG
//DEBUGLOG << "Community::outOccSuit(): i=" << i << " rep=" << rep
//	<< " occSuit[i][rep]=" << occSuit[i][rep]
//	<< " sum=" << sum << " ss=" << ss
//	<< endl;
#endif
	}
	mean = sum/(double)sim.reps;
	sd = (ss - (sum*sum/(double)sim.reps)) / (double)(sim.reps-1);
#if RSDEBUG
//DEBUGLOG << "Community::outOccSuit(): i=" << i
//	<< " mean=" << mean << " sd=" << sd << endl;
#endif
	if (sd > 0.0) sd = sqrt(sd);
	else sd = 0.0;
	se = sd / sqrt((double)(sim.reps));
#if RSDEBUG
//DEBUGLOG << "Community::outOccSuit(): i=" << i
//	<< " sd=" << sd << " se=" << se << endl;
#endif

//	outsuit << i*sim.outInt << "\t" << mean << "\t" << se << endl;
//	if (view) viewOccSuit(i*sim.outInt,mean,se);
	outsuit << i*sim.outIntOcc << "\t" << mean << "\t" << se << endl;
	if (view) viewOccSuit(i*sim.outIntOcc,mean,se);
}

}

// Open traits file and write header record
bool Community::outTraitsHeaders(Species *pSpecies,int landNr) {
return subComms[0]->outTraitsHeaders(pLandscape,pSpecies,landNr);
}

// Write records to traits file
/* NOTE: for summary traits by rows, which is permissible for a cell-based landscape
only, this function relies on the fact that subcommunities are created in the same
sequence as patches, which is in asecending order of x nested within descending
order of y
*/
//void Community::outTraits(emigCanvas ecanv,trfrCanvas tcanv,Species *pSpecies,
//	int rep,int yr,int gen)
void Community::outTraits(traitCanvas tcanv,Species *pSpecies, 
	int rep,int yr,int gen)
{
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();
landParams land = pLandscape->getLandParams();
//demogrParams dem = pSpecies->getDemogr();
//emigRules emig = pSpecies->getEmig();
//trfrRules trfr = pSpecies->getTrfr();
//locn loc;
traitsums *ts = 0;
traitsums sctraits;
if (sim.outTraitsRows && yr >= sim.outStartTraitRow && yr%sim.outIntTraitRow == 0) {
	// create array of traits means, etc., one for each row
	ts = new traitsums[land.dimY];
	for (int y = 0; y < land.dimY; y++) {
		for (int i = 0; i < NSEXES; i++) {
			ts[y].ninds[i] = 0;
			ts[y].sumD0[i]     = ts[y].ssqD0[i]     = 0.0;
			ts[y].sumAlpha[i]  = ts[y].ssqAlpha[i]  = 0.0;
			ts[y].sumBeta[i]   = ts[y].ssqBeta[i]   = 0.0;
			ts[y].sumDist1[i]  = ts[y].ssqDist1[i]  = 0.0;
			ts[y].sumDist2[i]  = ts[y].ssqDist2[i]  = 0.0;
			ts[y].sumProp1[i]  = ts[y].ssqProp1[i]  = 0.0;
			ts[y].sumStepL[i]  = ts[y].ssqStepL[i]  = 0.0;
			ts[y].sumRho[i]    = ts[y].ssqRho[i]    = 0.0;
			ts[y].sumS0[i]     = ts[y].ssqS0[i]     = 0.0;
			ts[y].sumAlphaS[i] = ts[y].ssqAlphaS[i] = 0.0;
			ts[y].sumBetaS[i]  = ts[y].ssqBetaS[i]  = 0.0;
		}
	}
}
if (v.viewTraits
|| ((sim.outTraitsCells && yr >= sim.outStartTraitCell && yr%sim.outIntTraitCell == 0) ||
		(sim.outTraitsRows && yr >= sim.outStartTraitRow && yr%sim.outIntTraitRow == 0)))
{
	// generate output for each sub-community (patch) in the community
	int nsubcomms = (int)subComms.size();
	for (int i = 1; i < nsubcomms; i++) { // // all except matrix sub-community
//		sctraits = subComms[i]->outTraits(ecanv,tcanv,pLandscape,rep,yr,gen);
//		sctraits = subComms[i]->outTraits(tcanv,pLandscape,rep,yr,gen,v.viewGrad);
		sctraits = subComms[i]->outTraits(tcanv,pLandscape,rep,yr,gen,false); 
		locn loc = subComms[i]->getLocn();
		int y = loc.y;
		if (sim.outTraitsRows && yr >= sim.outStartTraitRow && yr%sim.outIntTraitRow == 0)
		{
			for (int s = 0; s < NSEXES; s++) {
				ts[y].ninds[s]     += sctraits.ninds[s];
				ts[y].sumD0[s]     += sctraits.sumD0[s];     ts[y].ssqD0[s]     += sctraits.ssqD0[s];
				ts[y].sumAlpha[s]  += sctraits.sumAlpha[s];  ts[y].ssqAlpha[s]  += sctraits.ssqAlpha[s];
				ts[y].sumBeta[s]   += sctraits.sumBeta[s];   ts[y].ssqBeta[s]   += sctraits.ssqBeta[s];
				ts[y].sumDist1[s]  += sctraits.sumDist1[s];  ts[y].ssqDist1[s]  += sctraits.ssqDist1[s];
				ts[y].sumDist2[s]  += sctraits.sumDist2[s];  ts[y].ssqDist2[s]  += sctraits.ssqDist2[s];
				ts[y].sumProp1[s]  += sctraits.sumProp1[s];  ts[y].ssqProp1[s]  += sctraits.ssqProp1[s];
				ts[y].sumStepL[s]  += sctraits.sumStepL[s];  ts[y].ssqStepL[s]  += sctraits.ssqStepL[s];
				ts[y].sumRho[s]    += sctraits.sumRho[s];    ts[y].ssqRho[s]    += sctraits.ssqRho[s];
				ts[y].sumS0[s]     += sctraits.sumS0[s];     ts[y].ssqS0[s]     += sctraits.ssqS0[s];
				ts[y].sumAlphaS[s] += sctraits.sumAlphaS[s]; ts[y].ssqAlphaS[s] += sctraits.ssqAlphaS[s];
				ts[y].sumBetaS[s]  += sctraits.sumBetaS[s];  ts[y].ssqBetaS[s]  += sctraits.ssqBetaS[s];
			}
		}
	}
	if (nsubcomms > 0 && sim.outTraitsRows
	&& yr >= sim.outStartTraitRow && yr%sim.outIntTraitRow == 0) {
		for (int y = 0; y < land.dimY; y++) {
			if ((ts[y].ninds[0]+ts[y].ninds[1]) > 0) {
				writeTraitsRows(pSpecies,rep,yr,gen,y,ts[y]);
			}
		}
	}
}
//if (sim.outTraitsRows && yr >= sim.outStartTraitRow && yr%sim.outIntTraitRow == 0)
//{
//	if (ts != 0) delete[] ts;
//}
if (ts != 0) { delete[] ts; ts = 0; }
}

// Write records to trait rows file
void Community::writeTraitsRows(Species *pSpecies,int rep,int yr,int gen,int y,
	traitsums ts)
{
//simParams sim = paramsSim->getSim();
//simView v = paramsSim->getViews();
//landData land = pLandscape->getLandData();
//landOrigin origin = pLandscape->getOrigin();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
double mn,sd;

// calculate population size in case one phase is sex-dependent and the other is not
// (in which case numbers of individuals are recorded by sex)
int popsize = ts.ninds[0] + ts.ninds[1];
outtraitsrows << rep << "\t" << yr << "\t" << gen
	<< "\t" << y;
//	<< "\t"	<< y*land.resol + origin.minNorth;
if ((emig.indVar && emig.sexDep) || (trfr.indVar && trfr.sexDep))
	outtraitsrows << "\t" << ts.ninds[0] << "\t" << ts.ninds[1];
else
	outtraitsrows << "\t" << popsize;

if (emig.indVar) {
	if (emig.sexDep) {
		if (ts.ninds[0] > 0) mn = ts.sumD0[0]/(double)ts.ninds[0]; else mn = 0.0;
		if (ts.ninds[0] > 1) sd = ts.ssqD0[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		outtraitsrows << "\t" << mn << "\t" << sd;
		if (ts.ninds[1] > 0) mn = ts.sumD0[1]/(double)ts.ninds[1]; else mn = 0.0;
		if (ts.ninds[1] > 1) sd = ts.ssqD0[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		outtraitsrows << "\t" << mn << "\t" << sd;
		if (emig.densDep) {
			if (ts.ninds[0] > 0) mn = ts.sumAlpha[0]/(double)ts.ninds[0]; else mn = 0.0;
			if (ts.ninds[0] > 1) sd = ts.ssqAlpha[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (ts.ninds[1] > 0) mn = ts.sumAlpha[1]/(double)ts.ninds[1]; else mn = 0.0;
			if (ts.ninds[1] > 1) sd = ts.ssqAlpha[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (ts.ninds[0] > 0) mn = ts.sumBeta[0]/(double)ts.ninds[0]; else mn = 0.0;
			if (ts.ninds[0] > 1) sd = ts.ssqBeta[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (ts.ninds[1] > 0) mn = ts.sumBeta[1]/(double)ts.ninds[1]; else mn = 0.0;
			if (ts.ninds[1] > 1) sd = ts.ssqBeta[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
		}
	}
	else { // no sex dependence in emigration
		if (popsize > 0) mn = ts.sumD0[0]/(double)popsize; else mn = 0.0;
		if (popsize > 1) sd = ts.ssqD0[0]/(double)popsize - mn*mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		outtraitsrows << "\t" << mn << "\t" << sd;
		if (emig.densDep) {
			if (popsize > 0) mn = ts.sumAlpha[0]/(double)popsize; else mn = 0.0;
			if (popsize > 1) sd = ts.ssqAlpha[0]/(double)popsize - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (popsize > 0) mn = ts.sumBeta[0]/(double)popsize; else mn = 0.0;
			if (popsize > 1) sd = ts.ssqBeta[0]/(double)popsize - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
		}
	}
}

if (trfr.indVar) {
	if (trfr.moveModel) {
		if (trfr.moveType == 2) { // CRW
			// NB - CURRENTLY CANNOT BE SEX-DEPENDENT...
			if (popsize > 0) mn = ts.sumStepL[0]/(double)popsize; else mn = 0.0;
			if (popsize > 1) sd = ts.ssqStepL[0]/(double)popsize - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (popsize > 0) mn = ts.sumRho[0]/(double)popsize; else mn = 0.0;
			if (popsize > 1) sd = ts.ssqRho[0]/(double)popsize - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
//			if (ts.ninds[0] > 0) mn = ts.sumStepL[0]/(double)ts.ninds[0]; else mn = 0.0;
//			if (ts.ninds[0] > 1) sd = ts.ssqStepL[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
//			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//			outtraitsrows << "\t" << mn << "\t" << sd;
//			if (ts.ninds[0] > 0) mn = ts.sumRho[0]/(double)ts.ninds[0]; else mn = 0.0;
//			if (ts.ninds[0] > 1) sd = ts.ssqRho[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
//			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//			outtraitsrows << "\t" << mn << "\t" << sd;
		}
	}
	else { // dispersal kernel
		if (trfr.sexDep) {
			if (ts.ninds[0] > 0) mn = ts.sumDist1[0]/(double)ts.ninds[0]; else mn = 0.0;
			if (ts.ninds[0] > 1) sd = ts.ssqDist1[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (ts.ninds[1] > 0) mn = ts.sumDist1[1]/(double)ts.ninds[1]; else mn = 0.0;
			if (ts.ninds[1] > 1) sd = ts.ssqDist1[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (trfr.twinKern) 
			{
				if (ts.ninds[0] > 0) mn = ts.sumDist2[0]/(double)ts.ninds[0]; else mn = 0.0;
				if (ts.ninds[0] > 1) sd = ts.ssqDist2[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				outtraitsrows << "\t" << mn << "\t" << sd;
				if (ts.ninds[1] > 0) mn = ts.sumDist2[1]/(double)ts.ninds[1]; else mn = 0.0;
				if (ts.ninds[1] > 1) sd = ts.ssqDist2[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				outtraitsrows << "\t" << mn << "\t" << sd;
				if (ts.ninds[0] > 0) mn = ts.sumProp1[0]/(double)ts.ninds[0]; else mn = 0.0;
				if (ts.ninds[0] > 1) sd = ts.ssqProp1[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				outtraitsrows << "\t" << mn << "\t" << sd;
				if (ts.ninds[1] > 0) mn = ts.sumProp1[1]/(double)ts.ninds[1]; else mn = 0.0;
				if (ts.ninds[1] > 1) sd = ts.ssqProp1[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				outtraitsrows << "\t" << mn << "\t" << sd;
			}
		}
		else { // sex-independent
			if (popsize > 0) mn = ts.sumDist1[0]/(double)popsize; else mn = 0.0;
			if (popsize > 1) sd = ts.ssqDist1[0]/(double)popsize - mn*mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			outtraitsrows << "\t" << mn << "\t" << sd;
			if (trfr.twinKern) 
			{
				if (popsize > 0) mn = ts.sumDist2[0]/(double)popsize; else mn = 0.0;
				if (popsize > 1) sd = ts.ssqDist2[0]/(double)popsize - mn*mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				outtraitsrows << "\t" << mn << "\t" << sd;
				if (popsize > 0) mn = ts.sumProp1[0]/(double)popsize; else mn = 0.0;
				if (popsize > 1) sd = ts.ssqProp1[0]/(double)popsize - mn*mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				outtraitsrows << "\t" << mn << "\t" << sd;
			}
		}
	}
}

if (sett.indVar) {
	// NB - CURRENTLY CANNOT BE SEX-DEPENDENT...
//	if (sett.sexDep) {
//		if (ts.ninds[0] > 0) mn = ts.sumS0[0]/(double)ts.ninds[0]; else mn = 0.0;
//		if (ts.ninds[0] > 1) sd = ts.ssqS0[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
//		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//		outtraitsrows << "\t" << mn << "\t" << sd;
//		if (ts.ninds[1] > 0) mn = ts.sumS0[1]/(double)ts.ninds[1]; else mn = 0.0;
//		if (ts.ninds[1] > 1) sd = ts.ssqS0[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
//		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//		outtraitsrows << "\t" << mn << "\t" << sd;
//		if (ts.ninds[0] > 0) mn = ts.sumAlphaS[0]/(double)ts.ninds[0]; else mn = 0.0;
//		if (ts.ninds[0] > 1) sd = ts.ssqAlphaS[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
//		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//		outtraitsrows << "\t" << mn << "\t" << sd;
//		if (ts.ninds[1] > 0) mn = ts.sumAlphaS[1]/(double)ts.ninds[1]; else mn = 0.0;
//		if (ts.ninds[1] > 1) sd = ts.ssqAlphaS[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
//		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//		outtraitsrows << "\t" << mn << "\t" << sd;
//		if (ts.ninds[0] > 0) mn = ts.sumBetaS[0]/(double)ts.ninds[0]; else mn = 0.0;
//		if (ts.ninds[0] > 1) sd = ts.ssqBetaS[0]/(double)ts.ninds[0] - mn*mn; else sd = 0.0;
//		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//		outtraitsrows << "\t" << mn << "\t" << sd;
//		if (ts.ninds[1] > 0) mn = ts.sumBetaS[1]/(double)ts.ninds[1]; else mn = 0.0;
//		if (ts.ninds[1] > 1) sd = ts.ssqBetaS[1]/(double)ts.ninds[1] - mn*mn; else sd = 0.0;
//		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
//		outtraitsrows << "\t" << mn << "\t" << sd;
//	}
//	else { // no sex dependence in settlement
		if (popsize > 0) mn = ts.sumS0[0]/(double)popsize; else mn = 0.0;
		if (popsize > 1) sd = ts.ssqS0[0]/(double)popsize - mn*mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		outtraitsrows << "\t" << mn << "\t" << sd;
		if (popsize > 0) mn = ts.sumAlphaS[0]/(double)popsize; else mn = 0.0;
		if (popsize > 1) sd = ts.ssqAlphaS[0]/(double)popsize - mn*mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		outtraitsrows << "\t" << mn << "\t" << sd;
		if (popsize > 0) mn = ts.sumBetaS[0]/(double)popsize; else mn = 0.0;
		if (popsize > 1) sd = ts.ssqBetaS[0]/(double)popsize - mn*mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		outtraitsrows << "\t" << mn << "\t" << sd;
//	}
}

outtraitsrows << endl;
}

// Open trait rows file and write header record
bool Community::outTraitsRowsHeaders(Species *pSpecies,int landNr) {

if (landNr == -999) { // close file
	if (outtraitsrows.is_open()) outtraitsrows.close();
	outtraitsrows.clear();
	return true;
}

string name;
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
simParams sim = paramsSim->getSim();

string DirOut = paramsSim->getDir(2);
if (sim.batchMode) {
	name = DirOut
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(landNr) + "_TraitsXrow.txt";
}
else {
	name = DirOut + "Sim" + Int2Str(sim.simulation) + "_TraitsXrow.txt";
}
outtraitsrows.open(name.c_str());

outtraitsrows << "Rep\tYear\tRepSeason\ty";
if ((emig.indVar && emig.sexDep) || (trfr.indVar && trfr.sexDep))
	outtraitsrows << "\tN_females\tN_males";
else
	outtraitsrows << "\tN";

if (emig.indVar) {
	if (emig.sexDep) {
		if (emig.densDep) {
			outtraitsrows << "\tF_meanD0\tF_stdD0\tM_meanD0\tM_stdD0";
			outtraitsrows << "\tF_meanAlpha\tF_stdAlpha\tM_meanAlpha\tM_stdAlpha";
			outtraitsrows << "\tF_meanBeta\tF_stdBeta\tM_meanBeta\tM_stdBeta";
		}
		else {
			outtraitsrows << "\tF_meanEP\tF_stdEP\tM_meanEP\tM_stdEP";
		}
	}
	else {
		if (emig.densDep) {
			outtraitsrows << "\tmeanD0\tstdD0\tmeanAlpha\tstdAlpha";
			outtraitsrows << "\tmeanBeta\tstdBeta";
		}
		else {
			outtraitsrows << "\tmeanEP\tstdEP";
		}
	}
}
if (trfr.indVar) {
	if (trfr.moveModel) {
		if (trfr.moveType == 2) {
			outtraitsrows << "\tmeanStepLength\tstdStepLength\tmeanRho\tstdRho";
		}
	}
	else { // dispersal kernel
		if (trfr.sexDep) {
			outtraitsrows << "\tF_mean_distI\tF_std_distI\tM_mean_distI\tM_std_distI";
			if (trfr.twinKern)
				outtraitsrows << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
		}
		else {
			outtraitsrows << "\tmean_distI\tstd_distI";
			if (trfr.twinKern)
				outtraitsrows << "\tmean_distII\tstd_distII\tmeanPfirstKernel\tstdPfirstKernel";
		}
	}
}

if (sett.indVar) {
//	if (sett.sexDep) {
//		outtraitsrows << "\tF_meanS0\tF_stdS0\tM_meanS0\tM_stdS0";
//		outtraitsrows << "\tF_meanAlphaS\tF_stdAlphaS\tM_meanAlphaS\tM_stdAlphaS";
//		outtraitsrows << "\tF_meanBetaS\tF_stdBetaS\tM_meanBetaS\tM_stdBetaS";
//	}
//	else {
		outtraitsrows << "\tmeanS0\tstdS0";
		outtraitsrows << "\tmeanAlphaS\tstdAlphaS";
		outtraitsrows << "\tmeanBetaS\tstdBetaS";
//	}
}
outtraitsrows << endl;

return outtraitsrows.is_open();

}

Rcpp::IntegerMatrix Community::addYearToPopList(int rep, int yr) {  // TODO: define new simparams to control start and interval of output

	landParams ppLand = pLandscape->getLandParams();
	Rcpp::IntegerMatrix pop_map_year(ppLand.dimY,ppLand.dimX);
	intptr patch = 0;
	Patch* pPatch = 0;
	intptr subcomm = 0;
	SubCommunity *pSubComm = 0;
	popStats pop;
	//pop.breeding = false;
	pop.nInds = pop.nAdults = pop.nNonJuvs = 0;

	for (int y = 0; y < ppLand.dimY; y++) {
		for (int x = 0; x < ppLand.dimX; x++) {
			Cell *pCell = pLandscape->findCell(x,y); //if (pLandscape->cells[y][x] == 0) {
			if (pCell == 0) { // no-data cell
				pop_map_year(ppLand.dimY-1-y,x) = NA_INTEGER;
			} else {
				patch = pCell->getPatch();
				if (patch == 0) { // matrix cell
					pop_map_year(ppLand.dimY-1-y,x) = 0;
				}
				else{
					pPatch = (Patch*)patch;
					subcomm = pPatch->getSubComm();
					if (subcomm == 0) { // check if sub-community exists
						pop_map_year(ppLand.dimY-1-y,x) = 0;
					} else {
						pSubComm = (SubCommunity*)subcomm;
						pop = pSubComm->getPopStats();
						pop_map_year(ppLand.dimY-1-y,x) = pop.nInds; // use indices like this because matrix gets transposed upon casting it into a raster on R-level
						//pop_map_year(ppLand.dimY-1-y,x) = pop.nAdults;
					}
				}
			}
		}
	}
	//list_outPop.push_back(pop_map_year, "rep" + std::to_string(rep) + "_year" + std::to_string(yr));
    return pop_map_year;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
