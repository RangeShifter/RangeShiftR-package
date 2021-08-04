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

#include "SubCommunity.h"
//---------------------------------------------------------------------------

ofstream outtraits;

//---------------------------------------------------------------------------

SubCommunity::SubCommunity(Patch* pPch, int num) {
subCommNum = num;
pPatch = pPch;
// record the new sub-community no. in the patch
pPatch->setSubComm((intptr)this);
initial = false;
occupancy = 0;
}

SubCommunity::~SubCommunity() {
pPatch->setSubComm(0);
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	delete popns[i];
}
popns.clear();
if (occupancy != 0) delete[] occupancy;
}

intptr SubCommunity::getNum(void) { return subCommNum; }

Patch* SubCommunity::getPatch(void) { return pPatch; }

locn SubCommunity::getLocn(void) {
locn loc = pPatch->getCellLocn(0);
return loc;
}

void SubCommunity::setInitial(bool b) { initial = b; }

void SubCommunity::initialise(Landscape *pLandscape,Species *pSpecies) 
{
//patchLimits limits;
//locn loc;
int ncells;
landParams ppLand = pLandscape->getLandParams();
initParams init = paramsInit->getInit();
#if RSDEBUG
//DEBUGLOG << "SubCommunity::initialise(): subCommNum=" << subCommNum
//	<< " seedType="<< init.seedType
//	<< " popns.size()="<< popns.size()
//	<< endl;
#endif
// determine size of initial population
//int hx,nInds;
int nInds = 0;
if (subCommNum == 0 // matrix patch
||  !initial)   		// not in initial region or distribution
	nInds = 0;
else {
	float k = pPatch->getK();
	if (k > 0.0) { // patch is currently suitable for this species
		switch (init.initDens) {
		case 0: // at carrying capacity
			nInds = (int)k;
			break;
		case 1: // at half carrying capacity
			nInds = (int)(k/2.0);
			break;
		case 2: // specified no. per cell or density
			ncells = pPatch->getNCells();
			if (ppLand.patchModel) {
				nInds = (int)(init.indsHa * (float)(ncells*ppLand.resol*ppLand.resol) / 10000.0);      
			}
			else {
				nInds = init.indsCell * ncells;
			}
			break;
		}
	}
	else nInds = 0;
}

// create new population (even if it has no individuals)
//popns.push_back(new Population(pSpecies,pPatch,nInds));
//newPopn(pSpecies,pPatch,nInds);

// create new population only if it is non-zero or the matrix popn
if (subCommNum == 0 || nInds > 0) {
	newPopn(pLandscape,pSpecies,pPatch,nInds);
}

}

// initialise a specified individual
void SubCommunity::initialInd(Landscape *pLandscape,Species *pSpecies,
	Patch *pPatch,Cell *pCell,int ix) 
{

demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
genomeData gen = pSpecies->getGenomeData();
short stg,age,repInt;
Individual *pInd;
float probmale;
//	 bool movt;
//	 short moveType;

// create new population if not already in existence
int npopns = (int)popns.size();
if (npopns < 1) {
	newPopn(pLandscape,pSpecies,pPatch,0);
}

// create new individual
initInd iind = paramsInit->getInitInd(ix);
if (dem.stageStruct) {
	stg = iind.stage; age = iind.age;   repInt = sstruct.repInterval;
}
else {
	age = stg = 1;   repInt = 0;
}
if (dem.repType == 0) {
	probmale = 0.0;
}
else {
	if (iind.sex == 1) probmale = 1.0; else probmale = 0.0;
}
//if (trfr.moveModel) {
//	movt = true;
//	if (trfr.moveType) {
//
//	}
//}
pInd = new Individual(pCell,pPatch,stg,age,repInt,probmale,trfr.moveModel,trfr.moveType);

// add new individual to the population
// NB THIS WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES...
popns[0]->recruit(pInd);

if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers)
{
	// individual variation - set up genetics
	landData land = pLandscape->getLandData();
	pInd->setGenes(pSpecies,land.resol);
}

}

// Create a new population, and return its address
Population* SubCommunity::newPopn(Landscape *pLandscape,Species *pSpecies,
	Patch *pPatch,int nInds) 
{
#if RSDEBUG
//DEBUGLOG << "SubCommunity::newPopn(): subCommNum = " << subCommNum
//	<< " pPatch = " << pPatch << " nInds = "<< nInds << endl;
#endif
landParams land = pLandscape->getLandParams();
int npopns = (int)popns.size();
popns.push_back(new Population(pSpecies,pPatch,nInds,land.resol));
#if RSDEBUG
//DEBUGLOG << "SubCommunity::newPopn(): subCommNum = " << subCommNum
//	<< " npopns = " << npopns << " popns[npopns] = " << popns[npopns]
//	<< endl;
#endif
return popns[npopns];
}

popStats SubCommunity::getPopStats(void) {
popStats p,pop;               
p.pSpecies = 0; p.spNum = 0; p.nInds = p.nAdults = p.nNonJuvs = 0; p.breeding = false;
p.pPatch = pPatch;
// FOR SINGLE SPECIES IMPLEMENTATION, THERE IS ONLY ONE POPULATION IN THE PATCH
//p = popns[0]->getStats();
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
#if RSDEBUG
//DEBUGLOG << "SubCommunity::getPopStats(): npops = " << npops
//	<< " i = " << i
//	<< " popns[i] = " << popns[i] << endl;
#endif
	pop = popns[i]->getStats();
	p.pSpecies = pop.pSpecies;
	p.spNum = pop.spNum;
	p.nInds += pop.nInds;
	p.nNonJuvs += pop.nNonJuvs;
	p.nAdults += pop.nAdults;
	p.breeding = pop.breeding;
#if RSDEBUG
//DEBUGLOG << "SubCommunity::getPopStats():"
//	<< " p.pSpecies = " << p.pSpecies
//	<< " p.pPatch = " << p.pPatch
//	<< " p.spNum = " << p.spNum
//	<< " p.nInds = " << p.nInds
//	<< endl;
#endif
}
return p;
}

void SubCommunity::resetPopns(void) {
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	delete popns[i];
}
popns.clear();
// clear the list of populations in the corresponding patch
pPatch->resetPopn();
}

void SubCommunity::resetPossSettlers(void) {
if (subCommNum == 0) return; // not applicable in the matrix
//int npops = (int)popns.size();
//for (int i = 0; i < npops; i++) { // all populations
//	popns[i]->resetPossSettlers();
//}
pPatch->resetPossSettlers();
}

// Extirpate all populations according to
// option 0 - random local extinction probability
// option 1 - local extinction probability gradient
// NB only applied for cell-based model
void SubCommunity::localExtinction(int option) {
double pExtinct = 0.0;
if (option == 0) {
	envStochParams env = paramsStoch->getStoch();
	if (env.localExt) pExtinct = env.locExtProb;
}
else {
	envGradParams grad = paramsGrad->getGradient();
	Cell *pCell = pPatch->getRandomCell(); // get only cell in the patch
	// extinction prob is complement of cell gradient value plus any non-zero prob at the optimum
	pExtinct = 1.0 - pCell->getEnvVal() + grad.extProbOpt;
	if (pExtinct > 1.0) pExtinct = 1.0;
}
if (pRandom->Bernoulli(pExtinct)) {
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->extirpate();
	}
}
}

// Action in event of patch becoming unsuitable owing to landscape change
void SubCommunity::patchChange(void) {
if (subCommNum == 0) return; // no reproduction in the matrix
Species *pSpecies;
float localK = 0.0;
int npops = (int)popns.size();
// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
if (npops < 1) return;
localK = pPatch->getK();
if (localK <= 0.0) { // patch in dynamic landscape has become unsuitable
	for (int i = 0; i < npops; i++) { // all populations
		pSpecies = popns[i]->getSpecies();
		demogrParams dem = pSpecies->getDemogr();
		if (dem.stageStruct) {
			stageParams sstruct = pSpecies->getStage();
			if (sstruct.disperseOnLoss) popns[i]->allEmigrate();
			else popns[i]->extirpate();
		}
		else { // non-stage-structured species is destroyed
			popns[i]->extirpate();
		}
	}
}
}

void SubCommunity::reproduction(int resol,float epsGlobal,short rasterType,bool patchModel)
{
if (subCommNum == 0) return; // no reproduction in the matrix
float localK,envval;
//Species *pSpecies;
Cell *pCell;
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();

int npops = (int)popns.size();
// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
if (npops < 1) return;

localK = pPatch->getK();
if (localK > 0.0) {
	if (patchModel) {
		envval = 1.0; // environmental gradient is currently not applied for patch-based model
	}
	else { // cell-based model
		if (grad.gradient && grad.gradType == 2)
		{ // gradient in fecundity
			Cell *pCell = pPatch->getRandomCell(); // locate the only cell in the patch
			envval = pCell->getEnvVal();
		}
		else envval = 1.0;
	}
	if (env.stoch && !env.inK) { // stochasticity in fecundity
		if (env.local) {
			if (!patchModel) { // only permitted for cell-based model
				pCell = pPatch->getRandomCell();
				if (pCell != 0) envval += pCell->getEps();
			}
		}
		else { // global stochasticity
			envval += epsGlobal;
		}
	}
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->reproduction(localK,envval,resol);
		popns[i]->fledge();
	}
}
/*
else { // patch in dynamic landscape has become unsuitable
	// NB - THIS WILL NEED TO BE MADE SPECIES-SPECIFIC...
	Species *pSpecies;
	for (int i = 0; i < npops; i++) { // all populations
		pSpecies = popns[i]->getSpecies();
		demogrParams dem = pSpecies->getDemogr();
		if (dem.stageStruct) {
			stageParams sstruct = pSpecies->getStage();
			if (sstruct.disperseOnLoss) popns[i]->allEmigrate();
			else popns[i]->extirpate();
		}
		else { // non-stage-structured species is destroyed
			popns[i]->extirpate();
		}
	}
}
*/

/*
#if RSDEBUG
//if (npops > 0) {
//	DEBUGLOG << "SubCommunity::reproduction(): this = " << this
//		<< " npops = " << npops << endl;
//}
#endif
#if RSDEBUG
//DEBUGLOG << "SubCommunity::reproduction(): patchNum = " << pPatch->getPatchNum()
//	<< " nCells " << pPatch->getNCells()
//	<< " localK = " << localK
//	<< endl;
#endif
*/
}

void SubCommunity::emigration(void) 
{
if (subCommNum == 0) return; // no emigration from the matrix
float localK;
int npops = (int)popns.size();
// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
if (npops < 1) return;
localK = pPatch->getK();
// NOTE that even if K is zero, it could have been >0 in previous time-step, and there
// might be emigrants if there is non-juvenile emigration
for (int i = 0; i < npops; i++) { // all populations
//	localK = pPatch->getK();
	popns[i]->emigration(localK);
}
}

// Remove emigrants from their natal patch and add to patch 0 (matrix)
void SubCommunity::initiateDispersal(SubCommunity* matrix) {
if (subCommNum == 0) return; // no dispersal initiation in the matrix
#if RSDEBUG
//DEBUGLOG << "SubCommunity::initiateDispersal(): this=" << this
//	<< " subCommNum=" << subCommNum
//	<< endl;
#endif
popStats pop;
disperser disp;

int npops = (int)popns.size();
#if RSDEBUG
//DEBUGLOG << "SubCommunity::initiateDispersal(): patchNum = " << patchNum
//	<< " npops " << npops
//	<< endl;
#endif
for (int i = 0; i < npops; i++) { // all populations
	pop = popns[i]->getStats();
	for (int j = 0; j < pop.nInds; j++) {
#if RSDEBUG
//DEBUGLOG << "SubCommunity::initiateDispersal(): i = " << i
//	<< " j " << j
//	<< endl;
#endif
		disp = popns[i]->extractDisperser(j);
		if (disp.yes) { // disperser - has already been removed from natal population
			// add to matrix population
			matrix->recruit(disp.pInd,pop.pSpecies);
		}
	}
	// remove pointers to emigrants
	popns[i]->clean();
}

}

// Add an individual into the local population of its species in the patch
void SubCommunity::recruit(Individual *pInd,Species *pSpecies) {
#if RSDEBUG
//DEBUGLOG << "SubCommunity::recruit(): this=" << this
//	<< endl;
#endif
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	if (pSpecies == popns[i]->getSpecies()) {
		popns[i]->recruit(pInd);
	}
}
}

// Transfer through the matrix - run for the matrix sub-community only
int SubCommunity::transfer(Landscape *pLandscape,short landIx,short nextseason) 
{
#if RSDEBUG
//DEBUGLOG << "SubCommunity::transfer(): this=" << this
//	<< endl;
#endif
int ndispersers = 0;
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	ndispersers += popns[i]->transfer(pLandscape,landIx,nextseason);
#if RSDEBUG
//DEBUGLOG << "SubCommunity::transfer(): i = " << i
//	<< " this = " << this
//	<< " subCommNum = " << subCommNum
//	<< " npops = " << npops
//	<< " popns[i] = " << popns[i]
//	<< " ndispersers = " << ndispersers
//	<< endl;
#endif
}
return ndispersers;
}

//---------------------------------------------------------------------------

// Remove emigrants from patch 0 (matrix) and transfer to sub-community
// in which their destination co-ordinates fall
// This function is executed for the matrix patch only

void SubCommunity::completeDispersal(Landscape *pLandscape,bool connect) 
{
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): this=" << this
//	<< endl;
#endif
//popStats pop;
//int popsize,subcomm;
int popsize;
disperser settler;
Species *pSpecies;
Population *pPop;
Patch *pPrevPatch;
Patch *pNewPatch;
Cell *pPrevCell;
SubCommunity *pSubComm;

int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	pSpecies = popns[i]->getSpecies();
	popsize = popns[i]->getNInds();
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): i=" << i
//	<< " subCommNum=" << subCommNum
//	<< " pPatch=" << pPatch << " patchNum=" << pPatch->getPatchNum()
//	<< " npops=" << npops
//	<< " popns[i]=" << popns[i]
//	<< " popsize=" << popsize
//	<< endl;
#endif
	for (int j = 0; j < popsize; j++) {
		bool settled;
		settler = popns[i]->extractSettler(j);
		settled = settler.yes;
		if (settled) {
			// settler - has already been removed from matrix population
			// in popns[i]->extractSettler()
#if RSDEBUG
//locn loc = pNewCell->getLocn();
//DEBUGLOG << "SubCommunity::completeDispersal(): j=" << j << " settled=" << settled;
//#if GROUPDISP
//if (emig.groupdisp) {
//DEBUGLOG	<< " groupID=" << groupsettler.pGroup->getGroupID()
//	<< " groupsize=" << groupsettler.groupsize
//	<< " status=" << groupsettler.pGroup->getStatus();
//}
//else {
//DEBUGLOG	<< " settler.pInd=" << settler.pInd;
//}
//#else
//DEBUGLOG	<< " settler.pInd=" << settler.pInd;
//#endif // GROUPDISP
//DEBUGLOG	<< " loc.x=" << loc.x << " loc.y=" << loc.y << endl;
#endif // RSDEBUG
			// find new patch
			pNewPatch = (Patch*)settler.pCell->getPatch();
			// find population within the patch (if there is one)
//			subcomm = ppatch->getSubComm();
//			if (subcomm == 0) { // no sub-community has yet been set up
//				// CANNOT SET UP A NEW ONE FROM WITHIN AN EXISTING ONE!!!!!
//			}
			pPop = (Population*)pNewPatch->getPopn((intptr)pSpecies);
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): j = " << j
//	<< " pCell = " << pCell << " ppatch = " << ppatch << " pPop = " << pPop
//	<< endl;
#endif
			if (pPop == 0) { // settler is the first in a previously uninhabited patch
				// create a new population in the corresponding sub-community
				pSubComm = (SubCommunity*)pNewPatch->getSubComm();
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): j = " << j
//	<< " pSubComm = " << pSubComm << endl;
#endif
				pPop = pSubComm->newPopn(pLandscape,pSpecies,pNewPatch,0);
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): j=" << j
//	<< " pPop=" << pPop << endl;
#endif
			}
			pPop->recruit(settler.pInd);   
			if (connect) { // increment connectivity totals
				int newpatch = pNewPatch->getSeqNum();
				pPrevCell = settler.pInd->getLocn(0); // previous cell
				intptr patch = pPrevCell->getPatch();
				if (patch != 0) {
					pPrevPatch = (Patch*)patch;
					int prevpatch = pPrevPatch->getSeqNum();
					pLandscape->incrConnectMatrix(prevpatch,newpatch);
				}
			}
		}
		else { // for group dispersal only
		}
	}
	// remove pointers in the matrix popn to settlers
	popns[i]->clean();
#if RSDEBUG
//pop = popns[i]->getStats();
popsize = popns[i]->getNInds();
//DEBUGLOG << "SubCommunity::completeDispersal(): i=" << i
//	<< " popns[i]=" << popns[i]
////	<< " pop.pPatch = " << pop.pPatch
////	<< " pop.nInds = " << pop.nInds
//	<< " popsize=" << popsize
//	<< endl;
#endif
}

}

//---------------------------------------------------------------------------

void SubCommunity::survival(short part,short option0,short option1)
{
int npops = (int)popns.size();
if (npops < 1) return;
if (part == 0) {
	float localK = pPatch->getK();
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->survival0(localK,option0,option1);
	}
}
else {
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->survival1();
	}
}
}

void SubCommunity::ageIncrement(void) {
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popns[i]->ageIncrement();
}
}

// Find the population of a given species in a given patch
Population* SubCommunity::findPop(Species *pSp,Patch *pPch) {
#if RSDEBUG
DEBUGLOG << "SubCommunity::findPop(): this=" << this
	<< endl;
#endif

Population *pPop = 0;
popStats pop;
int npops = (int)popns.size();

for (int i = 0; i < npops; i++) { // all populations
	pop = popns[i]->getStats();
	if (pop.pSpecies == pSp && pop.pPatch == pPch) { // population located
		pPop = popns[i];
		break;
	}
	else pPop = 0;
}
return pPop;
}

//---------------------------------------------------------------------------

void SubCommunity::createOccupancy(int nrows) {
if (occupancy != 0) deleteOccupancy();
if (nrows > 0) {
	occupancy = new int[nrows];
	for (int i = 0; i < nrows; i++) occupancy[i] = 0;
}
}

void SubCommunity::updateOccupancy(int row) {
#if RSDEBUG
//DEBUGLOG << "SubCommunity::updateOccupancy(): this=" << this
//	<< endl;
#endif
popStats pop;
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) {
	pop = popns[i]->getStats();
	if (pop.nInds > 0 && pop.breeding) {
		occupancy[row]++;
		i = npops;
	}
}
}

int SubCommunity::getOccupancy(int row) {
if (row >= 0) return occupancy[row];
else return 0;
}

void SubCommunity::deleteOccupancy(void) {
	delete[] occupancy;
	occupancy = 0;
}

//---------------------------------------------------------------------------
// Open population file and write header record
bool SubCommunity::outPopHeaders(Landscape *pLandscape,Species *pSpecies,int option)
{
bool fileOK;
Population *pPop;
landParams land = pLandscape->getLandParams();

if (option == -999) { // close the file
	// as all populations may have been deleted, set up a dummy one
	// species is not necessary
	pPop = new Population();
	fileOK = pPop->outPopHeaders(-999,land.patchModel);
	delete pPop;
}
else { // open the file
	// as no population has yet been created, set up a dummy one
	// species is necessary, as columns depend on stage and sex structure
	pPop = new Population(pSpecies,pPatch,0,land.resol);
	fileOK = pPop->outPopHeaders(land.landNum,land.patchModel);
	delete pPop;
}
return fileOK;
}

// Write records to population file
void SubCommunity::outPop(Landscape *pLandscape,int rep,int yr,int gen)
{
landParams land = pLandscape->getLandParams();
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
bool writeEnv = false;
bool gradK = false;
if (grad.gradient) {
	writeEnv = true;
	if (grad.gradType == 1) gradK = true; // ... carrying capacity
}
if (env.stoch) writeEnv = true;

// generate output for each population within the sub-community (patch)
// provided that the patch is suitable (i.e. non-zero carrying capacity)
// or the population is above zero (possible if there is stochasticity or a moving gradient)
// or it is the matrix patch in a patch-based model
int npops = (int)popns.size();
int patchnum;
//Species* pSpecies;
Cell *pCell;
float localK;
float eps = 0.0;
if (env.stoch) {
	if (env.local) {
		pCell = pPatch->getRandomCell();
		if (pCell != 0) eps = pCell->getEps();
	}
	else {
		eps = pLandscape->getGlobalStoch(yr);
	}
}

patchnum = pPatch->getPatchNum();
for (int i = 0; i < npops; i++) { // all populations
//	pSpecies = popns[i]->getSpecies();
	localK = pPatch->getK();
	if (localK > 0.0 || (land.patchModel && patchnum == 0)) {
		popns[i]->outPopulation(rep,yr,gen,eps,land.patchModel,writeEnv,gradK);
	}
	else {
		if (popns[i]->totalPop() > 0) {
			popns[i]->outPopulation(rep,yr,gen,eps,land.patchModel,writeEnv,gradK);
		}
	}
}
}

// Write records to individuals file
void SubCommunity::outInds(Landscape *pLandscape,int rep,int yr,int gen,int landNr) {
landParams ppLand = pLandscape->getLandParams();
if (landNr >= 0) { // open the file
	popns[0]->outIndsHeaders(rep,landNr,ppLand.patchModel);
	return;
}
if (landNr == -999) { // close the file
	popns[0]->outIndsHeaders(rep,-999,ppLand.patchModel);
	return;
}
// generate output for each population within the sub-community (patch)
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popns[i]->outIndividual(pLandscape,rep,yr,gen,pPatch->getPatchNum());
}
}

// Write records to individuals file
void SubCommunity::outGenetics(int rep,int yr,int gen,int landNr)
{
//landParams ppLand = pLandscape->getLandParams();
if (landNr >= 0) { // open the file
	popns[0]->outGenetics(rep,yr,landNr);
	return;
}
if (landNr == -999) { // close the file
	popns[0]->outGenetics(rep,yr,landNr);
	return;
}
// generate output for each population within the sub-community (patch)
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popns[i]->outGenetics(rep,yr,landNr);
}
}

// Population size of a specified stage
int SubCommunity::stagePop(int stage) {
int popsize = 0;
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popsize += popns[i]->stagePop(stage);
}
return popsize;
}

// Open traits file and write header record
bool SubCommunity::outTraitsHeaders(Landscape *pLandscape,Species *pSpecies,int landNr)
{
//Population *pPop;
landParams land = pLandscape->getLandParams();
if (landNr == -999) { // close file
	if (outtraits.is_open()) outtraits.close();
	outtraits.clear();
	return true;
}

string name;
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
simParams sim = paramsSim->getSim();

string DirOut = paramsSim->getDir(2);
if (sim.batchMode) {
	if (land.patchModel){
		name = DirOut
			+ "Batch" + Int2Str(sim.batchNum) + "_"
			+ "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(landNr)
			+ "_TraitsXpatch.txt";
	}
	else{
		name = DirOut
			+ "Batch" + Int2Str(sim.batchNum) + "_"
			+ "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(landNr)
			+ "_TraitsXcell.txt";
	}
}
else {
	if (land.patchModel){
		name = DirOut + "Sim" + Int2Str(sim.simulation) + "_TraitsXpatch.txt";
	}
	else{
		name = DirOut + "Sim" + Int2Str(sim.simulation) + "_TraitsXcell.txt";
	}
}
outtraits.open(name.c_str());

outtraits << "Rep\tYear\tRepSeason";
if (land.patchModel) outtraits << "\tPatchID";
else
	outtraits << "\tx\ty";

if (emig.indVar) {
	if (emig.sexDep) {
		if (emig.densDep) {
			outtraits << "\tF_meanD0\tF_stdD0\tM_meanD0\tM_stdD0";
			outtraits << "\tF_meanAlpha\tF_stdAlpha\tM_meanAlpha\tM_stdAlpha";
			outtraits << "\tF_meanBeta\tF_stdBeta\tM_meanBeta\tM_stdBeta";
		}
		else {
			outtraits << "\tF_meanEP\tF_stdEP\tM_meanEP\tM_stdEP";
		}
	}
	else {
		if (emig.densDep) {
			outtraits << "\tmeanD0\tstdD0\tmeanAlpha\tstdAlpha";
			outtraits << "\tmeanBeta\tstdBeta";
		}
		else {
			outtraits << "\tmeanEP\tstdEP";
		}
	}
}
if (trfr.indVar) {
	if (trfr.moveModel) {
		if (trfr.moveType == 1) {
			outtraits << "\tmeanDP\tstdDP\tmeanGB\tstdGB";
			outtraits << "\tmeanAlphaDB\tstdAlphaDB\tmeanBetaDB\tstdBetaDB";
		}
		if (trfr.moveType == 2) {
			outtraits << "\tmeanStepLength\tstdStepLength\tmeanRho\tstdRho";
		}
	}
	else {
		if (trfr.sexDep) {
			outtraits << "\tF_mean_distI\tF_std_distI\tM_mean_distI\tM_std_distI";
			if (trfr.twinKern)
				outtraits << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
		}
		else {
			outtraits << "\tmean_distI\tstd_distI";
			if (trfr.twinKern)
				outtraits << "\tmean_distII\tstd_distII\tmeanPfirstKernel\tstdPfirstKernel";
		}
	}
}
if (sett.indVar) {
	if (sett.sexDep) {
		outtraits << "\tF_meanS0\tF_stdS0\tM_meanS0\tM_stdS0";
		outtraits << "\tF_meanAlphaS\tF_stdAlphaS\tM_meanAlphaS\tM_stdAlphaS";
		outtraits << "\tF_meanBetaS\tF_stdBetaS\tM_meanBetaS\tM_stdBetaS";
	}
	else {
		outtraits << "\tmeanS0\tstdS0";
		outtraits << "\tmeanAlphaS\tstdAlphaS";
		outtraits << "\tmeanBetaS\tstdBetaS";
	}
}
outtraits << endl;

return outtraits.is_open();
}

// Write records to traits file and return aggregated sums
traitsums SubCommunity::outTraits(traitCanvas tcanv,
	Landscape *pLandscape,int rep,int yr,int gen,bool commlevel)
{
int popsize,ngenes;
landParams land = pLandscape->getLandParams();
simParams sim = paramsSim->getSim();
bool writefile = false;
if (sim.outTraitsCells && yr%sim.outIntTraitCell == 0 && !commlevel)
	writefile = true;
traitsums ts,poptraits;
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

// generate output for each population within the sub-community (patch)
// provided that the patch is suitable (i.e. non-zero carrying capacity)
int npops = (int)popns.size();
Species* pSpecies;
float localK;

for (int i = 0; i < npops; i++) { // all populations
	localK = pPatch->getK();
	if (localK > 0.0 && popns[i]->getNInds() > 0) {
		pSpecies = popns[i]->getSpecies();
		demogrParams dem = pSpecies->getDemogr();
		emigRules emig = pSpecies->getEmig();
		trfrRules trfr = pSpecies->getTrfr();
		settleType sett = pSpecies->getSettle();
		poptraits = popns[i]->getTraits(pSpecies);

		if (writefile) {
			outtraits << rep << "\t" << yr << "\t" << gen;
			if (land.patchModel) {
				outtraits << "\t" << pPatch->getPatchNum();
			}
			else {
				locn loc = pPatch->getCellLocn(0);
				outtraits << "\t" << loc.x << "\t" << loc.y;
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
				if (ngenes == 2) popsize = poptraits.ninds[g];
				else popsize = poptraits.ninds[0] + poptraits.ninds[1];
				if (popsize > 0) {
					mnD0[g] = poptraits.sumD0[g] / (double)popsize;
					mnAlpha[g] = poptraits.sumAlpha[g] / (double)popsize;
					mnBeta[g] = poptraits.sumBeta[g] / (double)popsize;
					if (popsize > 1) {
						sdD0[g] = poptraits.ssqD0[g]/(double)popsize	- mnD0[g]*mnD0[g];
						if (sdD0[g] > 0.0) sdD0[g] = sqrt(sdD0[g]); else sdD0[g] = 0.0;
						sdAlpha[g] = poptraits.ssqAlpha[g]/(double)popsize	- mnAlpha[g]*mnAlpha[g];
						if (sdAlpha[g] > 0.0) sdAlpha[g] = sqrt(sdAlpha[g]); else sdAlpha[g] = 0.0;
						sdBeta[g] = poptraits.ssqBeta[g]/(double)popsize	- mnBeta[g]*mnBeta[g];
						if (sdBeta[g] > 0.0) sdBeta[g] = sqrt(sdBeta[g]); else sdBeta[g] = 0.0;
					}
					else {
						sdD0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
					}
				}
			}
			if (writefile) {
				if (emig.sexDep) {
					outtraits << "\t" << mnD0[0] << "\t" << sdD0[0];
					outtraits << "\t" << mnD0[1] << "\t" << sdD0[1];
					if (emig.densDep) {
						outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
						outtraits << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
						outtraits << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
						outtraits << "\t" << mnBeta[1]  << "\t" << sdBeta[1];
					}
				}
				else { // sex-independent
					outtraits << "\t" << mnD0[0] << "\t" << sdD0[0];
					if (emig.densDep) {
						outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
						outtraits << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
					}
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
			double mnDist1[2], mnDist2[2], mnProp1[2], mnStepL[2], mnRho[2];
			double sdDist1[2],sdDist2[2],sdProp1[2],sdStepL[2],sdRho[2];
			double mnDP[2], mnGB[2], mnAlphaDB[2], mnBetaDB[2];
			double sdDP[2],sdGB[2],sdAlphaDB[2],sdBetaDB[2];
			for (int g = 0; g < ngenes; g++) {
				mnDist1[g] = mnDist2[g] = mnProp1[g] = mnStepL[g] = mnRho[g] = 0.0;
				sdDist1[g] = sdDist2[g] = sdProp1[g] = sdStepL[g] = sdRho[g] = 0.0;
				mnDP[g] = mnGB[g] = mnAlphaDB[g] = mnBetaDB[g] = 0.0;
				sdDP[g] = sdGB[g] = sdAlphaDB[g] = sdBetaDB[g] = 0.0;
				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ngenes == 2) popsize = poptraits.ninds[g];
				else popsize = poptraits.ninds[0] + poptraits.ninds[1];
				if (popsize > 0) {
					mnDist1[g] = poptraits.sumDist1[g] / (double)popsize;
					mnDist2[g] = poptraits.sumDist2[g] / (double)popsize;
					mnProp1[g] = poptraits.sumProp1[g] / (double)popsize;
					mnStepL[g] = poptraits.sumStepL[g] / (double)popsize;
					mnRho[g] =   poptraits.sumRho[g]   / (double)popsize;
					mnDP[g] = poptraits.sumDP[g] / (double)popsize;
					mnGB[g] = poptraits.sumGB[g] / (double)popsize;
					mnAlphaDB[g] = poptraits.sumAlphaDB[g] / (double)popsize;
					mnBetaDB[g]  = poptraits.sumBetaDB[g]  / (double)popsize;
					if (popsize > 1) {
						sdDist1[g] = poptraits.ssqDist1[g]/(double)popsize	- mnDist1[g]*mnDist1[g];
						if (sdDist1[g] > 0.0) sdDist1[g] = sqrt(sdDist1[g]); else sdDist1[g] = 0.0;
						sdDist2[g] = poptraits.ssqDist2[g]/(double)popsize	- mnDist2[g]*mnDist2[g];
						if (sdDist2[g] > 0.0) sdDist2[g] = sqrt(sdDist2[g]); else sdDist2[g] = 0.0;
						sdProp1[g] = poptraits.ssqProp1[g]/(double)popsize	- mnProp1[g]*mnProp1[g];
						if (sdProp1[g] > 0.0) sdProp1[g] = sqrt(sdProp1[g]); else sdProp1[g] = 0.0;
						sdStepL[g] = poptraits.ssqStepL[g]/(double)popsize	- mnStepL[g]*mnStepL[g];
						if (sdStepL[g] > 0.0) sdStepL[g] = sqrt(sdStepL[g]); else sdStepL[g] = 0.0;
						sdRho[g] = poptraits.ssqRho[g]/(double)popsize	- mnRho[g]*mnRho[g];
						if (sdRho[g] > 0.0) sdRho[g] = sqrt(sdRho[g]); else sdRho[g] = 0.0;
						sdDP[g] = poptraits.ssqDP[g]/(double)popsize	- mnDP[g]*mnDP[g];
						if (sdDP[g] > 0.0) sdDP[g] = sqrt(sdDP[g]); else sdDP[g] = 0.0;
						sdGB[g] = poptraits.ssqGB[g]/(double)popsize	- mnGB[g]*mnGB[g];
						if (sdGB[g] > 0.0) sdGB[g] = sqrt(sdGB[g]); else sdGB[g] = 0.0;
						sdAlphaDB[g] = poptraits.ssqAlphaDB[g]/(double)popsize	- mnAlphaDB[g]*mnAlphaDB[g];
						if (sdAlphaDB[g] > 0.0) sdAlphaDB[g] = sqrt(sdAlphaDB[g]); else sdAlphaDB[g] = 0.0;
						sdBetaDB[g]  = poptraits.ssqBetaDB[g]/(double)popsize	- mnBetaDB[g]*mnBetaDB[g];
						if (sdBetaDB[g] > 0.0) sdBetaDB[g] = sqrt(sdBetaDB[g]); else sdBetaDB[g] = 0.0;
					}
				}
			}
			if (writefile) {
				if (trfr.moveModel) {
					if (trfr.moveType == 1) {
						outtraits << "\t" << mnDP[0] << "\t" << sdDP[0];
						outtraits << "\t" << mnGB[0] << "\t" << sdGB[0];
						outtraits << "\t" << mnAlphaDB[0] << "\t" << sdAlphaDB[0];
						outtraits << "\t" << mnBetaDB[0] << "\t" << sdBetaDB[0];
					}
					if (trfr.moveType == 2) {
						outtraits << "\t" << mnStepL[0] << "\t" << sdStepL[0];
						outtraits << "\t" << mnRho[0] << "\t" << sdRho[0];
					}
				}
				else {
					if (trfr.sexDep) {
						outtraits << "\t" << mnDist1[0] << "\t" << sdDist1[0];
						outtraits << "\t" << mnDist1[1] << "\t" << sdDist1[1];
						if (trfr.twinKern) 
						{
							outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
							outtraits << "\t" << mnDist2[1] << "\t" << sdDist2[1];
							outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
							outtraits << "\t" << mnProp1[1] << "\t" << sdProp1[1];
						}
					}
					else { // sex-independent
						outtraits << "\t" << mnDist1[0] << "\t" << sdDist1[0];
						if (trfr.twinKern) 
						{
							outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
							outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
						}
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
//			ngenes = 1;
			double mnS0[2],mnAlpha[2],mnBeta[2],sdS0[2],sdAlpha[2],sdBeta[2];
			for (int g = 0; g < ngenes; g++) {
				mnS0[g] = mnAlpha[g] = mnBeta[g] = sdS0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ngenes == 2) popsize = poptraits.ninds[g];
				else popsize = poptraits.ninds[0] + poptraits.ninds[1];
				if (popsize > 0) {
					mnS0[g] = poptraits.sumS0[g] / (double)popsize;
					mnAlpha[g] = poptraits.sumAlphaS[g] / (double)popsize;
					mnBeta[g] = poptraits.sumBetaS[g] / (double)popsize;
					if (popsize > 1) {
						sdS0[g] = poptraits.ssqS0[g]/(double)popsize	- mnS0[g]*mnS0[g];
						if (sdS0[g] > 0.0) sdS0[g] = sqrt(sdS0[g]); else sdS0[g] = 0.0;
						sdAlpha[g] = poptraits.ssqAlphaS[g]/(double)popsize	- mnAlpha[g]*mnAlpha[g];
						if (sdAlpha[g] > 0.0) sdAlpha[g] = sqrt(sdAlpha[g]); else sdAlpha[g] = 0.0;
						sdBeta[g] = poptraits.ssqBetaS[g]/(double)popsize	- mnBeta[g]*mnBeta[g];
						if (sdBeta[g] > 0.0) sdBeta[g] = sqrt(sdBeta[g]); else sdBeta[g] = 0.0;
					}
					else {
						sdS0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
					}
				}
			}
			if (writefile) {
				if (sett.sexDep) {
					outtraits << "\t" << mnS0[0] << "\t" << sdS0[0];
					outtraits << "\t" << mnS0[1] << "\t" << sdS0[1];
					outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
					outtraits << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
					outtraits << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
					outtraits << "\t" << mnBeta[1]  << "\t" << sdBeta[1];
				}
				else { // sex-independent
					outtraits << "\t" << mnS0[0] << "\t" << sdS0[0];
					outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
					outtraits << "\t" << mnBeta[0]  << "\t" << sdBeta[0];
				}
			}
		}

		if (writefile) outtraits << endl;

		for (int s = 0; s < NSEXES; s++) {
			ts.ninds[s] 	  += poptraits.ninds[s];
			ts.sumD0[s]     += poptraits.sumD0[s];     ts.ssqD0[s]     += poptraits.ssqD0[s];
			ts.sumAlpha[s]  += poptraits.sumAlpha[s];  ts.ssqAlpha[s]  += poptraits.ssqAlpha[s];
			ts.sumBeta[s]   += poptraits.sumBeta[s];   ts.ssqBeta[s]   += poptraits.ssqBeta[s];
			ts.sumDist1[s]  += poptraits.sumDist1[s];  ts.ssqDist1[s]  += poptraits.ssqDist1[s];
			ts.sumDist2[s]  += poptraits.sumDist2[s];  ts.ssqDist2[s]  += poptraits.ssqDist2[s];
			ts.sumProp1[s]  += poptraits.sumProp1[s];  ts.ssqProp1[s]  += poptraits.ssqProp1[s];
			ts.sumDP[s]     += poptraits.sumDP[s];     ts.ssqDP[s]     += poptraits.ssqDP[s];
			ts.sumGB[s]     += poptraits.sumGB[s];     ts.ssqGB[s]     += poptraits.ssqGB[s];
			ts.sumAlphaDB[s] += poptraits.sumAlphaDB[s]; ts.ssqAlphaDB[s] += poptraits.ssqAlphaDB[s];
			ts.sumBetaDB[s]  += poptraits.sumBetaDB[s];  ts.ssqBetaDB[s]  += poptraits.ssqBetaDB[s];
			ts.sumStepL[s]  += poptraits.sumStepL[s];  ts.ssqStepL[s]  += poptraits.ssqStepL[s];
			ts.sumRho[s]    += poptraits.sumRho[s];    ts.ssqRho[s]    += poptraits.ssqRho[s];
			ts.sumS0[s]     += poptraits.sumS0[s];     ts.ssqS0[s]     += poptraits.ssqS0[s];
			ts.sumAlphaS[s] += poptraits.sumAlphaS[s]; ts.ssqAlphaS[s] += poptraits.ssqAlphaS[s];
			ts.sumBetaS[s]  += poptraits.sumBetaS[s];  ts.ssqBetaS[s]  += poptraits.ssqBetaS[s];
#if RSDEBUG
//DEBUGLOG << "SubCommunity::outTraits(): i=" << i << " popns[i]=" << popns[i]
//	<< " s=" << s
////	<< " poptraits.sumRho[s]= " << poptraits.sumRho[s]
////	<< " ts.sumRho[s]= " << ts.sumRho[s]
//	<< " poptraits.sumDP[s]= " << poptraits.sumDP[s] << " poptraits.ssqDP[s]= " << poptraits.ssqDP[s]
//	<< " ts.sumDP[s]= " << ts.sumDP[s] << " ts.ssqDP[s]= " << ts.ssqDP[s]
//	<< " poptraits.sumGB[s]= " << poptraits.sumGB[s] << " poptraits.ssqGB[s]= " << poptraits.ssqGB[s]
//	<< " ts.sumGB[s]= " << ts.sumGB[s] << " ts.ssqGB[s]= " << ts.ssqGB[s]
//	<< endl;
#endif
		}
	}
}
return ts;
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


