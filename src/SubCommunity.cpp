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
#if RS_EMBARCADERO
#pragma hdrstop
#endif

#include "SubCommunity.h"
//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma package(smart_init) 
#endif

ofstream outtraits;

//---------------------------------------------------------------------------

SubCommunity::SubCommunity(Patch* pPch, int num) {
subCommNum = num;
pPatch = pPch;
// record the new sub-community no. in the patch
pPatch->setSubComm((intptr)this);
initial = false;
occupancy = 0;
#if RS_CONTAIN
habIndex = -1;
#endif // RS_CONTAIN 
#if RS_CONTAIN
cullTarget = false;
firstYear = -1;
cullCount = 0;
#endif // RS_CONTAIN 
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

#if RS_CONTAIN
void SubCommunity::setHabIndex(Species *pSpecies,short rastertype,short landIx) {
demogrParams dem = pSpecies->getDemogr();
if (dem.habDepDem) {   
	if (rastertype == 0) {
		habIndex = -1;
		Cell *pCell = pPatch->getRandomCell();
		if (pCell != 0) {
			habIndex = pCell->getHabIndex(landIx);
		}
	}
	else habIndex = 0;	
}
else habIndex = 0;
}
#endif // RS_CONTAIN 

void SubCommunity::setInitial(bool b) { initial = b; }

#if PEDIGREE
void SubCommunity::initialise(Landscape *pLandscape,Species *pSpecies,Pedigree *pPed) 
#else
void SubCommunity::initialise(Landscape *pLandscape,Species *pSpecies) 
#endif
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
#if SEASONAL
	float k = pPatch->getK(0);
#else
	float k = pPatch->getK();
#endif // SEASONAL 
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
#if PEDIGREE
	newPopn(pLandscape,pSpecies,pPed,pPatch,nInds);
#else
	newPopn(pLandscape,pSpecies,pPatch,nInds);
#endif
}

}

// initialise a specified individual
#if PEDIGREE
void SubCommunity::initialInd(Landscape *pLandscape,Species *pSpecies,
	Pedigree *pPed,Patch *pPatch,Cell *pCell,int ix) 
#else
void SubCommunity::initialInd(Landscape *pLandscape,Species *pSpecies,
	Patch *pPatch,Cell *pCell,int ix) 
#endif
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
#if PEDIGREE
	newPopn(pLandscape,pSpecies,pPed,pPatch,0);
#else
	newPopn(pLandscape,pSpecies,pPatch,0);
#endif
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
#if RS_CONTAIN
pInd = new Individual(pCell,pPatch,stg,age,repInt,1,probmale,trfr.moveModel,trfr.moveType);
#else
#if PARTMIGRN
pInd = new Individual(pSpecies,pCell,pPatch,stg,age,repInt,probmale,trfr.moveModel,trfr.moveType);
#else
pInd = new Individual(pCell,pPatch,stg,age,repInt,probmale,trfr.moveModel,trfr.moveType);
#endif // PARTMIGRN 
#endif // RS_CONTAIN 

// add new individual to the population
// NB THIS WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES...
popns[0]->recruit(pInd);

#if GOBYMODEL
if (true)
#else
if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers)
#endif
{
	// individual variation - set up genetics
	landData land = pLandscape->getLandData();
	pInd->setGenes(pSpecies,land.resol);
}

}

// Create a new population, and return its address
#if PEDIGREE
Population* SubCommunity::newPopn(Landscape *pLandscape,Species *pSpecies,
	Pedigree *pPed,Patch *pPatch,int nInds) 
#else
Population* SubCommunity::newPopn(Landscape *pLandscape,Species *pSpecies,
	Patch *pPatch,int nInds) 
#endif
{
#if RSDEBUG
//DEBUGLOG << "SubCommunity::newPopn(): subCommNum = " << subCommNum
//	<< " pPatch = " << pPatch << " nInds = "<< nInds << endl;
#endif
landParams land = pLandscape->getLandParams();
int npopns = (int)popns.size();
#if PEDIGREE
popns.push_back(new Population(pSpecies,pPed,pPatch,nInds,land.resol));
#else
popns.push_back(new Population(pSpecies,pPatch,nInds,land.resol));
#endif
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
#if GOBYMODEL
p.nSocial = p.nAsocial = 0;
#endif
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
#if RS_CONTAIN
	pop = popns[i]->getStats(habIndex);
#else
#if SPATIALDEMOG
	pop = popns[i]->getStats(pPatch->getDemoScaling());
#else
	pop = popns[i]->getStats();
#endif // SPATIALDEMOG
#endif // RS_CONTAIN 
	p.pSpecies = pop.pSpecies;
	p.spNum = pop.spNum;
	p.nInds += pop.nInds;
	p.nNonJuvs += pop.nNonJuvs;
	p.nAdults += pop.nAdults;
	p.breeding = pop.breeding;
#if GOBYMODEL
	p.nSocial += pop.nSocial;
	p.nAsocial += pop.nAsocial;
#endif
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
#if SEASONAL
demogrParams dem = pSpecies->getDemogr();
for (int i = 0; i < dem.nSeasons; i++) {
	if (pPatch->getK(i) > localK) localK = pPatch->getK(i);
}
#else
localK = pPatch->getK();
#endif // SEASONAL 
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

#if GROUPDISP
Individual* SubCommunity::getFather(int minbrdstage,int ix) {
int npops = (int)popns.size();
if (npops < 1) return 0;
return popns[0]->getFather(minbrdstage,ix);
}
#endif

#if SEASONAL
void SubCommunity::reproduction(int resol,float epsGlobal,short season,short rasterType,bool patchModel)
#else
#if GROUPDISP
void SubCommunity::reproduction(Landscape* pLandscape,Species* pSpecies,int minbrdstage,
	const std::vector <Individual*> *pfglobal,const int nfglobal,
	const int resol,const locn min,const locn max,
	const float epsGlobal,const short rasterType,const bool patchModel)
#else
#if BUTTERFLYDISP
void SubCommunity::reproduction(int resol,float epsGlobal,short dispersal,short option,
	short rasterType,bool patchModel)
#else
void SubCommunity::reproduction(int resol,float epsGlobal,short rasterType,bool patchModel)
#endif // BUTTERFLYDISP
#endif // GROUPDISP
#endif // SEASONAL
{
if (subCommNum == 0) return; // no reproduction in the matrix
float localK,envval;
#if SPATIALDEMOG
std::vector <float> localDemoScaling;
#endif //SPATIALDEMOG
//Species *pSpecies;
Cell *pCell;
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
#if BUTTERFLYDISP
//demogrParams dem = pSpecies->getDemogr();
#endif
#if GROUPDISP
demogrParams dem = pSpecies->getDemogr();
std::vector <Individual*> fnbrhd;
std::vector <Individual*> *pfnbrhd;
pfnbrhd = &fnbrhd;
int nfnbrhd = 0;

if (dem.repType == 3 && popns.size() > 0) { // hermaphrodite species (and present)
	if (dem.paternity == 2) { // pollen kernel
		// populate list of potential 'fathers' in neighbourhood cells
		popStats pop = popns[0]->getStats();
		Individual* pInd;
		if (pop.breeding) {
			locn loc = getLocn();
			int nbrcells = 0;
			for (int xx = -1; xx < 2; xx++) {
				for (int yy = -1; yy < 2; yy++) {
					if (xx != 0 || yy != 0) { // not current cell
						if (loc.x+xx >= min.x && loc.x+xx <= max.x
						&& 	loc.y+yy >= min.y && loc.y+yy <= max.y ) {
							nbrcells++;
							pCell = pLandscape->findCell(loc.x+xx,loc.y+yy);
							if (pCell != 0) {
								intptr pnbrpatch = pCell->getPatch();
								if (pnbrpatch != 0) {
									Patch* pNbrPatch = (Patch*)pnbrpatch;
									intptr pnbrpopn = pNbrPatch->getPopn((intptr)pSpecies);
									if (pnbrpopn != 0) {
										Population* pNbrPopn = (Population*)pnbrpopn;
										popStats nbrpop = pNbrPopn->getStats();
										if (nbrpop.breeding) {
											for (int j = 0; j < nbrpop.nInds; j++) {
												pInd = pNbrPopn->getFather(minbrdstage,j);
												if (pInd != 0) fnbrhd.push_back(pInd);
											}
										}
									}
								}
							}
						}
					} // end of not current cell
				} // end of for yy
			} // end of for xx
			nfnbrhd = (int)fnbrhd.size();
#if RSDEBUG
//DEBUGLOG << "SubCommunity::reproduction(): this=" << this
//	<< " loc.x=" << loc.x << " loc.y=" << loc.y << " nbrcells=" << nbrcells
//	<< " nfnbrhd=" << nfnbrhd
////	<< " min.x=" << min.x << " min.y=" << min.y
////	<< " max.x=" << max.x << " max.y=" << max.y
//	<< endl;
#endif
		} // end of population is breeding
	}
}

#if RSDEBUG
//DEBUGLOG << "SubCommunity::reproduction(): this=" << this
//	<< " pfnbrhd=" << pfnbrhd << " nfnbrhd=" << nfnbrhd
//	<< endl;
//DEBUGLOG << "SubCommunity::reproduction(): this=" << this
//	<< " pfglobal=" << pfglobal << " nglobal=" << nglobal
//	<< " pfglobal[0]=" << (*pfglobal)[0] << " pfglobal[1]=" << (*pfglobal)[1]
//	<< " pfglobal[nfathers-1]=" << (*pfglobal)[nfathers-1]
//	<< endl;
#endif
#endif

int npops = (int)popns.size();
// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
if (npops < 1) return;

#if SEASONAL
localK = pPatch->getK(season);
#else
localK = pPatch->getK();
#endif // SEASONAL 
#if SPATIALDEMOG
localDemoScaling = pPatch->getDemoScaling();
#endif //SPATIALDEMOG
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
#if RS_CONTAIN
#if SEASONAL
		popns[i]->reproduction(habIndex,season,localK,envval,resol);
#else
		popns[i]->reproduction(habIndex,localK,envval,resol);
#endif // SEASONAL
#else
#if SEASONAL
		popns[i]->reproduction(season,localK,envval,resol);
#else
#if GROUPDISP
		popns[i]->reproduction(pfglobal,nfglobal,pfnbrhd,nfnbrhd,localK,envval,resol);
#else
#if BUTTERFLYDISP
		popns[i]->reproduction(localK,envval,resol,option);
#else
#if SPATIALDEMOG
		popns[i]->reproduction(localK,envval,resol,localDemoScaling);
#else
		popns[i]->reproduction(localK,envval,resol);
#endif // SPATIALDEMOG
#endif // BUTTERFLYDISP 
#endif // GROUPDISP 
#endif // SEASONAL
#endif // RS_CONTAIN 
#if BUTTERFLYDISP
		if (option == 0) // complete classical reproduction
			popns[i]->fledge();
#else
		popns[i]->fledge();
#endif
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

#if BUTTERFLYDISP
// Complete reproduction when dispersal has occurred during reproduction
void SubCommunity::fledge(void) {
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popns[i]->fledge();
}
}
#endif

#if RS_DISEASE
void SubCommunity::emigration(Species *pSpecies,short season) 
#else
#if SEASONAL
void SubCommunity::emigration(short season) 
#else
void SubCommunity::emigration(void) 
#endif // SEASONAL 
#endif // RS_DISEASE  
{
if (subCommNum == 0) return; // no emigration from the matrix
float localK;
int npops = (int)popns.size();
// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
if (npops < 1) return;
#if RS_DISEASE
// use local K for the NEXT season in emigration decision
demogrParams dem = pSpecies->getDemogr();
if (season+1 < dem.nSeasons) localK = pPatch->getK(season+1);
else localK = pPatch->getK(0);
#else
#if SEASONAL
localK = pPatch->getK(season);
#else
localK = pPatch->getK();
#endif // SEASONAL 
#endif // RS_DISEASE  
// NOTE that even if K is zero, it could have been >0 in previous time-step, and there
// might be emigrants if there is non-juvenile emigration
for (int i = 0; i < npops; i++) { // all populations
//	localK = pPatch->getK();
#if SEASONAL
	popns[i]->emigration(localK,season);
#else
	popns[i]->emigration(localK);
#endif
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
#if RS_CONTAIN
	pop = popns[i]->getStats(habIndex);
#else
#if SPATIALDEMOG
	pop = popns[i]->getStats(pPatch->getDemoScaling());
#else
	pop = popns[i]->getStats();
#endif // SPATIALDEMOG
#endif // RS_CONTAIN 
#if GROUPDISP
	bool newgroup = true;
	int groupsize,currentsize;
	emigRules emig = pop.pSpecies->getEmig();
#endif // GROUPDISP
	for (int j = 0; j < pop.nInds; j++) {
#if RSDEBUG
//DEBUGLOG << "SubCommunity::initiateDispersal(): i = " << i
//	<< " j " << j
//	<< endl;
#endif
		disp = popns[i]->extractDisperser(j);
		if (disp.yes) { // disperser - has already been removed from natal population
#if GROUPDISP
			if (emig.groupdisp) {
				// add to matrix group
				if (newgroup) {
					groupsize = 1 + pRandom->Poisson(emig.groupmean-1);
					currentsize = 0;
					matrix->recruit(pPatch,disp.pInd,pop.pSpecies,true);
					currentsize++;
					if (currentsize < groupsize) newgroup = false;
				}
				else {
					matrix->recruit(pPatch,disp.pInd,pop.pSpecies,false);
					currentsize++;
					if (currentsize >= groupsize) newgroup = true;
				}
			}
			else {
				// add to matrix population
				matrix->recruit(disp.pInd,pop.pSpecies);
			}
#else
			// add to matrix population
			matrix->recruit(disp.pInd,pop.pSpecies);
#endif // GROUPDISP
		}
	}
	// remove pointers to emigrants
	popns[i]->clean();
}

}

#if GROUPDISP

// Add an individual into the local population of its species in the patch
void SubCommunity::recruit(Patch *pPch,Individual *pInd,Species *pSpecies,
	bool newgroup) {
#if RSDEBUG
//DEBUGLOG << "SubCommunity::recruit(): GROUPDISP this=" << this
//	<< " pPch=" << pPch << " pInd=" << pInd << " indID=" << pInd->getId() << endl;
#endif
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	if (pSpecies == popns[i]->getSpecies()) {
		trfrRules trfr = pSpecies->getTrfr();
		popns[i]->recruit(pPch,pInd,trfr.moveModel,trfr.moveType,newgroup);
	}
}
}

#if PEDIGREE
void SubCommunity::outGroups(Pedigree *pPed,int rep,int yr,int gen,bool patchmodel) {   
int npops = (int)popns.size();
int patchnum = pPatch->getPatchNum();
#if RSDEBUG
DEBUGLOG << "SubCommunity::outGroups(): this=" << this
	<< " patchnum=" << patchnum << " npops=" << npops << endl;
#endif
for (int i = 0; i < npops; i++) { // all populations
		popns[i]->outGroups(pPed,rep,yr,gen,patchmodel);
}
}
#endif // PEDIGREE

#endif // GROUPDISP 

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
#if SEASONAL || RS_RCPP
int SubCommunity::transfer(Landscape *pLandscape,short landIx,short nextseason) 
#else
int SubCommunity::transfer(Landscape *pLandscape,short landIx) 
#endif // SEASONAL || RS_RCPP
{
#if RSDEBUG
//DEBUGLOG << "SubCommunity::transfer(): this=" << this
//	<< endl;
#endif
int ndispersers = 0;
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
#if GROUPDISP
	Species *pSpecies = popns[i]->getSpecies();
	emigRules emig = pSpecies->getEmig();
	if (emig.groupdisp)  
		ndispersers += popns[i]->grouptransfer(pLandscape,landIx);
	else
		ndispersers += popns[i]->transfer(pLandscape,landIx);
#else
#if SEASONAL || RS_RCPP
	ndispersers += popns[i]->transfer(pLandscape,landIx,nextseason);
#else
	ndispersers += popns[i]->transfer(pLandscape,landIx);
#endif // SEASONAL || RS_RCPP
#endif // GROUPDISP
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

#if PEDIGREE
void SubCommunity::completeDispersal(Landscape *pLandscape,Pedigree *pPed,bool connect) 
#else
void SubCommunity::completeDispersal(Landscape *pLandscape,bool connect) 
#endif
{
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): this=" << this
//	<< endl;
#endif
//popStats pop;
//int popsize,subcomm;
int popsize;
#if GROUPDISP
dispgroup groupsettler;
Individual *pInd;
#endif
disperser settler;
Species *pSpecies;
Population *pPop;
Patch *pPrevPatch;
Patch *pNewPatch;
Cell *pPrevCell;
#if GROUPDISP
Cell *pNewCell;
#endif
SubCommunity *pSubComm;

int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	pSpecies = popns[i]->getSpecies();
#if GROUPDISP
	emigRules emig = pSpecies->getEmig();
	if (emig.groupdisp) popsize = popns[i]->getNGroups();
	else popsize = popns[i]->getNInds();
#else
	popsize = popns[i]->getNInds();
#endif
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
#if GROUPDISP
		if (emig.groupdisp) {
			groupsettler = popns[i]->extractGroupSettler(j);
			settled = groupsettler.yes;
			pNewCell = groupsettler.pGroup->getLocn(1); // current location of group
		}
		else {
			settler = popns[i]->extractSettler(j);
			settled = settler.yes;
			pNewCell = settler.pInd->getLocn(1); // current location of individual
		}
#else
		settler = popns[i]->extractSettler(j);
		settled = settler.yes;
#endif
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
#if GROUPDISP
			if (emig.groupdisp) {
				pNewPatch = (Patch*)groupsettler.pCell->getPatch();
			}
			else {
				pNewPatch = (Patch*)settler.pCell->getPatch();
			}
#else
			pNewPatch = (Patch*)settler.pCell->getPatch();
#endif
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
#if PEDIGREE
				pPop = pSubComm->newPopn(pLandscape,pSpecies,pPed,pNewPatch,0);
#else
				pPop = pSubComm->newPopn(pLandscape,pSpecies,pNewPatch,0);
#endif
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): j=" << j
//	<< " pPop=" << pPop << endl;
#endif
			}
#if GROUPDISP
			// recruit each individual within the group to the population in which
			// the group has settled
			// NB THIS CODE SHOULD BE RATIONALISED IF GROUP DISPERSAL BECOMES STANDARD
			if (emig.groupdisp) {
				short status = groupsettler.pGroup->getStatus();
				for (int k = 0; k < groupsettler.groupsize; k++) {
					pInd = groupsettler.pGroup->getMember(k);
					pPop->recruit(pInd);
					pInd->setStatus(status);
					pInd->moveTo(pNewCell);
					pathSteps p = groupsettler.pGroup->getSteps();
					pInd->setYearSteps(p.year);
					if (connect) { // increment connectivity totals
						int newpatch = pNewPatch->getSeqNum();
						pPrevCell = pInd->getLocn(0); // previous cell
						intptr patch = pPrevCell->getPatch();
						if (patch != 0) {
							pPrevPatch = (Patch*)patch;
							int prevpatch = pPrevPatch->getSeqNum();
							pLandscape->incrConnectMatrix(prevpatch,newpatch);
						}
					}
				}
				// all group processing is now finished
				delete groupsettler.pGroup;
			}
			else {
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
#else
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
#endif
		}
		else { // for group dispersal only
#if GROUPDISP
			if (emig.groupdisp) {
				// deal with groups which have died during dispersal
				// in order to be reported in individual-level output, each Individual
				// within the group must be recruited to the matrix Population itself
				short status = groupsettler.pGroup->getStatus();
#if RSDEBUG
//locn xloc = pNewCell->getLocn();
//DEBUGLOG << "SubCommunity::completeDispersal(): j=" << j
//	<< " settled=" << settled
//	<< " groupID=" << groupsettler.pGroup->getGroupID()
//	<< " status=" << status
//	<< " x=" << xloc.x << " y=" << xloc.y
//	<< endl;
#endif
				if (status > 5 && status < 9) {
					for (int k = 0; k < groupsettler.groupsize; k++) {
						pInd = groupsettler.pGroup->getMember(k);
#if RSDEBUG
//locn xloc = pNewCell->getLocn();
//DEBUGLOG << "SubCommunity::completeDispersal(): 1111 k=" << k
//	<< " pInd=" << pInd << " indID=" << pInd->getId() << " status=" << pInd->getStatus()
//	<< endl;
#endif
						pInd->setStatus(status);
						pInd->moveTo(pNewCell);
						pathSteps p = groupsettler.pGroup->getSteps();
						pInd->setYearSteps(p.year);
#if RSDEBUG
//DEBUGLOG << "SubCommunity::completeDispersal(): 2222 k=" << k
//	<< " pInd=" << pInd << " indID=" << pInd->getId() << " status=" << pInd->getStatus()
//	<< " p.year=" << p.year
//	<< endl;
#endif
						popns[i]->recruit(pInd);
					}
				// all group processing is now finished
				// delete the group so that Individuals are not repeatedly recruited
				// to the matrix Population
				popns[i]->deleteGroup(j);
				}
			}
#endif // GROUPDISP
		}
	}
#if GROUPDISP
	if (emig.groupdisp) {
		// remove pointers in the matrix popn to settler groups
		popns[i]->groupclean();
	}
	else {
		popns[i]->clean();
	}
#else
	// remove pointers in the matrix popn to settlers
	popns[i]->clean();
#endif
#if RSDEBUG
//pop = popns[i]->getStats();
#if GROUPDISP
if (emig.groupdisp) popsize = popns[i]->getNGroups();
else popsize = popns[i]->getNInds();
#else
popsize = popns[i]->getNInds();
#endif
//DEBUGLOG << "SubCommunity::completeDispersal(): i=" << i
//	<< " popns[i]=" << popns[i]
////	<< " pop.pPatch = " << pop.pPatch
////	<< " pop.nInds = " << pop.nInds
//	<< " popsize=" << popsize
//	<< endl;
#endif
}

}

#if GROUPDISP
// Delete dispersal groups once dispersal has finished
// This function is executed for the matrix patch only
void SubCommunity::deleteGroups(void) {
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popns[i]->deleteGroups();
}
}
#endif

//---------------------------------------------------------------------------

#if SEASONAL
void SubCommunity::survival(short season,short part,short option0,short option1)
#else
#if SPATIALMORT
void SubCommunity::survival(short part,short period,short option0,short option1)
#else
#if PEDIGREE
void SubCommunity::survival(Pedigree *pPed,short part,short option0,short option1)
#else
void SubCommunity::survival(short part,short option0,short option1)
#endif // PEDIGREE 
#endif // SPATIALMORT
#endif // SEASONAL
{
int npops = (int)popns.size();
if (npops < 1) return;
#if SPATIALMORT
// determine spatial mortality for the patch
// NOTE: if model is cell-based, then correct mortality for the cell is applied
// if model is patch-based, the mortality for a random cell is applied, and therefore,
// if mortalities differ within the patch, they will be averaged for a large population
// over the course of many years
float mort = 0.0;
simParams sim = paramsSim->getSim();
if (sim.mortMapLoaded) {
	Cell *pCell = pPatch->getRandomCell();
	if (pCell != 0) {
		mort = pCell->getMort(period);
	}
}
#endif
#if SPATIALDEMOG
std::vector <float> localDemoScaling;
#endif //SPATIALDEMOG
if (part == 0) {
#if SEASONAL
	float localK = pPatch->getK(season);
#else
#if SPATIALDEMOG
	localDemoScaling = pPatch->getDemoScaling();
#endif //SPATIALDEMOG
	float localK = pPatch->getK();
#endif // SEASONAL 
	for (int i = 0; i < npops; i++) { // all populations
#if RS_CONTAIN
#if SEASONAL
		popns[i]->survival0(localK,habIndex,season,option0,option1);
#else
		popns[i]->survival0(localK,habIndex,option0,option1);
#endif // SEASONAL
#else
#if SEASONAL
		popns[i]->survival0(localK,season,option0,option1);
#else
#if SPATIALMORT
		popns[i]->survival0(localK,mort,option0,option1); // SPATIALLY VARYING MORTALITY
#else
#if PEDIGREE
		popns[i]->survival0(pPed,localK,option0,option1);
#else
#if SPATIALDEMOG
		popns[i]->survival0(localK,option0,option1,localDemoScaling);
#else
		popns[i]->survival0(localK,option0,option1);
#endif // SPATIALDEMOG
#endif // PEDIGREE 
#endif // SPATIALMORT 
#endif // SEASONAL
#endif // RS_CONTAIN 
	}
}
else {
	for (int i = 0; i < npops; i++) { // all populations
#if PEDIGREE
		popns[i]->survival1(pPed);
#else
		popns[i]->survival1();
#endif 
	}
}
}

#if RS_CONTAIN

short SubCommunity::findCullTarget(Cull *pCull,int year,int nstages,int resol)
{
int ncells = pPatch->getNCells();
culldata c = pCull->getCullData();      
int nthreshold = (int)(c.densThreshold * (double)ncells * (double)resol * (double)resol / 10000.0);
#if RSDEBUG
//DEBUGLOG << "SubCommunity::findCullTarget(): year=" << year
//	<< " nstages=" << nstages << " densThreshold=" << c.densThreshold 
//	<< " PatchNum=" << pPatch->getPatchNum() << " ncells=" << ncells 
//	<< " nthreshold=" << nthreshold
//	<< endl;
#endif
short target = 0;
cullCount = 0;
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	if (nstages < 2) { // non-structured
//		if (popns[i]->getNInds() >= c.popnThreshold) 
		if (popns[i]->getNInds() > nthreshold)
		{
			cullTarget = true; target = 1;
			cullCount = popns[i]->getNInds();
			if (firstYear < 0) firstYear = year;
		}	
	}
	else { // stage-structured
		int nInds = 0;
		for (short stg = 0; stg < nstages; stg++) {
			if (pCull->getCullStage(stg)) {
				nInds += popns[i]->stagePop(stg);
			}
#if RSDEBUG
//DEBUGLOG << "SubCommunity::findCullTarget(): stg=" << stg
//	<< " nInds=" << nInds << " getCullStage=" << pCull->getCullStage(stg)
//	<< endl;
#endif
		}
//		if (nInds >= c.popnThreshold) 
		if (nInds > nthreshold) 
		{
			cullTarget = true; target = 1;
			cullCount = nInds;
			if (firstYear < 0) firstYear = year;
		}			
	}
}
// now convert the true cull count to an estimate allowing for the variance 
// in counting accuracy implied by the count c.v.
#if RSDEBUG
//DEBUGLOG << "SubCommunity::findCullTarget(): BEFORE cullCount=" << cullCount
//	<< " countCV=" << c.countCV 
//	<< endl;
#endif
if (cullTarget) {
	double countSD = (double)cullCount * c.countCV / 100.0;
	cullCount = (int)(0.5 + pRandom->Normal((double)cullCount,countSD));
	if (cullCount < 0) cullCount = 0;
}
#if RSDEBUG
//DEBUGLOG << "SubCommunity::findCullTarget(): AFTER cullCount=" << cullCount
//	<< endl;
#endif
#if RSDEBUG
//DEBUGLOG << "SubCommunity::findCullTarget(): target=" << target
//	<< " cullTarget=" << cullTarget << " firstYear=" << firstYear
//	<< endl;
#endif
return target;
}

bool SubCommunity::isCullTarget(void) { return cullTarget; }

int SubCommunity::initialYear(void) { return firstYear; }

double SubCommunity::damageIndex(void) { return pPatch->getDamageIndex(); }

void SubCommunity::resetCullTarget(void) { cullTarget = false; }

void SubCommunity::resetCull(void) {
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) popns[i]->resetCull();
} 

//void SubCommunity::cullPatch(Cull *pCull,int pop,float cullrate) 
void SubCommunity::cullPatch(Cull *pCull,int pop,int resol)
{ 
int ncells = pPatch->getNCells();
//popns[pop]->cull(pCull,cullrate);
popns[pop]->cull(pCull,(double)(ncells * resol * resol)/10000.0);
cullTarget = false;
}

// Record any damage caused by a resident population within the patch
void SubCommunity::updateDamage(Landscape *pLandscape,Species *pSpecies,Cull *pCull) {
int ncells = pPatch->getNCells();
int npops = (int)popns.size();
#if RSDEBUG
DEBUGLOG << "SubCommunity::updateDamage(): PatchNum=" << pPatch->getPatchNum()
	<< " ncells=" << ncells << " npops=" << npops
	<< " hasDamageLocns=" << (int)pPatch->hasDamageLocns() 
	<< endl;
#endif
if (!pPatch->hasDamageLocns()) return;
landParams ppLand = pLandscape->getLandParams();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
damageparams d = pDamageParams->getDamageParams();
popStats p;
Cell *pCell;
DamageLocn *pDamageLocn;
double densityfactor = 10000.0 / (double)(ncells * ppLand.resol * ppLand.resol);
// find damage locations within the patch
for (int i = 0; i < ncells; i++) {
	pCell = pPatch->getCell(i);
	if (pCell != 0) {
		pDamageLocn = pCell->getDamage();
		if (pDamageLocn != 0) {
			for (int j = 0; j < npops; j++) { // all populations
				p = popns[j]->getStats(0); 
				int ninds = 0;     
				switch (d.occOption) {				 
				case 0: // total population size
					pDamageLocn->updateOccupancyDamage((double)p.nInds);
					break;
				case 1: // total population density
#if RSDEBUG
DEBUGLOG << "SubCommunity::updateDamage(): occOption=1: nInds=" << p.nInds
	<< " ncells=" << ncells << " resol=" << ppLand.resol
	<< " density=" << (double)p.nInds * densityfactor 
	<< " pDamageLocn=" << pDamageLocn 
	<< endl;
#endif
					pDamageLocn->updateOccupancyDamage((double)p.nInds * densityfactor);
					break;
				case 2: // density of culled stages
					if (!dem.stageStruct) break; // condition should not occur
					for (int s = 0; s < sstruct.nStages; s++) {
						if (pCull->getCullStage(s)) { // stage is to be culled
							ninds += popns[j]->stagePop(s);
						}
					}		
#if RSDEBUG
DEBUGLOG << "SubCommunity::updateDamage(): occOption=2: ninds=" << ninds
	<< " ncells=" << ncells << " resol=" << ppLand.resol
	<< " density=" << (double)ninds * densityfactor 
	<< " pDamageLocn=" << pDamageLocn 
	<< endl;
#endif
					pDamageLocn->updateOccupancyDamage((double)ninds * densityfactor);				
					break;
				case 3: // stage-specific density
					if (!dem.stageStruct) break; // condition should not occur
#if RSDEBUG
DEBUGLOG << "SubCommunity::updateDamage(): occOption=3: nInds=" << p.nInds
	<< " stage=" << d.stage << " N=" << popns[j]->stagePop(d.stage)
	<< " ncells=" << ncells << " resol=" << ppLand.resol
	<< " density=" << (double)popns[j]->stagePop(d.stage) * densityfactor 
	<< " pDamageLocn=" << pDamageLocn 
	<< endl;
#endif
					pDamageLocn->updateOccupancyDamage((double)popns[j]->stagePop(d.stage) * densityfactor);				
					break;
				}
			}
		}     
	}
}
}

int SubCommunity::getCullCount(void) { return cullCount; }
	
double SubCommunity::prevDamage(void) { return pPatch->getPrevDamage(); }

#endif // RS_CONTAIN 

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
#if RS_CONTAIN
	pop = popns[i]->getStats(habIndex);
#else
#if SPATIALDEMOG
	pop = popns[i]->getStats(pPatch->getDemoScaling());
#else
	pop = popns[i]->getStats();
#endif // SPATIALDEMOG
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
	pop = popns[i]->getStats(habIndex);
#else
#if SPATIALDEMOG
	pop = popns[i]->getStats(pPatch->getDemoScaling());
#else
	pop = popns[i]->getStats();
#endif // SPATIALDEMOG
#endif // RS_CONTAIN 
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
#if PEDIGREE
	pPop = new Population(pSpecies,0,pPatch,0,land.resol);
#else
	pPop = new Population(pSpecies,pPatch,0,land.resol);
#endif // PEDIGREE
	fileOK = pPop->outPopHeaders(land.landNum,land.patchModel);
	delete pPop;
}
return fileOK;
}

// Write records to population file
#if RS_ABC
void SubCommunity::outPop(Landscape *pLandscape,int rep,int yr,int gen,
	ABCmaster *pABCmaster,bool abcYear,bool popOutputYear)
#else
void SubCommunity::outPop(Landscape *pLandscape,int rep,int yr,int gen)
#endif
{
#if RS_ABC
#if RSDEBUG
//DEBUGLOG << "SubCommunity::outPop(): subCommNum=" << subCommNum
//	<< " yr=" << yr << " pPatch="<< pPatch
//	<< " abcYear=" << abcYear << " popOutputYear="<< popOutputYear
//	<< endl;
#endif
#endif
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
#if SEASONAL
	localK = pPatch->getK(gen);
#else
	localK = pPatch->getK();
#endif // SEASONAL 
#if RS_ABC
	if (popOutputYear) {
		if (localK > 0.0 || (land.patchModel && patchnum == 0)) {
			popns[i]->outPopulation(rep,yr,gen,eps,land.patchModel,writeEnv,gradK,
				abcYear,pABCmaster);
		}
		else {
			if (popns[i]->totalPop() > 0) {
				popns[i]->outPopulation(rep,yr,gen,eps,land.patchModel,writeEnv,gradK,
					abcYear,pABCmaster);
			}
		}
	}
	else {
		if (abcYear) {
			if (pABCmaster->sampledPatch(pPatch)) {
				popns[i]->outPopulation(rep,yr,gen,eps,land.patchModel,writeEnv,gradK,
					abcYear,pABCmaster);
			}
		}
	}
#else
	if (localK > 0.0 || (land.patchModel && patchnum == 0)) {
		popns[i]->outPopulation(rep,yr,gen,eps,land.patchModel,writeEnv,gradK);
	}
	else {
		if (popns[i]->totalPop() > 0) {
			popns[i]->outPopulation(rep,yr,gen,eps,land.patchModel,writeEnv,gradK);
		}
	}
#endif // RS_ABC 
}
}

#if RS_CONTAIN

// Open cull file and write header record
bool SubCommunity::outCullHeaders(Landscape *pLandscape,Species *pSpecies,int option)
{
bool fileOK;
Population *pPop;
landParams land = pLandscape->getLandParams();

if (option == -999) { // close the file
	// as all populations may have been deleted, set up a dummy one
	// species is not necessary
	pPop = new Population();
	fileOK = pPop->outCullHeaders(pLandscape,-999,land.patchModel);
	delete pPop;
}
else { // open the file
	// as no population has yet been created, set up a dummy one
	// species is necessary, as columns depend on stage and sex structure
#if PEDIGREE
	pPop = new Population(pSpecies,0,pPatch,0,land.resol);
#else
	pPop = new Population(pSpecies,pPatch,0,land.resol);
#endif // PEDIGREE
	fileOK = pPop->outCullHeaders(pLandscape,land.landNum,land.patchModel);
	delete pPop;
}
return fileOK;
}

// Write records to cull file
void SubCommunity::outCull(Landscape *pLandscape,int rep,int yr,int gen)
{
#if RSDEBUG
//DEBUGLOG << "SubCommunity::outCull(): subCommNum=" << subCommNum
//	<< " yr=" << yr << " PatchNum="<< pPatch->getPatchNum()
//	<< " cullTarget=" << cullTarget << " firstYear="<< firstYear
//	<< endl;
#endif
landParams land = pLandscape->getLandParams();
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();

// generate output for each population within the sub-community (patch)
// provided that the patch is suitable (i.e. non-zero carrying capacity)
// or the population is above zero (possible if there is stochasticity or a moving gradient)
// or it is the matrix patch in a patch-based model
int npops = (int)popns.size();
int patchnum;
//Species* pSpecies;
Cell *pCell;
float localK;

patchnum = pPatch->getPatchNum();
for (int i = 0; i < npops; i++) { // all populations
//	pSpecies = popns[i]->getSpecies();
#if SEASONAL
	localK = pPatch->getK(gen);
#else
	localK = pPatch->getK();
#endif // SEASONAL 
	if (localK > 0.0 || (land.patchModel && patchnum == 0)) {
		popns[i]->outCullData(pLandscape,rep,yr,gen,land.patchModel);
	}
	else {
		if (popns[i]->totalPop() > 0) {
			popns[i]->outCullData(pLandscape,rep,yr,gen,land.patchModel);
		}
	}
}
}

#endif // RS_CONTAIN 

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
#if GROUPDISP || ROBFITT
void SubCommunity::outGenetics(int rep,int yr,int gen,int landNr,bool patchmodel)
{
//landParams ppLand = pLandscape->getLandParams();
if (landNr >= 0) { // open the file
	popns[0]->outGenetics(rep,yr,landNr,patchmodel);
	return;
}
if (landNr == -999) { // close the file
	popns[0]->outGenetics(rep,yr,landNr,patchmodel);
	return;
}
// generate output for each population within the sub-community (patch)
#if ROBFITT
if (patchRequired(pPatch->getPatchNum())) {
#endif
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popns[i]->outGenetics(rep,yr,landNr,patchmodel);
}
#if ROBFITT
}
#endif
}
#else
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
#endif

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
#if RS_CONTAIN
			if (trfr.kernType == 1)
#else
			if (trfr.twinKern)
#endif // RS_CONTAIN 
				outtraits << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
		}
		else {
			outtraits << "\tmean_distI\tstd_distI";
#if RS_CONTAIN
			if (trfr.kernType == 1)
#else
			if (trfr.twinKern)
#endif // RS_CONTAIN 
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
#if VCL
landPix p = pLandscape->getLandPix();
simView v = paramsSim->getViews();
rgb colour;
int nFactor;
double nsd = NSD; // no. of s.d. to use to control range for display
int ixt = 0; // indexes trait images
#endif

// generate output for each population within the sub-community (patch)
// provided that the patch is suitable (i.e. non-zero carrying capacity)
int npops = (int)popns.size();
Species* pSpecies;
float localK;

for (int i = 0; i < npops; i++) { // all populations
#if SEASONAL
	localK = pPatch->getK(gen);
#else
	localK = pPatch->getK();
#endif // SEASONAL 
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
#if VCL

			if (v.viewTraits && !commlevel) {
				emigParams elim0,elim1;
				double min,max;
				elim0 = pSpecies->getEmigParams(0,0);
				min = elim0.d0Mean - nsd * elim0.d0Scale; if (min < 0.0) min = 0.0;
				max = elim0.d0Mean + nsd * elim0.d0Scale; if (max > 1.0) max = 1.0;
				nFactor = (int)((mnD0[0] - min)*768/(max - min));
				colour = draw_wheel(nFactor);
				pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
				if (emig.sexDep) {
					elim1 = pSpecies->getEmigParams(0,1);
					min = elim1.d0Mean - nsd * elim1.d0Scale; if (min < 0.0) min = 0.0;
					max = elim1.d0Mean + nsd * elim1.d0Scale; if (max > 1.0) max = 1.0;
					nFactor = (int)((mnD0[1] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					if (emig.densDep) {
						min = elim0.alphaMean - nsd * elim0.alphaScale;
						max = elim0.alphaMean + nsd * elim0.alphaScale;
						nFactor = (int)((mnAlpha[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						min = elim1.alphaMean - nsd * elim1.alphaScale;
						max = elim1.alphaMean + nsd * elim1.alphaScale;
						nFactor = (int)((mnAlpha[1] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						min = elim0.betaMean - nsd * elim0.betaScale;
						max = elim0.betaMean + nsd * elim0.betaScale;
						nFactor = (int)((mnBeta[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						min = elim1.betaMean - nsd * elim1.betaScale;
						max = elim1.betaMean + nsd * elim1.betaScale;
						nFactor = (int)((mnBeta[1] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					}
				}
				else {
					if (emig.densDep) {
						min = elim0.alphaMean - nsd * elim0.alphaScale;
						max = elim0.alphaMean + nsd * elim0.alphaScale;
						nFactor = (int)((mnAlpha[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						min = elim0.betaMean - nsd * elim0.betaScale;
						max = elim0.betaMean + nsd * elim0.betaScale;
						nFactor = (int)((mnBeta[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					}
				}
			}

#endif
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
#if RS_CONTAIN
						if (trfr.kernType == 1) 
#else
						if (trfr.twinKern) 
#endif // RS_CONTAIN 
						{
							outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
							outtraits << "\t" << mnDist2[1] << "\t" << sdDist2[1];
							outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
							outtraits << "\t" << mnProp1[1] << "\t" << sdProp1[1];
						}
					}
					else { // sex-independent
						outtraits << "\t" << mnDist1[0] << "\t" << sdDist1[0];
#if RS_CONTAIN
						if (trfr.kernType == 1) 
#else
						if (trfr.twinKern) 
#endif // RS_CONTAIN 
						{
							outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
							outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
						}
					}
				}
			}
#if VCL
			if (v.viewTraits && !commlevel) {
				if (trfr.moveModel) {
					if (trfr.moveType == 1) {
						trfrSMSParams slim0 = pSpecies->getSMSParams(0,0);
						double min,max;
						min = slim0.dpMean - nsd * slim0.dpScale; if (min < 1.0) min = 1.0;
						max = slim0.dpMean + nsd * slim0.dpScale;
						nFactor = (int)((mnDP[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						min = slim0.gbMean - nsd * slim0.gbScale; if (min < 1.0) min = 1.0;
						max = slim0.gbMean + nsd * slim0.gbScale;
						nFactor = (int)((mnGB[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						trfrSMSTraits sms = pSpecies->getSMSTraits();
						if (sms.goalType == 2) {
							min = slim0.alphaDBMean - nsd * slim0.alphaDBScale; if (min < 1.0) min = 1.0;
							max = slim0.alphaDBMean + nsd * slim0.alphaDBScale;
							nFactor = (int)((mnAlphaDB[0] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
							min = slim0.betaDBMean - nsd * slim0.betaDBScale; if (min < 1.0) min = 1.0;
							max = slim0.betaDBMean + nsd * slim0.betaDBScale;
							nFactor = (int)((mnBetaDB[0] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						}
					}
					if (trfr.moveType == 2) {
						trfrCRWParams clim0 = pSpecies->getCRWParams(0,0);
						double min,max;
						min = clim0.stepLgthMean - nsd * clim0.stepLScale; if (min < 1.0) min = 1.0;
						max = clim0.stepLgthMean + nsd * clim0.stepLScale;
						nFactor = (int)((mnStepL[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						min = clim0.rhoMean - nsd * clim0.rhoScale; if (min < 0.0) min = 0.0;
						max = clim0.rhoMean + nsd * clim0.rhoScale; if (max > 1.0) max = 1.0;
						nFactor = (int)((mnRho[0] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					}
				}
				else { // dispersal kernel
					float minDist; // min distance limit to trait maps
					minDist = (float)land.resol;
					if (pSpecies->useFullKernel()) minDist = 0.0;
					trfrKernParams klim0,klim1;
					klim0 = pSpecies->getKernParams(0,0);
					double min,max;
					min = klim0.dist1Mean - nsd * klim0.dist1Scale; if (min < minDist) min = minDist;
					max = klim0.dist1Mean + nsd * klim0.dist1Scale;
					nFactor = (int)((mnDist1[0] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					if (trfr.sexDep) {
						klim1 = pSpecies->getKernParams(0,1);
						min = klim1.dist1Mean - nsd * klim1.dist1Scale; if (min < minDist) min = minDist;
						max = klim1.dist1Mean + nsd * klim1.dist1Scale;
						nFactor = (int)((mnDist1[1] - min)*768/(max - min));
						colour = draw_wheel(nFactor);
						pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
#if RS_CONTAIN
						if (trfr.kernType == 1) 
#else
						if (trfr.twinKern) 
#endif // RS_CONTAIN 
						{
							min = klim0.dist2Mean - nsd * klim0.dist2Scale; if (min < minDist) min = minDist;
							max = klim0.dist2Mean + nsd * klim0.dist2Scale;
							nFactor = (int)((mnDist2[0] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
							min = klim1.dist2Mean - nsd * klim1.dist2Scale; if (min < minDist) min = minDist;
							max = klim1.dist2Mean + nsd * klim1.dist2Scale;
							nFactor = (int)((mnDist2[1] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
							min = klim0.PKern1Mean - nsd * klim0.PKern1Scale; if (min < 0.0) min = 0.0;
							max = klim0.PKern1Mean + nsd * klim0.PKern1Scale; if (max > 1.0) max = 1.0;
							nFactor = (int)((mnProp1[0] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
							min = klim1.PKern1Mean - nsd * klim1.PKern1Scale; if (min < 0.0) min = 0.0;
							max = klim1.PKern1Mean + nsd * klim1.PKern1Scale; if (max > 1.0) max = 1.0;
							nFactor = (int)((mnProp1[1] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						}
					}
					else { // sex-independent
#if RS_CONTAIN
						if (trfr.kernType == 1) 
#else
						if (trfr.twinKern) 
#endif // RS_CONTAIN 
						{
							min = klim0.dist2Mean - nsd * klim0.dist2Scale; if (min < minDist) min = minDist;
							max = klim0.dist2Mean + nsd * klim0.dist2Scale;
							nFactor = (int)((mnDist2[0] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
							min = klim0.PKern1Mean - nsd * klim0.PKern1Scale; if (min < 0.0) min = 0.0;
							max = klim0.PKern1Mean + nsd * klim0.PKern1Scale; if (max > 1.0) max = 1.0;
							nFactor = (int)((mnProp1[0] - min)*768/(max - min));
							colour = draw_wheel(nFactor);
							pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
						}
					}
				}
			}
#endif
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
#if VCL
			if (v.viewTraits && !commlevel) {
				settParams slim0,slim1;
				double min,max;
				slim0 = pSpecies->getSettParams(0,0);
				min = slim0.s0Mean - nsd * slim0.s0Scale; if (min < 0.0) min = 0.0;
				max = slim0.s0Mean + nsd * slim0.s0Scale; if (max > 1.0) max = 1.0;
				nFactor = (int)((mnS0[0] - min)*768/(max - min));
				colour = draw_wheel(nFactor);
				pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
				if (sett.sexDep) {
					slim1 = pSpecies->getSettParams(0,1);
					min = slim1.s0Mean - nsd * slim1.s0Scale; if (min < 0.0) min = 0.0;
					max = slim1.s0Mean + nsd * slim1.s0Scale; if (max > 1.0) max = 1.0;
					nFactor = (int)((mnS0[1] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					min = slim0.alphaSMean - nsd * slim0.alphaSScale;
					max = slim0.alphaSMean + nsd * slim0.alphaSScale;
					nFactor = (int)((mnAlpha[0] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					min = slim1.alphaSMean - nsd * slim1.alphaSScale;
					max = slim1.alphaSMean + nsd * slim1.alphaSScale;
					nFactor = (int)((mnAlpha[1] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					min = slim0.betaSMean - nsd * slim0.betaSScale;
					max = slim0.betaSMean + nsd * slim0.betaSScale;
					nFactor = (int)((mnBeta[0] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					min = slim1.betaSMean - nsd * slim1.betaSScale;
					max = slim1.betaSMean + nsd * slim1.betaSScale;
					nFactor = (int)((mnBeta[1] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
				}
				else {
					min = slim0.alphaSMean - nsd * slim0.alphaSScale;
					max = slim0.alphaSMean + nsd * slim0.alphaSScale;
					nFactor = (int)((mnAlpha[0] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
					min = slim0.betaSMean - nsd * slim0.betaSScale;
					max = slim0.betaSMean + nsd * slim0.betaSScale;
					nFactor = (int)((mnBeta[0] - min)*768/(max - min));
					colour = draw_wheel(nFactor);
					pPatch->drawCells(tcanv.pcanvas[ixt++],p.gpix,land.dimY,colour);
				}
			}
#endif
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

#if RS_ABC
dispstats SubCommunity::getDispStats(float resol) {
dispstats d,popnd;
d.nPhilo = d.nDisp = d.nSuccess = 0; d.sumDist = d.sumDist2 = 0.0;

int npops = (int)popns.size();
for (int i = 0; i < npops; i++) { // all populations
	popnd = popns[i]->getDispStats(resol);
	d.nPhilo += popnd.nPhilo;
	d.nDisp += popnd.nDisp;
	d.nSuccess += popnd.nSuccess;
	d.sumDist += popnd.sumDist;
	d.sumDist2 += popnd.sumDist2;
}

#if RSDEBUG
//if (npops > 0) {
//DEBUGLOG << "SubCommunity::getDispStats(): this=" << this
//	<< " nPhilo=" << d.nPhilo << " nDisp=" << d.nDisp << " nSuccess=" << d.nSuccess
//	<< " sumDist=" << d.sumDist << " sumDist2=" << d.sumDist2
//	<< endl;
//}
#endif

return d;
}
#endif

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


