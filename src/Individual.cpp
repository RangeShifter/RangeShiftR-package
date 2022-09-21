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

#include "Individual.h"
//---------------------------------------------------------------------------
#if RS_EMBARCADERO 
#pragma package(smart_init)
#endif

int Individual::indCounter = 0;

//---------------------------------------------------------------------------

#if GROUPDISP
Individual::Individual() // Default constructor
{
pPrevCell = 0; pCurrCell = 0; pNatalPatch = 0;
path = 0; crw = 0; smsData = 0;
emigtraits = 0; kerntraits = 0; setttraits = 0;
pGenome = 0;
}
#endif

// Individual constructor
#if RS_CONTAIN
Individual::Individual(Cell *pCell,Patch *pPatch,short stg,short a,short repInt,
	short mstg,float probmale,bool movt,short moveType)
#else
#if PARTMIGRN
Individual::Individual(Species *pSpecies,Cell *pCell,Patch *pPatch,short stg,short a, 
	short repInt,float probmale,bool movt,short moveType)
#else
Individual::Individual(Cell *pCell,Patch *pPatch,short stg,short a,short repInt,
	float probmale,bool movt,short moveType)
#endif // PARTMIGRN 
#endif // RS_CONTAIN 
{
indId = indCounter; indCounter++; // unique identifier for each individual
#if RSDEBUG
//DEBUGLOG << "Individual::Individual(): indId=" << indId
//	<< " stg=" << stg << " a=" << a << " probmale=" << probmale
//	<< endl;
#endif
#if GROUPDISP
parentId[0] = parentId[1] = groupId = -1;
//groupId = -1;
#if PEDIGREE
pParent[0] = pParent[1] = 0;
matPosn = -1;	
#endif // PEDIGREE
#endif // GROUPDISP
#if BUTTERFLYDISP
//nJuvs = 0;
pMate = 0;
#endif

stage = stg;
#if RS_CONTAIN
motherstage = mstg;	
#endif // RS_CONTAIN 
if (probmale <= 0.0) sex = 0;
else sex = pRandom->Bernoulli(probmale);
age = a;
status = 0;
#if SEASONAL
npatches = 0;
#if PARTMIGRN
migrnstatus = 0;
npatches = 6; // <======= TO BE SET AS A PARAMETER =============================
// set dispersal/migration status
double cumprop = 0.0;
//cumprop[0] = 0.0;
double r = pRandom->Random();
for (int i = 1; i < 7; i++) {
	cumprop += pSpecies->getPropDispMigrn(i);
	if (r < cumprop) {
		setMigrnStatus(i);
#if RSDEBUG
//DEBUGLOG << "Individual::Individual(): indId=" << indId
//	<< " i=" << i << " cumprop=" << cumprop << " r=" << r
//	<< endl;
#endif
		i = 7;
	}
}
#endif // PARTMIGRN 
#endif // SEASONAL

if (sex == 0 && repInt > 0) { // set no. of fallow seasons for female
	fallow = pRandom->IRandom(0,repInt);
}
else fallow = 9999;
isDeveloping = false;
#if GOBYMODEL
asocial = false;
#endif
#if SOCIALMODEL
asocial = false;
#endif
pPrevCell = pCurrCell = pCell;
pNatalPatch = pPatch;
#if SEASONAL
pPrevPatch = pPatch;
#endif
if (movt) {
	locn loc = pCell->getLocn();
	path = new pathData;
	path->year = 0; path->total = 0; path->out = 0;
#if SEASONAL
	path->season= 0;
#endif
	path->pSettPatch = 0; path->settleStatus = 0;
//	path->leftNatalPatch = false;
#if RS_RCPP
	path->pathoutput = 1;
#endif
	if (moveType == 1) { // SMS
		// set up location data for SMS
		smsData = new smsdata;
		smsData->dp = smsData->gb = smsData->alphaDB = 1.0;
		smsData->betaDB = 1; 
		smsData->prev.x = loc.x; smsData->prev.y = loc.y; // previous location
		smsData->goal.x = loc.x; smsData->goal.y = loc.y; // goal location - initialised for dispersal bias
#if PARTMIGRN
		smsData->goalType = 0;
#endif  // PARTMIGRN 
	}
	else smsData = 0;
	if (moveType == 2) { // CRW
		// set up continuous co-ordinates etc. for CRW movement
		crw = new crwParams;
		crw->xc = ((float)pRandom->Random()*0.999f) + (float)loc.x;
		crw->yc = (float)(pRandom->Random()*0.999f) + (float)loc.y;
		crw->prevdrn = (float)(pRandom->Random()*2.0 * PI);
		crw->stepL = crw->rho = 0.0;
	}
	else crw = 0;
}
else {
	path = 0; crw = 0; smsData = 0;
}
emigtraits = 0;
kerntraits = 0;
setttraits = 0;
pGenome = 0;
#if RSDEBUG
//locn currloc = pCurrCell->getLocn();
//DEBUGLOG << "Individual::Individual(): indId=" << indId
//	<< " x=" << currloc.x << " y=" << currloc.y
////	<< " smsData=" << smsData << " dp=" << smsData->dp
//	<< endl;
#endif
}

Individual::~Individual(void) {
if (path != 0) delete path;
if (crw != 0) delete crw;
if (smsData != 0) delete smsData;
if (emigtraits != 0) delete emigtraits;
if (kerntraits != 0) delete kerntraits;
if (setttraits != 0) delete setttraits;

if (pGenome != 0) delete pGenome;

#if SEASONAL
patches.clear();
#endif
}

//---------------------------------------------------------------------------

#if BUTTERFLYDISP
void Individual::setMate(Individual *pmate) {
if (pmate >= 0) pMate = pmate;
else pMate = 0;
}
//void Individual::setMated(short njuvs,Individual *pmate) {
//if (njuvs >= 0) nJuvs = njuvs;
//pMate = pmate;
//}
//int Individual::getNJuvs(void) { return (int)nJuvs; }
Individual* Individual::getMate(void) { return pMate; }
#endif

//---------------------------------------------------------------------------

// Set genes for individual variation from species initialisation parameters
void Individual::setGenes(Species *pSpecies,int resol) {
demogrParams dem = pSpecies->getDemogr();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
genomeData gen = pSpecies->getGenomeData();
simParams sim = paramsSim->getSim();
int ntraits;	// first trait for all/female expression, second for male expression

if (gen.trait1Chromosome) {
	pGenome = new Genome(pSpecies->getNChromosomes(),pSpecies->getNLoci(0),
		pSpecies->isDiploid());
}
else {
	pGenome = new Genome(pSpecies);
}
#if RSDEBUG
//DEBUGLOG << endl;
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " sex=" << sex
//	<< " trait1Chromosome=" << gen.trait1Chromosome << " pGenome=" << pGenome
//	<< endl;
#endif

int gposn = 0;	// current position on genome
int expr = 0;		// gene expression type - NOT CURRENTLY USED

#if GOBYMODEL
// NB for consistency, and possible future inclusion in general RS model,
// some redundant coding here mimics that for dispersal traits
int socialposn = 0;
socialposn = gposn;
socialParams socparams;
ntraits = 1;
double socval;
for (int g = 0; g < ntraits; g++) { // ONLY trait for females/all
	socparams = pSpecies->getSocialParams();
	socval = pRandom->Normal(0.0,socparams.socSD) / socparams.socScale;
	if (gen.trait1Chromosome) {
		pGenome->setGene(gposn++,expr,socval,gen.alleleSD);
	}
	else {
		pGenome->setTrait(pSpecies,gposn++,socval,gen.alleleSD);
	}
}
// record phenotypic trait
double genval = 0.0;
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		genval = pGenome->express(0,0,0);
	}
	else {
		genval = pGenome->express(pSpecies,0);
	}
}
socparams = pSpecies->getSocialParams();
double phenval = genval*socparams.socScale + socparams.socMean;
if (phenval < 0.0) asocial = true; else asocial = false;
#endif // GOBYMODEL

#if SOCIALMODEL
// NB for consistency, and possible future inclusion in general RS model,
// some redundant coding here mimics that for dispersal traits
int socialposn = 0;
socialposn = gposn;
socialParams socparams;
ntraits = 1;
double socval;
for (int g = 0; g < ntraits; g++) { // ONLY trait for females/all
	socparams = pSpecies->getSocialParams();
	socval = pRandom->Normal(0.0,socparams.socSD) / socparams.socScale;
	if (gen.trait1Chromosome) {
		pGenome->setGene(gposn++,expr,socval,gen.alleleSD);
	}
	else {
		pGenome->setTrait(pSpecies,gposn++,socval,gen.alleleSD);
	}
}
// record phenotypic trait
double genval = 0.0;
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		genval = pGenome->express(0,0,0);
	}
	else {
		genval = pGenome->express(pSpecies,0);
	}
}
socparams = pSpecies->getSocialParams();
double phenval = genval*socparams.socScale + socparams.socMean;
if (phenval < 0.0) asocial = true; else asocial = false;
#endif // SOCIALMODEL

//int emigposn = 0;
#if RSDEBUG
//DEBUGLOG << "Individual::setGenes(): emigration genes" << endl;
#endif
if (emig.indVar) { // set emigration genes
	int emigposn = gposn;
	double d0,alpha,beta;
	emigParams eparams;
//	emigScales scale = pSpecies->getEmigScales();
	if (emig.sexDep) { // must be a sexual species
		ntraits = 2;
	}
	else {
		if (dem.repType == 0) { // asexual reproduction (haploid)
			ntraits = 1;
		}
		else { // sexual reproduction
			ntraits = 1;
		}
	}
	for (int g = 0; g < ntraits; g++) { // first trait for females/all, second for males
		eparams = pSpecies->getEmigParams(0,g);
		d0 = pRandom->Normal(0.0,eparams.d0SD) / eparams.d0Scale;
		if (emig.densDep) {
			alpha = pRandom->Normal(0.0,eparams.alphaSD) / eparams.alphaScale;
			beta  = pRandom->Normal(0.0,eparams.betaSD) / eparams.betaScale;
		}
#if RSDEBUG
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " g=" << g
//	<< " eparams.d0Mean=" << eparams.d0Mean << " eparams.d0SD=" << eparams.d0SD
//	<< " eparams.d0Scale=" << eparams.d0Scale << " d0=" << d0
////	<< " log(d0/(1.0-d0))=" << log(d0/(1.0-d0))
//	<< endl;
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " g=" << g
//	<< " eparams.alphaMean=" << eparams.alphaMean << " eparams.alphaSD=" << eparams.alphaSD
//	<< " eparams.alphaScale=" << eparams.alphaScale << " alpha=" << alpha
//	<< endl;
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " g=" << g
//	<< " eparams.betaMean=" << eparams.betaMean << " eparams.betaSD=" << eparams.betaSD
//	<< " eparams.betaScale=" << eparams.betaScale << " beta=" << beta
//	<< endl;
#endif
		if (gen.trait1Chromosome) {
			pGenome->setGene(gposn++,expr,d0,gen.alleleSD);
			if (emig.densDep) {
				pGenome->setGene(gposn++,expr,alpha,gen.alleleSD);
				pGenome->setGene(gposn++,expr,beta,gen.alleleSD);
			}
		}
		else {
			pGenome->setTrait(pSpecies,gposn++,d0,gen.alleleSD);
			if (emig.densDep) {
				pGenome->setTrait(pSpecies,gposn++,alpha,gen.alleleSD);
				pGenome->setTrait(pSpecies,gposn++,beta,gen.alleleSD);
			}
		}
	}
	// record phenotypic traits
	if (emig.densDep) {
		setEmigTraits(pSpecies,emigposn,3,emig.sexDep);
	}
	else {
		setEmigTraits(pSpecies,emigposn,1,emig.sexDep);
	}
}

//int trfrposn = 0;
if (trfr.indVar) { // set transfer genes
	int trfrposn = gposn;
	if (trfr.sexDep) { // must be a sexual species
		ntraits = 2;
	}
	else {
		if (dem.repType == 0) { // asexual reproduction
			ntraits = 1;
		}
		else { // sexual reproduction
			ntraits = 1;
		}
	}
//	trfrScales scale = pSpecies->getTrfrScales();
	if (trfr.moveModel) {
		if (trfr.moveType == 1) { // set SMS genes
			double dp,gb,alphaDB,betaDB;
			trfrSMSParams smsparams = pSpecies->getSMSParams(0,0); // only traits for females/all
			trfrSMSTraits smstraits = pSpecies->getSMSTraits();
			dp = pRandom->Normal(0.0,smsparams.dpSD) / smsparams.dpScale;
			gb = pRandom->Normal(0.0,smsparams.gbSD) / smsparams.gbScale;
			if (smstraits.goalType == 2) {
				alphaDB = pRandom->Normal(0.0,smsparams.alphaDBSD) / smsparams.alphaDBScale;
				betaDB  = pRandom->Normal(0.0,smsparams.betaDBSD)  / smsparams.betaDBScale;
			}
			if (gen.trait1Chromosome) {
				pGenome->setGene(gposn++,expr,dp,gen.alleleSD);
				pGenome->setGene(gposn++,expr,gb,gen.alleleSD);
				if (smstraits.goalType == 2) {
					pGenome->setGene(gposn++,expr,alphaDB,gen.alleleSD);
					pGenome->setGene(gposn++,expr,betaDB,gen.alleleSD);
				}
			}
			else {
				pGenome->setTrait(pSpecies,gposn++,dp,gen.alleleSD);
				pGenome->setTrait(pSpecies,gposn++,gb,gen.alleleSD);
				if (smstraits.goalType == 2) {
					pGenome->setTrait(pSpecies,gposn++,alphaDB,gen.alleleSD);
					pGenome->setTrait(pSpecies,gposn++,betaDB,gen.alleleSD);
				}
			}
			// record phenotypic traits
			if (smstraits.goalType == 2)
				setSMSTraits(pSpecies,trfrposn,4,false);
			else
				setSMSTraits(pSpecies,trfrposn,2,false);
		}
		if (trfr.moveType == 2) { // set CRW genes
			double stepL,rho;
			trfrCRWParams m = pSpecies->getCRWParams(0,0); // only traits for females/all
			stepL = pRandom->Normal(0.0,m.stepLgthSD) / m.stepLScale;
			rho   = pRandom->Normal(0.0,m.rhoSD) / m.rhoScale;
			if (gen.trait1Chromosome) {
				pGenome->setGene(gposn++,expr,stepL,gen.alleleSD);
				pGenome->setGene(gposn++,expr,rho,gen.alleleSD);
			}
			else {
				pGenome->setTrait(pSpecies,gposn++,stepL,gen.alleleSD);
				pGenome->setTrait(pSpecies,gposn++,rho,gen.alleleSD);
			}
			// record phenotypic traits
			setCRWTraits(pSpecies,trfrposn,2,false);
		}
	}
	else { // set kernel genes
		double dist1,dist2,prob1;
		trfrKernParams k;
		for (int g = 0; g < ntraits; g++) { // first traits for females/all, second for males
			k = pSpecies->getKernParams(0,g);
			dist1 = pRandom->Normal(0.0,k.dist1SD) / k.dist1Scale;
#if RS_CONTAIN
			if (trfr.kernType == 1) 
#else
			if (trfr.twinKern) 
#endif // RS_CONTAIN 
			{
				dist2 = pRandom->Normal(0.0,k.dist2SD) / k.dist2Scale;
				prob1 = pRandom->Normal(0.0,k.PKern1SD) / k.PKern1Scale;
			}
			if (gen.trait1Chromosome) {
				pGenome->setGene(gposn++,expr,dist1,gen.alleleSD);
#if RS_CONTAIN
				if (trfr.kernType == 1) 
#else
				if (trfr.twinKern) 
#endif // RS_CONTAIN 
				{
					pGenome->setGene(gposn++,expr,dist2,gen.alleleSD);
					pGenome->setGene(gposn++,expr,prob1,gen.alleleSD);
				}
			}
			else {
				pGenome->setTrait(pSpecies,gposn++,dist1,gen.alleleSD);
#if RS_CONTAIN
				if (trfr.kernType == 1) 
#else
				if (trfr.twinKern) 
#endif // RS_CONTAIN 
				{
					pGenome->setTrait(pSpecies,gposn++,dist2,gen.alleleSD);
					pGenome->setTrait(pSpecies,gposn++,prob1,gen.alleleSD);
				}
      }
		}
		// record phenotypic traits
#if RS_CONTAIN
		if (trfr.kernType == 1) 
#else
		if (trfr.twinKern) 
#endif // RS_CONTAIN 
		{
			setKernTraits(pSpecies,trfrposn,3,resol,trfr.sexDep);
		}
		else {
			setKernTraits(pSpecies,trfrposn,1,resol,trfr.sexDep);
		}
	}
}

//int settposn = 0;
#if RSDEBUG
//DEBUGLOG << "Individual::setGenes(): settlement genes" << endl;
#endif
if (sett.indVar) {
	int settposn = gposn;
	double s0,alpha,beta;
	settParams sparams;
//	settScales scale = pSpecies->getSettScales();
	if (sett.sexDep) { // must be a sexual species
		ntraits = 2;
	}
	else {
		if (dem.repType == 0) { // asexual reproduction
			ntraits = 1;
		}
		else { // sexual reproduction
			ntraits = 1;
		}
	}
	for (int g = 0; g < ntraits; g++) { // first trait for females/all, second for males
		if (sim.batchMode) {
			sparams = pSpecies->getSettParams(0,g);
		}
		else { // individual variability not (yet) implemented as sex-dependent in GUI
			sparams = pSpecies->getSettParams(0,0);        			
		}
		s0 = pRandom->Normal(0.0,sparams.s0SD) / sparams.s0Scale;
		alpha = pRandom->Normal(0.0,sparams.alphaSSD) / sparams.alphaSScale;
		beta  = pRandom->Normal(0.0,sparams.betaSSD) / sparams.betaSScale;
#if RSDEBUG
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " g=" << g
//	<< " sparams.s0Mean=" << sparams.s0Mean
//	<< " sparams.s0SD=" << sparams.s0SD
//	<< " sparams.s0Scale=" << sparams.s0Scale
//	<< " s0=" << s0
//	<< endl;
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " g=" << g
//	<< " sparams.alphaSMean=" << sparams.alphaSMean
//	<< " sparams.alphaSSD=" << sparams.alphaSSD
//	<< " sparams.alphaSScale=" << sparams.alphaSScale
//	<< " alpha=" << alpha
//	<< endl;
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " g=" << g
//	<< " sparams.betaSMean=" << sparams.betaSMean
//	<< " sparams.betaSSD=" << sparams.betaSSD
//	<< " sparams.betaSScale=" << sparams.betaSScale
//	<< " beta=" << beta
//	<< endl;
#endif
		if (gen.trait1Chromosome) {
			pGenome->setGene(gposn++,expr,s0,gen.alleleSD);
			pGenome->setGene(gposn++,expr,alpha,gen.alleleSD);
			pGenome->setGene(gposn++,expr,beta,gen.alleleSD);
		}
		else {
			pGenome->setTrait(pSpecies,gposn++,s0,gen.alleleSD);
			pGenome->setTrait(pSpecies,gposn++,alpha,gen.alleleSD);
			pGenome->setTrait(pSpecies,gposn++,beta,gen.alleleSD);
		}
	}
	// record phenotypic traits
	setSettTraits(pSpecies,settposn,3,sett.sexDep);
}

if (!gen.trait1Chromosome) {
	if (gen.neutralMarkers || pSpecies->getNNeutralLoci() > 0) {
		pGenome->setNeutralLoci(pSpecies,gen.alleleSD);
	}
}
#if RSDEBUG
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " finished"
//	<< endl;
#endif
}

// Inherit genome from parent(s)
void Individual::setGenes(Species *pSpecies,Individual *mother,Individual *father,
	int resol)
{
#if RSDEBUG
//locn currloc = pCurrCell->getLocn();
//DEBUGLOG << "Individual::setGenes(): indId=" << indId
//	<< " x=" << currloc.x << " y=" << currloc.y
////	<< " pSpecies=" << pSpecies
//	<< " mother=" << mother
//	<< " motherID=" << mother->getId()
//	<< " father=" << father;
//if (father != 0) DEBUGLOG << " fatherID=" << father->getId();
//DEBUGLOG << endl;
#endif
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();

#if GROUPDISP
parentId[0] = mother->getId();
#if PEDIGREE
pParent[0] = mother;
#endif // PEDIGREE
#endif // GROUPDISP 

Genome *pFatherGenome;
#if GROUPDISP
if (father == 0) pFatherGenome = 0;
else {
	pFatherGenome = father->pGenome;
	parentId[1] = father->getId();
#if PEDIGREE
	pParent[1] = father;
#endif // PEDIGREE
}
#else
if (father == 0) pFatherGenome = 0; else pFatherGenome = father->pGenome;
#endif // GROUPDISP 

pGenome = new Genome(pSpecies,mother->pGenome,pFatherGenome);

#if GOBYMODEL
// record phenotypic trait
double genval = 0.0;
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		genval = pGenome->express(0,0,0);
	}
	else {
		genval = pGenome->express(pSpecies,0);
	}
}
socialParams socparams = pSpecies->getSocialParams();
double phenval = genval*socparams.socScale + socparams.socMean;
if (phenval < 0.0) asocial = true; else asocial = false;
#endif
#if SOCIALMODEL
// record phenotypic trait
double genval = 0.0;
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		genval = pGenome->express(0,0,0);
	}
	else {
		genval = pGenome->express(pSpecies,0);
	}
}
socialParams socparams = pSpecies->getSocialParams();
double phenval = genval*socparams.socScale + socparams.socMean;
if (phenval < 0.0) asocial = true; else asocial = false;
#endif

if (emig.indVar) {
	// record emigration traits
	if (father == 0) { // haploid
		if (emig.densDep) {
			setEmigTraits(pSpecies,0,3,0);
		}
		else {
			setEmigTraits(pSpecies,0,1,0);
		}
	}
	else { // diploid
		if (emig.densDep) {
			setEmigTraits(pSpecies,0,3,emig.sexDep);
		}
		else {
			setEmigTraits(pSpecies,0,1,emig.sexDep);
		}
	}
}

if (trfr.indVar) {
	// record movement model traits
	if (trfr.moveModel) {
		if (trfr.moveType == 1) { // SMS
			trfrSMSTraits s = pSpecies->getSMSTraits();
			if (s.goalType == 2)
				setSMSTraits(pSpecies,trfr.movtTrait[0],4,0);
			else
				setSMSTraits(pSpecies,trfr.movtTrait[0],2,0);
		}
		if (trfr.moveType == 2) { // CRW
			setCRWTraits(pSpecies,trfr.movtTrait[0],2,0);
		}
	}
	else { // kernel
		if (father == 0) { // haploid
#if RS_CONTAIN
			if (trfr.kernType == 1) 
#else
			if (trfr.twinKern) 
#endif // RS_CONTAIN 
			{
				setKernTraits(pSpecies,trfr.movtTrait[0],3,resol,0);
			}
			else {
				setKernTraits(pSpecies,trfr.movtTrait[0],1,resol,0);
			}
		}
		else { // diploid
#if RS_CONTAIN
			if (trfr.kernType == 1) 
#else
			if (trfr.twinKern) 
#endif // RS_CONTAIN 
			{
				setKernTraits(pSpecies,trfr.movtTrait[0],3,resol,trfr.sexDep);
			}
			else {
				setKernTraits(pSpecies,trfr.movtTrait[0],1,resol,trfr.sexDep);
			}
		}
	}
}

if (sett.indVar) {
	// record settlement traits
	if (father == 0) { // haploid
		setSettTraits(pSpecies,sett.settTrait[0],3,0);
	}
	else { // diploid
		setSettTraits(pSpecies,sett.settTrait[0],3,sett.sexDep);
//		setSettTraits(pSpecies,sett.settTrait[0],3,0);
	}
}

#if RSDEBUG
//emigParams e = getEmigTraits(0,1,0);
//DEBUGLOG << "Individual::setGenes(): indId=" << indId << " finished "
//	<< " d0=" << e.d0
////	<< " alpha=" << e.alpha << " beta=" << e.beta
//	<< endl;
#endif
}

#if VIRTUALECOLOGIST
Genome* Individual::getGenome(void) { return pGenome; }
#endif

//---------------------------------------------------------------------------

// Identify whether an individual is a potentially breeding female -
// if so, return her stage, otherwise return 0
int Individual::breedingFem(void) {
if (sex == 0) {
	if (status == 0 || status == 4 || status == 5) return stage;
	else return 0;
}
else return 0;
}

int Individual::getId(void) { return indId; }

#if GROUPDISP
int Individual::getParentId(short i) {
if (i >= 0 && i <= 1) return parentId[i];
else return -1;
}
#if PEDIGREE
Individual* Individual::getParent(short i) { 
if (i >= 0 && i <= 1) return pParent[i];
else return 0;
}
#endif // PEDIGREE
void Individual::setGroupId(int g) { if (g >= 0) groupId = g; }
int Individual::getGroupId(void) { return groupId; }
#endif // GROUPDISP

int Individual::getSex(void) { return sex; }

int Individual::getStatus(void) { return status; }
#if SEASONAL
#if PARTMIGRN
int Individual::getMigrnStatus(void) { return migrnstatus; }
#endif // PARTMIGRN 
#endif // SEASONAL

indStats Individual::getStats(void) {
indStats s;
s.stage = stage; s.sex = sex; s.age = age; s.status = status; s.fallow = fallow;
s.isDeveloping = isDeveloping;
#if SEASONAL
#if PARTMIGRN
s.migrnstatus = migrnstatus;
#endif // PARTMIGRN 
#endif // SEASONAL
#if GOBYMODEL
s.asocial = asocial;
#endif
#if SOCIALMODEL
s.asocial = asocial;
#endif
return s;
}
#if GOBYMODEL
bool Individual::isAsocial(void) { return asocial; }
#endif
#if SOCIALMODEL
bool Individual::isAsocial(void) { return asocial; }
#endif

Cell* Individual::getLocn(const short option) {
#if SEASONAL
if (option == 2) { // return a random location in the previous patch
//	Cell *pCell;
//	pCell = pPrevPatch->getRandomCell();
	return pPrevPatch->getRandomCell();

}
#endif
if (option == 0) { // return previous location
	return pPrevCell;
}
else { // return current location
	return pCurrCell;
}
}

Patch* Individual::getNatalPatch(void) { return pNatalPatch; }

void Individual::setYearSteps(int t) {
if (path != 0 && t >= 0) {                     
	if (t >= 0) path->year = t;
	else path->year = 666;
}
#if RSDEBUG
//DEBUGLOG << "Individual::setYearSteps(): indId=" << indId
//	<< " t=" << t << " path->year=" << path->year
//	<< endl;
#endif
}

pathSteps Individual::getSteps(void) {
pathSteps s;
if (path == 0) {
	s.year = 0; s.total = 0; s.out = 0;
#if SEASONAL
	s.season = 0;
#endif
}
else {
	s.year = path->year; s.total = path->total; s.out = path->out;
#if SEASONAL
	s.season = path->season;
#endif
}
return s;
}

settlePatch Individual::getSettPatch(void) {
settlePatch s;
if (path == 0) {
	s.pSettPatch = 0; s.settleStatus = 0;
}
else {
	s.pSettPatch = path->pSettPatch; s.settleStatus = path->settleStatus;
}
return s;
}

void Individual::setSettPatch(const settlePatch s) {
if (path == 0) {
	path = new pathData;
	path->year = 0; path->total = 0; path->out = 0; path->settleStatus = 0;
#if SEASONAL
	path->season = 0;
#endif
#if RS_RCPP
	path->pathoutput = 1;
#endif
}
if (s.settleStatus >= 0 && s.settleStatus <= 2) path->settleStatus = s.settleStatus;
path->pSettPatch = s.pSettPatch;             
}

#if SEASONAL
void Individual::resetPathSeason(void) {
if (path != 0) path->season = 0;
}
#endif

// Set phenotypic emigration traits
void Individual::setEmigTraits(Species *pSpecies,short emiggenelocn,short nemigtraits,
	bool sexdep) {
#if RSDEBUG
//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
//	<< " emiggenelocn=" << emiggenelocn << " nemigtraits=" << nemigtraits << " sexdep=" << sexdep
//	<< endl;
#endif
emigTraits e; e.d0 = e.alpha = e.beta = 0.0;
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		if (sexdep) {
			if (nemigtraits == 3) { // emigration is density-dependent
				e.d0    = (float)pGenome->express(emiggenelocn+3*sex,0,0);
				e.alpha = (float)pGenome->express(emiggenelocn+3*sex+1,0,0);
				e.beta  = (float)pGenome->express(emiggenelocn+3*sex+2,0,0);
			}
			else {
				e.d0    = (float)pGenome->express(emiggenelocn+sex,0,0);
			}
		}
		else {
			e.d0 = (float)pGenome->express(emiggenelocn,0,0);
			if (nemigtraits == 3) { // emigration is density-dependent
				e.alpha = (float)pGenome->express(emiggenelocn+1,0,0);
				e.beta  = (float)pGenome->express(emiggenelocn+2,0,0);
			}
		}
	}
	else {
		if (sexdep) {              
			if (nemigtraits == 3) { // emigration is density-dependent
				e.d0 = (float)pGenome->express(pSpecies,emiggenelocn+3*sex);
				e.alpha = (float)pGenome->express(pSpecies,emiggenelocn+3*sex+1);
				e.beta  = (float)pGenome->express(pSpecies,emiggenelocn+3*sex+2);
			}
			else {
				e.d0 = (float)pGenome->express(pSpecies,emiggenelocn+sex);
			}
		}
		else {
			e.d0 = (float)pGenome->express(pSpecies,emiggenelocn);
			if (nemigtraits == 3) { // emigration is density-dependent
				e.alpha = (float)pGenome->express(pSpecies,emiggenelocn+1);
				e.beta  = (float)pGenome->express(pSpecies,emiggenelocn+2);
			}
		}
	}
}
#if RSDEBUG
//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
//	<< " e.d0=" << e.d0 << " e.alpha=" << e.alpha << " e.beta=" << e.beta
//	<< endl;
#endif

emigParams eparams;
if (sexdep) {
	eparams = pSpecies->getEmigParams(0,sex);
}
else {
	eparams = pSpecies->getEmigParams(0,0);
}
#if RSDEBUG
//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
//	<< " eparams.betaMean=" << eparams.betaMean << " eparams.betaSD=" << eparams.betaSD 
//	<< " eparams.betaScale=" << eparams.betaScale
//	<< endl;
#endif
emigtraits = new emigTraits;
emigtraits->d0 = (float)(e.d0*eparams.d0Scale + eparams.d0Mean);
emigtraits->alpha = (float)(e.alpha*eparams.alphaScale + eparams.alphaMean);
emigtraits->beta = (float)(e.beta*eparams.betaScale + eparams.betaMean);
#if RSDEBUG
//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
//	<< " emigtraits->d0=" << emigtraits->d0
//	<< " emigtraits->alpha=" << emigtraits->alpha << " emigtraits->beta=" << emigtraits->beta
//	<< endl;
#endif
if (emigtraits->d0 < 0.0) emigtraits->d0 = 0.0;
if (emigtraits->d0 > 1.0) emigtraits->d0 = 1.0;
#if RSDEBUG
//DEBUGLOG << "Individual::setEmigTraits(): indId=" << indId
//	<< " emigtraits->d0=" << emigtraits->d0
//	<< " emigtraits->alpha=" << emigtraits->alpha << " emigtraits->beta=" << emigtraits->beta
//	<< endl;
#endif
return;
}

// Get phenotypic emigration traits
emigTraits Individual::getEmigTraits(void) {
#if RSDEBUG
//DEBUGLOG << "Individual::getEmigTraits(): indId=" << indId
//	<< endl;
#endif
emigTraits e; e.d0 = e.alpha = e.beta = 0.0;
if (emigtraits != 0) {
	e.d0 = emigtraits->d0;
	e.alpha = emigtraits->alpha;
	e.beta = emigtraits->beta;
}
#if RSDEBUG
//DEBUGLOG << "Individual::getEmigTraits(): indId=" << indId
//	<< " e.d0=" << e.d0 << " e.alpha=" << e.alpha << " e.beta=" << e.beta
//	<< endl;
#endif

return e;
}

// Set phenotypic transfer by kernel traits
void Individual::setKernTraits(Species *pSpecies,short kerngenelocn,short nkerntraits,
	int resol,bool sexdep) {
#if RSDEBUG
//DEBUGLOG << "Individual::setKernTraits(): indId=" << indId
//	<< " kerngenelocn=" << kerngenelocn << " nkerntraits=" << nkerntraits << " sexdep=" << sexdep
//	<< endl;
#endif
trfrKernTraits k; k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		if (sexdep) {
			if (nkerntraits == 3) { // twin kernel
				k.meanDist1 = (float)pGenome->express(kerngenelocn+3*sex,0,sex);
				k.meanDist2 = (float)pGenome->express(kerngenelocn+3*sex+1,0,sex);
				k.probKern1 = (float)pGenome->express(kerngenelocn+3*sex+2,0,sex);
			}
			else {
				k.meanDist1 = (float)pGenome->express(kerngenelocn+sex,0,sex);
			}
		}
		else {
			k.meanDist1 = (float)pGenome->express(kerngenelocn,0,0);
			if (nkerntraits == 3) { // twin kernel
				k.meanDist2 = (float)pGenome->express(kerngenelocn+1,0,0);
				k.probKern1 = (float)pGenome->express(kerngenelocn+2,0,0);
			}
		}
	}
	else {
		if (sexdep) {
			if (nkerntraits == 3) { // twin kernel
				k.meanDist1 = (float)pGenome->express(pSpecies,kerngenelocn+3*sex);
				k.meanDist2 = (float)pGenome->express(pSpecies,kerngenelocn+3*sex+1);
				k.probKern1 = (float)pGenome->express(pSpecies,kerngenelocn+3*sex+2);
			}
			else {
				k.meanDist1 = (float)pGenome->express(pSpecies,kerngenelocn+sex);
			}
		}
		else {
			k.meanDist1 = (float)pGenome->express(pSpecies,kerngenelocn);
			if (nkerntraits == 3) { // twin kernel
				k.meanDist2 = (float)pGenome->express(pSpecies,kerngenelocn+1);
				k.probKern1 = (float)pGenome->express(pSpecies,kerngenelocn+2);
			}
		}
  }
}
#if RSDEBUG
//DEBUGLOG << "Individual::setKernTraits(): indId=" << indId
//	<< " k.meanDist1=" << k.meanDist1 << " k.meanDist2=" << k.meanDist2
//	<< " k.probKern1=" << k.probKern1
//	<< endl;
#endif

trfrKernParams kparams;
if (sexdep) {
	kparams = pSpecies->getKernParams(0,sex);
}
else {
	kparams = pSpecies->getKernParams(0,0);
}
kerntraits = new trfrKernTraits;
kerntraits->meanDist1 = (float)(k.meanDist1*kparams.dist1Scale + kparams.dist1Mean);
kerntraits->meanDist2 = (float)(k.meanDist2*kparams.dist2Scale + kparams.dist2Mean);
kerntraits->probKern1 = (float)(k.probKern1*kparams.PKern1Scale + kparams.PKern1Mean);
#if RSDEBUG
//DEBUGLOG << "Individual::setKernTraits(): indId=" << indId
//	<< " kerntraits->meanDist1=" << kerntraits->meanDist1
//	<< " kerntraits->meanDist2=" << kerntraits->meanDist2
//	<< " kerntraits->probKern1=" << kerntraits->probKern1
//	<< endl;
#endif
if (!pSpecies->useFullKernel()) {
	// kernel mean(s) may not be less than landscape resolution
	if (kerntraits->meanDist1 < resol) kerntraits->meanDist1 = (float)resol;
	if (kerntraits->meanDist2 < resol) kerntraits->meanDist2 = (float)resol;
}
if (kerntraits->probKern1 < 0.0) kerntraits->probKern1 = 0.0;
if (kerntraits->probKern1 > 1.0) kerntraits->probKern1 = 1.0;
#if RSDEBUG
//DEBUGLOG << "Individual::setKernTraits(): indId=" << indId
//	<< " kerntraits->meanDist1=" << kerntraits->meanDist1
//	<< " kerntraits->meanDist2=" << kerntraits->meanDist2
//	<< " kerntraits->probKern1=" << kerntraits->probKern1
//	<< endl;
#endif
return;
}

// Get phenotypic emigration traits
trfrKernTraits Individual::getKernTraits(void) {
#if RSDEBUG
//DEBUGLOG << "Individual::getKernTraits(): indId=" << indId
//	<< endl;
#endif
trfrKernTraits k; k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;
if (kerntraits != 0) {
	k.meanDist1 = kerntraits->meanDist1;
	k.meanDist2 = kerntraits->meanDist2;
	k.probKern1 = kerntraits->probKern1;
}
#if RSDEBUG
//DEBUGLOG << "Individual::getKernTraits(): indId=" << indId
//	<< " k.meanDist1=" << k.meanDist1 << " k.meanDist2=" << k.meanDist1
//	<< " k.probKern1=" << k.probKern1
//	<< endl;
#endif

return k;
}

// Set phenotypic transfer by SMS traits
void Individual::setSMSTraits(Species *pSpecies,short SMSgenelocn,short nSMStraits,
		bool sexdep) {
#if RSDEBUG
//DEBUGLOG << "Individual::setSMSTraits(): indId=" << indId
//	<< " SMSgenelocn=" << SMSgenelocn << " nSMStraits=" << nSMStraits << " sexdep=" << sexdep
//	<< endl;
#endif
trfrSMSTraits s = pSpecies->getSMSTraits();
double dp,gb,alphaDB,betaDB;
dp = gb = alphaDB = betaDB = 0.0;
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		if (sexdep) {
			dp = pGenome->express(SMSgenelocn,0,0);
			gb = pGenome->express(SMSgenelocn+1,0,0);
			if (nSMStraits == 4) {
				alphaDB = pGenome->express(SMSgenelocn+2,0,0);
				betaDB  = pGenome->express(SMSgenelocn+3,0,0);
			}
		}
		else {
			dp = pGenome->express(SMSgenelocn,0,0);
			gb = pGenome->express(SMSgenelocn+1,0,0);
			if (nSMStraits == 4) {
				alphaDB = pGenome->express(SMSgenelocn+2,0,0);
				betaDB  = pGenome->express(SMSgenelocn+3,0,0);
			}
		}
	}
	else {
		if (sexdep) {
			dp = pGenome->express(pSpecies,SMSgenelocn);
			gb = pGenome->express(pSpecies,SMSgenelocn+1);
			if (nSMStraits == 4) {
				alphaDB = pGenome->express(pSpecies,SMSgenelocn+2);
				betaDB  = pGenome->express(pSpecies,SMSgenelocn+3);
			}
		}
		else {
			dp = pGenome->express(pSpecies,SMSgenelocn);
			gb = pGenome->express(pSpecies,SMSgenelocn+1);
			if (nSMStraits == 4) {
				alphaDB = pGenome->express(pSpecies,SMSgenelocn+2);
				betaDB  = pGenome->express(pSpecies,SMSgenelocn+3);
			}
		}
	}
}
#if RSDEBUG
//DEBUGLOG << "Individual::setSMSTraits(): indId=" << indId
//	<< " dp=" << dp << " gb=" << gb
//	<< " alphaDB=" << alphaDB << " betaDB=" << betaDB
//	<< endl;
#endif

trfrSMSParams smsparams;
if (sexdep) {
	smsparams = pSpecies->getSMSParams(0,0);
}
else {
	smsparams = pSpecies->getSMSParams(0,0);
}
smsData->dp = (float)(dp*smsparams.dpScale + smsparams.dpMean);        
smsData->gb = (float)(gb*smsparams.gbScale + smsparams.gbMean);
if (s.goalType == 2) {
	smsData->alphaDB = (float)(alphaDB*smsparams.alphaDBScale + smsparams.alphaDBMean);
	smsData->betaDB  = (int)(betaDB*smsparams.betaDBScale + smsparams.betaDBMean + 0.5);
}
else {
	smsData->alphaDB = s.alphaDB;
	smsData->betaDB  = s.betaDB;
}
#if RSDEBUG
//DEBUGLOG << "Individual::setSMSTraits() 1111: indId=" << indId
//	<< " smsData->dp=" << smsData->dp	<< " smsData->gb=" << smsData->gb
//	<< " smsData->alphaDB=" << smsData->alphaDB	<< " smsData->betaDB=" << smsData->betaDB
//	<< endl;
#endif
if (smsData->dp < 1.0) smsData->dp = 1.0;
if (smsData->gb < 1.0) smsData->gb = 1.0;
if (smsData->alphaDB <= 0.0) smsData->alphaDB = 0.000001f;
if (smsData->betaDB < 1) smsData->betaDB = 1;
#if RSDEBUG
//DEBUGLOG << "Individual::setSMSTraits() 2222: indId=" << indId
//	<< " smsData->dp=" << smsData->dp	<< " smsData->gb=" << smsData->gb
//	<< " smsData->alphaDB=" << smsData->alphaDB	<< " smsData->betaDB=" << smsData->betaDB
//	<< endl;
#endif
return;
}

// Get phenotypic transfer by SMS traits
trfrSMSTraits Individual::getSMSTraits(void) {
#if RSDEBUG
//DEBUGLOG << "Individual::getSMSTraits(): indId=" << indId << " smsData=" << smsData
//	<< endl;
#endif
trfrSMSTraits s; s.dp = s.gb = s.alphaDB = 1.0; s.betaDB = 1;
if (smsData != 0) {
	s.dp = smsData->dp; s.gb = smsData->gb;
	s.alphaDB = smsData->alphaDB; s.betaDB = smsData->betaDB;
}
#if RSDEBUG
//DEBUGLOG << "Individual::getSMSTraits(): indId=" << indId
//	<< " s.dp=" << s.dp << " s.gb=" << s.gb
//	<< " s.alphaDB=" << s.alphaDB << " s.betaDB=" << s.betaDB
//	<< endl;
#endif
return s;
}

// Set phenotypic transfer by CRW traits
void Individual::setCRWTraits(Species *pSpecies,short CRWgenelocn,short nCRWtraits,
	bool sexdep) {
#if RSDEBUG
//DEBUGLOG << "Individual::setCRWTraits(): indId=" << indId
//	<< " CRWgenelocn=" << CRWgenelocn << " nCRWtraits=" << nCRWtraits << " sexdep=" << sexdep
//	<< endl;
#endif
trfrCRWTraits c; c.stepLength = c.rho = 0.0;          
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		if (sexdep) {
			c.stepLength = (float)pGenome->express(CRWgenelocn+sex,0,sex);
			c.rho = (float)pGenome->express(CRWgenelocn+2+sex,0,sex);
		}
		else {
			c.stepLength = (float)pGenome->express(CRWgenelocn,0,0);
			c.rho = (float)pGenome->express(CRWgenelocn+1,0,0);
		}
	}
	else {
		if (sexdep) {
			c.stepLength = (float)pGenome->express(pSpecies,CRWgenelocn+sex);
			c.rho = (float)pGenome->express(pSpecies,CRWgenelocn+2+sex);
		}
		else {
			c.stepLength = (float)pGenome->express(pSpecies,CRWgenelocn);
			c.rho = (float)pGenome->express(pSpecies,CRWgenelocn+1);
		}
	}
}
#if RSDEBUG
//DEBUGLOG << "Individual::setCRWTraits(): indId=" << indId
//	<< " c.stepLength=" << c.stepLength << " c.rho=" << c.rho
//	<< endl;
#endif

trfrCRWParams cparams;
if (sexdep) {
	cparams = pSpecies->getCRWParams(0,sex);
}
else {
	cparams = pSpecies->getCRWParams(0,0);
}
crw->stepL = (float)(c.stepLength*cparams.stepLScale + cparams.stepLgthMean);        
crw->rho   = (float)(c.rho*cparams.rhoScale + cparams.rhoMean);
#if RSDEBUG
//DEBUGLOG << "Individual::setCRWTraits(): indId=" << indId
//	<< " crw->stepL=" << crw->stepL	<< " crw->rho=" << crw->rho
//	<< endl;
#endif
if (crw->stepL < 1.0) crw->stepL = 1.0;
if (crw->rho < 0.0) crw->rho = 0.0;
if (crw->rho > 0.999) crw->rho = 0.999f;
#if RSDEBUG
//DEBUGLOG << "Individual::setCRWTraits(): indId=" << indId
//	<< " crw->stepL=" << crw->stepL	<< " crw->rho=" << crw->rho
//	<< endl;
#endif
return;
}

// Get phenotypic transfer by CRW traits
trfrCRWTraits Individual::getCRWTraits(void) {
#if RSDEBUG
//DEBUGLOG << "Individual::getCRWTraits(): indId=" << indId
//	<< endl;
#endif
trfrCRWTraits c; c.stepLength = c.rho = 0.0;
if (crw != 0) {
	c.stepLength = crw->stepL;
	c.rho = crw->rho;
}
#if RSDEBUG
//DEBUGLOG << "Individual::getCRWTraits(): indId=" << indId
//	<< " c.stepLength=" << c.stepLength << " c.rho=" << c.rho
//	<< endl;
#endif

return c;

}

// Set phenotypic settlement traits
void Individual::setSettTraits(Species *pSpecies,short settgenelocn,short nsetttraits,
	bool sexdep) {
#if RSDEBUG
//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId << " sex=" << sex
//	<< " settgenelocn=" << settgenelocn << " nsetttraits=" << nsetttraits << " sexdep=" << sexdep
//	<< endl;
#endif
//simParams sim = paramsSim->getSim();
settleTraits s; s.s0 = s.alpha = s.beta = 0.0;            
if (pGenome != 0) {
	if (pSpecies->has1ChromPerTrait()) {
		if (sexdep) {
			s.s0    = (float)pGenome->express(settgenelocn+3*sex,0,0);
			s.alpha = (float)pGenome->express(settgenelocn+3*sex+1,0,0);
			s.beta  = (float)pGenome->express(settgenelocn+3*sex+2,0,0);
		}
		else {
			s.s0    = (float)pGenome->express(settgenelocn,0,0);
			s.alpha = (float)pGenome->express(settgenelocn+1,0,0);
			s.beta  = (float)pGenome->express(settgenelocn+2,0,0);
		}
	}
	else {
		if (sexdep) {
			s.s0    = (float)pGenome->express(pSpecies,settgenelocn+3*sex);
			s.alpha = (float)pGenome->express(pSpecies,settgenelocn+3*sex+1);
			s.beta  = (float)pGenome->express(pSpecies,settgenelocn+3*sex+2);
		}
		else {
			s.s0    = (float)pGenome->express(pSpecies,settgenelocn);
			s.alpha = (float)pGenome->express(pSpecies,settgenelocn+1);
			s.beta  = (float)pGenome->express(pSpecies,settgenelocn+2);
		}

  }
}
#if RSDEBUG
//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
//	<< " s.s0=" << s.s0 << " s.alpha=" << s.alpha << " s.beta=" << s.beta
//	<< endl;
#endif

settParams sparams;
if (sexdep) {
	sparams = pSpecies->getSettParams(0,sex);
}
else {
	sparams = pSpecies->getSettParams(0,0);     
}
#if RSDEBUG
//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
//	<< " sparams.s0Mean=" << sparams.s0Mean << " sparams.s0SD=" << sparams.s0SD 
//	<< " sparams.s0Scale=" << sparams.s0Scale
//	<< endl;
#endif
setttraits = new settleTraits;
setttraits->s0    = (float)(s.s0*sparams.s0Scale + sparams.s0Mean);
setttraits->alpha = (float)(s.alpha*sparams.alphaSScale + sparams.alphaSMean);
setttraits->beta  = (float)(s.beta*sparams.betaSScale + sparams.betaSMean);
#if RSDEBUG
//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
//	<< " setttraits->s0=" << setttraits->s0
//	<< " setttraits->alpha=" << setttraits->alpha << " setttraits->beta=" << setttraits->beta
//	<< endl;
#endif
if (setttraits->s0 < 0.0) setttraits->s0 = 0.0;
if (setttraits->s0 > 1.0) setttraits->s0 = 1.0;
#if RSDEBUG
//DEBUGLOG << "Individual::setSettTraits(): indId=" << indId
//	<< " setttraits->s0=" << setttraits->s0
//	<< " setttraits->alpha=" << setttraits->alpha << " setttraits->beta=" << setttraits->beta
//	<< endl;
#endif
return;
}

// Get phenotypic settlement traits
settleTraits Individual::getSettTraits(void) {
#if RSDEBUG
//DEBUGLOG << "Individual::getSettTraits(): indId=" << indId
//	<< endl;
#endif
settleTraits s; s.s0 = s.alpha = s.beta = 0.0;
if (setttraits != 0) {
	s.s0    = setttraits->s0;
	s.alpha = setttraits->alpha;
	s.beta  = setttraits->beta;
}
#if RSDEBUG
//DEBUGLOG << "Individual::getSettTraits(): indId=" << indId
//	<< " s.s0=" << s.s0 << " s.alpha=" << s.alpha << " s.beta=" << s.beta
//	<< endl;
#endif

return s;
}

/*
locus Individual::getAlleles(int g) {
locus l; l.allele[0] = l.allele[1] = 0.0;
if (pGenome != 0) l = pGenome->getAlleles(g);
return l;
}
*/

void Individual::setStatus(short s) {
if (s >= 0 && s <= 9) status = s;
status = s;
}
#if SEASONAL
#if PARTMIGRN
void Individual::setMigrnStatus(short m) {
if (m >= 0 && m <= 6) migrnstatus = m;
}
void Individual::setPrevPatch(Patch *pPatch) {
pPrevPatch = pPatch;
}
#endif // PARTMIGRN
#endif // SEASONAL

#if SEASONAL
#if PARTMIGRN

//void Individual::setNpatches(const short n) {
//if (n > 0) npatches = n;
//}

void Individual::addPatch(patchlist p) {
int size = (int)patches.size();
#if RSDEBUG
//DEBUGLOG << "Individual::addPatch(): indId=" << indId
//	<< " size=" << size << " Patch=" << p.pPatch->getPatchNum()
//	<< " p.season=" << p.season << " p.breeding=" << p.breeding
//	<< endl;
#endif
if (migrnstatus >= 4 && p.pPatch == pNatalPatch ) {
	// a disperser may not remember its natal patch
	return;
}
bool present = false;
//bool fixed;
if (p.breeding) {
	p.fixed = true;
}
else {
	if (migrnstatus == 3 || migrnstatus == 6) p.fixed = false;
	else p.fixed = true;
}
for (int i = 0; i < size; i++) {                
	if (patches[i].pPatch == p.pPatch && patches[i].breeding == p.breeding) {
		present = true;
//		patches[i].fixed = p.fixed;
		patches[i].fixed = p.fixed;
	}
	if (patches[i].pPatch != p.pPatch && patches[i].breeding == p.breeding) {
		// un-fix another patch previous fixed for the same breeding status
		patches[i].fixed = false;     
	}
//	if (present) break;
}
if (!present) patches.push_back(p);
size = (int)patches.size();
#if RSDEBUG
//DEBUGLOG << "Individual::addPatch(): indId=" << indId
//	<< " size=" << size 
//	<< endl;
#endif
}

patchlist Individual::getPatch(const int n) {
patchlist p;      
int size = (int)patches.size();
if (n >= 0 && n < size) {
	p = patches[n];
}
else {
	p.pPatch = 0; p.season = 0; p.breeding = p.fixed = false;
}
return p;
}

void Individual::setGoal(const locn loc,const short	gtype,const bool breeding) {
#if RSDEBUG
//DEBUGLOG << "Individual::setGoal(): indId=" << indId
//	<< " gtype=" << gtype 
//	<< " breeding=" << breeding 
//	<< " smsData=" << smsData 
//	<< endl;
#endif
if (smsData != 0) {
	if (gtype >= 0 && gtype <= 2) {
		if (gtype == 1) {
			// set goal to fixed patch (if any) appropriate to season
			smsData->goalType = 0;
			int size = (int)patches.size();
//			patchlist p; p.pPatch = 0; p.fixed = false;
			for (int i = 0; i < size; i++) { 
#if RSDEBUG
//DEBUGLOG << "Individual::setGoal(): indId=" << indId
//	<< " size=" << size 
//	<< " i=" << i 
//	<< " patches[i].Patch=" << patches[i].pPatch->getPatchNum() 
//	<< " patches[i].fixed=" << patches[i].fixed 
//	<< " patches[i].breeding=" << patches[i].breeding 
//	<< endl;
#endif
				if (patches[i].fixed && patches[i].breeding == breeding) {   
					smsData->goalType = 1; 
					smsData->goal = patches[i].pPatch->getRandomCell()->getLocn();
#if RSDEBUG
//DEBUGLOG << "Individual::setGoal(): indId=" << indId
//	<< " goal.x=" << smsData->goal.x 
//	<< " goal.y=" << smsData->goal.y 
//	<< endl;
#endif
					i = size;
				}
			}
		}
		else {
			smsData->goalType = gtype; smsData->goal = loc;
		}
	}
}
}

#endif // PARTMIGRN 
#endif // SEASONAL

void Individual::developing(void) {
isDeveloping = true;
}

void Individual::develop(void) {
stage++; isDeveloping = false;
}

void Individual::ageIncrement(short maxage) {
if (status < 6) { // alive
	age++;
	if (age > maxage) status = 9;			// exceeds max. age - dies
	else {
		if (path != 0) path->year = 0;	// reset annual step count for movement models
		if (status == 3) // waiting to continue dispersal
			status = 1;
	}
}
}

void Individual::incFallow(void) { fallow++; }

void Individual::resetFallow(void) { fallow = 0; }

//---------------------------------------------------------------------------
// Move to a specified neighbouring cell
void Individual::moveto(Cell *newCell) {
// check that location is indeed a neighbour of the current cell
locn currloc = pCurrCell->getLocn();
locn newloc = newCell->getLocn();
double d = sqrt(((double)currloc.x-(double)newloc.x)*((double)currloc.x-(double)newloc.x)
	+ ((double)currloc.y-(double)newloc.y)*((double)currloc.y-(double)newloc.y));
if (d >= 1.0 && d < 1.5) { // ok
	pCurrCell = newCell; status = 5;
}
}

#if GROUPDISP
//---------------------------------------------------------------------------
// Move to any specified cell
void Individual::moveTo(Cell *newCell) {
if (newCell != 0) {
	pCurrCell = newCell;
}
}
#endif

//---------------------------------------------------------------------------
// Move to a new cell by sampling a dispersal distance from a single or double
// negative exponential kernel
#if RS_CONTAIN
// or the 2Dt kernel or the WALD kernel
#endif // RS_CONTAIN 
// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
#if SEASONAL
int Individual::moveKernel(Landscape *pLandscape,Species *pSpecies,
	const short repType,const short nextseason,const bool absorbing)
#else
int Individual::moveKernel(Landscape *pLandscape,Species *pSpecies,
	const short repType,const bool absorbing)
#endif // SEASONAL 
{

intptr patch;
int patchNum = 0;
int newX = 0,newY = 0;
int dispersing = 1;
double xrand,yrand,meandist,dist,r1,rndangle,nx,ny;
float localK;
trfrKernTraits kern;      
Cell* pCell;
Patch* pPatch;
locn loc = pCurrCell->getLocn();

landData land = pLandscape->getLandData();

bool usefullkernel = pSpecies->useFullKernel();
trfrRules trfr = pSpecies->getTrfr();
settleRules sett = pSpecies->getSettRules(stage,sex);

pCell = NULL; 
pPatch = NULL;

if (trfr.indVar) { // get individual's kernel parameters
	kern.meanDist1 = kern.meanDist2 = kern.probKern1 = 0.0;
	if (pGenome != 0) {
		kern.meanDist1 = kerntraits->meanDist1;
#if RS_CONTAIN
		if (trfr.kernType == 1) 
#else
		if (trfr.twinKern) 
#endif // RS_CONTAIN 
		{
			kern.meanDist2 = kerntraits->meanDist2;
			kern.probKern1 = kerntraits->probKern1;
		}
	}
}
else { // get kernel parameters for the species
	if (trfr.sexDep) {
		if (trfr.stgDep) {
			kern = pSpecies->getKernTraits(stage,sex);
		}
		else {
			kern = pSpecies->getKernTraits(0,sex);
		}
	}
	else {
		if (trfr.stgDep) {
			kern = pSpecies->getKernTraits(stage,0);
		}
		else {
			kern = pSpecies->getKernTraits(0,0);
		}
	}
}
#if RSDEBUG
//Patch *startPatch = (Patch*)startpatch;
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " x=" << loc.x << " y=" << loc.y
////	<< " natalPatch = " << natalPatch
////	<< " startpatch = " << startpatch << " patchNum = " << startPatch->getPatchNum()
//	<< " kern.meanDist1=" << kern.meanDist1;
//if (trfr.twinKern) {
//	DEBUGLOG << " meanDist2=" << kern.meanDist2 << " probKern1=" << kern.probKern1;
//}
//DEBUGLOG << endl;
#endif

// scale the appropriate kernel mean to the cell size
#if RS_CONTAIN
if (trfr.kernType == 1) 
#else
if (trfr.twinKern) 
#endif // RS_CONTAIN 
{
	if (pRandom->Bernoulli(kern.probKern1))
		meandist = kern.meanDist1 / (float)land.resol;
	else
		meandist = kern.meanDist2 / (float)land.resol;
}
else
	meandist = kern.meanDist1 / (float)land.resol;
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " meandist=" << meandist << endl;
#endif
// scaled mean may not be less than 1 unless emigration derives from the kernel
// (i.e. the 'use full kernel' option is applied)
if (!usefullkernel && meandist < 1.0) meandist = 1.0;
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " meandist=" << meandist << endl;
#endif

#if RSDEBUG
//Patch *startPatch = (Patch*)startpatch;
//DEBUGLOG << "Individual::moveKernel(): indId = " << indId << " x = " << x << " y = " << y
//	<< " natalPatch = " << natalPatch
////	<< " startpatch = " << startpatch << " patchNum = " << startPatch->getPatchNum()
//	<< " meanDist1 = " << kern.meanDist1;
//if (trfr.twinKern) {
//	DEBUGLOG << " probKern1 = " << kern.probKern1 << " meanDist2 = " << kern.meanDist2;
//}
//DEBUGLOG << " meandist = " << meandist << endl;
#endif

#if RS_CONTAIN

// simple 2Dt model (no environmental effects) 
trfr2Dt t2 = pSpecies->getTrfr2Dt();
double propkern1 = t2.propKernel1;
double p,u,f,f0;
bool reject;

// WALD model
trfrWald w = pSpecies->getTrfrWald(); 
double hr = pSpecies->getTrfrHr(motherstage);
double mu = w.meanU * hr / w.vt;
double gamma = (w.meanU * hr * hr) / (2.0 * w.kappa * w.hc * w.sigma_w);   

// select kernel to sample for this individual (if 2Dt) and
// find a suitable maximum x-value for the range of distances to sample
int maxdim = max(land.dimX,land.dimY) * land.resol;  
double maxx = (double)maxdim;
bool ok = false;
double fdim;
if (trfr.kernType == 2) {
	if (pRandom->Bernoulli(propkern1)) { // sample from kernel 1
		u = exp(t2.u0Kernel1); p = exp(t2.p0Kernel1); f0 = p / (PI * u);
		while (!ok) {
			fdim = p / (PI * u * pow((1.0 + (maxx*maxx/u)),(p + 1.0)));
			if (fdim >= f0/1000.0) ok = true; else maxx /= 1.25;
		}
	}
	else { // sample from kernel 2
		u = exp(t2.u0Kernel2); p = exp(t2.p0Kernel2); f0 = p / (PI * u);
		while (!ok) {
			fdim = p / (PI * u * pow((1.0 + (maxx*maxx/u)),(p + 1.0)));
			if (fdim >= f0/1000.0) ok = true; else maxx /= 1.25;
		}
	}
}
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " x=" << loc.x << " y=" << loc.y
//	<< " maxx=" << maxx 
//	<< " u=" << u << " p=" << p << " f0=" << f0 
//	<< endl;
#endif
if (trfr.kernType == 3) {
	f0 = 1.0;
	while (!ok) {
		fdim = sqrt(gamma / (2.0 * PI * maxx * maxx * maxx)) 
						* exp(-1.0 * gamma * (maxx-mu) * (maxx-mu) / (2.0 * maxx * mu * mu) );
		if (fdim >= f0/10000.0) ok = true; else maxx /= 1.25;
	}	
}
#if RSDEBUG
DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " x=" << loc.x << " y=" << loc.y
	<< " stage=" << stage << " hr=" << hr << " mu=" << mu << " gamma=" << gamma 
	<< " maxx=" << maxx 
	<< endl;
#endif

#endif // RS_CONTAIN 

int loopsteps = 0; // new counter to prevent infinite loop added 14/8/15
do {
	do {
		do {
			// randomise the cell within the patch, provided that the individual is still in
			// its natal cell (i.e. not waiting in the matrix)
			// this is because, if the patch is very large, the individual is near the centre
			// and the (single) kernel mean is (not much more than) the cell size, an infinite
			// loop could otherwise result, as the individual never reaches the patch edge
			// (in a cell-based model, this has no effect, other than as a processing overhead)
			if (status == 1) {
				pCell = pNatalPatch->getRandomCell();
				if (pCell != 0) {
					loc = pCell->getLocn();
				}
			}
			// randomise the position of the individual inside the cell
			xrand = (double)loc.x + pRandom->Random()*0.999;
			yrand = (double)loc.y + pRandom->Random()*0.999;

#if RS_CONTAIN

			switch (trfr.kernType) {
				
			case 0: // single negative exponential
			case 1: // single negative exponential
				r1 = 0.0000001 + pRandom->Random()*(1.0-0.0000001);
				dist = (-1.0*meandist)*log(r1);  // for LINUX_CLUSTER
				break;
				
			case 2: // 2Dt
				
			// sample distance from 2Dt kernel by method of REJECTION SAMPLING 

			// NOTE: sampling must be in real-world co-ordinates (not cell co-ordinates)
			// as kernel units are metres

			reject = true;
			while (reject) {
				// sample a random distance along the x-axis
				dist = pRandom->Random() * maxx;
				// sample a random y-axis variate between zero and max. possible 
				r1 = pRandom->Random() * f0;
				// calculate value of kernel at dist;
				f = p / (PI * u * pow((1.0 + (dist*dist/u)),(p+1.0)));
				if (r1 <= f) reject = false;
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " dist=" << dist << " r1=" << r1
//	<< " f=" << f << " reject=" << reject << endl;
#endif
				}
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " SAMPLED dist=" << dist 
//	<< endl;
#endif
			// convert sampled distance to cell co-ordinates 
			dist /= (double)land.resol;      
//			rndangle = pRandom->Random() * 2.0 * PI;
//			nx = (xrand + dist * cos(rndangle)) / land.resol;
//			ny = (yrand + dist * sin(rndangle)) / land.resol;

				break;
				
			case 3: // Wald

//			dist = 2 * land.resol;
				
			// sample distance from 2Dt kernel by method of REJECTION SAMPLING 

			// NOTE: sampling must be in real-world co-ordinates (not cell co-ordinates)
			// as kernel units are metres

			reject = true;
			while (reject) {
				// sample a random distance along the x-axis
				dist = pRandom->Random() * maxx;
				// sample a random y-axis variate between zero and max. possible 
				r1 = pRandom->Random() * f0;
				// calculate value of kernel at dist;
				f = sqrt(gamma / (2.0 * PI * dist * dist * dist)) 
						* exp(-1.0 * gamma * (dist-mu) * (dist-mu) / (2.0 * dist * mu * mu) );
				if (r1 <= f) reject = false;
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " dist=" << dist << " r1=" << r1
//	<< " f=" << f << " reject=" << reject << endl;
#endif
				}
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " SAMPLED dist=" << dist 
//	<< endl;
#endif
			// convert sampled distance to cell co-ordinates 
			dist /= (double)land.resol;      

				break;
				
			}

#else
			
			r1 = 0.0000001 + pRandom->Random()*(1.0-0.0000001);
//			dist = (-1.0*meandist)*std::log(r1);
			dist = (-1.0*meandist)*log(r1);  // for LINUX_CLUSTER
			
#endif // RS_CONTAIN 

#if RS_CONTAIN
			rndangle = pRandom->Normal(w.meanDirn,w.sdDirn);
			if (rndangle < 0.0) rndangle += 360.0;
			if (rndangle >= 360.0) rndangle -= 360.0;
			rndangle *= 2.0 * PI / 360.00;
#if RSDEBUG
DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " status=" << status
	<< " meanDirn=" << w.meanDirn << " sdDirn=" << w.sdDirn << " rndangle=" << rndangle
	<< " loopsteps=" << loopsteps 
	<< endl;
#endif
#else
			rndangle = pRandom->Random() * 2.0 * PI;
#endif // RS_CONTAIN 
			nx = xrand + dist * sin(rndangle);
			ny = yrand + dist * cos(rndangle);
			if (nx < 0.0) newX = -1; else newX = (int)nx;
			if (ny < 0.0) newY = -1; else newY = (int)ny;
#if RSDEBUG
			if (path != 0) (path->year)++;
#endif
			loopsteps++;
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " status=" << status
//	<< " loopsteps=" << loopsteps << " newX=" << newX << " newY=" << newY
//	<< " loc.x=" << loc.x << " loc.y=" << loc.y
//	<< endl;
#endif
		} while (loopsteps < 1000 &&
				((!absorbing && (newX < land.minX || newX > land.maxX
												|| newY < land.minY || newY > land.maxY))
				 || (!usefullkernel && newX == loc.x && newY == loc.y))
				);
		if (loopsteps < 1000) {
			if (newX < land.minX || newX > land.maxX
					|| newY < land.minY || newY > land.maxY) { // beyond absorbing boundary
				pCell = 0;
				patch = 0;
				patchNum = -1;
			}
			else {
				pCell = pLandscape->findCell(newX,newY);
				if (pCell == 0) { // no-data cell
					patch = 0;
					patchNum = -1;
				}
				else {
					patch = pCell->getPatch();
					if (patch == 0) { // matrix
						pPatch = 0;
						patchNum = 0;
					}
					else {
						pPatch = (Patch*)patch;
						patchNum = pPatch->getPatchNum();
					}
				}
			}
		}
		else {
			patch = 0;
			patchNum = -1;
		}
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId << " status=" << status
//	<< " loopsteps=" << loopsteps << " newX=" << newX << " newY=" << newY
//	<< " pCell=" << pCell << " patch=" << patch << " patchNum=" << patchNum
//	<< endl;
#endif
	} while (!absorbing && patchNum < 0 && loopsteps < 1000); 			 // in a no-data region
}
while (!usefullkernel && pPatch == pNatalPatch && loopsteps < 1000); 	// still in the original (natal) patch

if (loopsteps < 1000) {
	if (pCell == 0) { // beyond absorbing boundary or in no-data cell
		pCurrCell = 0;
		status = 6;
		dispersing = 0;
	}
	else {
		pCurrCell = pCell;
		if (pPatch == 0) localK = 0.0; // matrix
#if SEASONAL
		else localK = pPatch->getK(nextseason);
#else
		else localK = pPatch->getK();
#endif // SEASONAL 
		if (patchNum > 0 && localK > 0.0) { // found a new patch
			status = 2; // record as potential settler
		}
		else {
			dispersing = 0;
			// can wait in matrix if population is stage structured ...
			if (pSpecies->stageStructured()) {
				// ... and wait option is applied ...
				if (sett.wait) { // ... it is
					status = 3; // waiting
				}
				else // ... it is not
					status = 6; // dies (unless there is a suitable neighbouring cell)
			}
			else
				status = 6; // dies (unless there is a suitable neighbouring cell)
		}
	}
}
else {
	status = 6;
	dispersing = 0;
}
#if RSDEBUG
//DEBUGLOG << "Individual::moveKernel(): indId=" << indId
//	<< " newX=" << newX << " newY=" << newY
//	<< " patch=" << patch
//	<< " patchNum=" << patchNum << " status=" << status;
//DEBUGLOG << endl;
#endif

// apply dispersal-related mortality, which may be distance-dependent
dist *= (float)land.resol; // re-scale distance moved to landscape scale
if (status < 7) {
	double dispmort;
	trfrMortParams mort = pSpecies->getMortParams();
	if (trfr.distMort) {
		dispmort = 1.0 / (1.0 + exp(-(dist - mort.mortBeta)*mort.mortAlpha));  
	}
	else {
		dispmort = mort.fixedMort;
	}
	if (pRandom->Bernoulli(dispmort)) {
		status = 7; // dies
		dispersing = 0;
	}
}

return dispersing;
}

//---------------------------------------------------------------------------
// Make a single movement step according to a mechanistic movement model
// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
#if SEASONAL
int Individual::moveStep(Landscape *pLandscape,Species *pSpecies,
	const short landIx,const short nextseason,const bool absorbing)
#else
int Individual::moveStep(Landscape *pLandscape,Species *pSpecies,
	const short landIx,const bool absorbing)
#endif // SEASONAL 
{

if (status != 1) return 0; // not currently dispersing

intptr patch;
int patchNum;
#if VCL
int oldX, oldY;
#endif
int newX,newY;
locn loc;
int dispersing = 1;
double xcnew,ycnew;
double angle;
double mortprob,rho,steplen;
movedata move;
Patch* pPatch = 0;
bool absorbed = false;
//int popsize;

landData land = pLandscape->getLandData();
#if VCL
simView v = paramsSim->getViews();
#endif
simParams sim = paramsSim->getSim();

trfrRules trfr = pSpecies->getTrfr();
trfrCRWTraits movt = pSpecies->getCRWTraits();
settleSteps settsteps = pSpecies->getSteps(stage,sex);

patch = pCurrCell->getPatch();
#if RSDEBUG
//DEBUGLOG << "Individual::moveStep() AAAA: indId=" << indId
//	<< " pCurrCell=" << pCurrCell << " patch=" << patch
//	<< endl;
#endif

if (patch == 0) { // matrix
	pPatch = 0;
	patchNum = 0;
}
else {
	pPatch = (Patch*)patch;
	patchNum = pPatch->getPatchNum();
}
// apply step-dependent mortality risk ...
#if TEMPMORT
int h;
switch (trfr.smType) {
case 0: // constant
	mortprob = movt.stepMort;
	break;
case 1: // habitat-dependent
	h = pCurrCell->getHabIndex(landIx);
	if (h < 0) { // no-data cell - should not occur, but if it does, individual dies
		mortprob = 1.0;
	}
	else mortprob = pSpecies->getHabMort(h);
	break;
case 2: // temporally variable
	mortprob = pSpecies->getMortality();
	break;
} 
#else
if (trfr.habMort) 
{ // habitat-dependent
	int h = pCurrCell->getHabIndex(landIx);
	if (h < 0) { // no-data cell - should not occur, but if it does, individual dies
		mortprob = 1.0;
	}
	else mortprob = pSpecies->getHabMort(h);
#if RSDEBUG
//locn temploc = pCurrCell->getLocn();
//DEBUGLOG << "Individual::moveStep(): x=" << temploc.x << " y=" << temploc.x
//	<< " landIx=" << landIx << " h=" << h << " mortprob=" << mortprob
//	<< endl;
#endif
}
else mortprob = movt.stepMort;
#endif // TEMPMORT 
// ... unless individual has not yet left natal patch in emigration year
if (pPatch == pNatalPatch && path->out == 0 && path->year == path->total) {
	mortprob = 0.0;
}
#if RSDEBUG
locn loc0,loc1,loc2;
//loc0 = pCurrCell->getLocn();
//DEBUGLOG << "Individual::moveStep() BBBB: indId=" << indId << " status=" << status
//	<< " path->year=" << path->year << " path->out=" << path->out
//	<< " settleStatus=" << path->settleStatus
//	<< " x=" << loc0.x << " y=" << loc0.y
////	<< " patch=" << patch
//	<< " pPatch=" << pPatch
//	<< " patchNum=" << patchNum;
////	<< " natalPatch=" << natalPatch;
////if (crw != 0) {
////	DEBUGLOG << " xc=" << crw->xc << " yc=" << crw->yc;
////	DEBUGLOG << " rho=" << movt.rho << " stepLength=" << movt.stepLength;
////}
//DEBUGLOG << endl;
#endif
if (pRandom->Bernoulli(mortprob)) { // individual dies
	status = 7;
	dispersing = 0;
}
else { // take a step
	(path->year)++;
#if SEASONAL
	(path->season)++;
#endif
	(path->total)++;
//	if (pPatch != pNatalPatch || path->out > 0) (path->out)++;
	if (patch == 0 || pPatch == 0 || patchNum == 0) { // not in a patch
		if (path != 0) path->settleStatus = 0; // reset path settlement status
		(path->out)++;
	}
	loc = pCurrCell->getLocn();
	newX = loc.x; newY = loc.y;
#if VCL
	oldX = loc.x; oldY = loc.y;
#endif


	switch (trfr.moveType) {

	case 1: // SMS
#if RSDEBUG                   
//loc1 = pCurrCell->getLocn();      
//DEBUGLOG << "Individual::moveStep() FFFF: indId=" << indId << " status=" << status
////	<< " path->year=" << path->year
//	<< " path->season=" << path->season
//	<< " x=" << loc1.x << " y=" << loc1.y
//	<< " smsData->goalType=" << smsData->goalType
//	<< " goal.x=" << smsData->goal.x
//	<< " goal.y=" << smsData->goal.y
//	<< endl;
#endif
#if PARTMIGRN
		move = smsMove(pLandscape,pSpecies,landIx,pPatch==pPrevPatch,trfr.indVar,absorbing);
#else
		move = smsMove(pLandscape,pSpecies,landIx,pPatch==pNatalPatch,trfr.indVar,absorbing);
#endif  // PARTMIGRN 
#if RSDEBUG
//DEBUGLOG << "Individual::moveStep() GGGG: indId=" << indId << " status=" << status
//	<< " move.dist=" << move.dist
//	<< endl;
#endif
		if (move.dist < 0.0) {
			// either INTERNAL ERROR CONDITION - INDIVIDUAL IS IN NO-DATA SQUARE
			// or individual has crossed absorbing boundary ...
			// ... individual dies
			status = 6;
			dispersing = 0;
		}
		else {
#if RSDEBUG
//loc1 = pCurrCell->getLocn();
//DEBUGLOG << "Individual::moveStep() HHHH: indId=" << indId << " status=" << status
//	<< " path->year=" << path->year
//	<< " x=" << loc1.x << " y=" << loc1.y
////	<< " smsData = " << smsData
//	<< endl;
#endif

		// WOULD IT BE MORE EFFICIENT FOR smsMove TO RETURN A POINTER TO THE NEW CELL? ...

			patch = pCurrCell->getPatch();
			//int patchnum;
			if (patch == 0) {
				pPatch = 0; 
				//patchnum = 0;
			}
			else {
				pPatch = (Patch*)patch; 
				//patchnum = pPatch->getPatchNum();
			}
			if (sim.saveVisits && pPatch != pNatalPatch) {
				pCurrCell->incrVisits();
			}
#if VCL
			if (v.viewPaths) {
				if ((Patch*)patch != pNatalPatch) {
					loc = pCurrCell->getLocn();
					drawMove((float)oldX+0.5,(float)oldY+0.5,(float)loc.x+0.5,(float)loc.y+0.5);
				}
			}
#endif
#if RS_CONTAIN
			if (status < 6) {
				DamageLocn *pDamageLocn = pCurrCell->getDamage();
				if (pDamageLocn != 0) pDamageLocn->updateTraversalDamage();
			}
#endif // RS_CONTAIN 
		}
		break;

	case 2: // CRW
		if (trfr.indVar) {
			if (crw != 0) {
				movt.stepLength = crw->stepL;
				movt.rho        = crw->rho;
			}
		}

		steplen = movt.stepLength; if (steplen < 0.2*land.resol) steplen = 0.2*land.resol;
		rho = movt.rho; if (rho > 0.99) rho = 0.99;
		if (pPatch == pNatalPatch) {
			rho = 0.99; // to promote leaving natal patch
			path->out = 0;
		}
		if (movt.straigtenPath && path->settleStatus > 0) {
			// individual is in a patch and has already determined whether to settle
			rho = 0.99; // to promote leaving the patch
			path->out = 0;
		}
		int loopsteps = 0; // new counter to prevent infinite loop added 14/8/15
		do {
			do {
				// new direction
				if (newX < land.minX || newX > land.maxX || newY < land.minY || newY > land.maxY
				|| pCurrCell == 0) {
					// individual has tried to go out-of-bounds or into no-data area
					// allow random move to prevent repeated similar move
					angle = wrpcauchy(crw->prevdrn,0.0);
				}
				else
					angle = wrpcauchy(crw->prevdrn,rho);
				// new continuous cell coordinates
				xcnew = crw->xc + sin(angle) * steplen/(float)land.resol;
				ycnew = crw->yc + cos(angle) * steplen/(float)land.resol;
				if (xcnew < 0.0) newX = -1; else newX = (int)xcnew;
				if (ycnew < 0.0) newY = -1; else newY = (int)ycnew;
				loopsteps++;
#if RSDEBUG
//DEBUGLOG << "Individual::moveStep(): indId=" << indId
//	<< " xc=" << crw->xc << " yc=" << crw->yc << " pCurrCell=" << pCurrCell
//	<< " steps=" << path->year << " loopsteps=" << loopsteps
//	<< " steplen=" << steplen << " rho=" << rho << " angle=" << angle
//	<< " xcnew=" << xcnew << " ycnew=" << ycnew << " newX=" << newX << " newY=" << newY << endl;
#endif
			}
			while (!absorbing && loopsteps < 1000 &&
				(newX < land.minX || newX > land.maxX || newY < land.minY || newY > land.maxY));
			if (newX < land.minX || newX > land.maxX || newY < land.minY || newY > land.maxY)
				pCurrCell = 0;
			else
				pCurrCell = pLandscape->findCell(newX,newY);
			if (pCurrCell == 0) { // no-data cell or beyond absorbing boundary
				patch = 0;
				if (absorbing) absorbed = true;
			}
			else
				patch = pCurrCell->getPatch();
#if RSDEBUG
//DEBUGLOG << "Individual::moveStep(): indId=" << indId
//	<< " loopsteps=" << loopsteps << " absorbed=" << absorbed
//	<< " pCurrCell=" << pCurrCell << " patch=" << patch << endl;
#endif
		} while (!absorbing && pCurrCell == 0 && loopsteps < 1000);
#if VCL
		if (v.viewPaths) {
			if (newX >= land.minX && newX <= land.maxX && newY >= land.minY && newY <= land.maxY) {
				if (patch > 0) {
					if ((Patch*)patch != pNatalPatch) drawMove(crw->xc,crw->yc,xcnew,ycnew);
				}
				else {
					drawMove(crw->xc,crw->yc,xcnew,ycnew);
				}
			}
		}
#endif
		crw->prevdrn = (float)angle;
		crw->xc = (float)xcnew; crw->yc = (float)ycnew;
		if (absorbed) { // beyond absorbing boundary or in no-data square
			status = 6;
			dispersing = 0;
			pCurrCell = 0;
		}
		else {
			if (loopsteps >= 1000) { // unable to make a move
				// INTERNAL ERROR CONDITION - INDIVIDUAL IS IN NO-DATA SQUARE
				// NEED TO TAKE SOME FORM OF INFORMATIVE ACTION ...
				// ... individual dies as it cannot move
				status = 6;
				dispersing = 0;
				// current cell will be invalid (zero), so set back to previous cell
				pCurrCell = pPrevCell;
			}
		}
#if RSDEBUG
//DEBUGLOG << "Individual::moveStep(): indId=" << indId
//	<< " status=" << status
//	<< " pCurrCell=" << pCurrCell << " patch=" << patch << endl;
#endif
		break;

	} // end of switch (trfr.moveType)

#if RSDEBUG
//locn loc2;
//if (pCurrCell > 0) {
//	loc2 = pCurrCell->getLocn();
//}
//else {
//	loc2.x = -9999; loc2.y = -9999;
//}
//DEBUGLOG << "Individual::moveStep() ZZZZ: indId=" << indId
//	<< " status=" << status
//	<< " path->total=" << path->total
//	<< " x=" << loc2.x << " y=" << loc2.y
//	<< " patch=" << patch;
//if (patch > 0) {
//	pPatch = (Patch*)patch;
//	DEBUGLOG << " patchNum=" << pPatch->getPatchNum()
//		<< " getK()=" << pPatch->getK()
//		<< " popn=" << pPatch->getPopn((int)pSpecies);
//}
//	DEBUGLOG << endl;
#endif
	if (patch > 0  // not no-data area or matrix
	&&  path->total >= settsteps.minSteps) {
		pPatch = (Patch*)patch;
#if PARTMIGRN
		bool ok = false;
		if (pPatch != pPrevPatch
		&& (migrnstatus < 4 || pPatch != pNatalPatch)) {
			if (smsData->goalType == 1) {
				// aiming for a goal - check if goal patch has been reached
				Cell *pCell = pLandscape->findCell(smsData->goal.x,smsData->goal.y);
				if (pCell != 0) {
					Patch *pGoalPatch = (Patch*)pCell->getPatch();
					if (pPatch == pGoalPatch) ok = true;
				}
			}
			else {
				if (migrnstatus == 4) { // disperser resident
					// may settle in a patch only if it has non-zero K in ALL seasons
					// otherwise it is unsuitable
					ok = pPatch->suitableInAllSeasons();
				}
				else ok = true;
			}
		}
		if (ok)
#else
		if (pPatch != pNatalPatch) 
#endif  // PARTMIGRN 
		{
			// determine whether the new patch is potentially suitable
#if SEASONAL
			if (pPatch->getK(nextseason) > 0.0) 
#else
			if (pPatch->getK() > 0.0) 
#endif // SEASONAL 
			{ // patch is suitable
					status = 2;
			}
		}
	}
	if (status != 2 && status != 6) { // suitable patch not found, not already dead
#if PARTMIGRN
		if (path->season >= settsteps.maxStepsYr || path->season >= settsteps.maxSteps) {
			status = 6; 
		}
#else
		if (path->year >= settsteps.maxStepsYr) {
			status = 3; // waits until next year
		}
		if (path->total >= settsteps.maxSteps) {
			status = 6; // dies
			dispersing = 0;
		}
#endif
	}
} // end of single movement step

return dispersing;

}

//---------------------------------------------------------------------------

// Functions to implement the SMS algorithm

// Move to a neighbouring cell according to the SMS algorithm
movedata Individual::smsMove(Landscape *pLand,Species *pSpecies,
	const short landIx,const bool natalPatch,const bool indvar,const bool absorbing)
{

array3x3d nbr; 	// to hold weights/costs/probs of moving to neighbouring cells
array3x3d goal;	// to hold weights for moving towards a goal location
array3x3f hab;	// to hold weights for habitat (includes percep range)
int x2,y2; 			// x index from 0=W to 2=E, y index from 0=N to 2=S
int newX = 0,newY = 0;
Cell *pCell;
Cell *pNewCell = NULL;
double sum_nbrs = 0.0;
movedata move;
int cellcost,newcellcost;
locn current;

//if (write_out) {
//	out<<endl<<"ind = "<<get_id()<<" pr = "<<get_pr()<<" step = "<<step;
//	if (step == 0) out<<" start at "<<x<<" "<<y;
//	out<<endl;
//}
//pCell = pLand->findCell(x,y);
if (pCurrCell == 0)
{
// x,y is a NODATA square - this should not occur here
// return a negative distance to indicate an error
  move.dist = -69.0; move.cost = 0.0;
  return move;
}

#if RSDEBUG
//DEBUGLOG << "Individual::smsMove(): this=" << this << endl;
#endif

landData land = pLand->getLandData();
trfrSMSTraits movt = pSpecies->getSMSTraits();
current = pCurrCell->getLocn();

//get weights for directional persistence....
//if ((path->out > 0 && path->out < 10 && path->out < 2*movt.pr)
if ((path->out > 0 && path->out <= (movt.pr+1))
|| 	natalPatch
|| (movt.straigtenPath && path->settleStatus > 0)) {
#if PARTMIGRN
	if (smsData->goalType == 1) 
		// do not inflate DP so as not to swamp GB 
		nbr = getSimDir(current.x,current.y,movt.dp);		
	else {
#endif  // PARTMIGRN 
	// inflate directional persistence to promote leaving the patch
	if (indvar) nbr = getSimDir(current.x,current.y,10.0f*smsData->dp);
	else nbr = getSimDir(current.x,current.y,10.0f*movt.dp);
#if PARTMIGRN
	}
#endif  // PARTMIGRN 
}
else {
	if (indvar) nbr = getSimDir(current.x,current.y,smsData->dp);
	else nbr = getSimDir(current.x,current.y,movt.dp);
}
if (natalPatch || path->settleStatus > 0) path->out = 0;
//if (natalPatch) path->out = 0;
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove() 0000: nbr matrix" << endl;
//for (y2 = 2; y2 > -1; y2--) {
//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << nbr.cell[x2][y2] << " ";
//	DEBUGLOG << endl;
//}
#endif
//if (write_out) {
//	out<<endl<<"directional persistence weights:"<<endl;
//	for (y2 = 2; y2 > -1; y2--) {
//		for (x2 = 0; x2 < 3; x2++) out<<nbr.cell[x2][y2]<<" ";
//		out<<endl;
//	}
//}
//if (write_out2)
//	out2<<endl<<"ind = "<<get_id()<<" pr = "<<get_pr()<<" step = "<<step<<endl<<endl;

//get weights for goal bias....
double gb;
if (movt.goalType == 2) { // dispersal bias
	int nsteps = 0;
	if (path->year == path->total) { // first year of dispersal - use no. of steps outside natal patch
		nsteps = path->out;
	}
	else { // use total no. of steps
		nsteps = path->total;
	}
	if (indvar) {
		double exp_arg = -((double)nsteps - (double)smsData->betaDB)*(-smsData->alphaDB);
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove(): exp_arg=" << exp_arg;
#endif
		if (exp_arg > 100.0) exp_arg = 100.0; // to prevent exp() overflow error
		gb = 1.0 + (smsData->gb - 1.0)/(1.0 + exp(exp_arg));
	}
	else {
		double exp_arg = -((double)nsteps - (double)movt.betaDB)*(-movt.alphaDB);
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove(): exp_arg=" << exp_arg;
#endif
		if (exp_arg > 100.0) exp_arg = 100.0; // to prevent exp() overflow error
		gb = 1.0 + (movt.gb - 1.0)/(1.0 + exp(exp_arg));
	}
}
else gb = movt.gb;
#if PARTMIGRN
if (smsData->goalType == 1) {
	if ((path->out > 0 && path->out <= (movt.pr+1))
	|| 	natalPatch
	|| (movt.straigtenPath && path->settleStatus > 0)) {
		// inflate goal bias to promote leaving the patch in the bias direction
		goal = getGoalBias(current.x,current.y,smsData->goalType,10.0*gb);  
	}
	else {
		goal = getGoalBias(current.x,current.y,smsData->goalType,gb);  
	}
}
else
	goal = getGoalBias(current.x,current.y,smsData->goalType,gb);  
#else
goal = getGoalBias(current.x,current.y,movt.goalType,(float)gb);
#endif  // PARTMIGRN 
//if (write_out) {
//	out<<"goal bias weights:"<<endl;
//	for (y2 = 2; y2 > -1; y2--) {
//		for (x2 = 0; x2 < 3; x2++) out<<goal.cell[x2][y2]<<" ";
//		out<<endl;
//	}
//}
//if (write_out2)
//  out2<<endl<<"ind = "<<get_id()<<" pr = "<<get_pr()<<" step = "<<step<<endl<<endl;

// get habitat-dependent weights (mean effective costs, given perceptual range)
// first check if costs have already been calculated

hab = pCurrCell->getEffCosts();
#if RSDEBUG
//if (hab.cell[0][0] >= 0) {
//	DEBUGLOG << "Individual::smsMove() 1111: x=" << current.x << " y=" << current.y << endl;
//	for (y2 = 2; y2 > -1; y2--) {
//		for (x2 = 0; x2 < 3; x2++) DEBUGLOG << hab.cell[x2][y2] << " ";
//		DEBUGLOG << endl;
//	}
//}
#endif
//if (write_out) {
//  out<<"stored effective costs:"<<endl;
//	for (y2 = 2; y2 > -1; y2--) {
//    for (x2 = 0; x2 < 3; x2++) out<<hab.cell[x2][y2]<<" ";
//	  out<<endl;
//	}
//}
if (hab.cell[0][0] < 0.0) { // costs have not already been calculated
	hab = getHabMatrix(pLand,pSpecies,current.x,current.y,movt.pr,movt.prMethod,
					landIx,absorbing);
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove() 2222: " << endl;
//for (y2 = 2; y2 > -1; y2--) {
//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << hab.cell[x2][y2] << " ";
//	DEBUGLOG << endl;
//}
#endif
	pCurrCell->setEffCosts(hab);
}
else { // they have already been calculated - no action required
//	if (write_out) {
//		out<<"*** using previous effective costs ***"<<endl;
//	}
}
//if (write_out) {
//	out<<"mean effective costs:"<<endl;
//	for (y2 = 2; y2 > -1; y2--) {
//		for (x2 = 0; x2 < 3; x2++) {
//			out<<hab.cell[x2][y2]<<" ";
//		}
//		out<<endl;
//	}
//	out<<"weighted effective costs:"<<endl;
//}

// determine weighted effective cost for the 8 neighbours
// multiply directional persistence, goal bias and habitat habitat-dependent weights
for (y2 = 2; y2 > -1; y2--) {
  for (x2 = 0; x2 < 3; x2++) {
		if(x2 == 1 && y2 == 1) nbr.cell[x2][y2] = 0.0;
    else {
			if(x2 == 1 || y2 == 1) //not diagonal
        nbr.cell[x2][y2] = nbr.cell[x2][y2]*goal.cell[x2][y2]*hab.cell[x2][y2];
      else // diagonal
				nbr.cell[x2][y2] = (float)SQRT2*nbr.cell[x2][y2]*goal.cell[x2][y2]*hab.cell[x2][y2];
    }
//		if (write_out) {
//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"         ";
//		}
	}
//	if (write_out) out<<endl;
}
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove() 3333: " << endl;
//for (y2 = 2; y2 > -1; y2--) {
//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << nbr.cell[x2][y2] << " ";
//	DEBUGLOG << endl;
//}
#endif

// determine reciprocal of effective cost for the 8 neighbours
//if (write_out) out<<"reciprocal weighted effective costs:"<<endl;
for (y2 = 2; y2 > -1; y2--) {
  for (x2 = 0; x2 < 3; x2++) {
    if (nbr.cell[x2][y2] > 0.0) nbr.cell[x2][y2] = 1.0f/nbr.cell[x2][y2];
//		if (write_out) {
//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"         ";
//		}
	}
//	if (write_out) out<<endl;
}

// set any cells beyond the current landscape limits and any no-data cells
// to have zero probability
// increment total for re-scaling to sum to unity

#if RSDEBUG
//array3x3d temp;
//for (y2 = 2; y2 > -1; y2--) {
//	for (x2 = 0; x2 < 3; x2++) {
//		temp.cell[x2][y2] = nbr.cell[x2][y2];
//		if (current.x == 488 && current.y == 422) {
//			pCell = pLand->findCell((current.x+x2-1),(current.y+y2-1));
//			DEBUGLOG << "Individual::smsMove(): this=" << this
//				<< " IN THE PROBLEM CELL"
//				<< " y=" << current.y << " x=" << current.x
//				<< " y2=" << y2 << " x2=" << x2
//				<< " pCell=" << pCell;
//			if (pCell != 0) DEBUGLOG << " pCell->getCost=" << pCell->getCost();
//			DEBUGLOG << endl;
//		}
//	}
//}
#endif

for (y2 = 2; y2 > -1; y2--) {
	for (x2 = 0; x2 < 3; x2++) {
		if (!absorbing) {
			if ((current.y+y2-1) < land.minY || (current.y+y2-1) > land.maxY
			||  (current.x+x2-1) < land.minX || (current.x+x2-1) > land.maxX)
			// cell is beyond current landscape limits
				nbr.cell[x2][y2] = 0.0;
			else { // check if no-data cell
				pCell = pLand->findCell((current.x+x2-1),(current.y+y2-1));
				if (pCell == 0) nbr.cell[x2][y2] = 0.0; // no-data cell
			}
		}
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove(): this=" << this
//	<< " y=" << current.y << " x=" << current.x
//	<< " y2=" << y2 << " x2=" << x2
//	<< " pCell=" << pCell
//	<< endl;
#endif
//		if (write_out) {
//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"      ";
//		}
		sum_nbrs += nbr.cell[x2][y2];
	}
//	if (write_out) out<<endl;
}

// scale effective costs as probabilities summing to 1
//if (write_out) out<<"probabilities:"<<endl;
if (sum_nbrs > 0.0) { // should always be the case, but safest to check...
	for (y2 = 2; y2 > -1; y2--) {
		for (x2 = 0; x2 < 3; x2++) {
			nbr.cell[x2][y2] = nbr.cell[x2][y2]/(float)sum_nbrs;
//		if (write_out) {
//			out<<nbr.cell[x2][y2]<<" "; if (x2==1 && y2==1) out<<"      ";
//		}
		}
//	if (write_out) out<<endl;
	}
}
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove() 4444: " << endl;
//for (y2 = 2; y2 > -1; y2--) {
//	for (x2 = 0; x2 < 3; x2++) DEBUGLOG << nbr.cell[x2][y2] << " ";
//	DEBUGLOG << endl;
//}
#endif

// set up cell selection probabilities
//if (write_out) out<<"rnd = "<<rnd<<endl;
double cumulative[9];
int j = 0;
cumulative[0] = nbr.cell[0][0];
for (y2 = 0; y2 < 3; y2++) {
	for (x2 = 0; x2 < 3; x2++) {
		if (j != 0) cumulative[j] = cumulative[j-1] + nbr.cell[x2][y2];
		j++;
//    if (write_out) out<<"dx = "<<x2-1<<" dy = "<<y2-1<<" sum_rnd = "<<sum_rnd<<endl;
	}
}

// select direction at random based on cell selection probabilities
// landscape boundaries and no-data cells may be reflective or absorbing
cellcost = pCurrCell->getCost();
int loopsteps = 0; // new counter to prevent infinite loop added 14/8/15
do {
	do {
		double rnd = pRandom->Random();
		j = 0;
		for (y2 = 0; y2 < 3; y2++) {
			for (x2 = 0; x2 < 3; x2++) {
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove() 7777: rnd=" << rnd
//	<< " j=" << j	<< " cumulative[j]=" << cumulative[j]
//	<< endl;
#endif
				if (rnd < cumulative[j]) {
					newX = current.x + x2 - 1;
					newY = current.y + y2 - 1;
					if (x2 == 1 || y2 == 1) move.dist = (float)(land.resol);
					else move.dist = (float)(land.resol)*(float)SQRT2;
//			if (write_out) {
//				out<<"relative x and y "<<x2-1<<" "<<y2-1<<endl;
//				out<<"cost of move "<<move.cost<<endl;
//				out<<"move to: x = "<<x<<" y = "<<y;
//				if (oob) out<<" ***** OUT OF BOUNDS *****";
//				out<<endl;
//			}
					y2=999; x2=999; //to break out of x2 and y2 loops.
				}
				j++;
			}
		}
		loopsteps++;
	} while (loopsteps < 1000
				&& (!absorbing && (newX < land.minX || newX > land.maxX
												|| newY < land.minY || newY > land.maxY)));
	if (loopsteps >= 1000) pNewCell = 0;
	else {
		if (newX < land.minX || newX > land.maxX
		||  newY < land.minY || newY > land.maxY) {
			pNewCell = 0;
		}
		pNewCell = pLand->findCell(newX,newY);
	}
}
while (!absorbing && pNewCell == 0 && loopsteps < 1000); // no-data cell
#if RSDEBUG
//DEBUGLOG << "Individual::smsMove() 8888: pNewCell=" << pNewCell
//	<< " loopsteps=" << loopsteps
//	<< " current.x=" << current.x << " current.y=" << current.y
//	<< " newX=" << newX << " newY=" << newY
//	<< " land.minX=" << land.minX << " land.minY=" << land.minY
//	<< " land.maxX=" << land.maxX << " land.maxY=" << land.maxY
//	<< endl;
#endif
if (loopsteps >= 1000 || pNewCell == 0) {
	// unable to make a move or crossed absorbing boundary
	// flag individual to die
	move.dist = -123.0;
	if (pNewCell == 0) pCurrCell = pNewCell;
}
else {
	newcellcost = pNewCell->getCost();
	move.cost = move.dist*0.5f*((float)cellcost + (float)newcellcost);  
	// make the selected move
	if ((short)memory.size() == movt.memSize) {
		memory.pop(); // remove oldest memory element
	}
	memory.push(current); // record previous location in memory
	//if (write_out) out << "queue length is " << memory.size() << endl;
	pCurrCell = pNewCell;
}
return move;
}

// Weight neighbouring cells on basis of current movement direction
array3x3d Individual::getSimDir(const int x, const int y, const float dp) 
{

array3x3d d;
locn prev;
double theta;
int xx,yy;

//if (write_out) out<<"step 0"<<endl;
if (memory.empty())
{ // no previous movement, set matrix to unity
  for (xx = 0; xx < 3; xx++) {
    for (yy = 0; yy < 3; yy++) {
       d.cell[xx][yy] = 1;
    }
  }
}
else { // set up the matrix dependent on relationship of previous location to current
//  if (write_out) out<<"step 1"<<endl;
  d.cell[1][1]=0;
  prev = memory.front();
//  if (write_out) out<<"step 2"<<endl;
	if ((x-prev.x) == 0 && (y-prev.y) == 0) {
  // back to 'square 1' (first memory location) - use previous step drn only
		prev = memory.back();
//    if (write_out) out<<"step 3"<<endl;
//		if (write_out) out<<"*** using last step only: x,y = "<<prev.x<<","<<prev.y<<endl;
		if ((x-prev.x) == 0 && (y-prev.y) == 0) { // STILL HAVE A PROBLEM!
      for (xx = 0; xx < 3; xx++) {
        for (yy = 0; yy < 3; yy++) {
          d.cell[xx][yy] = 1.0;
        }
			}
//      if (write_out) out<<"step 4"<<endl;
      return d;
    }
  }
  else {
//    if (write_out) out<<"step 5"<<endl;
  }
//  if (write_out) out<<"step 6"<<endl;
	theta = atan2(((double)x-(double)prev.x),((double)y-(double)prev.y));
//  if (write_out) out<<"prev.x,prev.y: "<<prev.x<<","<<prev.y<<" theta: "<<theta<<endl;
	d = calcWeightings(dp,(float)theta);

}
return d;
}

// Weight neighbouring cells on basis of goal bias
//array3x3d Individual::getGoalBias(const int x,const int y,
//	const int goaltype,const float gb)
array3x3d Individual::getGoalBias(const int x, const int y, 
		const int goaltype, const float gb)
{

array3x3d d;
double theta;
int xx,yy;

if (goaltype == 0) { // no goal set
  for (xx = 0; xx < 3; xx++) {
    for (yy = 0; yy < 3; yy++) {
      d.cell[xx][yy] = 1.0;
    }
  }
}
else {
  d.cell[1][1]=0;
	if ((x - smsData->goal.x) == 0 && (y - smsData->goal.y) == 0) {
		// at goal, set matrix to unity
//    if (write_out) out<<"*** at goal: x,y = "<<goalx<<","<<goaly<<endl;
		for (xx = 0; xx < 3; xx++) {
			for (yy = 0; yy < 3; yy++) {
				d.cell[xx][yy] = 1.0;
			}
		}
		return d;
	}
	if (goaltype == 1) {
#if PARTMIGRN
		theta = atan2((double)(smsData->goal.x - x),(double)(smsData->goal.y - y));
#else
		// TEMPORARY CODE - GOAL TYPE 1 NOT YET IMPLEMENTED, AS WE HAVE NO MEANS OF
		// CAPTURING THE GOAL LOCATION OF EACH INDIVIDUAL
		for (xx = 0; xx < 3; xx++) {
			for (yy = 0; yy < 3; yy++) {
				d.cell[xx][yy] = 1.0;
			}
		}
		return d;
#endif  // PARTMIGRN 
	}
	else // goaltype == 2
		theta = atan2(((double)x -(double)smsData->goal.x),((double)y-(double)smsData->goal.y));
//  if (write_out) out<<"goalx,goaly: "<<goalx<<","<<goaly<<" theta: "<<theta<<endl;
	d = calcWeightings(gb,(float)theta);
}

return d;
}

// Calculate weightings for neighbouring cells
array3x3d Individual::calcWeightings(const double base,const double theta) {

array3x3d d; // 3x3 array indexed from SW corner by xx and yy
int dx,dy,xx,yy;

double i0 = 1.0; 					// direction of theta - lowest cost bias
double i1 = base;
double i2 = base * base;
double i3 = i2 * base;
double i4 = i3 * base;		// opposite to theta - highest cost bias

if (fabs(theta) > 7.0 * PI / 8.0) { dx = 0; dy = -1; }
else {
	if (fabs(theta) > 5.0 * PI / 8.0) { dy = -1; if (theta > 0) dx = 1; else dx = -1; }
	else {
		if (fabs(theta) > 3.0 * PI / 8.0) { dy = 0; if (theta > 0) dx = 1; else dx = -1; }
		else {
			if (fabs(theta) > PI / 8.0) { dy = 1; if (theta > 0) dx = 1; else dx = -1; }
			else { dy = 1; dx = 0; }
		}
  }
}
//  if (write_out) out<<"goalx,goaly: "<<goalx<<","<<goaly<<" dx,dy: "<<dx<<","<<dy
//    <<" theta: "<<theta<<endl;
d.cell[1][1] = 0; // central cell has zero weighting
d.cell[dx+1][dy+1] = (float)i0;
d.cell[-dx+1][-dy+1] = (float)i4;
if (dx == 0 || dy ==0) { // theta points to a cardinal direction
  d.cell[dy+1][dx+1] = (float)i2; d.cell[-dy+1][-dx+1] = (float)i2;
  if (dx == 0) { // theta points N or S
    xx = dx+1; if (xx > 1) dx -= 2; yy = dy;
    d.cell[xx+1][yy+1] = (float)i1; d.cell[-xx+1][yy+1] = (float)i1;
    d.cell[xx+1][-yy+1] = (float)i3; d.cell[-xx+1][-yy+1] = (float)i3;
  }
  else { // theta points W or E
		yy = dy+1; if (yy > 1) dy -= 2; xx = dx;
    d.cell[xx+1][yy+1] = (float)i1; d.cell[xx+1][-yy+1] = (float)i1;
    d.cell[-xx+1][yy+1] = (float)i3; d.cell[-xx+1][-yy+1] = (float)i3;
  }
}
else { // theta points to an ordinal direction
	d.cell[dx+1][-dy+1] = (float)i2; d.cell[-dx+1][dy+1] = (float)i2;          
  xx = dx+1; if (xx > 1) xx -= 2; d.cell[xx+1][dy+1] = (float)i1;
  yy = dy+1; if (yy > 1) yy -= 2; d.cell[dx+1][yy+1] = (float)i1;
  d.cell[-xx+1][-dy+1] = (float)i3; d.cell[-dx+1][-yy+1] = (float)i3;
  }

return d;
}

// Weight neighbouring cells on basis of (habitat) costs
array3x3f Individual::getHabMatrix(Landscape *pLand,Species *pSpecies,
	const int x,const int y,const short pr,const short prmethod,const short landIx,
	const bool absorbing)
{

array3x3f w; // array of effective costs to be returned
int ncells,x4,y4;
double weight,sumweights;
// NW and SE corners of effective cost array relative to the current cell (x,y):
int xmin = 0,ymin = 0,xmax = 0,ymax = 0;
int cost,nodatacost,h;
Cell *pCell;

landData land = pLand->getLandData();
if (absorbing) nodatacost = ABSNODATACOST;
else nodatacost = NODATACOST;

for (int x2=-1; x2<2; x2++) {   // index of relative move in x direction
	for (int y2=-1; y2<2; y2++) { // index of relative move in x direction

		w.cell[x2+1][y2+1] = 0.0; // initialise costs array to zeroes

		// set up corners of perceptual range relative to current cell
		if (x2==0 && y2==0) { // current cell - do nothing
			xmin=0; ymin=0; xmax=0; ymax=0;
		}
		else {
			if (x2==0 || y2==0) { // not diagonal (rook move)
				if (x2==0){ // vertical (N-S) move
					//out<<"ROOK N-S: x2 = "<<x2<<" y2 = "<<y2<<endl;
					if(pr%2==0) { xmin=-pr/2; xmax=pr/2; ymin=y2; ymax=y2*pr; } // PR even
					else { xmin=-(pr-1)/2; xmax=(pr-1)/2; ymin=y2; ymax=y2*pr; } // PR odd
				}
				if (y2==0) { // horizontal (E-W) move
					//out<<"ROOK E-W: x2 = "<<x2<<" y2 = "<<y2<<endl;
					if(pr%2==0) { xmin=x2; xmax=x2*pr; ymin=-pr/2; ymax=pr/2; } // PR even
					else { xmin=x2; xmax=x2*pr; ymin=-(pr-1)/2; ymax=(pr-1)/2; } // PR odd
				}
			}
			else { // diagonal (bishop move)
				//out<<"BISHOP: x2 = "<<x2<<" y2 = "<<y2<<endl;
				xmin=x2; xmax=x2*pr; ymin=y2; ymax=y2*pr;
			}
		}
		//out<<"pre  swap: xmin = "<<xmin<<" ymin = "<<ymin<<" xmax = "<<xmax<<" ymax = "<<ymax<<endl;
		if (xmin > xmax) { int z=xmax; xmax=xmin; xmin=z; } // swap xmin and xmax
		if (ymin > ymax) { int z=ymax; ymax=ymin; ymin=z; } // swap ymin and ymax
		//out<<"post swap: xmin = "<<xmin<<" ymin = "<<ymin<<" xmax = "<<xmax<<" ymax = "<<ymax<<endl;
//		if (write_out2) {
//			out2<<"current x and y  "<<x<<" "<<y<<endl;
//			out2<<"x2 and y2  "<<x2<<" "<<y2<<endl;
//			out2<<"xmin,ymin "<<xmin<<","<<ymin<<" xmax,ymax "<<xmax<<","<<ymax<<endl;
//		}

		// calculate effective mean cost of cells in perceptual range
		ncells = 0; weight = 0.0; sumweights = 0.0;
//		targetseen = 0;
		if (x2 != 0 || y2 != 0) { // not central cell (i.e. current cell)
			for (int x3=xmin; x3<=xmax; x3++) {
				for (int y3=ymin; y3<=ymax; y3++) {
          // if cell is out of bounds, treat landscape as a torus
          // for purpose of obtaining a cost,
					if ((x+x3) < 0) x4 = x+x3+land.maxX+1;
					else { if ((x+x3) > land.maxX) x4 = x+x3-land.maxX-1; else x4 = x+x3; }
					if ((y+y3) < 0) y4 = y+y3+land.maxY+1;
					else { if ((y+y3) > land.maxY) y4 = y+y3-land.maxY-1; else y4 = y+y3; }
//					if (write_out && (x4 < 0 || y4 < 0)) {
//						out<<"ERROR: x "<<x<<" y "<<y<<" x3 "<<x3<<" y3 "<<y3
//							<<" xbound "<<xbound<<" ybound "<<ybound<<" x4 "<<x4<<" y4 "<<y4<<endl;
//					}
					if (x4 < 0 || x4 > land.maxX || y4 < 0 || y4 > land.maxY) {
						// unexpected problem - e.g. due to ridiculously large PR
						// treat as a no-data cell
						cost = nodatacost;
					}
					else {
						// add cost of cell to total PR cost
						pCell = pLand->findCell(x4,y4);
						if (pCell == 0) { // no-data cell
							cost = nodatacost;
						}
						else {
							cost = pCell->getCost();
							if (cost < 0) cost = nodatacost;
							else {
								if (cost == 0) { // cost not yet set for the cell
									h = pCell->getHabIndex(landIx);
									cost = pSpecies->getHabCost(h);
#if RSDEBUG
//DEBUGLOG << "Individual::getHabMatrix(): x4=" << x4 << " y4=" << y4
//	<< " landIx=" << landIx << " h=" << h << " cost=" << cost
//	<< endl;
#endif
									pCell->setCost(cost);
								}
								else {
#if RSDEBUG
//DEBUGLOG << "Individual::getHabMatrix(): x4=" << x4 << " y4=" << y4
//	<< " cost=" << cost
//	<< endl;
#endif

                }
							}
						}
					}
					if (prmethod==1) { // arithmetic mean
						w.cell[x2+1][y2+1] += cost;
						ncells++;
					}
					if (prmethod==2) { // harmonic mean
            if (cost > 0) {
							w.cell[x2+1][y2+1] += (1.0f/(float)cost);
              ncells++;
            }
          }
          if (prmethod==3) { // arithmetic mean weighted by inverse distance
            if (cost>0) {
              // NB distance is still given by (x3,y3)
							weight = 1.0f /(double)sqrt((pow((double)x3,2)+pow((double)y3,2)));
              w.cell[x2+1][y2+1] += (float)(weight*(double)cost);
              ncells++; sumweights += weight;
            }
          }
//          if (write_out2) out2<<x+x3<<","<<y+y3<<","<<cost<<" wt. "<<weight<<endl;
//#if GO2TARGET
//					if (sq != 0) {
//						if (sq->get_target() > 50) targetseen++;
//					}
//#endif
				} //end of y3 loop
			}  //end of x3 loop
//		if (write_out) out<<"ncells in PR = "<<ncells<<" tot.wt. = "<<sumweights
//      <<" w.cell = "<<w.cell[x2+1][y2+1]<<endl;
		  if (ncells > 0) {
			  if (prmethod == 1) w.cell[x2+1][y2+1] /= ncells; // arithmetic mean
			  if (prmethod == 2) w.cell[x2+1][y2+1] = ncells/w.cell[x2+1][y2+1]; // hyperbolic mean
			  if (prmethod == 3 && sumweights > 0)
          w.cell[x2+1][y2+1] /= (float)sumweights; // weighted arithmetic mean
      }
//#if GO2TARGET
//      if (targetseen > 0) // target is within PR - set to a very low score
//        w.cell[x2+1][y2+1] = (1/(1000000*(double)targetseen));
//#endif
		}
		else { // central cell
			// record cost if not already recorded
			// has effect of preparing for storing effective costs for the cell
			pCell = pLand->findCell(x,y);
			cost = pCell->getCost();
			if (cost < 0) cost = nodatacost;
			else {
				if (cost == 0) { // cost not yet set for the cell
					h = pCell->getHabIndex(landIx);
					cost = pSpecies->getHabCost(h);
					pCell->setCost(cost);
				}
			}
		}
//		if (write_out2) out2<<"effective mean cost "<<w.cell[x2+1][y2+1]<<endl;

	}//end of y2 loop
}//end of x2 loop

return w;

}

//---------------------------------------------------------------------------
// Write records to individuals file
#if GROUPDISP  || ROBFITT
void Individual::outGenetics(const int rep,const int year,const int spnum,
	const int landNr,const bool patchmodel,const bool xtab)
{
#if RSDEBUG
//DEBUGLOG << "Individual::outGenetics(): indId=" << indId
//	<< " rep=" << rep << " landNr=" << landNr
//	<< endl;
#endif
if (landNr == -1) {
	if (pGenome != 0) {
		int X = -1; int Y = -1;
		if (patchmodel) {
			intptr ppatch = pCurrCell->getPatch();
			if (ppatch != 0) {
				Patch *pPatch = (Patch*) ppatch;
				X = pPatch->getPatchNum();
			}
		}
		else {
			locn loc = pCurrCell->getLocn();
      X = loc.x; Y = loc.y;
		}
		pGenome->outGenetics(rep,year,spnum,indId,X,Y,patchmodel,xtab);
	}
}
else { // open/close file
	pGenome->outGenHeaders(rep,landNr,patchmodel,xtab);
}

}
#else
void Individual::outGenetics(const int rep,const int year,const int spnum,
	const int landNr,const bool xtab)
{
#if RSDEBUG
//DEBUGLOG << "Individual::outGenetics(): indId=" << indId
//	<< " rep=" << rep << " landNr=" << landNr
//	<< endl;
#endif
if (landNr == -1) {
	if (pGenome != 0) {
		pGenome->outGenetics(rep,year,spnum,indId,xtab);
	}
}
else { // open/close file
	pGenome->outGenHeaders(rep,landNr,xtab);
}

}
#endif

#if RS_RCPP
//---------------------------------------------------------------------------
// Write records to movement paths file
void Individual::outMovePath(const int year)
{
	locn loc, prev_loc;

	//if (pPatch != pNatalPatch) {
	loc = pCurrCell->getLocn();
	// if still dispersing...
	if(status == 1){
		// at first step, record start cell first
		if(path->total == 1){
			prev_loc = pPrevCell->getLocn();
			outMovePaths << year << "\t" << indId << "\t"
						 << "0\t" << prev_loc.x << "\t" << prev_loc.y << "\t"
						 << "0\t"	// status at start cell is 0
						<< endl;
		}
		// then record current step
		outMovePaths << year << "\t" << indId << "\t"
			 << path->total << "\t" << loc.x << "\t" << loc.y << "\t"
			 << status << "\t"
			 << endl;
	}
	// if not anymore dispersing...
	if(status > 1 && status < 10){
		prev_loc = pPrevCell->getLocn();
		// record only if this is the first step as non-disperser
		if (path->pathoutput) {
			// if this is also the first step taken at all, record the start cell first
			if(path->total == 1){
				outMovePaths << year << "\t" << indId << "\t"
							 << "0\t" << prev_loc.x << "\t" << prev_loc.y << "\t"
							 << "0\t"	// status at start cell is 0
							<< endl;
			}
			outMovePaths << year << "\t" << indId << "\t"
						 << path->total << "\t" << loc.x << "\t" << loc.y << "\t"
						 << status << "\t"
						 << endl;
			// current cell will be invalid (zero), so set back to previous cell
			//pPrevCell = pCurrCell;
			path->pathoutput = 0;
		}
	}
}
#endif

//---------------------------------------------------------------------------

#if PEDIGREE
// Set position in relatedness matrix
void Individual::setMatPosn(unsigned int pos) { matPosn = pos; }
	// Get position in relatedness matrix
unsigned int Individual::getMatPosn(void)	{ return matPosn; }
#endif

//---------------------------------------------------------------------------

double wrpcauchy (double location, double rho) {
double result;

if(rho < 0.0 || rho > 1.0) {
	result = location;
}

if(rho == 0)
	result = pRandom->Random() * M_2PI;
else
	if(rho == 1) result = location;
	else {
		result = fmod(cauchy(location, -log(rho)), M_2PI);
	}
return result;
}

double cauchy(double location, double scale) {
if (scale < 0) return location;

return location + scale * tan(PI * pRandom->Random());
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
