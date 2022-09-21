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

#include "Population.h"
//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma package(smart_init)
#endif

ofstream outPop;
#if RS_CONTAIN
ofstream outCull;
#endif // RS_CONTAIN 
ofstream outInds;

//---------------------------------------------------------------------------

Population::Population(void) { 
nSexes = nStages = 0;
pPatch = NULL;
pSpecies = NULL;
return;
}

#if PEDIGREE
Population::Population(Species *pSp,Pedigree *pPed,Patch *pPch,int ninds,int resol) 
#else
Population::Population(Species *pSp,Patch *pPch,int ninds,int resol) 
#endif
{
// constructor for a Population of a specified size
#if RSDEBUG
//DEBUGLOG << "Population::Population(): this=" << this
//	<< " pPch=" << pPch << " ninds="<< ninds << endl;
#endif

int n,nindivs,age = 0,minage,maxage,nAges = 0;
int cumtotal = 0;
float probmale;
double ageprob,ageprobsum;
std::vector <double> ageProb; // for quasi-equilibrium initial age distribution
Cell *pCell;

if (ninds > 0) {
	inds.reserve(ninds);
	juvs.reserve(ninds);
}

pSpecies = pSp;
pPatch = pPch;
// record the new population in the patch
patchPopn pp;
pp.pSp = (intptr)pSpecies; pp.pPop = (intptr)this;
pPatch->addPopn(pp);
#if RSDEBUG
//DEBUGLOG << "Population::Population(): this=" << this
//	<< " added population to patch " << endl;
#endif

demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
//trfrSMSTraits sms = pSpecies->getSMSTraits();
settleType sett = pSpecies->getSettle();
genomeData gen = pSpecies->getGenomeData();
initParams init = paramsInit->getInit();

// determine no. of stages and sexes of species to initialise
if (dem.stageStruct) {
	nStages = sstruct.nStages;
}
else // non-structured population has 2 stages, but user only ever sees stage 1
	nStages = 2;
#if GROUPDISP
if (dem.repType == 0 || dem.repType == 3) { nSexes = 1; probmale = 0.0; }
else { nSexes = 2; probmale = dem.propMales; }
#else
if (dem.repType == 0) { nSexes = 1; probmale = 0.0; }
else { nSexes = 2; probmale = dem.propMales; }
#endif

// set up population sub-totals
for (int stg = 0; stg < NSTAGES; stg++) {
	for (int sex = 0; sex < NSEXES; sex++) {
		nInds[stg][sex] = 0;
	}
}

// set up local copy of minimum age table
short minAge[NSTAGES][NSEXES];
for (int stg = 0; stg < nStages; stg++) {
	for (int sex = 0; sex < nSexes; sex++) {
		if (dem.stageStruct) {
			if (dem.repType == 1) { // simple sexual model
				// both sexes use minimum ages recorded for females
				minAge[stg][sex] = pSpecies->getMinAge(stg,0);
			}
			else {
				minAge[stg][sex] = pSpecies->getMinAge(stg,sex);
			}
		}
		else { // non-structured population
			minAge[stg][sex] = 0;
		}
#if RSDEBUG
//DEBUGLOG << "Population::Population(): 1111 "
//	<< " minAge[" << stg << "][" << sex << "]=" << minAge[stg][sex]
//	<< endl;
#endif
	}
}

// individuals of new population must be >= stage 1
for (int stg = 1; stg < nStages; stg++) {
	if (dem.stageStruct) { // allocate to stages according to initialisation conditions
		// final stage is treated separately to ensure that correct total
		// no. of individuals is created
		if (stg == nStages-1) {
			n = ninds - cumtotal;
		}
		else {
			n = (int)(ninds * paramsInit->getProp(stg) + 0.5);
			cumtotal += n;
		}
	}
	else { // non-structured - all individuals go into stage 1
		n = ninds;
	}
//	for (int sex = 0; sex < nSexes; sex++) {
//		if (n < nSexes) n = nSexes; // to ensure at least one individual of each age is created
//		subPops.push_back(new SubPop(loc,stg,sex,n/nSexes));
//	}
	// establish initial age distribution
	minage = maxage = stg;
	if (dem.stageStruct) {
		// allow for stage-dependent minimum ages (use whichever sex is greater)
		if (minAge[stg][0] > 0 && minage < minAge[stg][0]) minage = minAge[stg][0];
		if (nSexes == 2 && minAge[stg][1] > 0 && minage < minAge[stg][1]) minage = minAge[stg][1];
		// allow for specified age distribution
		if (init.initAge != 0) { // not lowest age
			if (stg == nStages-1) maxage = sstruct.maxAge; // final stage
			else { // all other stages - use female max age, as sex of individuals is not predetermined
				maxage = minAge[stg+1][0] - 1;
			}
			if (maxage < minage) maxage = minage;
			nAges = maxage - minage + 1;
			if (init.initAge == 2) { // quasi-equilibrium distribution
#if RS_CONTAIN
#if SEASONAL
				double psurv = 1.0; // use female survival for the stage
				for (int s = 0; s < dem.nSeasons; s++) {
					psurv *= (double)pSpecies->getSurv(0,s,stg,0); // use female survival for the stage					
				}
#else
				double psurv = (double)pSpecies->getSurv(-1,stg,0); // use highest female survival for the stage
#endif // SEASONAL 
#else
#if SEASONAL
				double psurv = 1.0; // use female survival for the stage
				for (int s = 0; s < dem.nSeasons; s++) {
					psurv *= (double)pSpecies->getSurv(s,stg,0); // use female survival for the stage					
				}
#else
				double psurv = (double)pSpecies->getSurv(stg,0); // use female survival for the stage
#endif // SEASONAL 
#endif // RS_CONTAIN 
				ageProb.clear();
				ageprobsum = 0.0;
				ageprob = 1.0;
				for (int i = 0; i < nAges; i++) {
					ageProb.push_back(ageprob); ageprobsum += ageprob; ageprob *= psurv;
				}
				for (int i = 0; i < nAges; i++) {
					ageProb[i] /= ageprobsum;
					if (i > 0) ageProb[i] += ageProb[i-1]; // to give cumulative probability
				}
			}
		}
	}
#if RSDEBUG
//DEBUGLOG << "Population::Population(): this=" << this
//	<< " n=" << n << " stg=" << stg << " minage=" << minage << " maxage=" << maxage
//	<< endl;
#endif
	// create individuals
	int sex;
	nindivs = (int)inds.size();
	for (int i = 0; i < n; i++) {
		pCell = pPatch->getRandomCell();
		if (dem.stageStruct) {
			switch (init.initAge) {
			case 0: // lowest possible age
				age = minage;
				break;
			case 1: // randomised
				if (maxage > minage) age = pRandom->IRandom(minage,maxage);
				else age = minage;
				break;
			case 2: // quasi-equilibrium
				if (nAges > 1) {
					double rrr = pRandom->Random();
					int ageclass = 0;
					while (rrr > ageProb[ageclass]) ageclass++;
					age = minage + ageclass;
				}
				else age = minage;
				break;
			}
		}
		else age = stg;
#if RS_CONTAIN
#if RSDEBUG
	// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
		inds.push_back(new Individual(pCell,pPatch,stg,age,sstruct.repInterval,
			1,probmale,true,trfr.moveType));
#else
		inds.push_back(new Individual(pCell,pPatch,stg,age,sstruct.repInterval,
			1,probmale,trfr.moveModel,trfr.moveType));
#endif
#else
#if PARTMIGRN
#if RSDEBUG
	// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
		inds.push_back(new Individual(pSpecies,pCell,pPatch,stg,age,sstruct.repInterval,
			probmale,true,trfr.moveType));
#else
		inds.push_back(new Individual(pSpecies,pCell,pPatch,stg,age,sstruct.repInterval,
			probmale,trfr.moveModel,trfr.moveType));
#endif
#else
#if RSDEBUG
	// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
		inds.push_back(new Individual(pCell,pPatch,stg,age,sstruct.repInterval,
			probmale,true,trfr.moveType));
#else
		inds.push_back(new Individual(pCell,pPatch,stg,age,sstruct.repInterval,
			probmale,trfr.moveModel,trfr.moveType));
#endif
#endif // PARTMIGRN 
#endif // RS_CONTAIN 
		sex = inds[nindivs+i]->getSex();
#if PARTMIGRN
//		// TEMPORARY CODE TO SET MIGRATION STATUS
//		if (pRandom->Bernoulli(0.5)) inds[nindivs+i]->setMigrnStatus(1);
//		else {
//			if (pRandom->Bernoulli(0.5)) inds[nindivs+i]->setMigrnStatus(2);
//			else inds[nindivs+i]->setMigrnStatus(3);
//			patchlist p; p.pPatch = pPatch; p.season = 0; p.breeding = p.fixed = true;
//			inds[nindivs+i]->addPatch(p);
//		}
#endif  // PARTMIGRN 
#if GOBYMODEL
		if (true)
#else
#if SOCIALMODEL
		if (true)
#else
		if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers)
#endif
#endif
		{
			// individual variation - set up genetics
			inds[nindivs+i]->setGenes(pSpecies,resol);
		}
		nInds[stg][sex]++;
#if GROUPDISP
#if PEDIGREE
		// add new Individual to relationship table if it has breeding status
		bool addind = false;
		if (dem.stageStruct) {
			if (stg == nStages-1) addind = true;
			// NOTE - FOR CONVENIENCE CURRENTLY ASSUMED TO BE ONLY THE FINAL STAGE,
			// BUT SHOULD CHECK FOR POSSIBLE NON-ZERO FECUNDITY OF EARLIER STAGE(S)
		}
		else addind = true;
		if (addind) {
			int posn = pPed->addInd(inds[nindivs+i]);
			inds[nindivs+i]->setMatPosn(posn);
			pPed->setRelMat(posn,posn,1.0);
			for (int j = 0; j < posn; j++) {
				pPed->setRelMat(j,posn,0.0);
				pPed->setRelMat(posn,j,0.0);				
			}
		}
#endif
#endif
	}
}
#if RS_CONTAIN
nCulled = 0;
selectedForCull = false;
#endif // RS_CONTAIN 
#if RSDEBUG
//DEBUGLOG << "Population::Population(): this=" << this
//	<< " finished " << endl;
#endif
}

Population::~Population(void) {
int ninds = (int)inds.size();
for (int i = 0; i < ninds; i++) {
	if (inds[i] != NULL) delete inds[i];
}
inds.clear();
int njuvs = (int)juvs.size();
for (int i = 0; i < njuvs; i++) {
	if (juvs[i] != NULL) delete juvs[i];
}
juvs.clear();
}

traitsums Population::getTraits(Species *pSpecies) {
int g;
traitsums ts;       
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
//locus loc;

emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();

int ninds = (int)inds.size();
#if RSDEBUG
//DEBUGLOG << "Population::getTraits(): ninds = " << ts.ninds[0]
////	<< " nalleles = "<< nalleles
////	<< " nemiggenes = " << nemiggenes << " ntrfrgenes = " << ntrfrgenes
//	<< endl;
#endif
for (int i = 0; i < ninds; i++) {
	int sex = inds[i]->getSex();
	if (emig.sexDep || trfr.sexDep || sett.sexDep) g = sex; else g = 0;
	ts.ninds[g] += 1;
#if SOCIALMODEL
	// social traits
#endif
	// emigration traits
	emigTraits e = inds[i]->getEmigTraits();
	if (emig.sexDep) g = sex; else g = 0;
	ts.sumD0[g]    += e.d0;    ts.ssqD0[g] 		+= e.d0 * e.d0;
	ts.sumAlpha[g] += e.alpha; ts.ssqAlpha[g] += e.alpha * e.alpha;
	ts.sumBeta[g]  += e.beta;  ts.ssqBeta[g]	+= e.beta * e.beta;
	// transfer traits
	trfrKernTraits k = inds[i]->getKernTraits();
	if (trfr.sexDep) g = sex; else g = 0;
	ts.sumDist1[g] += k.meanDist1; ts.ssqDist1[g] += k.meanDist1 * k.meanDist1;
	ts.sumDist2[g] += k.meanDist2; ts.ssqDist2[g] += k.meanDist2 * k.meanDist2;
	ts.sumProp1[g] += k.probKern1; ts.ssqProp1[g] += k.probKern1 * k.probKern1;
	trfrSMSTraits sms = inds[i]->getSMSTraits();
	g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
	ts.sumDP[g] += sms.dp; ts.ssqDP[g] += sms.dp * sms.dp;
	ts.sumGB[g] += sms.gb; ts.ssqGB[g] += sms.gb * sms.gb;
	ts.sumAlphaDB[g] += sms.alphaDB; ts.ssqAlphaDB[g] += sms.alphaDB * sms.alphaDB;
	ts.sumBetaDB[g]  += sms.betaDB;  ts.ssqBetaDB[g]  += sms.betaDB * sms.betaDB;
#if RSDEBUG
//DEBUGLOG << "Population::getTraits():"
//	<< " i=" << i << " g=" << g
//	<< " sms.dp= " << sms.dp << " sms.gb= " << sms.gb
//	<< " ts.sumDP[g]= " << ts.sumDP[g] << " ts.ssqDP[g]= " << ts.ssqDP[g]
//	<< " ts.sumGB[g]= " << ts.sumGB[g] << " ts.ssqGB[g]= " << ts.ssqGB[g]
//	<< endl;
#endif
	trfrCRWTraits c = inds[i]->getCRWTraits();
	g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
	ts.sumStepL[g] += c.stepLength; ts.ssqStepL[g] += c.stepLength * c.stepLength;
	ts.sumRho[g]   += c.rho;        ts.ssqRho[g] += c.rho * c.rho;
	// settlement traits
	settleTraits s = inds[i]->getSettTraits();
	if (sett.sexDep) g = sex; else g = 0;
//	g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
	ts.sumS0[g]    += s.s0;     ts.ssqS0[g] 		+= s.s0 * s.s0;
	ts.sumAlphaS[g] += s.alpha; ts.ssqAlphaS[g] += s.alpha * s.alpha;
	ts.sumBetaS[g] += s.beta;   ts.ssqBetaS[g]	+= s.beta * s.beta;
#if RSDEBUG
//DEBUGLOG << "Population::getTraits():"
//	<< " i=" << i << " g=" << g << " a=" << a
//	<< " e.d0= " << e.d0 << " e.alpha= " << e.alpha << " e.beta= " << e.beta
//	<< " mnd0= " << emigTraits[g]->mnD0 << " mnAlpha= " << emigTraits[g]->mnAlpha << " mnBeta= " << emigTraits[g]->mnBeta
//	<< " sqd0= " << emigTraits[g]->sqD0 << " sqAlpha= " << emigTraits[g]->sqAlpha << " sqBeta= " << emigTraits[g]->sqBeta
//	<< endl;
#endif
}

return ts;
}

#if GROUPDISP
int Population::getNGroups(void) { return (int)groups.size(); }
#endif
int Population::getNInds(void) { return (int)inds.size(); }

#if VIRTUALECOLOGIST
Individual* Population::getInd(int i) {
if (i >= 0 && i < inds.size()) return inds[i];
else return 0;
}
#endif

#if RS_CONTAIN
popStats Population::getStats(short hab) 
#else
#if SPATIALDEMOG
popStats Population::getStats(std::vector <float> localDemoScaling)
#else
popStats Population::getStats(void) 
#endif // SPATIALDEMOG
#endif // RS_CONTAIN 
{
popStats p;
int ninds;
float fec;
bool breeders[2]; breeders[0] = breeders[1] = false;
demogrParams dem = pSpecies->getDemogr();
p.pSpecies = pSpecies;
p.pPatch = pPatch;
p.spNum = pSpecies->getSpNum();
p.nInds = (int)inds.size();
p.nNonJuvs = p.nAdults = 0;
p.breeding = false;
#if RSDEBUG
//DEBUGLOG << "Population::getStats(): this=" << this
////	<< " p.pSpecies=" << p.pSpecies << " p.spNum=" << p.spNum
//	<< " p.pPatch=" << p.pPatch << " patchNum=" << p.pPatch->getPatchNum()
//	<< " nStages=" << nStages << " nSexes=" << nSexes << " p.nInds=" << p.nInds
//	<< endl;
#endif
for (int stg = 1; stg < nStages; stg++) {
	for (int sex = 0; sex < nSexes; sex++) {
		ninds = nInds[stg][sex];
		p.nNonJuvs += ninds;
#if RSDEBUG
//DEBUGLOG << "Population::getStats(): this=" << this
//	<< " stg=" << stg << " sex=" << sex
//	<< " nInds[stg][sex]=" << nInds[stg][sex] << " p.nNonJuvs=" << p.nNonJuvs
//	<< endl;
#endif
		if (ninds > 0) {
			if (pSpecies->stageStructured()) {
#if RS_CONTAIN
#if SEASONAL
				fec = 0.0;
				float seasonfec;
				for (int s = 0; s < dem.nSeasons; s++) {
					if (dem.repType == 2) seasonfec = pSpecies->getFec(hab,s,stg,sex);
					else seasonfec = pSpecies->getFec(hab,s,stg,0);
					if (seasonfec > fec) fec = seasonfec;
				}
#else
				if (dem.repType == 2) fec = pSpecies->getFec(hab,stg,sex);
				else fec = pSpecies->getFec(hab,stg,0);
#endif // SEASONAL 
#else
#if SEASONAL
				fec = 0.0;
				float seasonfec;
				for (int s = 0; s < dem.nSeasons; s++) {
					if (dem.repType == 2) seasonfec = pSpecies->getFec(s,stg,sex);
					else seasonfec = pSpecies->getFec(s,stg,0);
					if (seasonfec > fec) fec = seasonfec;
				}
#else
#if SPATIALDEMOG
				if (dem.repType == 2) {
					if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,sex)>=0){
						fec = pSpecies->getFec(stg,sex)*localDemoScaling[pSpecies->getFecLayer(stg,sex)];
					}
					else fec = pSpecies->getFec(stg,sex);
				}
				else {
					if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,0)>=0){
						fec = pSpecies->getFec(stg,0)*localDemoScaling[pSpecies->getFecLayer(stg,0)];
					}
					else fec = pSpecies->getFec(stg,0);
				}
#else

				if (dem.repType == 2) fec = pSpecies->getFec(stg,sex);
				else fec = pSpecies->getFec(stg,0);
#endif //SPATIALDEMOG
#endif // SEASONAL 
#endif // RS_CONTAIN 
#if GROUPDISP
				if (dem.repType == 3) { // hermaphrodite
					if (dem.selfing) {
						if (fec > 0.0) { breeders[sex] = true; p.nAdults += ninds; }
					}
					else {
						// there must be at least 2 individuals present
						if (ninds > 1 && fec > 0.0) { breeders[sex] = true; p.nAdults += ninds; }
					}
				}
				else {
					if (fec > 0.0) { breeders[sex] = true; p.nAdults += ninds; }
				}
#else
				if (fec > 0.0) { breeders[sex] = true; p.nAdults += ninds; }
#endif
			}
			else breeders[sex] = true;
		}
	}
}
// is there a breeding population present?
if (nSexes == 1) {
	p.breeding = breeders[0];
}
else {
	if (breeders[0] && breeders[1]) p.breeding = true;
}
#if GOBYMODEL
// count number of each sociality phenotype
// NOTE: THIS IS VERY INEFFICIENT - IF THIS FEATURE BECOMES PERMANENT,
// MAINTAIN RECORD OF PHENOTYPES WITHIN POPULATION OF NUMBERS
p.nAsocial = p.nSocial = 0;
ninds = (int)inds.size();
indStats ind;
for (int i = 0; i < ninds; i++) {
	ind = inds[i]->getStats();
	if (ind.stage > 0) { // not a juvenile
		if (inds[i]->isAsocial()) p.nAsocial++; else p.nSocial++;
	}
}
#endif
#if SOCIALMODEL
// count number of each sociality phenotype
// NOTE: THIS IS INEFFICIENT - IF THIS FEATURE BECOMES PERMANENT, MAINTAIN
// RECORD WITHIN POPULATION OF NUMBERS
p.nAsocial = p.nSocial = 0;
ninds = (int)inds.size();
for (int i = 0; i < ninds; i++) {
	if (inds[i]->isAsocial()) p.nAsocial++; else p.nSocial++;
}
#endif
#if RSDEBUG
//DEBUGLOG << "Population::getStats(): this=" << this
//	<< " p.nInds=" << p.nInds << " p.nAdults=" << p.nAdults << " p.nNonJuvs=" << p.nNonJuvs
//	<< " breeders[0]=" << breeders[0] << " breeders[1]=" << breeders[1]
//	<< " p.breeding=" << p.breeding
//	<< endl;
#endif
return p;
}

Species* Population::getSpecies(void) { return pSpecies; }

int Population::totalPop(void) {
int t = 0;
for (int stg = 0; stg < nStages; stg++) {
	for (int sex = 0; sex < nSexes; sex++) {
		t += nInds[stg][sex];
#if RSDEBUG
//DEBUGLOG << "Population::totalPop(): this=" << this
//	<< " stg=" << stg << " sex=" << sex
//	<< " nInds[stg][sex]=" << nInds[stg][sex] << " t=" << t
//	<< endl;
#endif
	}
}
return t;
}

int Population::stagePop(int stg) {
int t = 0;
if (stg < 0 || stg >= nStages) return t;
for (int sex = 0; sex < nSexes; sex++) {
	t += nInds[stg][sex];
}
return t;
}

#if RS_ABC

int Population::sexPop(short sex) {
int t = 0;
if (sex < 0 || sex >= nSexes) return t;
for (int stg = 0; stg < nStages; stg++) {
	t += nInds[stg][sex];
}
return t;
}

int Population::stageSexPop(short	stg,short sex) {
if (stg < 0 || stg >= nStages) return 0;
if (sex < 0 || sex >= nSexes) return 0;
return nInds[stg][sex];
}
	
#endif // RS_ABC

//---------------------------------------------------------------------------
// Remove all Individuals
void Population::extirpate(void) {
int ninds = (int)inds.size();
for (int i = 0; i < ninds; i++) {
	if (inds[i] != NULL) delete inds[i];
}
inds.clear();
int njuvs = (int)juvs.size();
for (int i = 0; i < njuvs; i++) {
	if (juvs[i] != NULL) delete juvs[i];
}
juvs.clear();
for (int sex = 0; sex < nSexes; sex++) {
	for (int stg = 0; stg < nStages; stg++) {
		nInds[stg][sex] = 0;
	}
}
}

#if RS_CONTAIN

// Remove individuals with a specified probability
//void Population::cull(Cull *pCull,double pculled) 
void Population::cull(Cull *pCull, double hectares)
{
culldata c = pCull->getCullData();
cullstagedata cstage;
double pculled[NSTAGES];
for (int i = 0; i < nStages; i++) pculled[i] = 0.0;
double density = 0.0;
int nINNS = 0;
switch (c.cullRate) {	
case 0: // constant
	pculled[0] = c.cullMaxRate;
	break;
case 1: // logistic function of density of stages to be culled
	if (pSpecies->stageStructured()) {
		for (int i = 0; i < nStages; i++) {
			if (pCull->getCullStage(i)) { // stage is to be culled
				nINNS += stagePop(i);
			}
		}		
	}
	else {
		nINNS = getNInds();
	}
//	double density = (double)nINNS / hectares;
	density = (double)nINNS / hectares;
	pculled[0] = c.cullMaxRate / (1.0 + exp(-(density-c.cullBeta)*c.cullAlpha));	 
	break;
case 2: // stage-specific logistic functions of density
	for (int i = 0; i < nStages; i++) {
		if (pCull->getCullStage(i)) { // stage is to be culled
			cstage = pCull->getCullStageData(i);
			density = (double)stagePop(i) / hectares;
			pculled[i] = cstage.cullMaxRate / (1.0 + exp(-(density-cstage.cullBeta)*cstage.cullAlpha));
		}
	}
	break;
}
int ninds = (int)inds.size();
int njuvs = (int)juvs.size();
#if RSDEBUG
DEBUGLOG << "Population::cull(): this=" << this << " patch=" << pPatch->getPatchNum()	
	<< " hectares=" << hectares 
	<< " cullRate=" << c.cullRate << " cullMaxRate=" << c.cullMaxRate 
	<< " cullAlpha=" << c.cullAlpha << " cullBeta=" << c.cullBeta 
	<< " stagePop(0)=" << stagePop(0) << " stagePop(1)=" << stagePop(1) << " stagePop(2)=" << stagePop(2)
//	<< " getNInds=" << getNInds()
//	<< " ninds=" << ninds << " njuvs=" << njuvs
	<< " density=" << density 
	<< " pculled[0]=" << pculled[0] 
//	<< " nCulled=" << nCulled
	<< endl;
#endif
demogrParams dem = pSpecies->getDemogr();
nCulled = 0;
selectedForCull = true;
indStats ind;
for (int i = 0; i < ninds; i++) {
	if (dem.stageStruct) {
		ind = inds[i]->getStats();  
		if (pCull->getCullStage(ind.stage)) {
			if (c.cullRate == 2) {
				if (pRandom->Bernoulli(pculled[ind.stage])) { inds[i]->setStatus(10); nCulled++; }						
			}
			else {
				if (pRandom->Bernoulli(pculled[0])) { inds[i]->setStatus(10); nCulled++; }						
			}
#if RSDEBUG
//DEBUGLOG << "Population::cull(): i=" << i << " ind.stage=" << ind.stage 
//	<< " nCulled=" << nCulled 
//	<< endl;
#endif
		}			
	}
	else {
		if (pRandom->Bernoulli(pculled[0])) { inds[i]->setStatus(10); nCulled++; }		
	}
}
/*
for (int i = 0; i < njuvs; i++) {
	if (dem.stageStruct) {
		ind = juvs[i]->getStats();  
		if (pCull->getCullStage(ind.stage)) {
			if (pRandom->Bernoulli(pculled)) { juvs[i]->setStatus(10); nCulled++; }		
		}	
	}	
	else {
		if (pRandom->Bernoulli(pculled)) { juvs[i]->setStatus(10); nCulled++; }		
	}
}
*/
#if RSDEBUG
DEBUGLOG << "Population::cull(): this=" << this << " patch=" << pPatch->getPatchNum() 
	<< " nCulled=" << nCulled 
	<< endl;
#endif
} 

void Population::resetCull(void) { nCulled = 0; selectedForCull = false; } 

#endif // RS_CONTAIN 

#if GROUPDISP
Individual* Population::getFather(int minbrdstage,int ix) {
indStats ind = inds[ix]->getStats();
if (ind.stage >= minbrdstage) return inds[ix];
else return 0;
}
#endif

//---------------------------------------------------------------------------
// Produce juveniles and hold them in the juvs vector
#if RS_CONTAIN
#if SEASONAL
void Population::reproduction(const int hab,const int season,const float localK,const float envval,const int resol)
#else
void Population::reproduction(const int hab,const float localK,const float envval,const int resol)
#endif // SEASONAL 
#else
#if SEASONAL
void Population::reproduction(const int season,const float localK,const float envval,const int resol)
#else
#if GROUPDISP
void Population::reproduction(const std::vector <Individual*> *pfglobal,const int nfglobal,
	const std::vector <Individual*> *pfnbrhd,const int nfnbrhd,
	const float localK,const float envval,const int resol)
#else
#if BUTTERFLYDISP
void Population::reproduction(const float localK,const float envval,const int resol,
	const short option)
#else
#if SPATIALDEMOG
void Population::reproduction(const float localK,const float envval,const int resol, std::vector <float> localDemoScaling)
#else
void Population::reproduction(const float localK,const float envval,const int resol)
#endif // SPATIALDEMOG
#endif // BUTTERFLYDISP 
#endif // GROUPDISP 
#endif // SEASONAL 
#endif // RS_CONTAIN  
{

// get population size at start of reproduction
int ninds = (int)inds.size();
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): this=" << this
//#if BUTTERFLYDISP
//	<< " option=" << option
//#endif // BUTTERFLYDISP 
//#if RS_CONTAIN
//	<< " hab=" << hab
//#endif // RS_CONTAIN 
//	<< " ninds=" << ninds
//	<< endl;
#endif // RSDEBUG 
if (ninds == 0) return;

int nsexes,stage,sex,njuvs,nj,nmales,nfemales;
Cell *pCell;
indStats ind;
double expected;
bool skipbreeding;

//envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
genomeData gen = pSpecies->getGenomeData();

#if GROUPDISP
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): this=" << this
//	<< " pfglobal=" << pfglobal << " nglobal=" << nglobal
//	<< " pfglobal[0]=" << (*pfglobal)[0] << " pfglobal[1]=" << (*pfglobal)[1]
//	<< " pfglobal[nglobal-1]=" << (*pfglobal)[nglobal-1]
//	<< endl;
#endif
#endif

#if SOCIALMODEL

// ADDITIONAL VARIABLES FOR PROBIS SOCIAL POLYMORPHISM MODEL
float fecAsocial;
double fasocial,fsocial; // Allee effect below threshold
socialParams soc = pSpecies->getSocialParams();
int ncells = pPatch->getNCells();
if ((float)ninds > soc.Ts*localK) fsocial = 1.0;
//else fsocial  = 1.0 / (1.0 + soc.Cs*(soc.Ts*localK-(float)ninds)*(soc.Ts*localK-(float)ninds));
//else fsocial  = 1.0 / (1.0 + exp(soc.cs + soc.bs*(soc.Ts*localK-(float)ninds)));
else fsocial  = 1.0 / (1.0 + exp(soc.cs + soc.bs*(soc.Ts*localK-(float)ninds)/(float)ncells));
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): this=" << this
//	<< " totalPop=" << totalPop() << " ncells=" << ncells << " ninds=" << ninds
//	<< " soc.cs=" << soc.cs << " soc.bs=" << soc.bs
//	<< " soc.Ts=" << soc.Ts << " localK=" << localK
//	<< " EXParg=" << soc.cs + soc.bs*(soc.Ts*localK-(float)ninds)/(float)ncells
//	<< " fsocial=" << fsocial
//	<< endl;
#endif
if ((float)ninds > soc.Ta*localK) fasocial = 1.0;
//else fasocial = 1.0 / (1.0 + soc.Ca*(soc.Ta*localK-(float)ninds)*(soc.Ta*localK-(float)ninds));
//else fasocial  = 1.0 / (1.0 + exp(soc.ca + soc.ba*(soc.Ta*localK-(float)ninds)));
else fasocial  = 1.0 / (1.0 + exp(soc.ca + soc.ba*(soc.Ta*localK-(float)ninds)/(float)ncells));
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): this=" << this
//	<< " totalPop=" << totalPop() << " ncells=" << ncells << " ninds=" << ninds
//	<< " soc.ca=" << soc.ca << " soc.ba=" << soc.ba
//	<< " soc.Ta=" << soc.Ta << " localK=" << localK
//	<< " EXParg=" << soc.ca + soc.ba*(soc.Ta*localK-(float)ninds)/(float)ncells
//	<< " fasocial=" << fasocial
//	<< endl;
#endif

#endif // SOCIALMODEL

#if GROUPDISP
if (dem.repType == 0 || dem.repType == 3) nsexes = 1; else nsexes = 2;
//int nglobal = fglobal.size();
//int nneighbourhood = neighbourhood.size();
//int nneighbourhood = 0;
// CHANGE OF METHOD IMPLEMENTED 18/10/17
// if there are no potential fathers in the neighbourhood, then select
// father at random from local/global pools in their specified ratio
float propLocal,propNghbr;
propLocal = propNghbr = 0.0;
if (dem.paternity == 2) {
#if RSDEBUG
//DEBUGLOG << "Population::reproduction():"
//	<< " nfnbrhd=" << nfnbrhd
//	<< " dem.propLocal=" << dem.propLocal << " dem.propNghbr=" << dem.propNghbr
//	<< endl;
#endif
	if (nfnbrhd > 0) {
		propLocal = dem.propLocal; propNghbr = dem.propNghbr;
	}
else {
		propLocal = dem.propLocal / (1.0 - dem.propNghbr); propNghbr = 0.0;
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): AMENDING POLLEN KERNEL"
//	<< " propLocal=" << propLocal << " propNghbr=" << propNghbr
//	<< endl;
#endif
}
}
#else
if (dem.repType == 0) nsexes = 1; else nsexes = 2;
#endif // GROUPDISP 

#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): this=" << this
//#if RS_CONTAIN
//	<< " hab=" << hab
//#else
//	<< " pSpecies=" << pSpecies
//	<< " localK=" << localK << " envval=" << envval << " resol=" << resol
//#endif // RS_CONTAIN 
//	<< " sstruct.nStages=" << sstruct.nStages << " nsexes=" << nsexes << " ninds=" << ninds
//	<< endl;
#endif

// set up local copy of species fecundity table
float fec[NSTAGES][NSEXES];
#if GOBYMODEL
// also set up corresponding table for density-dependent effects, which cannot
// be applied until an individual female's sociality phenotype is known
float ddeffect[NSTAGES];
for (int i = 0; i < NSTAGES; i++) { ddeffect[i] = 1.0; }
#endif
#if SOCIALMODEL
// also set up corresponding table for density-dependent effects, which cannot
// be applied until an individual female's sociality phenotype is known
//float ddeffect[NSTAGES];
//for (int i = 0; i < NSTAGES; i++) { ddeffect[i] = 1.0; }
#endif
for (int stg = 0; stg < sstruct.nStages; stg++) {
	for (int sex = 0; sex < nsexes; sex++) {
		if (dem.stageStruct) {
			if (dem.repType == 1) { // simple sexual model
				// both sexes use fecundity recorded for females
#if RS_CONTAIN
#if SEASONAL
				fec[stg][sex] = pSpecies->getFec(hab,season,stg,0);
#else
				fec[stg][sex] = pSpecies->getFec(hab,stg,0);
#endif // SEASONAL 
#else
#if SEASONAL
				fec[stg][sex] = pSpecies->getFec(season,stg,0);
#else
#if SPATIALDEMOG
				if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,0)>=0){
					fec[stg][sex] = pSpecies->getFec(stg,0)*localDemoScaling[pSpecies->getFecLayer(stg,0)];
				}
				else fec[stg][sex] = pSpecies->getFec(stg,0);
#else
				fec[stg][sex] = pSpecies->getFec(stg,0);
#endif //SPATIALDEMOG
#endif // SEASONAL 
#endif // RS_CONTAIN 
			}
			else
#if RS_CONTAIN
#if SEASONAL
				fec[stg][sex] = pSpecies->getFec(hab,season,stg,sex);
#else
				fec[stg][sex] = pSpecies->getFec(hab,stg,sex);
#endif // SEASONAL 
#else
#if SEASONAL
				fec[stg][sex] = pSpecies->getFec(season,stg,sex);
#else
#if SPATIALDEMOG
				if (pSpecies->getFecSpatial() && pSpecies->getFecLayer(stg,sex)>=0){
					fec[stg][sex] = pSpecies->getFec(stg,sex)*localDemoScaling[pSpecies->getFecLayer(stg,sex)];
				}
				else fec[stg][sex] = pSpecies->getFec(stg,sex);
#else
				fec[stg][sex] = pSpecies->getFec(stg,sex);
#endif // SPATIALDEMOG
#endif // SEASONAL 
#endif // RS_CONTAIN 
//			if (sex == 0 && fec[stg][sex] > dem.lambda) dem.lambda = fec[stg][sex];
		}
		else { // non-structured population
			if (stg == 1) fec[stg][sex] = dem.lambda; // adults
			else fec[stg][sex] = 0.0; // juveniles
		}
#if RSDEBUG
//if (ninds > 0) {
//DEBUGLOG << "Population::reproduction(): fec[" << stg << "][" << sex << "] = " << fec[stg][sex]
//	<< endl;
//}
#endif
	}
}

if (dem.stageStruct) {
#if RSDEBUG
//if (ninds > 0) {
//	DEBUGLOG << "Population::reproduction(): ninds=" << ninds << " localK=" << localK
//		<< " effect of density dependence:" << endl;
//}
#endif
	// apply environmental effects and density dependence
	// to all non-zero female non-juvenile stages
	for (int stg = 1; stg < nStages; stg++) {
		if (fec[stg][0] > 0.0) {
			// apply any effect of environmental gradient and/or stochasticty
			fec[stg][0] *= envval;
			if (env.stoch && !env.inK) {
				// fecundity (at low density) is constrained to lie between limits specified
				// for the species
				float limit;
				limit = pSpecies->getMinMax(0);
				if (fec[stg][0] < limit) fec[stg][0] = limit;
				limit = pSpecies->getMinMax(1);
				if (fec[stg][0] > limit) fec[stg][0] = limit;
			}
			if (sstruct.fecDens) { // apply density dependence
				float effect = 0.0;
				if (sstruct.fecStageDens) { // stage-specific density dependence
					// NOTE: matrix entries represent effect of ROW on COLUMN 
					// AND males precede females
					float weight = 0.0;                  
					for (int effstg = 0; effstg < nStages; effstg++) {
						for (int effsex = 0; effsex < nSexes; effsex++) {
							if (dem.repType == 2) {
								if (effsex == 0) weight = pSpecies->getDDwtFec(2*stg+1,2*effstg+1); 
								else weight = pSpecies->getDDwtFec(2*stg+1,2*effstg);
							}
							else {
								weight = pSpecies->getDDwtFec(stg,effstg);								
							}
							effect += (float)nInds[effstg][effsex] * weight;								
#if RSDEBUG
//if (ninds > 0) {
//	DEBUGLOG << " effstg=" << effstg << " effsex=" << effsex << " nInds=" << nInds[effstg][effsex];
//	DEBUGLOG << " weight=" << weight << " effect=" << effect
//		<< endl;
//}
#endif
						}
					}
				}
				else // not stage-specific
					effect = (float)totalPop();
#if GOBYMODEL
				if (localK > 0.0) ddeffect[stg] = effect/localK;
#else
				if (localK > 0.0) fec[stg][0] *= exp(-effect/localK);
#endif
#if RSDEBUG
//if (ninds > 0) {
//	DEBUGLOG << " eff popn=" << effect << " exponential=" << exp(-effect/localK);
//	DEBUGLOG << " fec[" << stg << "][0]=" << fec[stg][0]
//		<< endl;
//}
#endif
			}
		}
	}
}
else { // non-structured - set fecundity for adult females only
	// apply any effect of environmental gradient and/or stochasticty
	fec[1][0] *= envval;
#if SOCIALMODEL
	fecAsocial = soc.asocRmax * dem.lambda * envval;
#endif
	if (env.stoch && !env.inK) {
		// fecundity (at low density) is constrained to lie between limits specified
		// for the species
		float limit;
		limit = pSpecies->getMinMax(0);
		if (fec[1][0] < limit) fec[1][0] = limit;
#if SOCIALMODEL
		if (fecAsocial < limit) fecAsocial = limit;
#endif
		limit = pSpecies->getMinMax(1);
		if (fec[1][0] > limit) fec[1][0] = limit;
#if SOCIALMODEL
		if (fecAsocial > limit) fecAsocial = limit;
#endif
	}
	// apply density dependence
	if (localK > 0.0) {
//#if GOBYMODEL
//		ddeffect[1] = (float)ninds/localK;
//#else
		if (dem.repType == 1 || dem.repType == 2) { // sexual model
			// apply factor of 2 (as in manual, eqn. 6)
			fec[1][0] *= 2.0;
		}
#if SOCIALMODEL
		fec[1][0] *= fsocial / (1.0 + fabs(dem.lambda-1.0)*pow(((float)ninds/localK),dem.bc));
		fecAsocial *= fasocial / (1.0 + fabs((soc.asocRmax*dem.lambda)-1.0)
			*pow(((float)ninds/(soc.asocK*localK)),(soc.asocBc*dem.bc)));
#else
		fec[1][0] /= (1.0f + fabs(dem.lambda-1.0f)*pow(((float)ninds/localK),dem.bc));
#endif
//#endif
	}
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): dem.lambda=" << dem.lambda << " ninds=" << ninds
//	<< " localK=" << localK << " dem.bc=" << dem.bc << " fec[1][0]=" << fec[1][0]
//	<< endl;
#endif
}

double propBreed;
Individual *father;
std::vector <Individual*> fathers;

switch (dem.repType) {

case 0: // asexual model
	for (int i = 0; i < ninds; i++) {
		stage = inds[i]->breedingFem();
#if RSDEBUG
//DEBUGLOG << "Population::reproduction():"
//	<< " i=" << i << " ID=" << inds[i]->getId() << " stage=" << stage
//	<< endl;
#endif
		if (stage > 0) { // female of breeding age
			if (dem.stageStruct) {
				// determine whether she must miss current breeding attempt
				ind = inds[i]->getStats();
				if (ind.fallow >= sstruct.repInterval) {
					if (pRandom->Bernoulli(sstruct.probRep)) skipbreeding = false;
					else skipbreeding = true;
				}
				else skipbreeding = true; // cannot breed this time
			}
			else skipbreeding = false; // not structured - always breed
			if (skipbreeding) {
				inds[i]->incFallow();
			}
			else { // attempt to breed
				inds[i]->resetFallow();
#if GOBYMODEL
				if (dem.stageStruct) {
					if (inds[i]->isAsocial()) {
						expected = fec[stage][0] * sstruct.asocF * exp(-sstruct.asocF*ddeffect[stage]);
					}
					else {
						expected = fec[stage][0] * exp(-ddeffect[stage]);
					}
				}
				else { // non-structured population
					// effect of sociality has not been defined (as at 3/6/2016)
					expected = fec[stage][0];
				}
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): 1111 ID=" << inds[i]->getId()
//	<< " asocial=" << inds[i]->isAsocial() << " dem.asocF=" << dem.asocF
//	<< " fec[stage][0]=" << fec[stage][0] << " ddeffect=" << ddeffect[stage]
//	<< " expected=" << expected
//	<< endl;
#endif
#else
	#if SOCIALMODEL
				if (inds[i]->isAsocial())
					expected = fecAsocial;
				else
					expected = fec[stage][0];
	#else
				expected = fec[stage][0];
	#endif // SOCIALMODEL
#endif // GOBYMODEL
				if (expected <= 0.0) njuvs = 0;
				else njuvs = pRandom->Poisson(expected);
#if BUTTERFLYDISP
				if (dem.stageStruct
				|| (!dem.stageStruct && dem.dispersal == 1)
					// classical reproduction before dispersal
				|| (!dem.stageStruct && dem.dispersal == 0 && option ==2))
          // parturition occurring after dispersal
				{
#endif
					nj = (int)juvs.size();
					pCell = pPatch->getRandomCell();
					for (int j = 0; j < njuvs; j++) {
#if RS_CONTAIN
#if RSDEBUG
					// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
						juvs.push_back(new Individual(pCell,pPatch,0,0,0,stage,0.0,true,trfr.moveType));
#else
						juvs.push_back(new Individual(pCell,pPatch,0,0,0,stage,0.0,trfr.moveModel,trfr.moveType));
#endif
#else
#if PARTMIGRN
#if RSDEBUG
					// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
						juvs.push_back(new Individual(pSpecies,pCell,pPatch,0,0,0,0.0,true,trfr.moveType));
#else
						juvs.push_back(new Individual(pSpecies,pCell,pPatch,0,0,0,0.0,trfr.moveModel,trfr.moveType));
#endif
#else
#if RSDEBUG
					// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
						juvs.push_back(new Individual(pCell,pPatch,0,0,0,0.0,true,trfr.moveType));
#else
						juvs.push_back(new Individual(pCell,pPatch,0,0,0,0.0,trfr.moveModel,trfr.moveType));
#endif
#endif // PARTMIGRN 
#endif // RS_CONTAIN 
						nInds[0][0]++;
#if GOBYMODEL
						if (true)
#else
	#if SOCIALMODEL
						if (true)
	#else
						if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers)
	#endif
#endif
						{
						// juv inherits genome from parent (mother)
							juvs[nj+j]->setGenes(pSpecies,inds[i],0,resol);
						}
					}
#if BUTTERFLYDISP
				}
#endif
			}
		}
	}
	break;

case 1: // simple sexual model
case 2: // complex sexual model
	// count breeding females and males
	// add breeding males to list of potential fathers
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): case 1:"
//	<< " fec[1][0]=" << fec[1][0] << " fec[1][1]=" << fec[1][1]
//	<< endl;
#endif
	nfemales = nmales = 0;
	for (int i = 0; i < ninds; i++) {
		ind = inds[i]->getStats();
		if (ind.sex == 0 && fec[ind.stage][0] > 0.0) nfemales++;
		if (ind.sex == 1 && fec[ind.stage][1] > 0.0) {
			fathers.push_back(inds[i]);
			nmales++;
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): i=" << i << " nmales=" << nmales
//	<< " inds[i]=" << inds[i] << endl;
#endif
		}
	}
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): breeding nfemales=" << nfemales
//	<< " breeding nmales=" << nmales << endl;
#endif
#if BUTTERFLYDISP
	if (nfemales > 0 && (nmales > 0 || option == 2))
#else
	if (nfemales > 0 && nmales > 0)
#endif
	{ // population can breed
		if (dem.repType == 2) { // complex sexual model
			// calculate proportion of eligible females which breed
			propBreed = (2.0*dem.harem*nmales)/(nfemales+dem.harem*nmales);
			if (propBreed > 1.0) propBreed = 1.0;
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): harem=" << dem.harem
//	<< " nfemales=" << nfemales << " nmales=" << nmales << " propBreed=" << propBreed
//	<< endl;
#endif
		}
		else propBreed = 1.0;
		for (int i = 0; i < ninds; i++) {
			stage = inds[i]->breedingFem();
			if (stage > 0 && fec[stage][0] > 0.0) { // (potential) breeding female
				if (dem.stageStruct) {
					// determine whether she must miss current breeding attempt
					ind = inds[i]->getStats();
					if (ind.fallow >= sstruct.repInterval) {
						if (pRandom->Bernoulli(sstruct.probRep)) skipbreeding = false;
						else skipbreeding = true;
					}
					else skipbreeding = true; // cannot breed this time
				}
				else skipbreeding = false; // not structured - always breed
				if (skipbreeding) {
					inds[i]->incFallow();
				}
				else { // attempt to breed
					inds[i]->resetFallow();
					// NOTE: FOR COMPLEX SEXUAL MODEL, NO. OF FEMALES *ACTUALLY* BREEDING DOES NOT
					// NECESSARILY EQUAL THE EXPECTED NO. FROM EQN. 7 IN THE MANUAL...
					if (pRandom->Bernoulli(propBreed)) {
#if GOBYMODEL
						if (dem.stageStruct) {
							if (inds[i]->isAsocial()) {
								expected = fec[stage][0] * sstruct.asocF * exp(-sstruct.asocF*ddeffect[stage]);
							}
							else {
								expected = fec[stage][0] * exp(-ddeffect[stage]);
							}
						}
						else { // non-structured population
							// effect of sociality has not been defined (as at 3/6/2016)
							expected = fec[stage][0];
						}
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): 2222 ID=" << inds[i]->getId()
//	<< " stage=" << stage << " asocial=" << inds[i]->isAsocial()
//	<< " sstruct.asocF=" << sstruct.asocF
//	<< " fec[stage][0]=" << fec[stage][0] << " ddeffect=" << ddeffect[stage]
//	<< " expected=" << expected
//	<< endl;
#endif
#else
						expected = fec[stage][0]; // breeds
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): THIS LINE SHOULD NOT APPEAR FOR GOBY MODEL"
//	<< " expected=" << expected
//	<< endl;
#endif
#endif // GOBYMODEL
					}
					else expected = 0.0; // fails to breed
					if (expected <= 0.0) njuvs = 0;
					else njuvs = pRandom->Poisson(expected);
#if RSDEBUG
//DEBUGLOG << "Population::reproduction():"
//	<< " i " << i << " ID=" << inds[i]->getId() << " stage=" << stage
//	<< " expected=" << expected << " njuvs=" << njuvs << endl;
#endif
#if BUTTERFLYDISP
					if (njuvs > 0 || option == 1)
#else
					if (njuvs > 0)
#endif
					{
						nj = (int)juvs.size();
#if BUTTERFLYDISP
						if (!dem.stageStruct && dem.dispersal == 0 && option ==2) {
							// identity of father is already held by dispersed female
							father = inds[i]->getMate();
							if (father == 0) {
								// female was not already mated before dispersal
								// obtain mate locally if any males present
								if (nmales == 0) {
									njuvs = 0; // female fails to breed
								}
								else {
									int rrr = 0;
									if (nmales > 1) rrr = pRandom->IRandom(0,nmales-1);
									father = fathers[rrr];
								}
							}
						}
						else {
#endif
							// select father at random from breeding males ...
							int rrr = 0;
							if (nmales > 1) rrr = pRandom->IRandom(0,nmales-1);
							father = fathers[rrr];
#if BUTTERFLYDISP
						}
#endif
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): i = " << i << " rrr = " << rrr
//	<< " father = " << father << endl;
#endif
#if BUTTERFLYDISP
						if (dem.stageStruct
						|| (!dem.stageStruct && dem.dispersal == 1)
							// mating occurring before dispersal
						|| (!dem.stageStruct && dem.dispersal == 0 && option == 2))
							// parturition occurring after dispersal
						{
#endif
							pCell = pPatch->getRandomCell();
							for (int j = 0; j < njuvs; j++) {
#if RS_CONTAIN
#if RSDEBUG
								// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
								juvs.push_back(new Individual(pCell,pPatch,0,0,0,stage,dem.propMales,true,trfr.moveType));
#else
								juvs.push_back(new Individual(pCell,pPatch,0,0,0,stage,dem.propMales,trfr.moveModel,trfr.moveType));
#endif
#else
#if PARTMIGRN
#if RSDEBUG
								// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
								juvs.push_back(new Individual(pSpecies,pCell,pPatch,0,0,0,dem.propMales,true,trfr.moveType));
#else
								juvs.push_back(new Individual(pSpecies,pCell,pPatch,0,0,0,dem.propMales,trfr.moveModel,trfr.moveType));
#endif
#else
#if RSDEBUG
								// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
								juvs.push_back(new Individual(pCell,pPatch,0,0,0,dem.propMales,true,trfr.moveType));
#else
								juvs.push_back(new Individual(pCell,pPatch,0,0,0,dem.propMales,trfr.moveModel,trfr.moveType));
#endif
#endif // PARTMIGRN 
#endif // RS_CONTAIN 
								sex = juvs[nj+j]->getSex();
								nInds[0][sex]++;
#if GOBYMODEL
								if (true)
#else
#if SOCIALMODEL
								if (true)
#else
								if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers)
#endif // SOCIALMODEL
#endif // GOBYMODEL
								{
									// juv inherits genome from parents
									juvs[nj+j]->setGenes(pSpecies,inds[i],father,resol);
								}
							}
#if BUTTERFLYDISP
						}
						else { // !dem.stageStruct && dem.dispersal == 0 && option == 1
							// record female's mate prior to possible dispersal
							inds[i]->setMate(father);
						}
#endif
					}
				}
			}
		}
	}
	fathers.clear();
	break;

#if GROUPDISP
case 3: // hermaphrodite
	// add all breeding individuals to list of potential local 'fathers'
	nfemales = nmales = 0;
	for (int i = 0; i < ninds; i++) {
		ind = inds[i]->getStats();
		if (fec[ind.stage][0] > 0.0) {
			fathers.push_back(inds[i]);
			nmales++;
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): i=" << i << " nmales=" << nmales
//	<< " inds[i]=" << inds[i] << endl;
#endif
		}
	}
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): breeding 'males'=" << nmales
//	<< endl;
#endif
	if (nmales > 0) { // population can breed
		propBreed = 1.0; // NOT NEEDED HERE BUT RETAINED FOR CONSISTENCY
		for (int i = 0; i < ninds; i++) {
			stage = inds[i]->breedingFem();
			if (stage > 0 && fec[stage][0] > 0.0) { // (potential) breeding 'female'
				if (dem.stageStruct) {
					// determine whether she must miss current breeding attempt
					ind = inds[i]->getStats();
					if (ind.fallow >= sstruct.repInterval) {
						if (pRandom->Bernoulli(sstruct.probRep)) skipbreeding = false;
						else skipbreeding = true;
					}
					else skipbreeding = true; // cannot breed this time
				}
				else skipbreeding = false; // not structured - always breed
				if (skipbreeding) {
					inds[i]->incFallow();
				}
				else { // attempt to breed
					inds[i]->resetFallow();
					if (pRandom->Bernoulli(propBreed)) {
						expected = fec[stage][0]; // breeds
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): "
//	<< " expected=" << expected
//	<< endl;
#endif
					}
					else expected = 0.0; // fails to breed
					if (expected <= 0.0) njuvs = 0;
					else njuvs = pRandom->Poisson(expected);
#if RSDEBUG
//DEBUGLOG << "Population::reproduction():"
//	<< " i " << i << " ID=" << inds[i]->getId() << " stage=" << stage
//	<< " expected=" << expected << " njuvs=" << njuvs << endl;
#endif
					if (njuvs > 0) {
						nj = (int)juvs.size();
						int jj = 0; // index of viable juveniles in the brood
						bool selectnewfather = true;
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): i = " << i << " rrr = " << rrr
//	<< " father = " << father << endl;
#endif
						pCell = pPatch->getRandomCell();
						for (int j = 0; j < njuvs; j++) {
							bool viable = true;
							if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers) {
								if (selectnewfather) {
									// select 'father' at random from breeding 'males' ...
									int rrr = 0;
									if (dem.paternity == 2) { // pollen kernel method
										double randnum = pRandom->Random();
										if (randnum <= propLocal) { // pick local father
											if (nmales > 1) rrr = pRandom->IRandom(0,nmales-1);
											father = fathers[rrr];
										}
										else {
											if (randnum <= propLocal+propNghbr) {
												// pick neighbourhood father
												if (nfnbrhd > 0) {
													if (nfnbrhd > 1) rrr = pRandom->IRandom(0,nfnbrhd-1);
													else rrr = 0;
													father = (*pfnbrhd)[rrr];
												}
												else { // there are no neighbourhood fathers
													// treat as an Allee effect
													father = 0; viable = false;
												}
											}
											else { // pick global father
												if (nfglobal > 1) rrr = pRandom->IRandom(0,nfglobal-1);
												father = (*pfglobal)[rrr];
											}
										}
									}
									else { // pick local father
										if (nmales > 1) rrr = pRandom->IRandom(0,nmales-1);
										father = fathers[rrr];
									}
								}
								if (father == inds[i]) { // self-fertilization ...
									if (!dem.selfing) { // ... is not allowed
										// juvenile dies as zygote
										viable = false;
									}
								}
								if (dem.paternity == 0 || (dem.paternity == 1 && nmales == 1)) {
									// prevent further selection of new 'father'
									selectnewfather = false;
								}
							}
							if (viable) {
#if RSDEBUG
								// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
								juvs.push_back(new Individual(pCell,pPatch,0,0,0,0.0,true,trfr.moveType));
#else
								juvs.push_back(new Individual(pCell,pPatch,0,0,0,0.0,trfr.moveModel,trfr.moveType));
#endif
								sex = juvs[nj+jj]->getSex();
								nInds[0][sex]++;
								if (emig.indVar || trfr.indVar || sett.indVar || gen.neutralMarkers) {
									// juv inherits genome from parents
									juvs[nj+jj]->setGenes(pSpecies,inds[i],father,resol);
								}
								jj++;
							}
						}
					}
				}
			}
		}
	}
	fathers.clear();
	break;
#endif

} // end of switch (dem.repType)

#if RS_CONTAIN
nCulled = 0;
selectedForCull = false;
#endif // RS_CONTAIN 

#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): before reprodn. " << " inds.size() = " << inds.size()
//	<< endl;
#endif

// THIS MAY NOT BE CORRECT FOR MULTIPLE SPECIES IF THERE IS SOME FORM OF
// CROSS-SPECIES DENSITY-DEPENDENT FECUNDITY

#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): after reprodn. this = " << this
//	<< " juvs.size() = " << juvs.size() << " inds.size() = " << inds.size()
//	<< endl;
#endif

}

// Following reproduction of ALL species, add juveniles to the population prior to dispersal
void Population::fledge(void)
{
#if RSDEBUG
//DEBUGLOG << "Population::fledge(): this=" << this
//	<< " ninds=" << (int)inds.size()
//	<< " njuvs=" << (int)juvs.size()
//	<< endl;
#endif
demogrParams dem = pSpecies->getDemogr();

if (dem.stageStruct) { // juveniles are added to the individuals vector
	inds.insert(inds.end(),juvs.begin(),juvs.end());
//	nInds += nJuvs; nJuvs = 0;
}
else { // all adults die and juveniles replace adults
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		delete inds[i];
	}
	inds.clear();
	for (int sex = 0; sex < nSexes; sex++) {
		nInds[1][sex] = 0; // set count of adults to zero
	}
	inds = juvs;
}
juvs.clear();

}

// Determine which individuals will disperse
#if SEASONAL
void Population::emigration(float localK,short season)
#else
void Population::emigration(float localK)
#endif
{
int nsexes;
double disp,Pdisp,NK;
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
emigTraits eparams;

indStats ind;
#if RSDEBUG
//DEBUGLOG << "Population::emigration(): this=" << this
//	<< " nStages=" << sstruct.nStages
////	<< " emig.emigGenes[0]=" << emig.emigGenes[0]
////	<< " emig.emigGenes[1]=" << emig.emigGenes[1]
//	<< endl;
#endif

// to avoid division by zero, assume carrying capacity is at least one individual
// localK can be zero if there is a moving gradient or stochasticity in K
if (localK < 1.0) localK = 1.0;
NK = (float)totalPop() / localK;

int ninds = (int)inds.size();

#if SOCIALMODEL

// ADDITIONAL VARIABLES FOR PROBIS SOCIAL POLYMORPHISM MODEL
double emigAsocial,emigSocial;
//double Wasocial,Wsocial;
double Rasocial,Rsocial; // expected mean fecundity;
double fasocial,fsocial; // Allee effect below threshold
socialParams soc = pSpecies->getSocialParams();
/*
// calculate fitness of each morph as in Fogarty et al. (2011)
// NB Ts and Ta are expressed in terms of relative density (N/K) by Fogarty et al. (2011)
// and therefore must be rescaled to absolute density here
Wsocial = soc.rs * (1.0 - ninds/localK) * (ninds/(soc.Ts*localK) - 1.0);
Wasocial = soc.ra * (1.0 - ninds/(soc.asocK*localK)) * (ninds/(soc.Ta*localK) - 1.0);
emigSocial = 1.0 / (1.0 + exp(soc.alpha*Wsocial + log((1.0-soc.dK)/soc.dK)));
emigAsocial = 1.0 / (1.0 + exp(soc.alpha*Wasocial + log((1.0-soc.dK)/soc.dK)));
*/
// calculate expected reproduction of each morph
// NB Ts and Ta are expressed in terms of relative density (N/K) by Fogarty et al. (2011)
// and therefore must be rescaled to absolute density here

int ncells = pPatch->getNCells();
if ((float)ninds > soc.Ts*localK) fsocial = 1.0;
//else fsocial  = 1.0 / (1.0 + soc.Cs*(soc.Ts*localK-(float)ninds)*(soc.Ts*localK-(float)ninds));
//else fsocial  = 1.0 / (1.0 + exp(soc.cs + soc.bs*(soc.Ts*localK-(float)ninds)));
else fsocial  = 1.0 / (1.0 + exp(soc.cs + soc.bs*(soc.Ts*localK-(float)ninds)/(float)ncells));
if ((float)ninds > soc.Ta*localK) fasocial = 1.0;
//else fasocial = 1.0 / (1.0 + soc.Ca*(soc.Ta*localK-(float)ninds)*(soc.Ta*localK-(float)ninds));
//else fasocial  = 1.0 / (1.0 + exp(soc.ca + soc.ba*(soc.Ta*localK-(float)ninds)));
else fasocial  = 1.0 / (1.0 + exp(soc.ca + soc.ba*(soc.Ta*localK-(float)ninds)/(float)ncells));
Rsocial	= fsocial * dem.lambda / (1.0 + fabs(dem.lambda-1.0)*pow(((float)ninds/localK),dem.bc));
Rasocial	= fasocial * soc.asocRmax*dem.lambda
	/ (1.0 + fabs(soc.asocRmax*dem.lambda-1.0)*pow(((float)ninds/(soc.asocK*localK)),(soc.asocBc*dem.bc)));
emigSocial = 1.0 / (1.0 + exp(soc.alpha*(Rsocial-1.0) + log((1.0-soc.dK)/soc.dK)));
emigAsocial = 1.0 / (1.0 + exp(soc.alpha*(Rasocial-1.0) + log((1.0-soc.dK)/soc.dK)));

#if RSDEBUG
int nTotal,Nasocial,Nsocial;
//if (int(this)%50 == 0) {
	// looping through individuals twice is computationally inefficient but efficient in
	// terms of programming time!
	Nasocial = Nsocial = 0;
	for (int i = 0; i < ninds; i++) {
		if (inds[i]->isAsocial()) Nasocial++; else Nsocial++;
	}
//	DEBUGLOG << "Population::emigration(): ninds=" << ninds
//		<< " Nsocial=" << Nsocial << " Nasocial=" << Nasocial
//		<< " Wsocial=" << Wsocial << " Wasocial=" << Wasocial
//		<< " emigSocial=" << emigSocial << " emigAsocial=" << emigAsocial
//		<< endl;
	DEBUGLOG << "Population::emigration(): ninds=" << ninds
		<< " Nsocial=" << Nsocial << " Nasocial=" << Nasocial
		<< " fsocial=" << fsocial << " fasocial=" << fasocial
		<< " Rsocial=" << Rsocial << " Rasocial=" << Rasocial
		<< " emigSocial=" << emigSocial << " emigAsocial=" << emigAsocial
		<< endl;
//}
#endif

#endif // SOCIALMODEL

// set up local copy of emigration probability table
// used when there is no individual variability
// NB - IT IS DOUBTFUL THIS CONTRIBUTES ANY SUBSTANTIAL TIME SAVING
if (dem.repType == 0) nsexes = 1; else nsexes = 2;
double Pemig[NSTAGES][NSEXES];

#if PARTMIGRN
bool breeding = pSpecies->getBreeding(season);
bool nextbreeding;
if (season+1 < dem.nSeasons) nextbreeding = pSpecies->getBreeding(season+1);
else nextbreeding = pSpecies->getBreeding(0);
#endif // PARTMIGRN 

for (int stg = 0; stg < sstruct.nStages; stg++) {
	for (int sex = 0; sex < nsexes; sex++) {
		if (emig.indVar) Pemig[stg][sex] = 0.0;
		else {
			if (emig.densDep) {
				if (emig.sexDep) {
					if (emig.stgDep) {
						eparams = pSpecies->getEmigTraits(stg,sex);
					}
					else {
						eparams = pSpecies->getEmigTraits(0,sex);
					}
				}
				else { // !emig.sexDep
					if (emig.stgDep) {
						eparams = pSpecies->getEmigTraits(stg,0);
					}
					else {
						eparams = pSpecies->getEmigTraits(0,0);
					}
				}
				Pemig[stg][sex] = eparams.d0/(1.0+exp(-(NK - eparams.beta)*eparams.alpha));
#if RSDEBUG
//if (ppLand.patch_model) {
//	DEBUGLOG << "Population::emigration(): stg=" << stg << " sex=" << sex
//		<< " totalPop=" << totalPop() << " localK=" << localK << " NK=" << NK
//		<< " d0=" << eparams.d0 << " beta=" << eparams.beta << " alpha=" << eparams.alpha
//		<< " Pemig[stg][sex]=" << Pemig[stg][sex]
//		<< endl;
//}
#endif
			}
			else { // density-independent
				if (emig.sexDep) {
					if (emig.stgDep) {
						Pemig[stg][sex] = pSpecies->getEmigD0(stg,sex);
					}
					else { // !emig.stgDep
						Pemig[stg][sex] = pSpecies->getEmigD0(0,sex);
					}
				}
				else { // !emig.sexDep
					if (emig.stgDep) {
						Pemig[stg][sex] = pSpecies->getEmigD0(stg,0);
					}
					else { // !emig.stgDep
						Pemig[stg][sex] = pSpecies->getEmigD0(0,0);
					}
				}
			}
		} // end of !emig.indVar
#if RSDEBUG
//DEBUGLOG << "Population::emigration(): this=" << (int)this
//	<< " totalPop()=" << totalPop()
//	<< " Pemig[" << stg << "][" << sex << "]="<< Pemig[stg][sex] << endl;
#endif
	}
}

//#if GROUPDISP
//bool newgroup = true;
//int currentsize = 0;
//#endif // GROUPDISP

//#if PARTMIGRN
//double cumprop[7];
//cumprop[0] = 0.0;
//for (int i = 1; i < 7; i++) {
//	cumprop[i] = cumprop[i-1] + pSpecies->getPropDispMigrn(i);
//#if RSDEBUG
//DEBUGLOG << "Population::emigration(): i=" << i
//	<< " cumprop[i]=" << cumprop[i]
//	<< endl;
//#endif
//}
//#endif // PARTMIGRN
//

for (int i = 0; i < ninds; i++) {
	ind = inds[i]->getStats();          
#if SEASONAL
	inds[i]->resetPathSeason();
#if PARTMIGRN
//	if (ind.migrnstatus == 0) { // migration status not yet set
//		// TEMPORARY CODE TO SET MIGRATION STATUS
//		if (pRandom->Bernoulli(0.5)) {
//			if (pRandom->Bernoulli(0.5)) inds[i]->setMigrnStatus(1);
//			else {
//				if (pRandom->Bernoulli(0.5)) inds[i]->setMigrnStatus(2);
//				else inds[i]->setMigrnStatus(3);
//			}				
//		}
//		else {
//			if (pRandom->Bernoulli(0.5)) inds[i]->setMigrnStatus(4);
//			else {
//				if (pRandom->Bernoulli(0.5)) inds[i]->setMigrnStatus(5);
//				else inds[i]->setMigrnStatus(6);
//			}				
//		}
//		double r = pRandom->Random();
//		for (int m = 1; m < 7; m++) {
//			if (r < cumprop[m]) {
//				inds[i]->setMigrnStatus(m);
//				m = 7;
//			}
//		}
//	ind = inds[i]->getStats();
//	}
	if (ind.status < 6
	&& (ind.migrnstatus == 0 || ind.migrnstatus == 2 || ind.migrnstatus == 3
			|| (ind.migrnstatus == 4 && ind.status == 0)
			|| ind.migrnstatus == 5 || ind.migrnstatus == 6)
	)
#else
	if (ind.status < 6) 
#endif // PARTMIGRN
#else
	if (ind.status < 1) 
#endif // SEASONAL
	{
#if GOBYMODEL
		if (emig.densDep) {
			if (emig.sexDep) {
				if (emig.stgDep) {
					eparams = pSpecies->getEmigTraits(ind.stage,ind.sex);
				}
				else {
					eparams = pSpecies->getEmigTraits(0,ind.sex);
				}
			}
			else { // !emig.sexDep
				if (emig.stgDep) {
					eparams = pSpecies->getEmigTraits(ind.stage,0);
				}
				else {
					eparams = pSpecies->getEmigTraits(0,0);
				}
			}
			if (ind.asocial) {
				Pdisp = eparams.d0/(1.0+exp(-(NK - (eparams.beta/emig.asocD))*eparams.alpha));
			}
			else {
				Pdisp = eparams.d0/(1.0+exp(-(NK - eparams.beta)*eparams.alpha));
			}
#if RSDEBUG
//if (ppLand.patch_model) {
//	DEBUGLOG << "Population::emigration(): i=" << i << " sex=" << ind.sex << " stage=" << ind.stage
//		<< " totalPop=" << totalPop() << " localK=" << localK << " NK=" << NK
//		<< " d0=" << eparams.d0 << " beta=" << eparams.beta << " alpha=" << eparams.alpha
//		<< " Pdisp=" << Pdisp
//		<< endl;
//}
#endif
			}
			else { // density-independent
				if (emig.sexDep) {
					if (emig.stgDep) {
						Pdisp = Pemig[ind.stage][ind.sex];
					}
					else { // !emig.stgDep
						Pdisp = Pemig[0][ind.sex];
					}
				}
				else { // !emig.sexDep
					if (emig.stgDep) {
						Pdisp = Pemig[ind.stage][0];
					}
					else { // !emig.stgDep
						Pdisp = Pemig[0][0];
					}
				}
			}
#else // !GOBYMODEL
		if (emig.indVar) { // individual variability in emigration
			if (dem.stageStruct && ind.stage != emig.emigStage) {
				// emigration may not occur
				Pdisp = 0.0;
			}
			else { // non-structured or individual is in emigration stage
				eparams = inds[i]->getEmigTraits();
				if (emig.densDep) { // density-dependent
					NK = (float)totalPop() / localK;
					Pdisp = eparams.d0/(1.0+exp(-(NK - eparams.beta)*eparams.alpha));
				}
				else { // density-independent
					if (emig.sexDep) {
						Pdisp = Pemig[0][ind.sex] + eparams.d0;
					}
					else {
						Pdisp = Pemig[0][0] + eparams.d0;
					}
				}
			}
		} // end of individual variability
		else { // no individual variability
#if SOCIALMODEL
			if (ind.asocial) Pdisp = emigAsocial;
			else Pdisp = emigSocial;
#if RSDEBUG
DEBUGLOG << "Population::emigration(): i=" << i << " asocial=" << ind.asocial
	<< " emigAsocial=" << emigAsocial << " emigSocial=" << emigSocial
	<< " Pdisp=" << Pdisp << endl;
#endif
#else // !SOCIALMODEL

#if SEASONAL
#endif

			if (emig.densDep) {
				if (emig.sexDep) {
					if (emig.stgDep) {
						Pdisp = Pemig[ind.stage][ind.sex];
					}
					else {
						Pdisp = Pemig[0][ind.sex];
					}
				}
				else { // !emig.sexDep
					if (emig.stgDep) {
						Pdisp = Pemig[ind.stage][0];
					}
					else {
						Pdisp = Pemig[0][0];
					}
				}
#if RSDEBUG
//if (ppLand.patch_model) {
//	DEBUGLOG << "Population::emigration(): i=" << i << " sex=" << ind.sex << " stage=" << ind.stage
//		<< " totalPop=" << totalPop() << " localK=" << localK << " NK=" << NK
//		<< " Pdisp=" << Pdisp
//		<< endl;
//}
#endif
			}
			else { // density-independent
				if (emig.sexDep) {
					if (emig.stgDep) {
						Pdisp = Pemig[ind.stage][ind.sex];
					}
					else { // !emig.stgDep
						Pdisp = Pemig[0][ind.sex];
					}
				}
				else { // !emig.sexDep
					if (emig.stgDep) {
						Pdisp = Pemig[ind.stage][0];
					}
					else { // !emig.stgDep
						Pdisp = Pemig[0][0];
					}
				}
//#if GROUPDISP
//				if (emig.groupdisp) {
//					if (Pdisp > 0) {
//
//					}
//
//				}
//#endif // GROUPDISP
			}

			
#endif // SOCIALMODEL
		} // end of no individual variability
#endif // GOBYMODEL

#if PARTMIGRN
		disp = 1;
#else
		disp = pRandom->Bernoulli(Pdisp);
#endif  // PARTMIGRN 
#if RSDEBUG
//DEBUGLOG << "Population::emigration(): i=" << i << " sex=" << ind.sex << " stage=" << ind.stage
//	<< " Pdisp=" << Pdisp << " disp=" << disp << endl;
#endif

		if (disp == 1) { // emigrant
			inds[i]->setStatus(1);
#if PARTMIGRN
			inds[i]->setPrevPatch(pPatch);
			patchlist p;
			p.pPatch = pPatch; p.season = season; p.breeding = breeding; 
//			p.fixed = false;
//			p.fixed = true;
			inds[i]->addPatch(p);
			locn dummy;
			inds[i]->setGoal(dummy,1,nextbreeding);
#endif // PARTMIGRN 
		}
#if PARTMIGRN
//		else { // not an emigrant
//			if (ind.migrnstatus == 0) { // migration status not yet set
//				// TEMPORARY CODE TO SET MIGRATION STATUS
//				if (pRandom->Bernoulli(0.5)) inds[i]->setMigrnStatus(1);
//				else {
//					if (pRandom->Bernoulli(0.5)) inds[i]->setMigrnStatus(2);
//					else inds[i]->setMigrnStatus(3);
//					patchlist p; p.pPatch = pPatch; p.season = 0; p.breeding = p.fixed = true;
//					inds[i]->addPatch(p);
//				}
//			}
//		}
#endif  // PARTMIGRN 
	} // end of if (ind.status < 1) condition
} // end of for loop
}

// All individuals emigrate after patch destruction
void Population::allEmigrate(void) {
int ninds = (int)inds.size();
for (int i = 0; i < ninds; i++) {
	inds[i]->setStatus(1);
}
}

// If an Individual has been identified as an emigrant, remove it from the Population
disperser Population::extractDisperser(int ix) {
disperser d;
indStats ind = inds[ix]->getStats();
#if RSDEBUG
//if (ind.status > 0) {
//	DEBUGLOG << "Population::extractDisperser(): ix = " << ix << " inds[ix] = " << inds[ix]
//		<< " indId = " << inds[ix]->getId() << " ind.status = " << ind.status
//		<< endl;
//}
#endif
if (ind.status == 1) { // emigrant
	d.pInd = inds[ix]; d.yes = true;
	inds[ix] = 0;
	nInds[ind.stage][ind.sex]--;
}
else {
	d.pInd = NULL; d.yes = false;
}
return d;
}

#if GROUPDISP
// For a group identified as being in the matrix population:
// if it is a settler, return its new location and remove it from the current population
// otherwise, leave it in the matrix population for possible reporting before deletion
dispgroup Population::extractGroupSettler(int ix) {
dispgroup d;
Cell* pCell;
//Patch* pPatch;

indStats ind = groups[ix]->getStats();

pCell = groups[ix]->getLocn(1);
#if RSDEBUG
//DEBUGLOG << "Population::extractGroupSettler(): ix = " << ix << " inds[ix] = " << inds[ix]
//	<< " indId = " << inds[ix]->getId() << " ind.status = " << ind.status << endl;
#endif
d.pGroup = groups[ix];  d.pCell = pCell; d.groupsize = groups[ix]->getSize(); d.yes = false;
if (ind.status == 4 || ind.status == 5) { // settled
	d.yes = true;
//	nInds[ind.stage][0] -= groups[ix]->getSize();
	groups[ix] = 0;
}
return d;
}

//void Population::setGroupStatus(int j,short s) {
//if (j >= 0 && j < groups.size()) {
//	groups[j]->setStatus(s);
////	groups[j] = 0;
//}
//}

void Population::deleteGroup(int j) {
if (j >= 0 && j < groups.size()) {
	if (groups[j] != 0) {
		delete groups[j];
		groups[j] = 0;
	}
}
}

#endif

// For an individual identified as being in the matrix population:
// if it is a settler, return its new location and remove it from the current population
// otherwise, leave it in the matrix population for possible reporting before deletion
disperser Population::extractSettler(int ix) {
disperser d;
Cell* pCell;
//Patch* pPatch;

indStats ind = inds[ix]->getStats();

pCell = inds[ix]->getLocn(1);
#if RSDEBUG
//DEBUGLOG << "Population::extractSettler(): ix = " << ix << " inds[ix] = " << inds[ix]
//	<< " indId = " << inds[ix]->getId() << " ind.status = " << ind.status << endl;
#endif
d.pInd = inds[ix];  d.pCell = pCell; d.yes = false;
if (ind.status == 4 || ind.status == 5) { // settled
	d.yes = true;
	inds[ix] = 0;
	nInds[ind.stage][ind.sex]--;
}
return d;
}

// Add a specified individual to the new/current dispersal group
#if GROUPDISP
void Population::recruit(Patch *pPch,Individual *pInd,bool movt,
	short moveType,bool newgroup) {
Cell *pCell;
if (newgroup) {
	pCell = pInd->getLocn(1); // find current cell
	indStats ind = pInd->getStats();
	pGroup = new Group(pCell,pPch,ind.stage,movt,moveType);
	groups.push_back(pGroup);
}
pGroup->addMember(pInd);
pInd->setGroupId(pGroup->getGroupID());
indStats ind = pInd->getStats();
//nInds[ind.stage][0] += pGroup->getSize();
//nInds[ind.stage][0]++;
#if RSDEBUG
//DEBUGLOG << "Population::recruit(): GROUPDISP patchNum=" << pPatch->getPatchNum()
//	<< " newgroup=" << newgroup
//	<< " groupsize=" << pGroup->getSize()
//	<< " indId=" << pInd->getId()
//	<< " nInds[" << ind.stage << "][0]=" << nInds[ind.stage][0]
//	<< endl;
#endif
}
#endif
// Add a specified individual to the population
void Population::recruit(Individual *pInd) {
inds.push_back(pInd);
indStats ind = pInd->getStats();
nInds[ind.stage][ind.sex]++;
#if RSDEBUG
//DEBUGLOG << "Population::recruit(): patchNum=" << pPatch->getPatchNum()
//	<< " indId=" << pInd->getId()
//	<< " nInds[" << ind.stage << "][" << ind.sex << "]=" << nInds[ind.stage][ind.sex]
//	<< endl;
#endif
}

//---------------------------------------------------------------------------

#if GROUPDISP

// Transfer is run for populations in the matrix only
int Population::grouptransfer(Landscape *pLandscape,short landIx) {
int ndispersers = 0;
int disperser;
short othersex;
bool mateOK,densdepOK;
intptr patch,popn;
int patchnum;
double localK,popsize,settprob;
Patch *pPatch;
Cell *pCell;
indStats ind;
Population *pNewPopn;
locn newloc,nbrloc;

landData ppLand = pLandscape->getLandData();
short reptype = pSpecies->getRepType();
trfrRules trfr = pSpecies->getTrfr();
settleType settletype = pSpecies->getSettle();
settleRules sett;
settleTraits settDD;
settlePatch settle;
simParams sim = paramsSim->getSim();

// each individual takes one step
// for dispersal by kernel, this should be the only step taken
int ngroups = (int)groups.size();
#if RSDEBUG
DEBUGLOG << "Population::grouptransfer(): 0000: this=" << this
	<< " ngroups=" << ngroups << " ndispersers=" << ndispersers
	<< endl;
#endif
for (int i = 0; i < ngroups; i++) {
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): 1111: i = " << i << " ID = " << groups[i]->getId()
//	<< endl;
#endif
#if RSDEBUG
//pCell = groups[i]->getLocn(1);
//locn loc = pCell->getLocn();
//DEBUGLOG << "Population::grouptransfer(): 1112: i=" << i
//	<< " groupID=" << groups[i]->getGroupID() << " groupsize=" << groups[i]->getSize()
//	<< " before:" << " loc.x=" << loc.x << " loc.y=" << loc.y
//	<< endl;
#endif
	if (trfr.moveModel) {
		disperser = groups[i]->moveStep(pLandscape,pSpecies,landIx,sim.absorbing);
	}
	else {
		disperser = groups[i]->moveKernel(pLandscape,pSpecies,reptype,sim.absorbing);
	}
	ndispersers += disperser;
	if (disperser) {
		if (reptype > 0)
		{ // sexual species - record as potential settler in new patch
			if (groups[i]->getStatus() == 2)
			{ // group has found a patch
				pCell = groups[i]->getLocn(1);
				patch = pCell->getPatch();
				if (patch != 0) { // not no-data area
					pPatch = (Patch*)patch;
//					pPatch->incrPossSettler(pSpecies,inds[i]->getSex());
					pPatch->incrPossSettler(pSpecies,0);
				}
			}
		}
	}
#if RSDEBUG
//pCell = groups[i]->getLocn(1);
//locn nloc = pCell->getLocn();
//DEBUGLOG << "Population::grouptransfer(): 1113: i=" << i
//	<< " groupID=" << groups[i]->getGroupID() << " groupsize=" << groups[i]->getSize()
//	<< " status=" << groups[i]->getStatus() << " disperser=" << disperser
//	<< " after:" << " x=" << nloc.x << " y=" << nloc.y
//	<< endl;
#endif
}
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): 5555: ngroups=" << ngroups
//	<< " ndispersers=" << ndispersers << endl;
#endif

// each group which has reached a potential patch decides whether to settle
for (int i = 0; i < ngroups; i++) {
	ind = groups[i]->getStats();
	if (ind.sex == 0) othersex = 1; else othersex = 0;
	sett = pSpecies->getSettRules(ind.stage,ind.sex);
//	sett = pSpecies->getSettRules(0,0);
	if (ind.status == 2)
	{ // awaiting settlement
		pCell = groups[i]->getLocn(1);
#if RSDEBUG
//newloc = pCell->getLocn();
//DEBUGLOG << "Population::grouptransfer(): 6666: i=" << i << " ID=" << groups[i]->getGroupID()
//	<< " groupsize=" << groups[i]->getSize() << " status=" << ind.status
//	<< " stage=" << ind.stage << " sex=" << ind.sex
//	<< " pCell=" << pCell << " x=" << newloc.x << " y=" << newloc.y
////	<< " findMate=" << sett.findMate
////	<< " wait=" << sett.wait
////	<< " go2nbrLocn=" << sett.go2nbrLocn
//	<< endl;
#endif

		mateOK = false;
		if (sett.findMate) {
			// determine whether at least one individual of the opposite sex is present in the
			// new population
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): 7777: othersex = " << othersex
//	<< " this = " << this << " pNewPopn = " << pNewPopn << " popsize = " << popsize
//	<< endl;
#endif
			if (matePresent(pCell,othersex)) {
				mateOK = true;
			}
		}
		else { // no requirement to find a mate
			mateOK = true;
		}

		densdepOK = false;
		settle = groups[i]->getSettPatch();
		if (sett.densDep) {
			patch = pCell->getPatch();
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): 8880: i=" << i << " patch=" << patch
//	<< endl;
#endif
			if (patch != 0) { // not no-data area
				pPatch = (Patch*)patch;
				if (settle.settleStatus == 0
				||  settle.pSettPatch != pPatch)
				// note: second condition allows for having moved from one patch to another
				// adjacent one
				{
//					groups[i]->resetPathOut(); // reset steps out of patch to zero
					// determine whether settlement occurs in the (new) patch
					localK = (double)pPatch->getK();
					popn = pPatch->getPopn((intptr)pSpecies);
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): 8881: i=" << i << " pPatch=" << pPatch
//	<< " localK=" << localK << " popn=" << popn << endl;
#endif
					if (popn == 0) { // population has not been set up in the new patch
						popsize = 0.0;
					}
					else {
						pNewPopn = (Population*)popn;
						popsize = (double)pNewPopn->totalPop();
					}
					if (localK > 0.0) {
						// make settlement decision
						if (settletype.indVar) settDD = groups[i]->getSettTraits();
						else settDD = pSpecies->getSettTraits(ind.stage,ind.sex);
#if GOBYMODEL
						if (ind.asocial) {
							settprob = settDD.s0 / (1.0 + exp(-(popsize/localK
								- (double)(settDD.beta*settletype.betaSasoc)) * (double)(settDD.alpha*settletype.alphaSasoc)));
						}
						else {
							settprob = settDD.s0 /
								(1.0 + exp(-(popsize/localK - (double)settDD.beta) * (double)settDD.alpha));
						}
#else
						settprob = settDD.s0 /
							(1.0 + exp(-(popsize/localK - (double)settDD.beta) * (double)settDD.alpha));
#endif
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): 8888: i=" << i << " ind.stage=" << ind.stage
//	<< " this=" << this << " pNewPopn=" << pNewPopn << " popsize=" << popsize
//	<< " localK=" << localK << " alpha=" << settDD.alpha << " beta=" << settDD.beta
//	<< " settprob=" << settprob
//	<< endl;
#endif
						if (pRandom->Bernoulli(settprob)) { // settlement allowed
							densdepOK = true;
							settle.settleStatus = 2;
						}
						else { // settlement procluded
							settle.settleStatus = 1;
						}
						settle.pSettPatch = pPatch;
					}
					groups[i]->setSettPatch(settle);
				}
				else {
					if (settle.settleStatus == 2) { // previously allowed to settle
						densdepOK = true;
					}
				}
			}
		}
		else { // no density-dependent settlement
			densdepOK = true;
			settle.settleStatus = 2;
			settle.pSettPatch = pPatch;
			groups[i]->setSettPatch(settle);
		}

		if (mateOK && densdepOK) { // can recruit to patch
			ind.status = 4;
			ndispersers--;
		}
		else { // does not recruit
			if (trfr.moveModel) {
				ind.status = 1; // continue dispersing, unless ...
				// ... maximum steps has been exceeded
				pathSteps steps = groups[i]->getSteps();
				settleSteps settsteps = pSpecies->getSteps(ind.stage,ind.sex);
				if (steps.year >= settsteps.maxStepsYr) {
					ind.status = 3; // waits until next year
				}
				if (steps.total >= settsteps.maxSteps) {
					ind.status = 6; // dies
				}
			}
			else { // dispersal kernel
				if (sett.wait) {
					ind.status = 3; // wait until next dispersal event
				}
				else {
					ind.status = 6; // (dies unless a neighbouring cell is suitable)
				}
				ndispersers--;
			}
		}

	groups[i]->setStatus(ind.status);
	}

	if (!trfr.moveModel && sett.go2nbrLocn && (ind.status == 3 || ind.status == 6))
	{
		// for kernel-based transfer only ...
		// determine whether recruitment to a neighbouring cell is possible
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): neighbour cell search: sett.go2nbrLocn = " << sett.go2nbrLocn
//	<< " ind.status = " << ind.status
//	<< endl;
#endif
		pCell = groups[i]->getLocn(1);
		newloc = pCell->getLocn();
		vector <Cell*> nbrlist;
		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				if (dx !=0 || dy != 0) { //cell is not the current cell
					nbrloc.x = newloc.x + dx; nbrloc.y = newloc.y + dy;
					if (nbrloc.x >= 0 && nbrloc.x <= ppLand.maxX
					&&  nbrloc.y >= 0 && nbrloc.y <= ppLand.maxY) { // within landscape
						// add to list of potential neighbouring cells if suitable, etc.
						pCell = pLandscape->findCell(nbrloc.x,nbrloc.y);
						if (pCell != 0) { // not no-data area
							patch = pCell->getPatch();
							if (patch != 0) { // not no-data area
								pPatch = (Patch*)patch;
								patchnum = pPatch->getPatchNum();
								if (patchnum > 0 &&  pPatch != groups[i]->getNatalPatch())
								{ // not the matrix or natal patch
									if (pPatch->getK() > 0.0) { // suitable
										if (sett.findMate) {
											if (matePresent(pCell,othersex)) {
												nbrlist.push_back(pCell);
											}
										}
										else
											nbrlist.push_back(pCell);
									}
								}
							}
						}
					}
				}
			}
		}
		int listsize = (int)nbrlist.size();
		if (listsize > 0) { // there is at least one suitable neighbouring cell
			if (listsize == 1) {
				groups[i]->moveto(nbrlist[0]);
			}
			else { // select at random from the list
				int rrr = pRandom->IRandom(0,listsize-1);
				groups[i]->moveto(nbrlist[rrr]);
			}
		}
		// else list empty - do nothing - group retains its current location and status
	}
#if RSDEBUG
//newloc = pCell->getLocn();
//DEBUGLOG << "Population::grouptransfer(): 9990: i=" << i << " ID=" << groups[i]->getGroupID()
//	<< " groupsize=" << groups[i]->getSize() << " status=" << ind.status
//	<< " pCell=" << pCell << " x=" << newloc.x << " y=" << newloc.y
//	<< endl;
#endif
}
#if RSDEBUG
//DEBUGLOG << "Population::grouptransfer(): 9999: ninds=" << ninds
//	<< " ndispersers=" << ndispersers << endl;
#endif

return ndispersers;
}

// Delete dispersal groups once dispersal has finished
void Population::deleteGroups(void) {
int ngroups = (int)groups.size();
#if RSDEBUG
//DEBUGLOG << "Population::deleteGroups(): this=" << this
//	<< " ngroups=" << ngroups << endl;
#endif
for (int i = 0; i < ngroups; i++) { // all groups
	delete groups[i];
}
groups.clear();
}

#if PEDIGREE
void Population::outGroups(Pedigree *pPed,int rep,int yr,int gen,bool patchmodel) {
int ngroups = (int)groups.size();
#if RSDEBUG
DEBUGLOG << "Population::calcRelatedness(): this=" << this
	<< " ngroups=" << ngroups << endl;
#endif
for (int i = 0; i < ngroups; i++) { // all groups
	groups[i]->outGroup(pPed,rep,yr,gen,patchmodel);
}
}
#endif

#endif // GROUPDISP 

// Transfer is run for populations in the matrix only
#if SEASONAL || RS_RCPP
int Population::transfer(Landscape *pLandscape,short landIx,short nextseason) 
#else
int Population::transfer(Landscape *pLandscape,short landIx)
#endif
{
int ndispersers = 0;
int disperser;
short othersex;
bool mateOK,densdepOK;
intptr patch,popn;
int patchnum;
double localK,popsize,settprob;
Patch *pPatch = 0;
Cell *pCell = 0;
indStats ind;
Population *pNewPopn = 0;
locn newloc,nbrloc;

landData ppLand = pLandscape->getLandData();
short reptype = pSpecies->getRepType();
trfrRules trfr = pSpecies->getTrfr();
settleType settletype = pSpecies->getSettle();
settleRules sett;
settleTraits settDD;
settlePatch settle;
simParams sim = paramsSim->getSim();

#if SEASONAL
//bool breeding = pSpecies->getBreeding(season);
#endif // SEASONAL 

// each individual takes one step
// for dispersal by kernel, this should be the only step taken
int ninds = (int)inds.size();
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 0000: ninds = " << ninds
//	<< " ndispersers = " << ndispersers << endl;
#endif
for (int i = 0; i < ninds; i++) {
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 1111: i = " << i << " ID = " << inds[i]->getId()
//	<< endl;
#endif
	if (trfr.moveModel) {
#if RSDEBUG
//pCell = inds[i]->getLocn(1);
//locn loc = pCell->getLocn();
//DEBUGLOG << "Population::transfer(): 1112: i = " << i << " ID = " << inds[i]->getId()
//	<< " before:" << " x = " << loc.x << " y = " << loc.y
//	<< endl;
#endif
#if SEASONAL
		disperser = inds[i]->moveStep(pLandscape,pSpecies,landIx,nextseason,sim.absorbing);
#else
		disperser = inds[i]->moveStep(pLandscape,pSpecies,landIx,sim.absorbing);
#endif // SEASONAL 
#if RSDEBUG
//pCell = inds[i]->getLocn(1);
//newloc = pCell->getLocn();
//DEBUGLOG << "Population::transfer(): 1113: i = " << i << " ID = " << inds[i]->getId()
//	<< " after: " << " x = " << newloc.x << " y = " << newloc.y
//	<< endl;
#endif
	}
	else {
#if SEASONAL
		disperser = inds[i]->moveKernel(pLandscape,pSpecies,reptype,nextseason,sim.absorbing);
#else
		disperser = inds[i]->moveKernel(pLandscape,pSpecies,reptype,sim.absorbing);
#endif // SEASONAL 
	}
	ndispersers += disperser;
	if (disperser) {
		if (reptype > 0)
		{ // sexual species - record as potential settler in new patch
			if (inds[i]->getStatus() == 2)
			{ // disperser has found a patch
				pCell = inds[i]->getLocn(1);
				patch = pCell->getPatch();
				if (patch != 0) { // not no-data area
					pPatch = (Patch*)patch;
					pPatch->incrPossSettler(pSpecies,inds[i]->getSex());
				}
			}
		}
	}
}
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 5555: ninds=" << ninds
//	<< " ndispersers=" << ndispersers << endl;
#endif

// each individual which has reached a potential patch decides whether to settle
for (int i = 0; i < ninds; i++) {
	ind = inds[i]->getStats();
	if (ind.sex == 0) othersex = 1; else othersex = 0;
	if (settletype.stgDep) {
		if (settletype.sexDep) sett = pSpecies->getSettRules(ind.stage,ind.sex);
		else sett = pSpecies->getSettRules(ind.stage,0);
	}
	else {
		if (settletype.sexDep) sett = pSpecies->getSettRules(0,ind.sex);
		else sett = pSpecies->getSettRules(0,0);
	}
	if (ind.status == 2)
	{ // awaiting settlement
		pCell = inds[i]->getLocn(1);
		if (pCell == 0) {
			// this condition can occur in a patch-based model at the time of a dynamic landscape
			// change when there is a range restriction in place, since a patch can straddle the
			// range restriction and an individual forced to disperse upon patch removal could
			// start its trajectory beyond the boundary of the restrictyed range - such a model is 
			// not good practice, but the condition must be handled by killing the individual conceerned
			ind.status = 6;			
		}
		else {

#if RSDEBUG
//newloc = pCell->getLocn();
//DEBUGLOG << "Population::transfer(): 6666: i=" << i << " ID=" << inds[i]->getId()
//	<< " sex=" << ind.sex << " status=" << ind.status
//	<< " pCell=" << pCell << " x=" << newloc.x << " y=" << newloc.y
//	<< " findMate=" << sett.findMate
////	<< " wait=" << sett.wait
////	<< " go2nbrLocn=" << sett.go2nbrLocn
//	<< endl;
#endif

			mateOK = false;
			if (sett.findMate) {
				// determine whether at least one individual of the opposite sex is present in the
				// new population
#if SEASONAL
				if (matePresent(pCell,othersex,nextseason)) mateOK = true;
#else
				if (matePresent(pCell,othersex)) mateOK = true;
#endif // SEASONAL
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 7777: othersex=" << othersex
//	<< " this=" << this << " pNewPopn=" << pNewPopn << " popsize=" << popsize << " mateOK=" << mateOK
//	<< endl;
#endif
			}
			else { // no requirement to find a mate
				mateOK = true;
			}

			densdepOK = false;
			settle = inds[i]->getSettPatch();      
#if PARTMIGRN
			if (sett.densDep && nextseason == 1 && (ind.migrnstatus == 3 || ind.migrnstatus == 6)) 
#else
			if (sett.densDep) 
#endif // PARTMIGRN  
			{       
				patch = pCell->getPatch();
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 8880: i=" << i << " patch=" << patch
//	<< endl;
#endif
				if (patch != 0) { // not no-data area
					pPatch = (Patch*)patch;
					if (settle.settleStatus == 0
					||  settle.pSettPatch != pPatch)
					// note: second condition allows for having moved from one patch to another
					// adjacent one
					{
//						inds[i]->resetPathOut(); // reset steps out of patch to zero
						// determine whether settlement occurs in the (new) patch
#if SEASONAL
						localK = (double)pPatch->getK(nextseason);
#else
						localK = (double)pPatch->getK();
#endif // SEASONAL 
						popn = pPatch->getPopn((intptr)pSpecies);
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 8881: i=" << i << " patchNum=" << pPatch->getPatchNum()
//	<< " localK=" << localK << " popn=" << popn << endl;
#endif
						if (popn == 0) { // population has not been set up in the new patch
							popsize = 0.0;
						}
						else {
							pNewPopn = (Population*)popn;
							popsize = (double)pNewPopn->totalPop();
						}
						if (localK > 0.0) {
							// make settlement decision
							if (settletype.indVar) settDD = inds[i]->getSettTraits();
							else settDD = pSpecies->getSettTraits(ind.stage,ind.sex);
#if GOBYMODEL
							if (ind.asocial) {
								settprob = settDD.s0 / (1.0 + exp(-(popsize/localK
									- (double)(settDD.beta*settletype.betaSasoc)) * (double)(settDD.alpha*settletype.alphaSasoc)));
							}
							else {
								settprob = settDD.s0 /
									(1.0 + exp(-(popsize/localK - (double)settDD.beta) * (double)settDD.alpha));
							}
#else
							settprob = settDD.s0 /
								(1.0 + exp(-(popsize/localK - (double)settDD.beta) * (double)settDD.alpha));
#endif
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 8888: i=" << i << " ind.stage=" << ind.stage
//	<< " this=" << this << " pNewPopn=" << pNewPopn << " popsize=" << popsize
//	<< " localK=" << localK << " alpha=" << settDD.alpha << " beta=" << settDD.beta
//	<< " settprob=" << settprob
//	<< endl;
#endif
							if (pRandom->Bernoulli(settprob)) { // settlement allowed
								densdepOK = true;
								settle.settleStatus = 2;
							}
							else { // settlement procluded
								settle.settleStatus = 1;
							}
							settle.pSettPatch = pPatch;
						}
						inds[i]->setSettPatch(settle);
					}
					else {
						if (settle.settleStatus == 2) { // previously allowed to settle
							densdepOK = true;
						}
					}
				}
			}
			else { // no density-dependent settlement
				densdepOK = true;
				settle.settleStatus = 2;
				settle.pSettPatch = pPatch;
				inds[i]->setSettPatch(settle);
			}

			if (mateOK && densdepOK) { // can recruit to patch
				ind.status = 4;
#if SEASONAL
//			ind.migrnstatus = 6;
#endif
				ndispersers--;
			}
			else { // does not recruit
				if (trfr.moveModel) {
					ind.status = 1; // continue dispersing, unless ...
					// ... maximum steps has been exceeded
					pathSteps steps = inds[i]->getSteps();
					settleSteps settsteps = pSpecies->getSteps(ind.stage,ind.sex);
#if PARTMIGRN
					if (steps.season >= settsteps.maxStepsYr || steps.season >= settsteps.maxSteps) {        
						ind.status = 6; // dies
					}
#else
					if (steps.year >= settsteps.maxStepsYr) {
						ind.status = 3; // waits until next year
					}
					if (steps.total >= settsteps.maxSteps) {        
						ind.status = 6; // dies
					}
#endif
				}
				else { // dispersal kernel
					if (sett.wait) {
						ind.status = 3; // wait until next dispersal event
					}
					else {
						ind.status = 6; // (dies unless a neighbouring cell is suitable)
					}
					ndispersers--;
				}
			}
		}
	inds[i]->setStatus(ind.status);
#if SEASONAL
//	inds[i]->setMigrnStatus(ind.migrnstatus);
#endif
	}
#if RS_RCPP
	// write each individuals current movement step and status to paths file
	if (trfr.moveModel && sim.outPaths) {
		if(nextseason >= sim.outStartPaths && nextseason%sim.outIntPaths==0) {
				inds[i]->outMovePath(nextseason);
		}
	}
#endif

	if (!trfr.moveModel && sett.go2nbrLocn && (ind.status == 3 || ind.status == 6))
	{
		// for kernel-based transfer only ...
		// determine whether recruitment to a neighbouring cell is possible
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): neighbour cell search: sett.go2nbrLocn = " << sett.go2nbrLocn
//	<< " ind.status = " << ind.status
//	<< endl;
#endif
		pCell = inds[i]->getLocn(1);
		newloc = pCell->getLocn();
		vector <Cell*> nbrlist;
		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				if (dx !=0 || dy != 0) { //cell is not the current cell
					nbrloc.x = newloc.x + dx; nbrloc.y = newloc.y + dy;
					if (nbrloc.x >= 0 && nbrloc.x <= ppLand.maxX
					&&  nbrloc.y >= 0 && nbrloc.y <= ppLand.maxY) { // within landscape
						// add to list of potential neighbouring cells if suitable, etc.
						pCell = pLandscape->findCell(nbrloc.x,nbrloc.y);
						if (pCell != 0) { // not no-data area
							patch = pCell->getPatch();
							if (patch != 0) { // not no-data area
								pPatch = (Patch*)patch;
								patchnum = pPatch->getPatchNum();
								if (patchnum > 0 &&  pPatch != inds[i]->getNatalPatch())
								{ // not the matrix or natal patch
#if SEASONAL
									if (pPatch->getK(nextseason) > 0.0) 
#else
									if (pPatch->getK() > 0.0) 
#endif // SEASONAL 
									{ // suitable
										if (sett.findMate) {
#if SEASONAL
											if (matePresent(pCell,othersex,nextseason)) nbrlist.push_back(pCell);
#else
											if (matePresent(pCell,othersex)) nbrlist.push_back(pCell);
#endif // SEASONAL 
										}
										else
											nbrlist.push_back(pCell);
									}
								}
							}
						}
					}
				}
			}
		}
		int listsize = (int)nbrlist.size();
		if (listsize > 0) { // there is at least one suitable neighbouring cell
			if (listsize == 1) {
				inds[i]->moveto(nbrlist[0]);
			}
			else { // select at random from the list
				int rrr = pRandom->IRandom(0,listsize-1);
				inds[i]->moveto(nbrlist[rrr]);
			}
		}
		// else list empty - do nothing - individual retains its current location and status
	}
}
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): 9999: ninds = " << ninds
//	<< " ndispersers = " << ndispersers << endl;
#endif


return ndispersers;
}

// Determine whether there is a potential mate present in a patch which a potential
// settler has reached
#if SEASONAL
bool Population::matePresent(Cell *pCell,short othersex,short nextseason)
#else
bool Population::matePresent(Cell *pCell,short othersex)
#endif // SEASONAL 
{
int patch;
Patch *pPatch;
Population *pNewPopn;
int popsize = 0;
bool matefound = false;

patch = (int)pCell->getPatch();
if (patch != 0) {
	pPatch = (Patch*)pCell->getPatch();
	if (pPatch->getPatchNum() > 0) { // not the matrix patch
#if SEASONAL
		if (pPatch->getK(nextseason) > 0.0) 
#else
		if (pPatch->getK() > 0.0) 
#endif // SEASONAL 
		{ // suitable
			pNewPopn = (Population*)pPatch->getPopn((intptr)pSpecies);
			if (pNewPopn != 0) {
				// count members of other sex already resident in the patch
				for (int stg = 0; stg < nStages; stg++) {
					popsize += pNewPopn->nInds[stg][othersex];
				}
			}
			if (popsize < 1) {
				// add any potential settlers of the other sex
				popsize += pPatch->getPossSettlers(pSpecies,othersex);
			}
		}
	}
}
if (popsize > 0) matefound = true;
return matefound;
}

//---------------------------------------------------------------------------
// Determine survival and development and record in individual's status code
// Changes are NOT applied to the Population at this stage

// FOR MULTIPLE SPECIES, MAY NEED TO SEPARATE OUT THIS IDENTIFICATION STAGE,
// SO THAT IT CAN BE PERFORMED FOR ALL SPECIES BEFORE ANY UPDATING OF POPULATIONS

#if RS_CONTAIN
#if SEASONAL
void Population::survival0(float localK,short hab,short season,short option0,short option1)
#else
void Population::survival0(float localK,short hab,short option0,short option1)
#endif // SEASONAL
#else
#if SEASONAL
void Population::survival0(float localK,short season,short option0,short option1)
#else
#if SPATIALMORT
void Population::survival0(float localK,float mort,short option0,short option1)
#else
#if PEDIGREE
void Population::survival0(Pedigree *pPed,float localK,short option0,short option1)
#else
#if SPATIALDEMOG
void Population::survival0(float localK,short option0,short option1, std::vector <float> localDemoScaling)
#else
void Population::survival0(float localK,short option0,short option1)
#endif // SPATIALDEMOG
#endif // PEDIGREE
#endif // SPATIALMORT
#endif // SEASONAL
#endif // RS_CONTAIN 
{
// option0:	0 - stage 0 (juveniles) only
//					1 - all stages
//					2 - stage 1 and above (all non-juveniles)
// option1:	0 - development only (when survival is annual)
//	  	 		1 - development and survival
//	  	 		2 - survival only (when survival is annual)

#if RSDEBUG
//DEBUGLOG << "Population::survival0():"
//	<< " pSpecies=" << pSpecies << " this=" << this << " PatchNum=" << pPatch->getPatchNum()
//#if SEASONAL
//	<< " season=" << season 
//#endif // SEASONAL
//	<< " localK=" << localK << " option=" << option
//	<< endl;
#endif

densDepParams ddparams = pSpecies->getDensDep();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();

// get surrent population size
int ninds = (int)inds.size();
if (ninds == 0) return;

// set up local copies of species development and survival tables
int nsexes;
if (dem.repType == 0) nsexes = 1; else nsexes = 2;
float dev[NSTAGES][NSEXES];
float surv[NSTAGES][NSEXES];
short minAge[NSTAGES][NSEXES];
for (int stg = 0; stg < sstruct.nStages; stg++) {
	for (int sex = 0; sex < nsexes; sex++) {
		if (dem.stageStruct) {
			if (dem.repType == 1) { // simple sexual model
				// both sexes use development and survival recorded for females
#if RS_CONTAIN
#if SEASONAL
				dev[stg][sex]  = pSpecies->getDev(hab,season,stg,0);
				surv[stg][sex] = pSpecies->getSurv(hab,season,stg,0);
#else
				dev[stg][sex]  = pSpecies->getDev(hab,stg,0);
				surv[stg][sex] = pSpecies->getSurv(hab,stg,0);
#endif // SEASONAL 
#else
#if SEASONAL
				dev[stg][sex]  = pSpecies->getDev(season,stg,0);
				surv[stg][sex] = pSpecies->getSurv(season,stg,0);
#else
#if SPATIALDEMOG
				if (pSpecies->getDevSpatial() && pSpecies->getDevLayer(stg,0)>=0){
					dev[stg][sex] = pSpecies->getDev(stg,0)*localDemoScaling[pSpecies->getDevLayer(stg,0)];
				}
				else dev[stg][sex] = pSpecies->getDev(stg,0);
				if (pSpecies->getSurvSpatial() && pSpecies->getSurvLayer(stg,0)>=0){
					surv[stg][sex] = pSpecies->getSurv(stg,0)*localDemoScaling[pSpecies->getSurvLayer(stg,0)];
				}
				else surv[stg][sex] = pSpecies->getSurv(stg,0);
#else
				dev[stg][sex]  = pSpecies->getDev(stg,0);
				surv[stg][sex] = pSpecies->getSurv(stg,0);
#endif //SPATIALDEMOG
#endif // SEASONAL 
#endif // RS_CONTAIN 
				minAge[stg][sex] = pSpecies->getMinAge(stg,0);
			}
			else {
#if RS_CONTAIN
#if SEASONAL
				dev[stg][sex]  = pSpecies->getDev(hab,season,stg,sex);
				surv[stg][sex] = pSpecies->getSurv(hab,season,stg,sex);
#else
				dev[stg][sex]  = pSpecies->getDev(hab,stg,sex);
				surv[stg][sex] = pSpecies->getSurv(hab,stg,sex);
#endif // SEASONAL 
#else
#if SEASONAL
				dev[stg][sex]  = pSpecies->getDev(season,stg,sex);
				surv[stg][sex] = pSpecies->getSurv(season,stg,sex);
#else
#if SPATIALDEMOG
				if (pSpecies->getDevSpatial() && pSpecies->getDevLayer(stg,sex)>=0){
					dev[stg][sex] = pSpecies->getDev(stg,sex)*localDemoScaling[pSpecies->getDevLayer(stg,sex)];
				}
				else dev[stg][sex] = pSpecies->getDev(stg,sex);
				if (pSpecies->getSurvSpatial() && pSpecies->getSurvLayer(stg,sex)>=0){
									surv[stg][sex] = pSpecies->getSurv(stg,sex)*localDemoScaling[pSpecies->getSurvLayer(stg,sex)];
				}
				else surv[stg][sex] = pSpecies->getSurv(stg,sex);
#else
				dev[stg][sex]  = pSpecies->getDev(stg,sex);
				surv[stg][sex] = pSpecies->getSurv(stg,sex);
#endif // SPATIALDEMOG
#endif // SEASONAL 
#endif // RS_CONTAIN 
				minAge[stg][sex] = pSpecies->getMinAge(stg,sex);
			}
			if (option1 == 0) surv[stg][sex] = 1.0; // development only - all survive
			if (option1 == 2) dev[stg][sex] = 0.0;  // survival only - none develops
		}
		else { // non-structured population
			if (stg == 1) { // adults
				dev[stg][sex] = 0.0; surv[stg][sex] = 0.0; minAge[stg][sex] = 0;
			}
			else { // juveniles
				dev[stg][sex] = 1.0; surv[stg][sex] = 1.0; minAge[stg][sex] = 0;
			}
		}
#if RSDEBUG
//DEBUGLOG << "Population::survival0(): 1111 "
//	<< " dev[" << stg << "][" << sex << "] = " << dev[stg][sex]
//	<< " surv[" << stg << "][" << sex << "] = " << surv[stg][sex]
//	<< endl;
#endif
	}
}

if (dem.stageStruct) {
#if RSDEBUG
//	DEBUGLOG << "Population::survival0(): 2222 "
//		<< " ninds=" << ninds << " localK=" << localK
//		<< " effect of density dependence:" << endl;
#endif
	// apply density dependence in development and/or survival probabilities
	for (int stg = 0; stg < nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (option1 != 2 && sstruct.devDens && stg > 0) {
#if RSDEBUG
//	DEBUGLOG << "DD in DEVELOPMENT for stg=" << stg << " sex=" << sex << endl;
#endif
				// NB DD in development does NOT apply to juveniles,
				// which must develop to stage 1 if they survive
				float effect = 0.0;
				if (sstruct.devStageDens) { // stage-specific density dependence
					// NOTE: matrix entries represent effect of ROW on COLUMN 
					// AND males precede females
					float weight = 0.0;
					for (int effstg = 0; effstg < nStages; effstg++) {  
						for (int effsex = 0; effsex < nSexes; effsex++) {
							if (dem.repType == 2) {
								int rowincr,colincr;
								if (effsex == 0) rowincr = 1; else rowincr = 0; 
								if (sex == 0) colincr = 1; else colincr = 0; 
								weight = pSpecies->getDDwtDev(2*stg+colincr,2*effstg+rowincr);
							}
							else {
								weight = pSpecies->getDDwtDev(stg,effstg);								
							}
						effect += (float)nInds[effstg][effsex] * weight;
#if RSDEBUG
//	DEBUGLOG << " effstg=" << effstg << " effsex=" << effsex;
//	DEBUGLOG << " weight=" << weight << " effect=" << effect
//		<< endl;
#endif
						}
					}
				}
				else // not stage-specific
					effect = (float)totalPop();
				if (localK > 0.0)
					dev[stg][sex] *= exp(-(ddparams.devCoeff*effect)/localK);
#if RSDEBUG
//DEBUGLOG << "Population::survival0(): 2288 " << " effect=" << effect;
//if (localK > 0.0)
//	DEBUGLOG << " exp=" << exp(-(ddparams.devCoeff*effect)/localK);
//DEBUGLOG << " dev[" << stg << "][" << sex << "] = " << dev[stg][sex]
//	<< endl;
#endif
			} // end of if (sstruct.devDens && stg > 0)
			if (option1 != 0 && sstruct.survDens) {
#if RSDEBUG
//	DEBUGLOG << "DD in SURVIVAL for stg=" << stg << " sex=" << sex << endl;
#endif
				float effect = 0.0;
				if (sstruct.survStageDens) { // stage-specific density dependence
					// NOTE: matrix entries represent effect of ROW on COLUMN 
					// AND males precede females
					float weight = 0.0;
					for (int effstg = 0; effstg < nStages; effstg++) {
						for (int effsex = 0; effsex < nSexes; effsex++) {
							if (dem.repType == 2) {
								int rowincr,colincr;
								if (effsex == 0) rowincr = 1; else rowincr = 0; 
								if (sex == 0) colincr = 1; else colincr = 0; 
								weight = pSpecies->getDDwtSurv(2*stg+colincr,2*effstg+rowincr);
							}
							else {
								weight = pSpecies->getDDwtSurv(stg,effstg);								
							}
							effect += (float)nInds[effstg][effsex] * weight;
#if RSDEBUG
//	DEBUGLOG << " effstg=" << effstg << " effsex=" << effsex;
//	DEBUGLOG << " weight=" << weight << " effect=" << effect
//		<< endl;
#endif
						}
					}
				}
				else // not stage-specific
					effect = (float)totalPop();
				if (localK > 0.0)
					surv[stg][sex] *= exp(-(ddparams.survCoeff*effect)/localK);
#if RSDEBUG
//DEBUGLOG << "Population::survival0(): 3333 " << " effect=" << effect;
//if (localK > 0.0)
//	DEBUGLOG << " exp = " << exp(-(ddparams.survCoeff*effect)/localK);
//DEBUGLOG << " surv[" << stg << "][" << sex << "] = " << surv[stg][sex]
//	<< endl;
#endif
			} // end of if (sstruct.survDens)
		}
	}
}

// identify which individuals die or develop
#if RSDEBUG
//DEBUGLOG << "Population::survival0():"  << " ninds " << ninds
//	<< endl;
#endif
#if SEASONAL
extrmevent e = pSpecies->getExtreme(season);
bool affected = false;
if (e.prob > 0.0) {
	if (pRandom->Bernoulli((double)e.prob)) affected = true; // population is subject to extreme event	
}
#if RSDEBUG
//DEBUGLOG << "Population::survival0(): option=" << option 
//	<< " e.prob=" << e.prob << " affected=" << affected
//	<< endl;
#endif
#endif // SEASONAL 
for (int i = 0; i < ninds; i++) {
	indStats ind = inds[i]->getStats();       
#if RSDEBUG
//DEBUGLOG << "Population::survival0():"
//	<< " i=" << i << " indId=" << inds[i]->getId()
//	<< " stage=" << ind.stage << " status=" << ind.status << " sex=" << ind.sex
//#if PARTMIGRN
//	<< " migrnstatus=" << inds[i]->getMigrnStatus()
//#endif  // PARTMIGRN 
//	<< endl;
#endif
	if ((ind.stage == 0 && option0 < 2) || (ind.stage > 0 && option0 > 0)) {
		// condition for processing the stage is met...
		if (ind.status < 6) { // not already doomed
			double probsurv = surv[ind.stage][ind.sex];
			// does the individual survive?
			if (pRandom->Bernoulli(probsurv)) { // survives
				// does the individual develop?
				double probdev = dev[ind.stage][ind.sex];
				if (ind.stage < nStages-1) { // not final stage
#if RSDEBUG
//DEBUGLOG << "Population::survival0():"
//	<< " i=" << i << " indId=" << inds[i]->getId()
//	<< " age=" << ind.age << " minAge[stage+1]=" << minAge[ind.stage+1][ind.sex]
//	<< " probdev=" << probdev 
//	<< endl;
#endif
					if (ind.age >= minAge[ind.stage+1][ind.sex]) { // old enough to enter next stage
#if RSDEBUG
//DEBUGLOG << "Population::survival0():"
//	<< " i=" << i << " indId=" << inds[i]->getId() << " OLD ENOUGH"
//	<< endl;
#endif
						if (pRandom->Bernoulli(probdev)) {
							inds[i]->developing();
#if GROUPDISP
#if PEDIGREE
							// add new Individual to relationship table if it now has breeding status
							bool addind = false;
							if (dem.stageStruct) {
								stageParams sstruct = pSpecies->getStage();
								if (ind.stage+1 == sstruct.nStages-1) addind = true;
								// NOTE - FOR CONVENIENCE CURRENTLY ASSUMED TO BE ONLY THE FINAL STAGE,
								// BUT SHOULD CHECK FOR POSSIBLE NON-ZERO FECUNDITY OF EARLIER STAGE(S)
							}
							else addind = true;
							if (addind) {
								int posn = pPed->addInd(inds[i]);
								inds[i]->setMatPosn(posn);
								pPed->setRelMat(posn,posn,666.6);
//								for (i = 0; i < posn; i++) {
//									pSpecies->setRelMat(i,posn,0.0);
//									pSpecies->setRelMat(posn,i,0.0);				
//								}
							}
#endif
#endif
						}
					}
				}
			}
			else { // doomed to die
#if RSDEBUG
//DEBUGLOG << "Population::survival0():"
//	<< " i=" << i << " indId=" << inds[i]->getId() << " DIES"
//	<< endl;
#endif
				inds[i]->setStatus(8);
			}
		}
	}
#if SPATIALMORT
	// additional spatial mortality
	if (mort > 0.0 && ind.stage > 0 && inds[i]->getStatus() < 6) {
		if (pRandom->Bernoulli((double)mort)) inds[i]->setStatus(9); // dies
	}
#endif
#if SEASONAL
	// additional mortality due to extreme event
	if (affected && e.mort > 0.0 && inds[i]->getStatus() < 6) {
		if (pRandom->Bernoulli((double)e.mort)) inds[i]->setStatus(9); // dies		
	}
#endif // SEASONAL 
#if RSDEBUG
//ind = inds[i]->getStats();
//DEBUGLOG << "Population::survival0():"
//	<< " i = " << i << " ID = " << inds[i]->getId()
//	<< " stage = " << ind.stage << " status = " << ind.status
//	<< endl;
#endif
}
}

// Apply survival changes to the population
#if PEDIGREE
void Population::survival1(Pedigree *pPed)
#else
void Population::survival1(void)
#endif
{
#if PEDIGREE
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
#endif

int ninds = (int)inds.size();
#if RSDEBUG
//DEBUGLOG << "Population::survival1(): this=" << this
//	<< " patchNum=" << pPatch->getPatchNum() << " ninds=" << ninds
//	<< endl;
#endif
for (int i = 0; i < ninds; i++) {
	indStats ind = inds[i]->getStats();       
#if RSDEBUG
//DEBUGLOG << "Population::survival1(): i=" << i
//	<< " indId=" << inds[i]->getId() << " stage=" << ind.stage << " sex=" << ind.sex
//	<< " isDeveloping=" << ind.isDeveloping << " status=" << ind.status
//	<< endl;
#endif
	if (ind.status > 5) { // doomed to die
#if PEDIGREE
		// NB. FOR EXPANDTREE GROUP DISPERSAL, RETAIN ANY ADULTS WHICH HAVE DIED,
		// AS THEY MAY HAVE PRODUCED YOUNG WHICH MATURE SEVERAL YEARS
		// AFTER THE PARENT'S DEMISE 
//		if (ind.stage < sstruct.nStages-1 || ind.status < 8) 
		if (inds[i]->getMatPosn() < 0 || ind.status < 8) 
		{
			delete inds[i];
			inds[i] = NULL;
			nInds[ind.stage][ind.sex]--;				
		}
#else
		delete inds[i];
		inds[i] = NULL;
		nInds[ind.stage][ind.sex]--;
#endif
	}
	else {
		if (ind.isDeveloping) { // develops to next stage
			nInds[ind.stage][ind.sex]--;
			inds[i]->develop();
			nInds[ind.stage+1][ind.sex]++;
		}
	}
}
#if RSDEBUG
//DEBUGLOG << "Population::survival1(): this=" << this
//	<< " patchNum=" << pPatch->getPatchNum() << " completed individuals loop"
//	<< endl;
#endif

#if RSDEBUG
//for (int i = 0; i < inds.size(); i++) {
//DEBUGLOG << "Population::survival1():" << " inds[" << i << "] = " << inds[i] << endl;
//}
#endif

// remove pointers to dead individuals
clean();
#if RSDEBUG
//DEBUGLOG << "Population::survival1(): this=" << this
//	<< " patchNum=" << pPatch->getPatchNum() << " finished"
//	<< endl;
#endif

}

#if SEASONAL && PARTMIGRN
void Population::extremeEvent(float probMort) {
int ninds = (int)inds.size();
#if RSDEBUG
DEBUGLOG << "Population::extremeEvent(): this=" << this
	<< " patchNum=" << pPatch->getPatchNum() << " ninds=" << ninds << " probMort=" << probMort
	<< endl;
#endif
for (int i = 0; i < ninds; i++) {
	indStats ind = inds[i]->getStats(); 
	if (ind.status < 6) { // not already doomed
		if (pRandom->Bernoulli(probMort)) inds[i]->setStatus(8); // dies
	}
}	      
}
#endif // SEASONAL && PARTMIGRN 

void Population::ageIncrement(void) {
int ninds = (int)inds.size();
stageParams sstruct = pSpecies->getStage();
#if RSDEBUG
//DEBUGLOG << "Population::ageIncrement():" << " inds.size() = " << inds.size()
//	<< endl;
#endif
for (int i = 0; i < ninds; i++) {
	inds[i]->ageIncrement(sstruct.maxAge);
}
}

#if GROUPDISP
//---------------------------------------------------------------------------
// Remove zero pointers to dead or dispersed individuals
void Population::groupclean(void)
{
int ngroups =  (int)groups.size();
#if RSDEBUG
//DEBUGLOG << "Population::groupclean():" << " ngroups=" << ngroups
//	<< endl;
#endif
if (ngroups > 0) {
	std::vector <Group*> survivors; // all surviving individuals
	for (int i = 0; i < ngroups; i++) {
#if RSDEBUG
//DEBUGLOG << "Population::groupclean():" << " i=" << i << " groups[i]=" << groups[i]
//	<< endl;
#endif
		if (groups[i] != NULL) {
			survivors.push_back(groups[i]);
		}
	}
	groups.clear();
	groups = survivors;
#if RS_RCPP
	shuffle(groups.begin(), groups.end(), pRandom->getRNG() );
#else

#if !RSDEBUG
	// do not randomise individuals in RSDEBUG mode, as the function uses rand()
	// and therefore the randomisation will differ between identical runs of RS
#if RS_EMBARCADERO
	random_shuffle(groups.begin(), groups.end());
#else
	shuffle(groups.begin(), groups.end(), pRandom->getRNG() );
#endif
#endif // !RSDEBUG

#endif // RS_RCPP
}
}
#endif

//---------------------------------------------------------------------------
// Remove zero pointers to dead or dispersed individuals
void Population::clean(void)
{
int ninds =  (int)inds.size();
if (ninds > 0) {
//	sort (inds.begin(), inds.end());
//	reverse (inds.begin(), inds.end());
//
//	while (inds.size() > 0 && inds[inds.size()-1] == NULL ) {
//		inds.pop_back();
//	}
	// ALTERNATIVE METHOD: AVOIDS SLOW SORTING OF POPULATION
	std::vector <Individual*> survivors; // all surviving individuals
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != NULL) {
			survivors.push_back(inds[i]);
		}
	}
	inds.clear();
	inds = survivors;
#if RS_RCPP
	shuffle(inds.begin(), inds.end(), pRandom->getRNG() );
#else

#if !RSDEBUG
	// do not randomise individuals in RSDEBUG mode, as the function uses rand()
	// and therefore the randomisation will differ between identical runs of RS
#if RS_EMBARCADERO
	random_shuffle(inds.begin(), inds.end());
#else
	shuffle(inds.begin(), inds.end(), pRandom->getRNG() );
#endif
#endif // !RSDEBUG

#endif // RS_RCPP
}
}

//---------------------------------------------------------------------------
// Open population file and write header record
bool Population::outPopHeaders(int landNr,bool patchModel) {

if (landNr == -999) { // close file
	if (outPop.is_open()) outPop.close();
	outPop.clear();
	return true;
}

string name;
//landParams ppLand = pLandscape->getLandParams();
//envStochParams env = paramsStoch->getStoch();
simParams sim = paramsSim->getSim();
envGradParams grad = paramsGrad->getGradient();

// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
// ATTRIBUTES OF *ALL* SPECIES AS DETECTED AT MODEL LEVEL
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();

if (sim.batchMode)  {
	name = paramsSim->getDir(2)
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(landNr) + "_Pop.txt";
}
else{
	name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation) +"_Pop.txt";
}
outPop.open(name.c_str());
#if SEASONAL
outPop << "Rep\tYear\tSeason";
#else
outPop << "Rep\tYear\tRepSeason";
#endif
if (patchModel) outPop << "\tPatchID\tNcells";
else outPop << "\tx\ty";
// determine whether environmental data need be written for populations
bool writeEnv = false;
if (grad.gradient) writeEnv = true;
if (paramsStoch->envStoch()) writeEnv = true;
if (writeEnv) outPop << "\tEpsilon\tGradient\tLocal_K";
outPop << "\tSpecies\tNInd";
#if RSDEBUG
//DEBUGLOG << "Population::outPopHeaders(): this=" << this
//	<< " patchNum=" << pPatch->getPatchNum()
//	<< " totalPop()=" << totalPop()
//	<< " nStages=" << nStages << " nSexes=" << nSexes
//	<< endl;
#endif
if (dem.stageStruct) {
#if GROUPDISP
	if (dem.repType == 0 || dem.repType == 3)
#else
	if (dem.repType == 0)
#endif
	{
		for (int i = 1; i < sstruct.nStages; i++) outPop << "\tNInd_stage" << i ;
		outPop << "\tNJuvs";
	}
	else {
		for (int i = 1; i < sstruct.nStages; i++)
			outPop << "\tNfemales_stage" << i << "\tNmales_stage" << i ;
		outPop << "\tNJuvFemales\tNJuvMales";
	}
#if GOBYMODEL
	outPop << "\tNasocial\tNsocial";
#endif
}
else {
#if GROUPDISP
	if (dem.repType == 1 || dem.repType == 2) outPop << "\tNfemales\tNmales";
#else
	if (dem.repType != 0) outPop << "\tNfemales\tNmales";
#endif
#if SOCIALMODEL
	outPop << "\tNasocial\tNsocial";
#endif
}
outPop << endl;

return outPop.is_open();
}

//---------------------------------------------------------------------------
// Write record to population file
#if RS_ABC
void Population::outPopulation(int rep,int yr,int gen,float eps,
	bool patchModel,bool writeEnv,bool gradK,bool abcYear,ABCmaster *pABCmaster)
#else
void Population::outPopulation(int rep,int yr,int gen,float eps,
	bool patchModel,bool writeEnv,bool gradK)
#endif // RS_ABC 
{
Cell *pCell;

#if RSDEBUG
//DEBUGLOG << "Population::outPopulations(): this=" << this
//	<< " writeEnv " << (int)writeEnv
//	<< endl;
#endif

// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
// ATTRIBUTES OF *ALL* SPECIES AS DETECTED AT MODEL LEVEL
demogrParams dem = pSpecies->getDemogr();

#if RS_ABC
simParams sim = paramsSim->getSim();
#endif
popStats p;

outPop << rep << "\t" << yr << "\t" << gen;
if (patchModel) {
	outPop << "\t" << pPatch->getPatchNum();
	outPop << "\t" << pPatch->getNCells();
}
else {
	locn loc = pPatch->getCellLocn(0);
	outPop << "\t" << loc.x << "\t" << loc.y;
}
if (writeEnv) {
	if (pPatch->getPatchNum() == 0) { // matrix
		outPop << "\t0\t0\t0";
	}
	else {
#if SEASONAL
		float k = pPatch->getK(gen);
#else
		float k = pPatch->getK();
#endif // SEASONAL 
		float envval = 0.0;
		pCell = pPatch->getRandomCell();
		if (pCell != 0) envval = pCell->getEnvVal();
		outPop << "\t" << eps << "\t" << envval << "\t" << k;
	}
}
#if RSDEBUG
//DEBUGLOG << "Population::outPopulation(): this=" << this
//	<< " patchNum=" << pPatch->getPatchNum()
//	<< " totalPop()=" << totalPop()
//	<< " nStages=" << nStages << " nSexes=" << nSexes
//	<< endl;
#endif
outPop << "\t" << pSpecies->getSpNum();
if (dem.stageStruct) {
#if RS_CONTAIN
	// here we can set the habitat index parameter to 0 as we are not concerned 
	// whether or not the poplation includes breeders
	p = getStats(0);
#else
#if SPATIALDEMOG
	p = getStats(pPatch->getDemoScaling());
#else
	p = getStats();
#endif // SPATIALDEMOG
#endif // RS_CONTAIN 
	outPop << "\t" << p.nNonJuvs;
	// non-juvenile stage totals from permanent array
	for (int stg = 1; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			outPop << "\t" << nInds[stg][sex];
		}
	}
	// juveniles from permanent array
	for (int sex = 0; sex < nSexes; sex++) {
		outPop << "\t" << nInds[0][sex];
	}
#if GOBYMODEL
	outPop << "\t" << p.nAsocial << "\t" << p.nSocial;
#endif
}
else { // non-structured population
	outPop << "\t" << totalPop();
#if GROUPDISP
	if (dem.repType == 1 || dem.repType == 2)
#else
	if (dem.repType != 0)
#endif
	{ // sexual model
		outPop << "\t" << nInds[1][0] << "\t" << nInds[1][1];
	}
#if SOCIALMODEL
	popStats p = getStats();
	outPop << "\t" << p.nAsocial << "\t" << p.nSocial;
#endif
}
outPop << endl;

/*
#if RS_ABC
obsdata obs;
if (abcYear) {
	int nobs = (int)pABCmaster->NObs();
	for (int i = 0; i < nobs; i++) {
		obs = pABCmaster->getObsData(i);
#if RSDEBUG
//DEBUGLOG << "Population::outPopulation(): this=" << this << " i=" << i << " yr=" << yr
//	<< " obs.year=" << obs.year << " obs.type=" << obs.type << " obs.name=" << obs.name
//	<< " obs.x=" << obs.x << " obs.y=" << obs.y
//	<< endl;
#endif
		if (obs.year == yr && obs.type == 2) {
			if (obs.name == "NInds" || obs.name == "Occupied") {
				bool match = false;
				if (patchModel) {
					if (obs.x == pPatch->getPatchNum()) {
						match = true;
#if RSDEBUG
//DEBUGLOG << "Population::outPopulation(): i=" << i << " PROCESS Population NInds"
//	<< " obs.id=" << obs.id << " obs.value=" << obs.value << " obs.x=" << obs.x
//	<< " pPatch->PatchNum()=" << pPatch->getPatchNum()
//	<< " totalPop()=" << totalPop() << " p.nNonJuvs=" << p.nNonJuvs
//	<< endl;
#endif
					}
				}
				else {
					locn loc = pPatch->getCentroid();
					if (obs.x == loc.x && obs.y == loc.y) {
						match = true;
#if RSDEBUG
DEBUGLOG << "Population::outPopulation(): i=" << i << " PROCESS Population NInds"
	<< " obs.id=" << obs.id << " obs.value=" << obs.value << " obs.x="
	<< obs.x << " obs.y=" << obs.y << " loc.x=" << loc.x << " loc.y=" << loc.y
	<< " totalPop()=" << totalPop() << " p.nNonJuvs=" << p.nNonJuvs
	<< endl;
#endif
					}
				}
				if (match) {
					if (obs.name == "NInds") {
						if (dem.stageStruct)
							pABCmaster->AddNewPred(sim.simulation,obs.id,rep,obs.value,p.nNonJuvs,obs.weight);
						else
							pABCmaster->AddNewPred(sim.simulation,obs.id,rep,obs.value,totalPop(),obs.weight);
					}
					else { // obs.name == "Occupied"
						pABCmaster->AddNewPred(sim.simulation,obs.id,rep,obs.value,p.breeding,obs.weight);
					}
				}
			}
		}
	}
}
#endif // ABC
*/
}

//---------------------------------------------------------------------------

#if RS_CONTAIN

// Open population file and write header record
bool Population::outCullHeaders(Landscape *pLand,int landNr,bool patchModel) {

if (landNr == -999) { // close file
	if (outCull.is_open()) outCull.close();
	outCull.clear();
	return true;
}

string name;
landParams ppLand = pLand->getLandParams();
//envStochParams env = paramsStoch->getStoch();
simParams sim = paramsSim->getSim();
//envGradParams grad = paramsGrad->getGradient();

//demogrParams dem = pSpecies->getDemogr();
//stageParams sstruct = pSpecies->getStage();

if (sim.batchMode)  {
	name = paramsSim->getDir(2)
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(landNr) + "_Cull.txt";
}
else{
	name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation) +"_Cull.txt";
}
outCull.open(name.c_str());
#if SEASONAL
outCull << "Rep\tYear\tSeason";
#else
outCull << "Rep\tYear\tRepSeason";
#endif
if (patchModel) outCull << "\tPatchID";
else {
	outCull << "\tx\ty\tHabitat";
	if (ppLand.rasterType == 1) outCull << "1";
}	
outCull << "\tSpecies\tSelected\tNculled";
outCull << endl;

return outPop.is_open();
}

// Write record to cull file
void Population::outCullData(Landscape *pLand,int rep,int yr,int gen,
	bool patchModel)
{
Cell *pCell;
int habitat,ix;
float habprop;
landParams ppLand = pLand->getLandParams();  
//demogrParams dem = pSpecies->getDemogr();
//stageParams sstruct = pSpecies->getStage();
//popStats p;

#if RSDEBUG
//DEBUGLOG << "Population::outCullData(): this=" << this
//	<< endl;
#endif

outCull << rep << "\t" << yr << "\t" << gen;
if (patchModel) {
	outCull << "\t" << pPatch->getPatchNum();
}
else {
	locn loc = pPatch->getCellLocn(0);
	outCull << "\t" << loc.x << "\t" << loc.y;
	pCell = pPatch->getRandomCell();
#if RSDEBUG
//DEBUGLOG << "Population::outCullData(): this=" << this
//	<< " x=" << loc.x << " y=" << loc.y << " pCell=" << pCell
//	<< endl;
#endif
	habitat = -9;
	if (pCell == 0) { 
		outCull << "\t" << habitat;
	}
	else {
		if (ppLand.rasterType == 0) {
			ix = pCell->getHabIndex(0); // CURRENTLY INITIAL HABITAT ONLY 
			habitat = pLand->getHabCode(ix); 			
			outCull << "\t" << habitat;
		}
		else {
			habprop = pCell->getHabitat(0);
			outCull << "\t" << habprop;
		}  
	}
}
#if RSDEBUG
//DEBUGLOG << "Population::outCullData(): this=" << this
//	<< " patchNum=" << pPatch->getPatchNum()
//	<< " totalPop()=" << totalPop()
//	<< " nStages=" << nStages << " nSexes=" << nSexes
//	<< endl;
#endif
outCull << "\t" << pSpecies->getSpNum();
if (selectedForCull) outCull << "\t1"; else outCull << "\t0";
outCull << "\t" << nCulled;
outCull << endl;

}

#endif // RS_CONTAIN 

//---------------------------------------------------------------------------
// Open individuals file and write header record
void Population::outIndsHeaders(int rep,int landNr,bool patchModel)
{

if (landNr == -999) { // close file
	if (outInds.is_open()) {
		outInds.close(); outInds.clear();
	}
	return;
}

string name;
demogrParams dem = pSpecies->getDemogr();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
simParams sim = paramsSim->getSim();

if (sim.batchMode) {
	name = paramsSim->getDir(2)
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation)
		+ "_Land" + Int2Str(landNr) + "_Rep" + Int2Str(rep) + "_Inds.txt";
}
else {
	name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation)
		+ "_Rep" + Int2Str(rep) +"_Inds.txt";
}
outInds.open(name.c_str());

#if GROUPDISP
outInds << "Rep\tYear\tRepSeason\tSpecies\tIndID\tMotherID";
if (dem.repType > 0) outInds << "\tFatherID";
outInds << "\tGroupID\tStatus";
#else
#if SEASONAL
#if PARTMIGRN
outInds << "Rep\tYear\tSeason\tSpecies\tIndID\tStatus\tMigrnStatus";
#else
outInds << "Rep\tYear\tSeason\tSpecies\tIndID\tStatus";
#endif // PARTMIGRN 
#else
outInds << "Rep\tYear\tRepSeason\tSpecies\tIndID\tStatus";
#endif // SEASONAL 
#endif // GROUPDISP 
if (patchModel) outInds << "\tNatal_patch\tPatchID";
else outInds << "\tNatal_X\tNatal_Y\tX\tY";
#if GOBYMODEL
outInds <<"\tAsocial";
#endif
#if SOCIALMODEL
outInds <<"\tAsocial";
#endif
if (dem.repType != 0) outInds << "\tSex";
if (dem.stageStruct) outInds << "\tAge\tStage";
if (emig.indVar) {
	if (emig.densDep) outInds << "\tD0\tAlpha\tBeta";
	else outInds << "\tEP";
}
if (trfr.indVar) {
	if (trfr.moveModel) {
		if (trfr.moveType == 1) { // SMS
			outInds << "\tDP\tGB\tAlphaDB\tBetaDB";
		}
		if (trfr.moveType == 2) { // CRW
			outInds << "\tStepLength\tRho";
		}
	}
	else { // kernel
		outInds << "\tMeanDistI";
#if RS_CONTAIN
		if (trfr.kernType == 1) outInds << "\tMeanDistII\tPKernelI";
#else
		if (trfr.twinKern) outInds << "\tMeanDistII\tPKernelI";
#endif // RS_CONTAIN 
	}
}
if (sett.indVar) {
	outInds << "\tS0\tAlphaS\tBetaS";
}
outInds << "\tDistMoved";
#if SEASONAL
outInds << "\tPrevX\tPrevY\tX\tY";
#endif
#if RSDEBUG
// ALWAYS WRITE NO. OF STEPS
outInds << "\tNsteps";
#else
if (trfr.moveModel) outInds << "\tNsteps";
#endif
#if SEASONAL
#if PARTMIGRN
outInds << "\tP0\tP1\tP2\tP3\tP4\tP5";
#endif // PARTMIGRN 
#endif // SEASONAL
outInds << endl;
}

//---------------------------------------------------------------------------
// Write records to individuals file
void Population::outIndividual(Landscape *pLandscape,int rep,int yr,int gen,
	int patchNum)
{
//int x, y, p_id;
bool writeInd;
pathSteps steps;
Cell *pCell;

landParams ppLand = pLandscape->getLandParams();
//landOrigin lim = pLandscape->getOrigin();
demogrParams dem = pSpecies->getDemogr();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
short spNum = pSpecies->getSpNum();

int ninds = (int)inds.size();

for (int i = 0; i < ninds; i++) {
	indStats ind = inds[i]->getStats();
	if (yr == -1) { // write all initialised individuals
		writeInd = true;
#if SEASONAL
		outInds << rep << "\t" << yr << "\t" << dem.nSeasons-1;
#else
		outInds << rep << "\t" << yr << "\t" << dem.repSeasons-1;
#endif
	}
	else {
		if (dem.stageStruct && gen < 0) { // write status 9 individuals only
			if (ind.status == 9) {
				writeInd = true;
#if SEASONAL
				outInds << rep << "\t" << yr << "\t" << dem.nSeasons-1;
#else
				outInds << rep << "\t" << yr << "\t" << dem.repSeasons-1;
#endif
			}
			else writeInd = false;
		}
		else {
			writeInd = true;
			outInds << rep << "\t" << yr << "\t" << gen;
		}
	}
	if (writeInd) {
		outInds << "\t" << spNum << "\t" << inds[i]->getId();
#if GROUPDISP
		outInds << "\t" << inds[i]->getParentId(0);
		if (dem.repType > 0) outInds << "\t" << inds[i]->getParentId(1);
//		Individual *pParent = inds[i]->getMother();
//		if (pParent == 0) outInds << "\t-1";
//		else outInds << "\t" << pParent->getId();
//		if (dem.repType > 0) {
//			pParent = inds[i]->getFather();
//			if (pParent == 0) outInds << "\t-1";
//			else outInds << "\t" << pParent->getId();
//		}	
		outInds << "\t" << inds[i]->getGroupId();
#endif
#if SEASONAL
#if PARTMIGRN
		outInds << "\t" << ind.status << "\t" << ind.migrnstatus;
#else
		outInds << "\t" << ind.status;
#endif // PARTMIGRN 
#else
		if (dem.stageStruct) outInds << "\t" << ind.status;
		else { // non-structured population
			outInds << "\t" << ind.status;
		}
#endif
		pCell = inds[i]->getLocn(1);
		locn loc;
		if (pCell == 0) loc.x = loc.y = -1; // beyond boundary or in no-data cell
		else loc = pCell->getLocn();
		pCell = inds[i]->getLocn(0);
		locn natalloc = pCell->getLocn();
//#if SEASONAL
//		pCell = inds[i]->getLocn(2);
//		locn prevloc = pCell->getLocn();
//#endif
		if (ppLand.patchModel) {
			outInds << "\t" << inds[i]->getNatalPatch()->getPatchNum();
			if (loc.x == -1) outInds << "\t-1";
			else outInds << "\t" << patchNum;
		}
		else { // cell-based model
			// EITHER write co-ordinates in cell units ...
			outInds << "\t" << (float)natalloc.x << "\t" << natalloc.y;
			outInds  << "\t" << (float)loc.x << "\t" << (float)loc.y ;
			// ... OR write co-ordinates in real-world units
//		outInds << "\t" << (float)natalloc.x * (float)ppLand.resol + (float)lim.minEast
//			 << "\t" << natalloc.y * (float)ppLand.resol + (float)lim.minNorth;
//		outInds  << "\t" << (float)loc.x * (float)ppLand.resol + (float)lim.minEast
//			<< "\t" << (float)loc.y * (float)ppLand.resol + (float)lim.minNorth;
		}
#if GOBYMODEL
		outInds <<"\t" << (int)ind.asocial;
#endif
#if SOCIALMODEL
		outInds <<"\t" << (int)ind.asocial;
#endif
		if (dem.repType != 0) outInds <<"\t" << ind.sex;
		if (dem.stageStruct) outInds <<"\t" << ind.age <<"\t"<< ind.stage;

		if (emig.indVar) {
			emigTraits e = inds[i]->getEmigTraits();
			if (emig.densDep) {
				outInds << "\t" << e.d0 << "\t" << e.alpha << "\t" << e.beta;
			}
			else {
				outInds << "\t" << e.d0;
			}
		} // end of if (emig.indVar)

		if (trfr.indVar) {
			if (trfr.moveModel) {
				if (trfr.moveType == 1) { // SMS
					trfrSMSTraits s = inds[i]->getSMSTraits();
					outInds << "\t" << s.dp << "\t" << s.gb;
					outInds << "\t" << s.alphaDB << "\t" << s.betaDB;
				} // end of SMS
				if (trfr.moveType == 2) { // CRW
					trfrCRWTraits c = inds[i]->getCRWTraits();
					outInds << "\t" << c.stepLength << "\t" << c.rho;
#if RSDEBUG
//DEBUGLOG << "Population::outIndividual():"
//	<< " patchNum=" << patchNum << " i=" << i << " ID=" << inds[i]->getId()
//	<< " nTrfrGenes=" << nTrfrGenes << " loc[0][0].allele[0]=" << loc[0][0].allele[0]
//	<< endl;
#endif
				} // end of CRW
			}
			else { // kernel
				trfrKernTraits k = inds[i]->getKernTraits();
#if RS_CONTAIN
				if (trfr.kernType == 1) 
#else
				if (trfr.twinKern) 
#endif // RS_CONTAIN 
				{
					outInds << "\t" << k.meanDist1 << "\t" << k.meanDist2 << "\t" << k.probKern1;
				}
				else {
					outInds << "\t" << k.meanDist1;
				}
			}
		}

		if (sett.indVar) {
			settleTraits s = inds[i]->getSettTraits();
			outInds << "\t" << s.s0 << "\t" << s.alpha << "\t" << s.beta;
		}

		// distance moved (metres)
		if (loc.x == -1) outInds << "\t-1";
		else {
#if SEASONAL
			pCell = inds[i]->getLocn(2);
			locn prevloc = pCell->getLocn();
			steps = inds[i]->getSteps();
			float d = 0.0;
			if (steps.season > 0) {
				d = ppLand.resol * sqrt((float)((prevloc.x-loc.x)*(prevloc.x-loc.x)
																					+ (prevloc.y-loc.y)*(prevloc.y-loc.y)));
			}
			outInds << "\t" << d;
			outInds << "\t" << prevloc.x << "\t" << prevloc.y;
			outInds << "\t" << loc.x << "\t" << loc.y;
#else
			float d = ppLand.resol * sqrt((float)((natalloc.x-loc.x)*(natalloc.x-loc.x)
																					+ (natalloc.y-loc.y)*(natalloc.y-loc.y)));
			outInds << "\t" << d;
#endif
		}
#if SEASONAL
//		steps = inds[i]->getSteps();
		outInds << "\t" << steps.season;
#if PARTMIGRN
		for (int j = 0; j < 6; j++) {
			patchlist p = inds[i]->getPatch(j);
			if (p.pPatch == 0) outInds << "\t-9";
			else outInds << "\t" << p.pPatch->getPatchNum() << "(" << p.season << ")/" << p.breeding << "/" << p.fixed;  
		}
#endif // PARTMIGRN 
#else
#if RSDEBUG
		// ALWAYS WRITE NO. OF STEPS
		steps = inds[i]->getSteps();
		outInds << "\t" << steps.year;
#else
		if (trfr.moveModel) {
			steps = inds[i]->getSteps();
			outInds << "\t" << steps.year;
		}
#endif
#endif // SEASONAL
		outInds << endl;
	} // end of writeInd condition

}
}

//---------------------------------------------------------------------------

#if RS_ABC
dispstats Population::getDispStats(float resol) {
dispstats d; d.nPhilo = d.nDisp = d.nSuccess = 0; d.sumDist = d.sumDist2 = 0.0;

demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig;
emigTraits e;

// array to identify dispersing stage(s) and sex(es)
int nstages = 1;
int nsexes = 1;
if (dem.stageStruct) {
	nstages = sstruct.nStages;
	emig = pSpecies->getEmig();
	if (emig.sexDep) nsexes = 2;
}
bool **dispersing;
dispersing = new bool *[nstages];
for (int stg = 0; stg < nstages; stg++) {
	dispersing[stg] = new bool [nsexes];
	for (int sx = 0; sx < nsexes; sx++) {
		dispersing[stg][sx] = false;
		if (dem.stageStruct) {
			e = pSpecies->getEmigTraits(stg,sx);
			if (e.d0 > 0.0) dispersing[stg][sx] = true;
		}
		else dispersing[stg][sx] = true;
	}
}

int ninds = (int)inds.size();
#if RSDEBUG
//DEBUGLOG << "Population::getDispStats(): this=" << this
//	<< " resol=" << resol << " nstages=" << nstages << " nsexes=" << nsexes
//	<< " ninds=" << ninds
//	<< endl;
//for (int stg = 0; stg < nstages; stg++) {
//	for (int sx = 0; sx < nsexes; sx++) {
//DEBUGLOG << "Population::getDispStats(): stg=" << stg << " sx=" << sx
//	<< " dispersing=" << dispersing[stg][sx] << endl;
//	}
//}
#endif

for (int i = 0; i < ninds; i++) {
	indStats ind = inds[i]->getStats();
	if (!emig.sexDep) ind.sex = 0;
	if (dispersing[ind.stage][ind.sex]) {
		if (ind.status < 6) {
			if (ind.status == 0) d.nPhilo++;
			else {
				d.nDisp++; d.nSuccess++;
				Cell *pCell = inds[i]->getLocn(1);
				locn loc;
				if (pCell == 0) loc.x = loc.y = -1; // beyond boundary or in no-data cell
				else loc = pCell->getLocn();
				pCell = inds[i]->getLocn(0);
				locn natalloc = pCell->getLocn();
				double dist = resol * sqrt((float)((natalloc.x-loc.x)*(natalloc.x-loc.x)
																				 + (natalloc.y-loc.y)*(natalloc.y-loc.y)));
				d.sumDist += dist; d.sumDist2 += dist * dist;
			}
		}
		else d.nDisp++; // died during dispersal
	}
}

for (int stg = 0; stg < nstages; stg++) {
	delete[] dispersing[stg];
}
delete[] dispersing;

#if RSDEBUG
//DEBUGLOG << "Population::getDispStats(): this=" << this
//	<< " nPhilo=" << d.nPhilo << " nDisp=" << d.nDisp << " nSuccess=" << d.nSuccess
//	<< " sumDist=" << d.sumDist << " sumDist2=" << d.sumDist2
//	<< endl;
#endif

return d;
}
#endif // RS_ABC

//---------------------------------------------------------------------------
// Write records to genetics file
#if GROUPDISP || ROBFITT
void Population::outGenetics(const int rep,const int year,const int landNr,
	const bool patchmodel)
#else
void Population::outGenetics(const int rep,const int year,const int landNr)
#endif
{

simParams sim = paramsSim->getSim();

if (landNr >= 0) { // open file
	Genome *pGenome;
	genomeData gen = pSpecies->getGenomeData();
	if (gen.trait1Chromosome) {
		pGenome = new Genome(pSpecies->getNChromosomes(),pSpecies->getNLoci(0),
			pSpecies->isDiploid());
	}
	else {
		pGenome = new Genome(pSpecies);
	}
#if GROUPDISP || ROBFITT
	pGenome->outGenHeaders(rep,landNr,patchmodel,sim.outGenXtab);
#else
	pGenome->outGenHeaders(rep,landNr,sim.outGenXtab);
#endif
	delete pGenome;
	return;
}

if (landNr == -999) { // close file
	Genome *pGenome = new Genome();
#if GROUPDISP || ROBFITT
	pGenome->outGenHeaders(rep,landNr,patchmodel,sim.outGenXtab);
#else
	pGenome->outGenHeaders(rep,landNr,sim.outGenXtab);
#endif
	delete pGenome;
	return;
}

short spNum = pSpecies->getSpNum();
short nstages = 1;
if (pSpecies->stageStructured()) {
	stageParams sstruct = pSpecies->getStage();
	nstages = sstruct.nStages;
}


int ninds = (int)inds.size();
for (int i = 0; i < ninds; i++) {
	indStats ind = inds[i]->getStats();
	if (year == 0 || sim.outGenType == 1
	|| (sim.outGenType == 0 && ind.stage == 0)
	|| (sim.outGenType == 2 && ind.stage == nstages-1)) {
#if GROUPDISP || ROBFITT
		inds[i]->outGenetics(rep,year,spNum,landNr,patchmodel,sim.outGenXtab);
#else
		inds[i]->outGenetics(rep,year,spNum,landNr,sim.outGenXtab);
#endif
	}
}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


