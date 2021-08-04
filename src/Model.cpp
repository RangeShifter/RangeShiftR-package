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

#include "Model.h"

ofstream outPar;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Rcpp::List RunModel(Landscape *pLandscape,int seqsim)
{
//int Nsuit,yr,totalInds;
int yr,totalInds;
//float gradval,gradmin,gradmax;
bool filesOK;
//int t0,t1;

landParams ppLand = pLandscape->getLandParams();
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
//emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
initParams init = paramsInit->getInit();
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();

//t0 = time(0);

#if RSDEBUG
landPix p = pLandscape->getLandPix();
DEBUGLOG << "RunModel(): reps=" << sim.reps
	<< " ppLand.nHab=" << ppLand.nHab
	<< " p.pix=" << p.pix
	<< endl;
//DEBUGLOG << "RunModel(): random integers:";
//for (int i = 0; i < 5; i++) {
//	int rrrr = pRandom->IRandom(1000,2000); DEBUGLOG << " " << rrrr;
//}
DEBUGLOG << endl;
#endif

if (!ppLand.generated) {
	if (!ppLand.patchModel) { // cell-based landscape
		// create patches for suitable cells, adding unsuitable cells to the matrix
		// NB this is an overhead here, but is necessary in case the identity of
		// suitable habitats has been changed from one simulation to another (GUI or batch)
		// substantial time savings may result during simulation in certain landscapes
		pLandscape->allocatePatches(pSpecies);
	}
	pComm = new Community(pLandscape); // set up community
	// set up a sub-community associated with each patch (incl. the matrix)
	pLandscape->updateCarryingCapacity(pSpecies,0,0);
//	if (ppLand.rasterType <= 2 && ppLand.dmgLoaded) 
//		pLandscape->updateDamageIndices();
	patchData ppp;
	int npatches = pLandscape->patchCount();
	for (int i = 0; i < npatches; i++) {
		ppp = pLandscape->getPatchData(i);
#if RSDEBUG
//DEBUGLOG << "RunModel(): i = " << i
//	<< " ppp.pPatch = " << ppp.pPatch << " ppp.patchNum = " << ppp.patchNum
//	<< endl;
#endif
		pComm->addSubComm(ppp.pPatch,ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
//	if (i == 0 || ppp.pPatch->getK() > 0.0) {
//		// SET UP SUB-COMMUNITY FOR MATRIX PATCH AND ANY PATCH HAVING NON-ZERO CARRYING CAPACITY
//		pComm->addSubComm(ppp.pPatch,ppp.patchNum);
//	}
	}
	if (init.seedType == 0 && init.freeType < 2 && init.initFrzYr > 0) {
		// restrict available landscape to initialised region
		pLandscape->setLandLimits(init.minSeedX,init.minSeedY,
			init.maxSeedX,init.maxSeedY);
	}
	else {
		pLandscape->resetLandLimits();
	}
}

Rcpp::List list_outPop;

// Loop through replicates
for (int rep = 0; rep < sim.reps; rep++) {
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): starting simulation=" << sim.simulation << " rep=" << rep << endl;
#endif
		Rcpp::Rcout << endl << "starting replicate " << rep << endl;

	MemoLine(("Running replicate " + Int2Str(rep) + "...").c_str());

	if (sim.saveVisits && !ppLand.generated) {
		pLandscape->resetVisits();
	}

	patchChange patchchange;
	costChange costchange;
	int npatchchanges = pLandscape->numPatchChanges();
	int ncostchanges = pLandscape->numCostChanges();
	int ixpchchg = 0;
	int ixcostchg = 0;
#if RSDEBUG
DEBUGLOG << "RunModel(): npatchchanges=" << npatchchanges << " ncostchanges=" << ncostchanges << endl;
#endif

	if (ppLand.generated) {
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): generating new landscape ..." << endl;
#endif
		// delete previous community (if any)
		// Note: this must be BEFORE the landscape is reset (as a sub-community accesses
		// its corresponding patch upon deletion)
		if (pComm != 0) delete pComm;
		// generate new cell-based landscape
		MemoLine("...generating new landscape...");
		pLandscape->resetLand();
#if RSDEBUG
DEBUGLOG << "RunModel(): finished resetting landscape" << endl << endl;
#endif
		pLandscape->generatePatches();
//#if VCL
		if (v.viewLand || sim.saveMaps) {
			pLandscape->setLandMap();
			pLandscape->drawLandscape(rep,0,ppLand.landNum);
		}
//#endif
//#if BATCH
//		if (sim.saveMaps) {
//			pLandscape->drawLandscape(rep,0,ppLand.landNum);
//		}
//#endif
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): finished generating patches" << endl;
#endif
		pComm = new Community(pLandscape); // set up community
		// set up a sub-community associated with each patch (incl. the matrix)
//		pLandscape->updateCarryingCapacity(pSpecies,0);
		pLandscape->updateCarryingCapacity(pSpecies,0,0);
		patchData ppp;
		int npatches = pLandscape->patchCount();
#if RSDEBUG
DEBUGLOG << "RunModel(): patch count is " << npatches << endl;
#endif
		for (int i = 0; i < npatches; i++) {
			ppp = pLandscape->getPatchData(i);
#if RSDEBUG
//DEBUGLOG << "RunModel(): i = " << i
//	<< " ppp.pPatch = " << ppp.pPatch << " ppp.patchNum = " << ppp.patchNum
//	<< endl;
#endif
#if RSWIN64
			SubCommunity *pSubComm = pComm->addSubComm(ppp.pPatch,ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
//			if (ppp.y >= 9995) {
//				DEBUGLOG << "RunModel(): i=" << i << " pSubComm=" << pSubComm
//					<< endl;
//			}
#else
			pComm->addSubComm(ppp.pPatch,ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
#endif
//	if (i == 0 || ppp.pPatch->getK() > 0.0) {
//		// SET UP SUB-COMMUNITY FOR MATRIX PATCH AND ANY PATCH HAVING NON-ZERO CARRYING CAPACITY
//		pComm->addSubComm(ppp.pPatch,ppp.patchNum);
//	}
		}
		MemoLine("...completed...");
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): finished generating populations" << endl;
#endif
	}
	if (init.seedType == 0 && init.freeType < 2 && init.initFrzYr > 0) {
		// restrict available landscape to initialised region
		pLandscape->setLandLimits(init.minSeedX,init.minSeedY,
			init.maxSeedX,init.maxSeedY);
	}
	else {
		pLandscape->resetLandLimits();
	}

	filesOK = true;
	if (rep == 0) {
		// open output files
		if (sim.outRange) { // open Range file
			if (!pComm->outRangeHeaders(pSpecies,ppLand.landNum)) {
				MemoLine("UNABLE TO OPEN RANGE FILE");
				filesOK = false;
			}
		}
		if (sim.outOccup && sim.reps > 1)
			if (!pComm->outOccupancyHeaders(0)) {
				MemoLine("UNABLE TO OPEN OCCUPANCY FILE(S)");
				filesOK = false;
			}
		if (sim.outPop) {
			// open Population file
			if (!pComm->outPopHeaders(pSpecies,ppLand.landNum)) {
				MemoLine("UNABLE TO OPEN POPULATION FILE");
				filesOK = false;
			}
		}
		if (sim.outTraitsCells)
			if (!pComm->outTraitsHeaders(pSpecies,ppLand.landNum)) {
				MemoLine("UNABLE TO OPEN TRAITS FILE");
				filesOK = false;
			}
		if (sim.outTraitsRows)
			if (!pComm->outTraitsRowsHeaders(pSpecies,ppLand.landNum)) {
				MemoLine("UNABLE TO OPEN TRAITS ROWS FILE");
				filesOK = false;
			}
		if (sim.outConnect && ppLand.patchModel) // open Connectivity file
			if (!pLandscape->outConnectHeaders(0)) {
				MemoLine("UNABLE TO OPEN CONNECTIVITY FILE");
				filesOK = false;
			}
	}
#if RSDEBUG
DEBUGLOG << "RunModel(): completed opening output files" << endl;
#endif
	if (!filesOK) {
#if RSDEBUG
DEBUGLOG << "RunModel(): PROBLEM - closing output files" << endl;
#endif
		// close any files which may be open
		if (sim.outRange) {
			pComm->outRangeHeaders(pSpecies,-999);
		}
		if (sim.outOccup && sim.reps > 1)
			pComm->outOccupancyHeaders(-999);
		if (sim.outPop) {
			pComm->outPopHeaders(pSpecies,-999);
		}
		if (sim.outTraitsCells)
			pComm->outTraitsHeaders(pSpecies,-999);
		if (sim.outTraitsRows)
			pComm->outTraitsRowsHeaders(pSpecies,-999);
		if (sim.outConnect && ppLand.patchModel)
			pLandscape->outConnectHeaders(-999);
		return Rcpp::List::create(Rcpp::Named("Errors") = 666);
	}

	if (env.stoch && !env.local) {
		// create time series in case of global environmental stochasticity
		pLandscape->setGlobalStoch(sim.years+1);
	}

	if (grad.gradient) { // set up environmental gradient
		pLandscape->setEnvGradient(pSpecies,true);
	}

	if (sim.outConnect && ppLand.patchModel)
		pLandscape->createConnectMatrix();

	// variables to control dynamic landscape
	landChange landChg; landChg.chgnum = 0; landChg.chgyear = 999999;
	if (!ppLand.generated && ppLand.dynamic) {
		landChg = pLandscape->getLandChange(0); // get first change year
	}

	// set up populations in the community
	pLandscape->updateCarryingCapacity(pSpecies,0,0);
#if RSDEBUG
DEBUGLOG << "RunModel(): completed updating carrying capacity" << endl;
#endif
//	if (init.seedType != 2) {
		pComm->initialise(pSpecies,-1);
//	}
	bool updateland = false;
	int landIx = 0; // landscape change index

#if RSDEBUG
DEBUGLOG << "RunModel(): completed initialisation, rep=" << rep
	<< " pSpecies=" << pSpecies << endl;
#endif
	Rcpp::Rcout << "RunModel(): completed initialisation " << endl;

	// open a new individuals file for each replicate
	if (sim.outInds)
		pComm->outInds(rep,0,0,ppLand.landNum);
	// open a new genetics file for each replicate
	if (sim.outGenetics) {
		pComm->outGenetics(rep,0,0,ppLand.landNum);
		if (!dem.stageStruct && sim.outStartGenetic == 0) {
			// write genetic data for initialised individuals of non-strucutred population
			pComm->outGenetics(rep,0,0,-1);
		}
	}
#if RSDEBUG
	// output initialised Individuals
	if (sim.outInds)
		pComm->outInds(rep,-1,-1,-1);
#endif
		// open a new movement paths file for each replicate
		if (sim.outPaths)
			pLandscape->outPathsHeaders(rep,0);

	// years loop
	MemoLine("...running...");
	for (yr = 0; yr < sim.years; yr++) {
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): starting simulation=" << sim.simulation 
	<< " rep=" << rep << " yr=" << yr << endl;
#endif
		Rcpp::checkUserInterrupt();
		bool updateCC = false;
		if (yr < 4
		|| (yr < 31 && yr%10 == 0)
		|| (yr < 301 && yr%100 == 0)
		|| (yr < 3001 && yr%1000 == 0)
		|| (yr < 30001 && yr%10000 == 0)
		|| (yr < 300001 && yr%100000 == 0)
		|| (yr < 3000001 && yr%1000000 == 0)
		) {
			Rcpp::Rcout << "starting year " << yr << "..." << endl;
		}
		if (init.seedType == 0 && init.freeType < 2) {
			// apply any range restrictions
			if (yr == init.initFrzYr) {
				// release initial frozen range - reset landscape to its full extent
				pLandscape->resetLandLimits();
				updateCC = true;
			}
			if (init.restrictRange) {
				if (yr > init.initFrzYr && yr < init.finalFrzYr) {
					if ((yr-init.initFrzYr)%init.restrictFreq == 0) {
						// apply dynamic range restriction
						commStats s = pComm->getStats();
						int minY = s.maxY-init.restrictRows;
						if (minY < 0) minY = 0;
#if RSDEBUG
DEBUGLOG << "RunModel(): restriction yr=" << yr
	<< " s.minY=" << s.minY << " s.maxY=" << s.maxY
	<< " init.restrictRows=" << init.restrictRows
	<< " minY=" << minY
	<< endl;
#endif
						pLandscape->setLandLimits(ppLand.minX,minY,ppLand.maxX,ppLand.maxY);
						updateCC = true;
#if RSDEBUG
//landData d = pLandscape->getLandData();
//DEBUGLOG << "RunModel(): landscape yr=" << yr
//	<< " minX=" << d.minX << " minY=" << d.minY << " maxX=" << d.maxX << " maxY=" << d.maxY
//	<< endl;
#endif
					}
				}
				if (yr == init.finalFrzYr) {
					// apply final range restriction
					commStats s = pComm->getStats();
#if RSDEBUG
DEBUGLOG << "RunModel(): final restriction yr=" << yr
	<< " s.minY=" << s.minY << " s.maxY=" << s.maxY
	<< endl;
#endif
					pLandscape->setLandLimits(ppLand.minX,s.minY,ppLand.maxX,s.maxY);
					updateCC = true;
#if RSDEBUG
//landData d = pLandscape->getLandData();
//DEBUGLOG << "RunModel(): landscape yr=" << yr
//	<< " minX=" << d.minX << " minY=" << d.minY << " maxX=" << d.maxX << " maxY=" << d.maxY
//	<< endl;
#endif
				}
			}
		}
		// environmental gradient, stochasticity & local extinction
		// or dynamic landscape
		updateland = false;
		if (env.stoch || grad.gradient || ppLand.dynamic) {
			if (grad.shifting && yr > grad.shift_begin && yr < grad.shift_stop) {
				paramsGrad->incrOptY();
				pLandscape->setEnvGradient(pSpecies,false);
				updateCC = true;
			}
#if RSDEBUG
//DEBUGLOG << "RunModel(): yr=" << yr << " shift_begin=" << grad.shift_begin
//	<< " shift_stop=" << grad.shift_stop << " opt_y=" << grad.opt_y << endl;
#endif
			if (env.stoch) {
				if (env.local) pLandscape->updateLocalStoch();
				updateCC = true;
			}
			if (ppLand.dynamic) {
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " landChg.chgnum=" << landChg.chgnum
	<< " landChg.chgyear=" << landChg.chgyear
	<< " npatchchanges=" << npatchchanges << " ncostchanges=" << ncostchanges
	<< " ixpchchg=" << ixpchchg << " ixcostchg=" << ixcostchg
	<< endl;
#endif
				if (yr == landChg.chgyear) { // apply landscape change
					landIx = landChg.chgnum;       
					updateland = updateCC = true;
					if (ppLand.patchModel) { // apply any patch changes
						Patch *pPatch;
						Cell *pCell;
						patchchange = pLandscape->getPatchChange(ixpchchg++);
						while (patchchange.chgnum <= landIx && ixpchchg <= npatchchanges) {
#if RSDEBUG
//DEBUGLOG << "RunModel(): yr=" << yr << " landIx=" << landIx
//	<< " npatchchanges=" << npatchchanges << " ixpchchg=" << ixpchchg
//	<< " patchchange.chgnum=" << patchchange.chgnum
//	<< " .oldpatch=" << patchchange.oldpatch
//	<< " .newpatch=" << patchchange.newpatch
//	<< " .x=" << patchchange.x << " .y=" << patchchange.y
//	<< endl;
#endif
							// move cell from original patch to new patch
							pCell = pLandscape->findCell(patchchange.x,patchchange.y);
							if (patchchange.oldpatch != 0) { // not matrix
								pPatch = pLandscape->findPatch(patchchange.oldpatch);
								pPatch->removeCell(pCell);
							}
							if (patchchange.newpatch == 0) { // matrix
								pPatch = 0;
							}
							else {
								pPatch = pLandscape->findPatch(patchchange.newpatch);
								pPatch->addCell(pCell,patchchange.x,patchchange.y);
							}
							pCell->setPatch((intptr)pPatch);
							// get next patch change
							patchchange = pLandscape->getPatchChange(ixpchchg++);
						}
						ixpchchg--;
						pLandscape->resetPatches(); // reset patch limits
					}
					if (landChg.costfile != "NULL") { // apply any SMS cost changes
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " landChg.costfile=" << landChg.costfile << endl;
#endif
						Cell *pCell;
						costchange = pLandscape->getCostChange(ixcostchg++);
						while (costchange.chgnum <= landIx && ixcostchg <= ncostchanges) {
#if RSDEBUG
//DEBUGLOG << "RunModel(): yr=" << yr << " landIx=" << landIx
//	<< " ncostchanges=" << ncostchanges << " ixcostchg=" << ixcostchg
//	<< " costchange.chgnum=" << costchange.chgnum
//	<< " .x=" << costchange.x << " .y=" << costchange.y
//	<< " .oldcost=" << costchange.oldcost
//	<< " .newcost=" << costchange.newcost
//	<< endl;
#endif
							pCell = pLandscape->findCell(costchange.x,costchange.y);
							if (pCell != 0) {
								pCell->setCost(costchange.newcost);
							}
							costchange = pLandscape->getCostChange(ixcostchg++);
						}
						ixcostchg--;
						pLandscape->resetEffCosts();
					}
					if (landIx < pLandscape->numLandChanges()) { // get next change
						landChg = pLandscape->getLandChange(landIx);
					}
					else {
						landChg.chgyear = 9999999;
					}
				}
			}
		} // end of environmental gradient, etc.

		if (updateCC) {
			pLandscape->updateCarryingCapacity(pSpecies,yr,landIx);
		}

		
		if (sim.outConnect && ppLand.patchModel)
			pLandscape->resetConnectMatrix();

		if (ppLand.dynamic && updateland) {
//			trfrRules trfr = pSpecies->getTrfr();
			if (trfr.moveModel && trfr.moveType == 1) { // SMS
				if (!trfr.costMap) pLandscape->resetCosts(); // in case habitats have changed
			}
			// apply effects of landscape change to species present in changed patches
			pComm->patchChanges();
			pComm->dispersal(landIx,yr);
		}
		if (init.restrictRange) {
			// remove any population from region removed from restricted range
			if (yr > init.initFrzYr && yr < init.finalFrzYr) {
				if ((yr-init.initFrzYr)%init.restrictFreq == 0) {
					pComm->patchChanges();
				}
			}
		}

		if (init.seedType == 2) {
			// add any new initial individuals for the current year
			pComm->initialise(pSpecies,yr);
		}

		for(int gen = 0; gen < dem.repSeasons; gen++) // generation loop
		{
#if RSDEBUG
// TEMPORARY RANDOM STREAM CHECK
//if (yr%1 == 0 && gen == 0)
if (yr%1 == 0)
{
DEBUGLOG << endl << "RunModel(): start of gen " << gen << " in year " << yr
	<< " for rep " << rep << " (";
for (int i = 0; i < 5; i++) {
	int rrrr = pRandom->IRandom(1000,2000);
	DEBUGLOG << " " << rrrr;
}
DEBUGLOG << " )"	<< endl;
}
#endif

			if (v.viewPop || (sim.saveMaps && yr%sim.mapInt == 0)) {
				if (updateland && gen == 0) {
					pLandscape->drawLandscape(rep,landIx,ppLand.landNum);
				}
				pComm->draw(rep,yr,gen,ppLand.landNum);
			}
			// Output and pop. visualisation before reproduction
			if (v.viewPop || v.viewTraits || sim.outOccup
			|| 	sim.outTraitsCells || sim.outTraitsRows || sim.saveMaps)
				PreReproductionOutput(pLandscape,pComm,rep,yr,gen);
			// for non-structured population, also produce range and population output now
			if (!dem.stageStruct && (sim.outRange || sim.outPop))
				RangePopOutput(pComm,rep,yr,gen);
			if ( sim.ReturnPopRaster && sim.outPop && yr >= sim.outStartPop && yr%sim.outIntPop == 0) {
				list_outPop.push_back(pComm->addYearToPopList(rep,yr), "rep" + std::to_string(rep) + "_year" + std::to_string(yr));
			}

#if RSDEBUG
//DEBUGLOG << "RunModel(): completed RangePopOutput()"
//	<< " Total_Size = " << Total_Size << endl;
#endif

			// apply local extinction for generation 0 only
			// CHANGED TO *BEFORE* RANGE & POPN OUTPUT PRODUCTION IN v1.1,
			// SO THAT NOS. OF JUVENILES BORN CAN BE REPORTED
			if (!ppLand.patchModel && gen == 0) {
				if (env.localExt) pComm->localExtinction(0);
				if (grad.gradient && grad.gradType == 3) pComm->localExtinction(1);
			}

			// reproduction
			pComm->reproduction(yr);

			if (dem.stageStruct) {
				if (sstruct.survival == 0) { // at reproduction
					pComm->survival(0,2,1); // survival of all non-juvenile stages
				}
			}

			// Output and pop. visualisation AFTER reproduction
			if (dem.stageStruct && (sim.outRange || sim.outPop))
				RangePopOutput(pComm,rep,yr,gen);

#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed reproduction" << endl;
#endif

			// Dispersal

			pComm->emigration();
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed emigration" << endl;
#endif
			pComm->dispersal(landIx,yr);
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed dispersal" << endl;
#endif

			// survival part 0
			if (dem.stageStruct) {
				if (sstruct.survival == 0) { // at reproduction
					pComm->survival(0,0,1); // survival of juveniles only
				}
				if (sstruct.survival == 1) { // between reproduction events
					pComm->survival(0,1,1); // survival of all stages
				}
				if (sstruct.survival == 2) { // annually
					pComm->survival(0,1,0); // development only of all stages
//					pComm->survival(0,1,0); // development only of all stages
				}
			}
			else { // non-structured population
				pComm->survival(0,1,1);
			}
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed survival part 0" << endl;
#endif


				// output Individuals
				if (sim.outInds && yr >= sim.outStartInd && yr%sim.outIntInd == 0)
					pComm->outInds(rep,yr,gen,-1);
				// output Genetics
				if (sim.outGenetics && yr >= sim.outStartGenetic && yr%sim.outIntGenetic == 0)
					pComm->outGenetics(rep,yr,gen,-1);

			// survival part 1
			if (dem.stageStruct) {
//				if (sstruct.survival != 2) { // at reproduction or between reproduction events
					pComm->survival(1,0,1);
//				}
			}
			else { // non-structured population
				pComm->survival(1,0,1);
			}
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed survival part 1" << endl;
#endif

		} // end of the generation loop
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed generation loop" << endl;
#endif

		totalInds =  pComm->totalInds();
		if (totalInds <= 0) { yr++; break; }

		// Connectivity Matrix
		if (sim.outConnect && ppLand.patchModel
		&& yr >= sim.outStartConn && yr%sim.outIntConn == 0)
			pLandscape->outConnect(rep,yr);

		if (dem.stageStruct && sstruct.survival == 2) {  // annual survival - all stages
			pComm->survival(0,1,2);
			pComm->survival(1,0,1);
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed annual survival" << endl;
#endif
		}

		if (dem.stageStruct) {
			pComm->ageIncrement(); // increment age of all individuals
			if (sim.outInds && yr >= sim.outStartInd && yr%sim.outIntInd == 0)
				pComm->outInds(rep,yr,-1,-1); // list any individuals dying having reached maximum age
			pComm->survival(1,0,1);						// delete any such individuals
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed Age_increment and final survival" << endl;
#endif
		totalInds =  pComm->totalInds();
		if (totalInds <= 0) { yr++; break; }
		}

	} // end of the years loop

	// Final output and popn. visualisation
	if (sim.saveMaps && yr%sim.mapInt == 0) {
		if (updateland) {
			pLandscape->drawLandscape(rep,landIx,ppLand.landNum);
		}
		pComm->draw(rep,yr,0,ppLand.landNum);
	}
		// produce final summary output
		if (v.viewPop || v.viewTraits || sim.outOccup
		|| 	sim.outTraitsCells || sim.outTraitsRows || sim.saveMaps)
			PreReproductionOutput(pLandscape,pComm,rep,yr,0);
		if (sim.outRange || sim.outPop)
			RangePopOutput(pComm,rep,yr,0);
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed final summary output" << endl;
#endif

	pComm->resetPopns();

//	if (batchMode) {
//		// delete the community of species using the landscape
//		pComm->resetPopns();
//	}

	//Reset the gradient optimum
	if (grad.gradient) paramsGrad->resetOptY();

	pLandscape->resetLandLimits();
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " landIx=" << "reset"
	<< " npatchchanges=" << npatchchanges << " ncostchanges=" << ncostchanges
	<< " ixpchchg=" << ixpchchg << " ixcostchg=" << ixcostchg
	<< endl;
#endif
	if (ppLand.patchModel && ppLand.dynamic && ixpchchg > 0) {
		// apply any patch changes to reset landscape to original configuration
		// (provided that at least one has already occurred)
		patchChange patchchange;
		Patch *pPatch;
		Cell *pCell;
		patchchange = pLandscape->getPatchChange(ixpchchg++);
		while (patchchange.chgnum <= 666666 && ixpchchg <= npatchchanges) {
#if RSDEBUG
//DEBUGLOG << "RunModel(): yr=" << yr << " landIx=" << "reset"
//	<< " npatchchanges=" << npatchchanges << " ixpchchg=" << ixpchchg
//	<< " patchchange.chgnum=" << patchchange.chgnum
//	<< " .oldpatch=" << patchchange.oldpatch
//	<< " .newpatch=" << patchchange.newpatch
//	<< " .x=" << patchchange.x << " .y=" << patchchange.y
//	<< endl;
#endif
			// move cell from original patch to new patch
			pCell = pLandscape->findCell(patchchange.x,patchchange.y);
			if (patchchange.oldpatch != 0) { // not matrix
				pPatch = pLandscape->findPatch(patchchange.oldpatch);
				pPatch->removeCell(pCell);
			}
			if (patchchange.newpatch == 0) { // matrix
				pPatch = 0;
			}
			else {
				pPatch = pLandscape->findPatch(patchchange.newpatch);
				pPatch->addCell(pCell,patchchange.x,patchchange.y);
			}
			pCell->setPatch((intptr)pPatch);
			// get next patch change
			patchchange = pLandscape->getPatchChange(ixpchchg++);
		}
		ixpchchg--;
		pLandscape->resetPatches();
	}
	if (ppLand.dynamic) {
		trfrRules trfr = pSpecies->getTrfr();
		if (trfr.moveModel && trfr.moveType == 1) { // SMS
			if (ixcostchg > 0) {
				// apply any cost changes to reset landscape to original configuration
				// (provided that at least one has already occurred)
						Cell *pCell;
						costchange = pLandscape->getCostChange(ixcostchg++);
						while (costchange.chgnum <= 666666 && ixcostchg <= ncostchanges) {
#if RSDEBUG
//DEBUGLOG << "RunModel(): yr=" << yr << " landIx=" << landIx
//	<< " ncostchanges=" << ncostchanges << " ixcostchg=" << ixcostchg
//	<< " costchange.chgnum=" << costchange.chgnum
//	<< " .x=" << costchange.x << " .y=" << costchange.y
//	<< " .oldcost=" << costchange.oldcost
//	<< " .newcost=" << costchange.newcost
//	<< endl;
#endif
							pCell = pLandscape->findCell(costchange.x,costchange.y);
							if (pCell != 0) {
								pCell->setCost(costchange.newcost);
							}
							costchange = pLandscape->getCostChange(ixcostchg++);
						}
						ixcostchg--;
				pLandscape->resetEffCosts();
			}
			if (!trfr.costMap) pLandscape->resetCosts(); // in case habitats have changed
		}
	}
//					if (landIx < pLandscape->numLandChanges()) { // get next change
//						landChg = pLandscape->getLandChange(landIx);
//					}
//					else {
//						landChg.chgyear = 9999999;
//					}
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed reset"
	<< endl;
#endif

	if (sim.outConnect && ppLand.patchModel)
		pLandscape->resetConnectMatrix(); // set connectivity matrix to zeroes

	if (sim.outInds) // close Individuals output file
		pComm->outInds(rep,0,0,-999);
	if (sim.outGenetics) // close Genetics output file
		pComm->outGenetics(rep,0,0,-999);

	if (sim.saveVisits) {
		pLandscape->outVisits(rep,ppLand.landNum);
		pLandscape->resetVisits();
	}

	if (sim.outPaths)
		pLandscape->outPathsHeaders(rep,-999);
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): finished rep=" << rep << endl;
#endif

} // end of the replicates loop

if (sim.outConnect && ppLand.patchModel) {
	pLandscape->deleteConnectMatrix();
	pLandscape->outConnectHeaders(-999); // close Connectivity Matrix file
}

// Occupancy outputs
if (sim.outOccup && sim.reps > 1) {
	MemoLine("Writing final occupancy output...");
	pComm->outOccupancy();
	pComm->outOccSuit(v.viewGraph);
//	pComm->deleteOccupancy((sim.years/sim.outInt)+1);
	pComm->deleteOccupancy((sim.years/sim.outIntOcc)+1);
	pComm->outOccupancyHeaders(-999);
	MemoLine("...finished");
}

if (sim.outRange) {
	pComm->outRangeHeaders(pSpecies,-999); // close Range file
}
if (sim.outPop) {
	pComm->outPopHeaders(pSpecies,-999); // close Population file
}
if (sim.outTraitsCells)
	pComm->outTraitsHeaders(pSpecies,-999); // close Traits file
if (sim.outTraitsRows)
	pComm->outTraitsRowsHeaders(pSpecies,-999); // close Traits rows file
// close Individuals & Genetics output files if open
// they can still be open if the simulation was stopped by the user
if (sim.outInds) pComm->outInds(0,0,0,-999);
if (sim.outGenetics) pComm->outGenetics(0,0,0,-999);

MemoLine("Deleting community...");
delete pComm; pComm = 0;
MemoLine("...finished");

// Write performance data
//t1 = time(0);
//RSlog << "Simulation," << sim.simulation << "," << sim.reps << "," << sim.years
//	<< "," << t1-t0 << endl;

	return list_outPop;

}

// Check whether a specified directory path exists
bool is_directory(const char *pathname) {
struct stat info;
if (stat(pathname, &info) != 0) return false; // path does not exist
if (S_ISDIR(info.st_mode)) return true;
return false;
}

//---------------------------------------------------------------------------
bool CheckDirectory(void)
{
bool errorfolder = false;

string subfolder;

subfolder = paramsSim->getDir(0) + "Inputs";
const char *inputs = subfolder.c_str();
if (!is_directory(inputs)) errorfolder = true;
subfolder = paramsSim->getDir(0) + "Outputs";
const char *outputs = subfolder.c_str();
if (!is_directory(outputs)) errorfolder = true;
subfolder = paramsSim->getDir(0) + "Output_Maps";
const char *outputmaps = subfolder.c_str();
if (!is_directory(outputmaps)) errorfolder = true;

return errorfolder;
}

//---------------------------------------------------------------------------
//For outputs and population visualisations pre-reproduction
void PreReproductionOutput(Landscape *pLand,Community *pComm,int rep,int yr,int gen)
{
#if RSDEBUG || VCL
landParams ppLand = pLand->getLandParams();
#endif
simParams sim = paramsSim->getSim();
simView v = paramsSim->getViews();

#if RSDEBUG
DEBUGLOG << "PreReproductionOutput(): 11111 rep=" << rep << " yr=" << yr << " gen=" << gen
	<< " landNum=" << ppLand.landNum << " maxX=" << ppLand.maxX << " maxY=" << ppLand.maxY
	<< endl;
DEBUGLOG << "PreReproductionOutput(): 11112 outRange=" << sim.outRange
	<< " outIntRange=" << sim.outIntRange
	<< " outPop=" << sim.outPop << " outIntPop=" << sim.outIntPop
	<< endl;
#endif

#if RSDEBUG
//DEBUGLOG << "PreReproductionOutput(): 22222 " << endl;
#endif

traitCanvas tcanv;
for (int i = 0; i < NTRAITS; i++) {
		tcanv.pcanvas[i] = 0;
}

// trait outputs and visualisation

if (v.viewTraits) {
	tcanv = SetupTraitCanvas();
}

if (v.viewTraits
|| ((sim.outTraitsCells && yr >= sim.outStartTraitCell && yr%sim.outIntTraitCell == 0) ||
		(sim.outTraitsRows && yr >= sim.outStartTraitRow && yr%sim.outIntTraitRow == 0)))
{
	pComm->outTraits(tcanv,pSpecies,rep,yr,gen);
}

#if RSDEBUG
//DEBUGLOG << "PreReproductionOutput(): 33333 " << endl;
#endif

if (sim.outOccup && yr%sim.outIntOcc == 0 && gen == 0)
	pComm->updateOccupancy(yr/sim.outIntOcc,rep);

#if RSDEBUG
//DEBUGLOG << "PreReproductionOutput(): 88888 " << endl;
#endif

// Remaining graphical output actions are performed for GUI only

#if RSDEBUG
//DEBUGLOG << "PreReproductionOutput(): finished " << endl;
#endif
}

//For outputs and population visualisations pre-reproduction
void RangePopOutput(Community *pComm,int rep,int yr,int gen)
{
simParams sim = paramsSim->getSim();

if (sim.outRange && (yr%sim.outIntRange == 0 || pComm->totalInds() <= 0))
	pComm->outRange(pSpecies,rep,yr,gen);

if (sim.outPop && yr >= sim.outStartPop && yr%sim.outIntPop == 0)
	pComm->outPop(rep,yr,gen);

}

//---------------------------------------------------------------------------
void OutParameters(Landscape *pLandscape)
{
double k;
//int nrows,ncols,nsexes,nstages;
int nsexes,nstages;

landParams ppLand = pLandscape->getLandParams();
genLandParams ppGenLand = pLandscape->getGenLandParams();
envGradParams grad = paramsGrad->getGradient();
envStochParams env = paramsStoch->getStoch();
demogrParams dem = pSpecies->getDemogr();
stageParams sstruct = pSpecies->getStage();
emigRules emig = pSpecies->getEmig();
trfrRules trfr = pSpecies->getTrfr();
settleType sett = pSpecies->getSettle();
settleRules srules;
settleSteps ssteps;
settleTraits settleDD;
simParams sim = paramsSim->getSim();

string name;
if (sim.batchMode)
	name = paramsSim->getDir(2)
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation)
		+ "_Land" + Int2Str(ppLand.landNum) + "_Parameters.txt";
else
	name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation) + "_Parameters.txt";
outPar.open(name.c_str());

outPar << "RangeShifter 2.0 ";

outPar << endl;

outPar << "================ ";

outPar << "   =====================";
outPar << endl << endl;

outPar << "BATCH MODE \t";
if (sim.batchMode) outPar << "yes" << endl; else outPar << "no" << endl;
outPar << "SEED \t" << RS_random_seed << endl;
outPar << "REPLICATES \t" << sim.reps << endl;
outPar << "YEARS \t" << sim.years << endl;
outPar << "REPRODUCTIVE SEASONS / YEAR\t" << dem.repSeasons << endl;
if (ppLand.patchModel){
	outPar << "PATCH-BASED MODEL" << endl;
	outPar << "No. PATCHES \t" << pLandscape->patchCount()-1 << endl;
}
else
	outPar << "CELL-BASED MODEL" << endl;
outPar << "BOUNDARIES \t";
if (sim.absorbing) outPar << "absorbing" << endl;
else outPar << "reflective" << endl;
outPar << endl;

outPar << "LANDSCAPE:\t";
if (ppLand.generated) {
	outPar << "artificially generated map" << endl;
	outPar << "TYPE: \t";
	if (ppGenLand.continuous) outPar << "continuous \t";
	else outPar << "discrete \t";
	if (ppGenLand.fractal) outPar << "fractal";
	else outPar << "random";
	outPar << endl << "PROPORTION OF SUITABLE HABITAT (p)\t" << ppGenLand.propSuit <<endl;
	if (ppGenLand.fractal) outPar << "HURST EXPONENT\t" << ppGenLand.hurst <<endl;
}
else {
	outPar << "imported map" << endl;
	outPar << "TYPE: \t";
	switch (ppLand.rasterType) {
		case 0:
			outPar << "habitat codes" << endl;
			break;
		case 1:
			outPar << "habitat % cover" << endl;
			break;
		case 2:
			outPar << "habitat quality" << endl;
			break;
	}
	outPar << "FILE NAME: ";
	if (ppLand.dynamic) {
		outPar << name_landscape << endl;
	}
	else{
		outPar << name_landscape << endl;
	}
	if (ppLand.patchModel) {
		outPar << "PATCH FILE: " << name_patch << endl;
	}
	if (trfr.costMap) {
		outPar << "COSTS FILE: " << name_costfile << endl;
	}
	outPar << "No. HABITATS:\t" << ppLand.nHab << endl;
}
outPar << "RESOLUTION (m): \t" << ppLand.resol << endl;
outPar << "DIMENSIONS:  X " << ppLand.dimX << "  Y " << ppLand.dimY << endl;
outPar << "AVAILABLE:   min.X " << ppLand.minX << " min.Y " << ppLand.minY
	<< "  max.X " << ppLand.maxX << " max.Y " << ppLand.maxY << endl;
if (!ppLand.generated && ppLand.dynamic) {
	landChange chg;
	outPar << "DYNAMIC LANDSCAPE: " << endl;
	int nchanges = pLandscape->numLandChanges();
	for (int i = 0; i < nchanges; i++) {
		chg = pLandscape->getLandChange(i);
		outPar << "Change no. " << chg.chgnum << " in year " << chg.chgyear << endl;
		outPar << "Landscape: " << chg.habfile << endl;
		if (ppLand.patchModel) {
			outPar << "Patches  : " << chg.pchfile << endl;
		}
		if (chg.costfile != "none" && chg.costfile != "NULL") {
			outPar << "Costs    : " << chg.costfile << endl;			
		}
//		outPar << "Change no. " << chg.chgnum << " in year " << chg.chgyear
//			<< " habitat map: " << chg.habfile << endl;
	}
}

outPar << endl << "SPECIES DISTRIBUTION LOADED: \t";
//if (sim.initDistLoaded)
if (ppLand.spDist)
{
	outPar << "yes" << endl;
	outPar << "RESOLUTION (m)\t" << ppLand.spResol << endl;
	outPar << "FILE NAME: ";
	outPar << name_sp_dist << endl;
}
else outPar << "no" << endl;

outPar << endl << "ENVIRONMENTAL GRADIENT:\t ";
if (grad.gradient)
{
	switch (grad.gradType) {
	case 1:
		if (dem.stageStruct) outPar << "Density dependence strength (1/b)" << endl;
		else outPar << "Carrying capacity (K)" << endl;
		break;
	case 2:
		if (dem.stageStruct) outPar << "Fecundity" << endl;
		else outPar << "Intrinsic growth rate (r)" << endl;
		break;
	case 3:
		outPar << "Local extinction probability" << endl;
		break;
	default:
		outPar << "ERROR ERROR ERROR" << endl;
		;
	}
	outPar << "G:\t\t " << grad.grad_inc << endl;
	outPar << "optimum Y:\t " << grad.opt_y << endl;
	outPar << "f:\t\t " << grad.factor << endl;
	if (grad.gradType == 3) outPar << "Local extinction prob. at optimum:\t "
		<< grad.extProbOpt << endl;
	outPar << "GRADIENT SHIFTING:\t ";
	if (grad.shifting)
	{
		outPar << "yes" << endl;
		outPar << "SHIFTING RATE  (rows/year):\t " << grad.shift_rate << endl;
		outPar << "SHIFTING START (year):\t\t " << grad.shift_begin << endl;
		outPar << "SHIFTING STOP  (year):\t\t " << grad.shift_stop << endl;
	}
  else   outPar << "no" << endl;
}
else outPar << "no";
outPar << endl;
outPar << "ENVIRONMENTAL STOCHASTICITY:\t";
if (env.stoch) {
	outPar << "yes" << endl;
	outPar << "TYPE\t in ";
	if (dem.stageStruct) {
		if (env.inK) outPar << "1/b" << endl;
		else outPar << "fecundity" << endl;
	}
	else{
		if (env.inK) outPar << "K" << endl;
		else outPar << "R" << endl;
	}
	outPar << "SPATIAL AUTOCORRELATION\t ";
	if (env.local) outPar << "local" << endl;
	else outPar << "global" << endl;
	outPar << "TEMPORAL AUTOCORRELATION (ac)\t" << env.ac << endl;
	outPar << "AMPLITUDE (std)\t" << env.std << endl;
	if (dem.stageStruct) {
		if (env.inK) {
			outPar << "MIN. 1/b\t" << pSpecies->getMinMax(0)
				* (10000.0/(float)(ppLand.resol*ppLand.resol)) << endl;
			outPar << "MAX. 1/b\t" << pSpecies->getMinMax(1)
				* (10000.0/(float)(ppLand.resol*ppLand.resol)) << endl;
		}
		else {
			outPar << "MIN. fecundity\t" << pSpecies->getMinMax(0) << endl;
			outPar << "MAX. fecundity\t" << pSpecies->getMinMax(1) << endl;
		}
	}
	else {
		if (env.inK) {
			outPar << "MIN. K\t" << pSpecies->getMinMax(0)
				* (10000.0/(float)(ppLand.resol*ppLand.resol)) << endl;
			outPar << "MAX. K\t" << pSpecies->getMinMax(1)
				* (10000.0/(float)(ppLand.resol*ppLand.resol)) << endl;
		}
		else {
			outPar << "MIN. r\t" << pSpecies->getMinMax(0) << endl;
			outPar << "MAX. r\t" << pSpecies->getMinMax(1) << endl;
		}
	}
}
else outPar << "no" << endl;
outPar << "LOCAL EXTINCTION PROBABILITY:\t";
if (env.localExt) outPar << env.locExtProb << endl;
else outPar << "0.0" << endl;

outPar << endl << "SPECIES' PARAMETERS." << endl;
outPar << "REPRODUCTION:" << endl;
outPar << "TYPE: ";
switch (dem.repType) {
	case 0:
		outPar << "Asexual / Only female model" << endl;
		break;
	case 1:
		outPar << "Sexual model (simple)";
		outPar << endl;
		outPar << "PROP. of MALES\t" << dem.propMales << endl;
		break;
	case 2:
		outPar << "Sexual model (explicit mating system)" << endl;
		outPar << "PROP. of MALES\t" << dem.propMales << endl;
		outPar << "MAX. HAREM SIZE (h)\t" << dem.harem << endl;
	break;
}
outPar << "STAGE STRUCTURE:\t";
if (dem.stageStruct){
	outPar << "yes" << endl;
	outPar << "PROBABILITY OF REPRODUCING IN SUBSEQUENT SEASONS\t" << sstruct.probRep << endl;
	outPar << "No. OF REP. SEASONS BEFORE SUBSEQUENT REPRODUCTIONS\t" << sstruct.repInterval << endl;
	if (!ppLand.generated && ppLand.dynamic) {
		outPar << "ACTION AFTER POPULATION DESTRUCTION: all individuals ";
		if (sstruct.disperseOnLoss) outPar << "disperse" << endl;
		else outPar << "die" << endl;
	}
	outPar << "No. STAGES\t" << sstruct.nStages << endl;
	outPar << "MAX. AGE\t" << sstruct.maxAge << endl;
	// no sex-specific demographic parameters
	if (dem.repType != 2) {
		outPar << "MIN. AGES:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "stage\t" << i << ":\t" << pSpecies->getMinAge(i,0) << "\tyears"<< endl;
		}
		outPar << "FECUNDITIES:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "stage\t" << i << ":\t" << pSpecies->getFec(i,0) << endl;
		}
		outPar << "DEVELOPMENT PROB.:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "stage\t" << i << ":\t" << pSpecies->getDev(i,0) << endl;
		}
		outPar << "SURVIVAL PROB.:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(i,0) << endl;
		}
	}
	// sex-specific demographic parameters
	else {
		outPar << "MIN. AGES:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "males " << i << ":\t" << pSpecies->getMinAge(i,1) << " years;\t";
			outPar << "females " << i << ":\t" << pSpecies->getMinAge(i,0) << " years" << endl;
		}
		outPar << "FECUNDITIES:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "males   " << i << ":\t" << pSpecies->getFec(i,1) << endl;
			outPar << "females " << i << ":\t" << pSpecies->getFec(i,0) << endl;
		}
		outPar << "DEVELOPMENT PROB.:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "males   " << i << ":\t" << pSpecies->getDev(i,1) << endl;
			outPar << "females " << i << ":\t" << pSpecies->getDev(i,0) << endl;
		}
		outPar << "SURVIVAL PROB.:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "males   " << i << ":\t" << pSpecies->getSurv(i,1) << endl;
			outPar << "females " << i << ":\t" << pSpecies->getSurv(i,0) << endl;
		}
	}
/*
#if RSDEBUG
	outPar << endl << "TRANSITION MATRIX AS ENTERED:" << endl;
	if (batchMode) {
		DEBUGLOG << "outParameters(): matrix = " << matrix
			<< " matrix[1][1] = " << matrix[1][1]
			<< endl;
		if (dem.repType == 2) {
			nrows = sstruct.nStages*2-1; ncols = sstruct.nStages*2;
		}
		else {
			nrows = sstruct.nStages; ncols = sstruct.nStages;
		}
		for (int i = 0; i < nrows; i++) {
			for (int j = 0; j < ncols; j++) {
				outPar << matrix[j][i] << "\t";
			}
			outPar << endl;
		}
	}
#if VCL
	else {

	// NOTE: TO PREVENT COMPILING FOR BATCH MODE, THIS CODE NEEDS TO BE INCLUDED IN COMPILER
	// CONDITIONAL BLOCK  AS SHOWN

//	outPar << "Row count: " << frmSpecies->transMatrix->RowCount << endl;
//	outPar << "Col count: " << frmSpecies->transMatrix->ColCount << endl;
		for (int i = 1; i < frmSpecies->transMatrix->RowCount; i++) {
			for (int j = 1; j < frmSpecies->transMatrix->ColCount; j++) {
				outPar << frmSpecies->transMatrix->Cells[j][i].ToDouble() << "\t";
			}
			outPar << endl;
		}
	}
#endif
	outPar << endl;
#endif
*/

	outPar << "SCHEDULING OF SURVIVAL: ";
	switch (sstruct.survival) {
	case 0:
		outPar << "At reproduction" << endl;
		break;
	case 1:
		outPar << "Between reproductive events" << endl;
		break;
	case 2:
		outPar << "Annually" << endl;
		break;
	}

	int mSize; // index for weights matrices
	if (dem.repType == 2) mSize = sstruct.nStages * NSEXES;
	else mSize = sstruct.nStages;

	outPar << "DENSITY-DEPENDENCE IN FECUNDITY:\t";
	if (sstruct.fecDens) {
		outPar << "yes" << endl;
		if (sstruct.fecStageDens) {
			outPar << "STAGE'S WEIGHTS:" << endl;
			for (int i = 0; i < mSize; i++) {
				if (dem.repType == 2) {
					outPar << "stage " << i/NSEXES << " ";
					if (i%NSEXES == 0) outPar << "males  : \t";
					else outPar << "females: \t";
				}
				else outPar << "stage " << i << ": \t";
				for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtFec(j,i) << "\t";
				outPar << endl;
			}
		}
		else outPar << "not stage-dependent" << endl;
	}
	else outPar << "no" << endl;

	densDepParams ddparams = pSpecies->getDensDep();

	outPar << "DENSITY-DEPENDENCE IN DEVELOPMENT:\t";
	if (sstruct.devDens) {
		outPar << "yes - coefficient: " << ddparams.devCoeff << endl;
		if (sstruct.devStageDens) {
			outPar << "STAGE'S WEIGHTS:" << endl;
			for (int i = 0; i < mSize; i++) {
				if (dem.repType == 2) {
					outPar << "stage " << i/NSEXES << " ";
					if (i%NSEXES == 0) outPar << "males  : \t";
					else outPar << "females: \t";
				}
				else outPar << "stage " << i << ": \t";
				for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtDev(j,i) << "\t";
				outPar << endl;
			}
		}
		else outPar << "not stage-dependent" << endl;
	}
	else outPar << "no" << endl;

	outPar << "DENSITY-DEPENDENCE IN SURVIVAL:\t\t";
	if (sstruct.survDens) {
		outPar << "yes - coefficient: " << ddparams.survCoeff << endl;
		if (sstruct.survStageDens) {
			outPar << "STAGE'S WEIGHTS:" << endl;
			for (int i = 0; i < mSize; i++) {
				if (dem.repType == 2) {
					outPar << "stage " << i/NSEXES << " ";
					if (i%NSEXES == 0) outPar << "males  : \t";
					else outPar << "females: \t";
				}
				else outPar << "stage " << i << ": \t";
				for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtSurv(j,i) << "\t";
				outPar << endl;
			}
		}
		else outPar << "not stage-dependent" << endl;
	}
	else outPar << "no" << endl;
} // end of if (dem.stageStruct)
else { // not stage-strutured
  outPar << "no" << endl;
  outPar << "Rmax\t" << dem.lambda << endl;
	outPar << "bc\t" << dem.bc << endl;
}

if (dem.stageStruct) {
	outPar << endl << "HABITAT SPECIFIC 1/b:" << endl;
}
else {
	outPar << endl << "CARRYING CAPACITIES:" << endl;
}
int nhab = ppLand.nHab;             
if (ppLand.generated) { 
	if (ppGenLand.continuous) nhab = 1;
}
for (int i = 0; i < nhab; i++) {
	k = pSpecies->getHabK(i) * (10000.0/(float)(ppLand.resol*ppLand.resol));
	if (!ppLand.generated && ppLand.rasterType == 0) { // imported & habitat codes
		outPar << "Habitat " << pLandscape->getHabCode(i) << ": \t";
	}
	else {
		outPar << "Habitat " << i << ": ";
	}
	if (dem.stageStruct) outPar << "1/b ";
	else outPar << "K ";
	outPar << k << endl;
}

emigTraits ep0,ep1;
emigParams eparams0,eparams1;
string sexdept 		= "SEX-DEPENDENT:   ";
string stgdept  = "STAGE-DEPENDENT: ";
string indvar = "INDIVIDUAL VARIABILITY: ";
string emigstage = "EMIGRATION STAGE: ";

outPar << endl << "DISPERSAL - EMIGRATION:\t";
if (emig.densDep) {
	outPar << "density-dependent" << endl;

	if (emig.sexDep) {
		outPar << sexdept << "yes" << endl;
		if (emig.stgDep) {
			outPar << stgdept << "yes" << endl;
			outPar << indvar << "no" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage " << i << ":" << endl;
				ep0 = pSpecies->getEmigTraits(i,0);
				ep1 = pSpecies->getEmigTraits(i,1);
				outPar << "D0:    females "<< ep0.d0 		<< "  males " << ep1.d0    << endl;
				outPar << "alpha: females "<< ep0.alpha << "  males " << ep1.alpha << endl;
				outPar << "beta:  females "<< ep0.beta 	<< "  males " << ep1.beta  << endl;
			}
		}
		else { // !emig.stgDep
			outPar << stgdept << "no" << endl;
			outPar << indvar;
			if (emig.indVar) {
				eparams0 = pSpecies->getEmigParams(0,0);
				eparams1 = pSpecies->getEmigParams(0,1);
				outPar << "yes" << endl;
				if (dem.stageStruct) {
					outPar << emigstage << emig.emigStage << endl;
				}
				outPar << "D0 females:     mean " << eparams0.d0Mean << "  s.d. " << eparams0.d0SD
					<< "  scaling factor " << eparams0.d0Scale << endl;
				outPar << "D0 males:       mean " << eparams1.d0Mean << "  s.d. " << eparams1.d0SD
					<< "  scaling factor " << eparams1.d0Scale << endl;
				outPar << "Alpha females:  mean " << eparams0.alphaMean << "  s.d. " << eparams0.alphaSD
					<< "  scaling factor " << eparams0.alphaScale << endl;
				outPar << "Alpha males:    mean " << eparams1.alphaMean << "  s.d. " << eparams1.alphaSD
					<< "  scaling factor " << eparams1.alphaScale << endl;
				outPar << "Beta females:   mean " << eparams0.betaMean << "  s.d. " << eparams0.betaSD
					<< "  scaling factor " << eparams0.betaScale << endl;
				outPar << "Beta males:     mean " << eparams1.betaMean << "  s.d. " << eparams1.betaSD
					<< "  scaling factor " << eparams1.betaScale << endl;
			}
			else {
				outPar << "no" << endl;
				ep0 = pSpecies->getEmigTraits(0,0);
				ep1 = pSpecies->getEmigTraits(0,1);
				outPar << "D0:    females "<< ep0.d0		<< "  males " << ep1.d0    << endl;
				outPar << "alpha: females "<< ep0.alpha << "  males " << ep1.alpha << endl;
				outPar << "beta:  females "<< ep0.beta 	<< "  males " << ep1.beta  << endl;
			}
		}
	}
	else { // !emig.sexDep
		outPar << sexdept << "no" << endl;
		if (emig.stgDep) {
			outPar << stgdept << "yes" << endl;
			outPar << indvar << "no" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				ep0 = pSpecies->getEmigTraits(i,0);
				outPar << "stage " << i << ": \t"  << "D0: " << ep0.d0;
				outPar<< " \talpha: " << ep0.alpha << " \tbeta: " << ep0.beta << endl;
			}
		}
		else { // !emig.stgDep
			outPar << stgdept << "no" << endl;
			outPar << indvar;
			if (emig.indVar) {
				eparams0 = pSpecies->getEmigParams(0,0);
				emigScales scale = pSpecies->getEmigScales();
				outPar << "yes" << endl;
				if (dem.stageStruct) {
					outPar << emigstage << emig.emigStage << endl;
				}
				outPar << "D0 mean:    " << eparams0.d0Mean 		<< "  s.d.: " << eparams0.d0SD
					<< "  scaling factor: " << scale.d0Scale << endl;
				outPar << "Alpha mean: " << eparams0.alphaMean 	<< "  s.d.: " << eparams0.alphaSD
					<< "  scaling factor: " << scale.alphaScale << endl;
				outPar << "Beta mean:  " << eparams0.betaMean 	<< "  s.d.: " << eparams0.betaSD
					<< "  scaling factor: " << scale.betaScale << endl;
			}
			else{
				outPar << "no" << endl;
				ep0 = pSpecies->getEmigTraits(0,0);
				outPar << "D0:    " << ep0.d0 << endl;
				outPar << "alpha: " << ep0.alpha << endl;
				outPar << "beta:  " << ep0.beta << endl;
			}
		}
	}
}
else { // not density-dependent
	string initprob = "INITIAL EMIGRATION PROB. ";
	outPar << "density-independent" << endl;
	if (!trfr.moveModel) { // transfer by kernel
		outPar << "USE FULL KERNEL TO DETERMINE EMIGRATION: ";
		if (pSpecies->useFullKernel()) outPar << "yes";
		else outPar << "no";
		outPar << endl;
	}

	if (emig.sexDep) {
		outPar << sexdept << "yes" << endl;
		if (emig.stgDep) {
			outPar << stgdept << "yes" << endl;
			outPar << indvar << "no" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage " << i << ": \t" << "EMIGRATION PROB.: \tfemales "
					<< pSpecies->getEmigD0(i,0) <<" \tmales "<< pSpecies->getEmigD0(i,1) << endl;
			}
		}
		else { // !emig.stgDep
			outPar << stgdept << "no" << endl;
			outPar << indvar;
			if (emig.indVar) {
				eparams0 = pSpecies->getEmigParams(0,0);
				eparams1 = pSpecies->getEmigParams(0,1);
				emigScales scale = pSpecies->getEmigScales();
				outPar << "yes" << endl;
				if (dem.stageStruct) {
					outPar << emigstage << emig.emigStage << endl;
				}
				outPar << initprob << "mean: " << "females " << eparams0.d0Mean
					<< "  males " << eparams1.d0Mean << endl;
				outPar << initprob << "s.d.: " << "females " << eparams0.d0SD
					<< "  males " << eparams1.d0SD << endl;
				outPar << initprob << "scaling factor: " << scale.d0Scale
					<< endl;
			}
			else{
				outPar << "no" << endl;
				outPar << "EMIGRATION PROB.: \tfemales "<< pSpecies->getEmigD0(0,0)
					<<"\t males " << pSpecies->getEmigD0(0,1) << endl;
			}
		}
	}
	else { // !emig.sexDep
		outPar << sexdept << "no" << endl;
		if (emig.stgDep) {
			outPar << stgdept << "yes" << endl;
			outPar << indvar << "no" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage " << i << ": \t" << "EMIGRATION PROB.: "
					<< pSpecies->getEmigD0(i,0) << endl;
			}
		}
		else { // !emig.stgDep
			outPar << stgdept << "no" << endl;
			outPar << indvar;
			if (emig.indVar) {
				eparams0 = pSpecies->getEmigParams(0,0);
				emigScales scale = pSpecies->getEmigScales();
				outPar << "yes" << endl;
				if (dem.stageStruct) {
					outPar << emigstage << emig.emigStage << endl;
				}
				outPar << initprob << "mean: " << eparams0.d0Mean << endl;
				outPar << initprob << "s.d.: " << eparams0.d0SD << endl;
				outPar << initprob << "scaling factor: " << scale.d0Scale << endl;
			}
			else {
				outPar << "no" << endl;
				outPar << "EMIGRATION PROB.:\t" << pSpecies->getEmigD0(0,0) << endl;
			}
		}
	}
}

// Transfer

outPar << endl << "DISPERSAL - TRANSFER: \t";

if (trfr.moveModel) {
	bool straigtenPath;
	if (trfr.moveType == 1) { // SMS
		trfrSMSTraits move = pSpecies->getSMSTraits();
		straigtenPath = move.straigtenPath;
		if (trfr.costMap) {
			outPar << "SMS\tcosts from imported cost map" << endl;
		}
		else {
			outPar << "SMS\tcosts:" << endl;
			if (!ppLand.generated && ppLand.rasterType == 0) {
				for (int i = 0; i < ppLand.nHab; i++)
					outPar << "\thab. " << pLandscape->getHabCode(i) << "\t"
						<< pSpecies->getHabCost(i) << endl;
			}
			else {
				for (int i = 0; i < ppLand.nHab; i++)
					outPar << "\thab. " << i << "\t"
						<< pSpecies->getHabCost(i) << endl;
			}
		}
		string pr = "PERCEPTUAL RANGE";
		outPar << pr << ":        " << move.pr << endl;
		outPar << pr << " METHOD: " << move.prMethod << endl;
		if (!trfr.indVar) outPar << "DIRECTIONAL PERSISTENCE: " << move.dp << endl;
		outPar << "MEMORY SIZE: " << move.memSize << endl;
		outPar << "GOAL TYPE:   " << move.goalType << endl;
		if (!trfr.indVar) {
			if (move.goalType == 2) { //  dispersal bias
				outPar << "GOAL BIAS:   " << move.gb << endl;
				outPar << "ALPHA DB:    " << move.alphaDB << endl;
				outPar << "BETA DB:     " << move.betaDB << endl;
			}
		}
		if (trfr.indVar) {
			trfrSMSParams s = pSpecies->getSMSParams(0,0);
			outPar << indvar << "yes " << endl;
			outPar << "DP mean: " << s.dpMean << "  s.d.: " << s.dpSD
				<< "  scaling factor: " << s.dpScale << endl;
			outPar << "GB mean: " << s.gbMean << "  s.d.: " << s.gbSD
				<< "  scaling factor: " << s.gbScale << endl;
			if (move.goalType == 2) { //  dispersal bias
				outPar << "Alpha DB mean: " << s.alphaDBMean << "  s.d.: " << s.alphaDBSD
					<< "  scaling factor: " << s.alphaDBScale << endl;
				outPar << "Beta DB mean:  " << s.betaDBMean << "  s.d.: " << s.betaDBSD
					<< "  scaling factor: " << s.betaDBScale << endl;
			}
		}
		else {
			outPar << indvar << "no " << endl;
		}
	}
	else { // CRW
		trfrCRWTraits move = pSpecies->getCRWTraits();
		straigtenPath = move.straigtenPath;
		outPar << "CRW" << endl;
		string lgth = "STEP LENGTH (m) ";
		string corr = "STEP CORRELATION";
		if (trfr.indVar) {
			trfrCRWParams m = pSpecies->getCRWParams(0,0);
			outPar << indvar << "yes" << endl;
			outPar << lgth << " mean: " << m.stepLgthMean;
			outPar << "  s.d.: " << m.stepLgthSD;
			outPar << "  scaling factor: " << m.stepLScale << endl;
			outPar << corr << " mean: " << m.rhoMean;
			outPar << "  s.d.: " << m.rhoSD;
			outPar << "  scaling factor: " << m.rhoScale << endl;
		}
		else {
			outPar << indvar << "no" << endl;
			outPar << lgth << ": " << move.stepLength << endl;
			outPar << corr << ": " << move.rho << endl;
		}
	}
	outPar << "STRAIGHTEN PATH AFTER DECISION NOT TO SETTLE: ";
	if (straigtenPath) outPar << "yes" << endl;
	else outPar << "no" << endl;
	outPar << "STEP MORTALITY:\t" << endl;
	if (trfr.habMort) 
	{
		outPar << "habitat dependent:\t" << endl;
		if (!ppLand.generated && ppLand.rasterType == 0) {
			for (int i = 0; i < ppLand.nHab; i++)
				outPar << "\thab. " << pLandscape->getHabCode(i) << "\t"
					<< pSpecies->getHabMort(i) << endl;
		}
		else{
			for (int i = 0; i < ppLand.nHab; i++)
				outPar << "\thab. " << i << "\t"
					<< pSpecies->getHabMort(i) << endl;
		}
	}
	else 
	{
		trfrCRWTraits move = pSpecies->getCRWTraits();
		outPar << "constant " << move.stepMort << endl;
	}
} // end of movement process
else { // kernel
	string meandist = "MEAN DISTANCE";
	string probkern = "PROB. KERNEL I";
	trfrKernTraits kern0,kern1;
	trfrKernParams k0,k1;
	outPar << "dispersal kernel" << endl << "TYPE: \t";
	if (trfr.twinKern) outPar << "double ";
	outPar << "negative exponential" << endl;

	if (trfr.sexDep) {
		outPar << sexdept << "yes" << endl;
		if (trfr.stgDep) {
			outPar << stgdept << "yes" << endl;
			outPar << indvar << "no" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage " << i << ":" << endl;
				kern0 = pSpecies->getKernTraits(i,0);
				kern1 = pSpecies->getKernTraits(i,1);
				outPar << meandist << " I: \tfemales "<< kern0.meanDist1 << " \tmales " << kern1.meanDist1 << endl;
				if (trfr.twinKern)
				{
					outPar << meandist << " II: \tfemales "<< kern0.meanDist2 << " \tmales " << kern1.meanDist2 << endl;
					outPar << probkern << ": \tfemales "<< kern0.probKern1 << " \tmales " << kern1.probKern1 << endl;
				}
			}
		}
		else { // !trfr.stgDep
			outPar << stgdept << "no" << endl;
			outPar << indvar;
			if (trfr.indVar) {
				k0 = pSpecies->getKernParams(0,0);
				k1 = pSpecies->getKernParams(0,1);
				outPar << "yes" << endl;
				outPar << meandist << " I  (mean): \tfemales " << k0.dist1Mean
					<< " \tmales " << k1.dist1Mean << endl;
				outPar << meandist << " I  (s.d.): \tfemales " << k0.dist1SD
					<< " \tmales " << k1.dist1SD << endl;
				outPar << meandist << " I  (scaling factor): \tfemales " << k0.dist1Scale
					<< " \tmales " << k1.dist1Scale << endl;
				if (trfr.twinKern)
				{
					outPar << meandist << " II (mean): \tfemales " << k0.dist2Mean
						<< " \tmales " << k1.dist2Mean << endl;
					outPar << meandist << " II (s.d.): \tfemales " << k0.dist2SD
						<< " \tmales " << k1.dist2SD << endl;
					outPar << meandist << " II (scaling factor): \tfemales " << k0.dist2Scale
						<< " \tmales " << k1.dist2Scale << endl;
					outPar << probkern << "   (mean): \tfemales " << k0.PKern1Mean
						<< " \tmales " << k1.PKern1Mean << endl;
					outPar << probkern << "   (s.d.): \tfemales " << k0.PKern1SD
						<< " \tmales " << k1.PKern1SD << endl;
					outPar << probkern << "   (scaling factor): \tfemales " << k0.PKern1Scale
						<< " \tmales " << k1.PKern1Scale << endl;
				}
			}
			else {
				outPar << "no" << endl;
				kern0 = pSpecies->getKernTraits(0,0);
				kern1 = pSpecies->getKernTraits(0,1);
				outPar << meandist << " I: \tfemales "<< kern0.meanDist1 << " \tmales " << kern1.meanDist1 << endl;
				if (trfr.twinKern)
				{
					outPar << meandist << " II: \tfemales "<< kern0.meanDist2 << " \tmales " << kern1.meanDist2 << endl;
					outPar << probkern << ": \tfemales "<< kern0.probKern1 << " \tmales " << kern1.probKern1 << endl;
				}
			}
		}
	}
	else { // !trfr.sexDep
		outPar << sexdept << "no" << endl;
		if (trfr.stgDep) {
			outPar << stgdept << "yes" << endl;
			outPar << indvar << "no" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				kern0 = pSpecies->getKernTraits(i,0);
				outPar << "stage " << i << ": \t"  << meandist << " I: " << kern0.meanDist1;
				if (trfr.twinKern)
				{
					outPar << " \t" << meandist << " II: " << kern0.meanDist2;
					outPar << " \t" << probkern << ": " << kern0.probKern1;
				}
				outPar << endl;
			}
		}
		else { // !trfr.stgDep
			outPar << stgdept << "no" << endl;
			outPar << indvar;
			if (trfr.indVar) {
				k0 = pSpecies->getKernParams(0,0);
				outPar << "yes" << endl;
				outPar << meandist << " I  (mean): " << k0.dist1Mean
					<< " \t(s.d.): " << k0.dist1SD
					<< " \t(scaling factor): " << k0.dist1Scale << endl;
				if (trfr.twinKern) 
				{
					outPar << meandist << " II (mean): " << k0.dist2Mean
						<< " \t(s.d.): " << k0.dist2SD
						<< " \t(scaling factor): " << k0.dist2Scale << endl;
					outPar << probkern << "   (mean): " << k0.PKern1Mean
						<< " \t(s.d.): " << k0.PKern1SD
						<< " \t(scaling factor): " << k0.PKern1Scale << endl;
				}
			}
			else{
				outPar << "no" << endl;
				kern0 = pSpecies->getKernTraits(0,0);
				outPar << meandist << " I: \t" << kern0.meanDist1 <<endl;
				if (trfr.twinKern) 
				{
					outPar << meandist << " II: \t" << kern0.meanDist2 <<endl;
					outPar << probkern << ": \t" << kern0.probKern1 <<endl;
				}
			}
		}
	}

	outPar << "DISPERSAL MORTALITY:   ";
	trfrMortParams mort = pSpecies->getMortParams();
	if (trfr.distMort) {
		outPar << "distance-dependent" << endl;
		outPar << "SLOPE: " << mort.mortAlpha << " \tINFLECTION POINT: " << mort.mortBeta << endl;
	}
	else {
		outPar << "constant" << endl << "MORTALITY PROBABILITY: " << mort.fixedMort << endl;
	}
} // end of kernel transfer

// Settlement

outPar << endl << "DISPERSAL - SETTLEMENT:" << endl;

if (trfr.moveModel) {
	string plusmating = "+ mating requirements";
	ssteps = pSpecies->getSteps(0,0);

	outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
	outPar << "MAX. No. OF STEPS:\t ";
	if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
	else outPar << ssteps.maxSteps << endl;

	if (sett.sexDep) {
		nsexes = 2;
		outPar << sexdept << "yes" << endl;
		if (sett.stgDep) {
			nstages = sstruct.nStages;
			outPar << stgdept << "yes" << endl;
		}
		else { // !sett.stgDep
			nstages = 1;
			outPar << stgdept << "no" << endl;
		}
	}
	else { // !sett.sexDep
		nsexes = 1;
		outPar << sexdept << "no" << endl;
		if (sett.stgDep) {
			nstages = sstruct.nStages;
			outPar << stgdept << "yes" << endl;
		}
		else { // !sett.stgDep
			nstages = 1;
			outPar << stgdept << "no" << endl;
		}
	}
	for (int sx = 0; sx < nsexes; sx++) {
		if (sett.sexDep) {
			if (sx == 0) outPar << "FEMALES:" << endl;
			else outPar << "MALES:" << endl;
		}
		outPar << "SETTLE IF: ";
		for (int i = 0; i < nstages; i++) {
			if (dem.stageStruct && nstages > 1) outPar << "stage " << i << ": " << endl;
			outPar << "find a suitable cell/patch ";
			srules = pSpecies->getSettRules(i,sx);
			if (srules.densDep) {
				settleDD = pSpecies->getSettTraits(i,sx);
				outPar << "+ density dependence ";
				if (srules.findMate) outPar << plusmating;
				outPar << endl;
				if (!sett.indVar) {
					outPar << "S0: " << settleDD.s0 << "  AlphaS: " << settleDD.alpha
						<< "  BetaS: " << settleDD.beta << endl;
				}
			}
			else {
				if (srules.findMate) outPar << plusmating << endl;
				else outPar << "(not the natal one)" << endl;
			}
			if (dem.stageStruct) {
				ssteps = pSpecies->getSteps(i,sx);
				outPar << "MAX. No. OF STEPS/YEAR:\t ";
				if (ssteps.maxStepsYr == 99999999) outPar << "not applied" << endl;
				else outPar << ssteps.maxStepsYr << endl;
			}
		}
	}
	if (sett.indVar) {
		settParams sparams0;
		outPar << "DENSITY DEPENDENCE + " << indvar << "yes" << endl;
		for (int sex = 0; sex < nsexes; sex++) {
			if (sett.sexDep) {
				if (sex == 0) outPar << "FEMALES:" << endl;
				else outPar << "MALES:" << endl;
			}
			sparams0 = pSpecies->getSettParams(0,sex);				
			settScales scale = pSpecies->getSettScales();
			outPar << "S0     - mean: " << sparams0.s0Mean 		<< "  s.d.: " << sparams0.s0SD
				<< "  scaling factor: " << scale.s0Scale << endl;
			outPar << "AlphaS - mean: " << sparams0.alphaSMean 	<< "  s.d.: " << sparams0.alphaSSD
				<< "  scaling factor: " << scale.alphaSScale << endl;
			outPar << "BetaS  - mean: " << sparams0.betaSMean 	<< "  s.d.: " << sparams0.betaSSD
				<< "  scaling factor: " << scale.betaSScale << endl;
		}
	}
}
else { // kernel-based transfer
	string notsuit = "IF THE ARRIVAL CELL/PATCH IS UNSUITABLE: ";
	string rchoose = " randomly choose a suitable neighb. cell/patch or ";
	string matereq = "MATING REQUIREMENTS: ";
	if (sett.sexDep) {
		nsexes = 2;
		outPar << sexdept << "yes" << endl;
		if (sett.stgDep) {
			nstages = sstruct.nStages;
			outPar << stgdept << "yes" << endl;
			outPar << notsuit << endl;
		}
		else {
			nstages = 1;
			outPar << stgdept << "no" << endl;
		}
	}
	else {
		nsexes = 1;
		outPar << sexdept << "no" << endl;
		if (sett.stgDep) {
			nstages = sstruct.nStages;
			outPar << stgdept << "yes" << endl;
			outPar << notsuit << endl;
		}
		else {
			nstages = 1;
			outPar << stgdept << "no" << endl;
			outPar << notsuit;
		}
	}
	for (int i = 0; i < nstages; i++) {
		if (sett.stgDep) {
			outPar << "stage " << i << ":" << endl;
		}
		for (int sx = 0; sx < nsexes; sx++) {
			if (sett.sexDep) {
				if (sx == 0) outPar << "FEMALES: ";
				else outPar << "MALES:   ";
				if (!sett.stgDep) outPar << notsuit;
			}
			srules = pSpecies->getSettRules(i,sx);
			if (srules.go2nbrLocn) {
				outPar << rchoose;
				if (srules.wait) outPar << "wait" << endl;
				else outPar << "die" << endl;
			}
			else {
				if (srules.wait) outPar << "wait" << endl;
				else outPar << "die" << endl;
			}
			outPar << matereq;
			if (srules.findMate) outPar << "yes" <<endl;
			else outPar << "no" <<endl;
		}
	}
}

// Genetics

outPar << endl << "GENETICS:" << endl;
int nspptraits = pSpecies->getNTraits();
outPar << "No. of variable traits:  " << nspptraits << endl;

genomeData d = pSpecies->getGenomeData();
if (emig.indVar || trfr.indVar || sett.indVar || d.neutralMarkers)
{
	if (d.diploid) outPar << "DIPLOID" << endl; else outPar << "HAPLOID" << endl;
	int nchromosomes = pSpecies->getNChromosomes();
	outPar << "No. of chromosomes:      " << nchromosomes;
	if (d.trait1Chromosome) {
		outPar << endl << "No. of loci/chromosome:  " << d.nLoci << endl;
	}
	else {
		outPar << " (chrom:loci)";
		for (int i = 0; i < nchromosomes; i++) {
			outPar << "  " << i << ":" << pSpecies->getNLoci(i);
		}
		outPar << endl;
	}
	outPar << "Mutation probability:    " << d.probMutn << endl;
	outPar << "Crossover probability:   " << d.probCrossover << endl;
	outPar << "Initial allele s.d.:     " << d.alleleSD << endl;
	outPar << "Mutation s.d.:           " << d.mutationSD << endl;
	if (d.neutralMarkers) {
		outPar << "NEUTRAL MARKERS ONLY" << endl;
	}
	else {
		if (!d.trait1Chromosome) {
			traitAllele allele;
			outPar << "TRAIT MAPPING:" << endl;
			outPar << "Architecture file:     " << genfilename << endl;
			int ntraitmaps = pSpecies->getNTraitMaps();
			outPar << "No. of traits defined: " << ntraitmaps << endl;
			for (int i = 0; i < ntraitmaps; i++) {
				int nalleles = pSpecies->getNTraitAlleles(i);
				outPar << "Trait " << i << ": (" << pSpecies->getTraitName(i)
					<< ") alleles: " << nalleles << " (chrom:locus)";
				for (int j = 0; j < nalleles; j++) {
					allele = pSpecies->getTraitAllele(i,j);
					outPar << "  " << allele.chromo << ":" << allele.locus;
				}
				outPar << endl;
			}
			if (ntraitmaps < nspptraits) { // list undefined traits
				outPar << "WARNING - the following traits were not defined"
					<< " in the genetic architecture file:" << endl;
				for (int i = ntraitmaps; i < nspptraits; i++) {
					outPar << "Trait " << i << ": (" << pSpecies->getTraitName(i)
						<< ") all individuals have mean phenotype" << endl;
				}
			}
			int nneutral = pSpecies->getNNeutralLoci();
			if (nneutral > 0) {
				outPar << "Neutral loci: " << nneutral << " (chrom:locus)";
				for (int i = 0; i < nneutral; i++) {
					allele = pSpecies->getNeutralAllele(i);
					outPar << "  " << allele.chromo << ":" << allele.locus;
				}
				outPar << endl;
			}
			if (d.pleiotropic)
				outPar << "Genome exhibits pleiotropy" << endl;
		}
	}
}

// Initialisation

initParams init = paramsInit->getInit();
outPar << endl << "INITIALISATION CONDITIONS:" << endl;
switch (init.seedType) {
case 0:
	outPar << "Free initialisation: \t" ;
	switch (init.freeType) {
	case 0:
		outPar << "Random \t" ;
		outPar << "No. of cells/patches: " << init.nSeedPatches << endl;
		break;
	case 1:
		outPar << "all suitable cells/patches" << endl;
		break;
	case 2:
		outPar << "manually selected cells/patches" << endl;
		break;
	}
	break;
case 1:
	outPar << "From species distribution: \t" << endl;
	switch (init.spDistType) {
		case 0:
			outPar << "all presence cells/patches" << endl;
			break;
		case 1:
			outPar << "some random presence cells/patches" << endl;
			break;
		case 2:
			outPar << "all cells/patches within selected distribution cells" << endl;
			break;
	}
	break;
case 2:
	outPar << "From initial individuals file: " << paramsSim->getDir(1) + init.indsFile << endl;
	break;
case 3:
	outPar << "From file" << endl;
	break;
}
if (init.seedType != 2) {
	outPar << "INITIAL NO. OF INDIVIDUALS: \t";
	switch (init.initDens) {
	case 0:
		outPar << "at carrying capacity" << endl;
		break;
	case 1:
		outPar << "at half carrying capacity" << endl;
		break;
	case 2:
		if (ppLand.patchModel) {
			outPar << init.indsHa << " individuals per ha" << endl;
		}
		else {
			outPar << init.indsCell << " individuals per cell" << endl;
		}
		break;
	}
	if (dem.stageStruct) {
		outPar << "INITIAL STAGE PROPORTIONS:" << endl;
		for (int i = 1; i < sstruct.nStages; i++) {
			outPar << "stage " << i << ": " << paramsInit->getProp(i) << " \t";
		}
		outPar << endl;
		outPar << "Initial age distribution: ";
		switch (init.initAge) {
		case 0:
			outPar << "lowest possible age";
			break;
		case 1:
			outPar << "randomised";
			break;
		case 2:
			outPar << "quasi-equilibrium";
			break;
		}
		outPar << endl;
	}
	outPar << "GEOGRAPHICAL CONSTRAINTS (cell numbers): " << endl;
	outPar << "min X: " << init.minSeedX << " max X: " << init.maxSeedX << endl;
	outPar << "min Y: " << init.minSeedY << " max Y: " << init.maxSeedY << endl;
//	if (init.seedType != 1 && init.freeType < 2 && init.initFrzYr > 0) {
//		outPar << "Freeze initial range until year " << init.initFrzYr << endl;
//	}
	if (init.seedType == 0 && init.freeType < 2) {
		if (init.initFrzYr > 0) {
			outPar << "Freeze initial range until year " << init.initFrzYr << endl;
		}
		if (init.restrictRange) {
			outPar << "Restrict range to northern " << init.restrictRows
				<< " rows every " << init.restrictFreq << " years" << endl;
			if (init.finalFrzYr < sim.years) {
				outPar << "Freeze range at year " << init.finalFrzYr << endl;
			}
		}
	}
}

outPar << endl << "OUTPUTS:" << endl;
if (sim.outRange) {
	outPar << "Range - every " << sim.outIntRange << " year";
	if (sim.outIntRange > 1) outPar << "s";
//	if (sim.outStartRange > 0) outPar << " starting year " << sim.outStartRange;
	outPar << endl;
}
if (sim.outOccup) {
	outPar << "Occupancy - every " << sim.outIntOcc << " year";
	if (sim.outIntOcc > 1) outPar << "s";
//	if (sim.outStartOcc > 0) outPar << " starting year " << sim.outStartOcc;
	outPar << endl;
}
if (sim.outPop) {
	outPar << "Populations - every " << sim.outIntPop << " year";
	if (sim.outIntPop > 1) outPar << "s";
	if (sim.outStartPop > 0) outPar << " starting year " << sim.outStartPop;
	outPar << endl;
}
if (sim.outInds) {
	outPar << "Individuals - every " << sim.outIntInd << " year";
	if (sim.outIntInd > 1) outPar << "s";
	if (sim.outStartInd > 0) outPar << " starting year " << sim.outStartInd;
	outPar << endl;
}
if (sim.outGenetics) {
	outPar << "Genetics - every " << sim.outIntGenetic << " year";
	if (sim.outIntGenetic > 1) outPar << "s";
	if (sim.outStartGenetic > 0) outPar << " starting year " << sim.outStartGenetic;
	if (dem.stageStruct) {
		switch (sim.outGenType) {
		case 0:
			outPar << " - juveniles only";
			break;
		case 1:
			outPar << " - all individuals";
			break;
		case 2:
			outPar << " - adults only";
			break;
		}
	}
	if (sim.outGenXtab) outPar << " (as cross table)";
	outPar << endl;
}

if (sim.outTraitsCells) {
	outPar << "Traits per ";
	if (ppLand.patchModel) outPar << "patch"; else outPar << "cell";
	outPar << " - every " << sim.outIntTraitCell << " year";
	if (sim.outIntTraitCell > 1) outPar << "s";
	if (sim.outStartTraitCell > 0) outPar << " starting year " << sim.outStartTraitCell;
	outPar << endl;
}
if (sim.outTraitsRows) {
	outPar << "Traits per row - every " << sim.outIntTraitRow << " year";
	if (sim.outIntTraitRow > 1) outPar << "s";
	if (sim.outStartTraitRow > 0) outPar << " starting year " << sim.outStartTraitRow;
	outPar << endl;
}
if (sim.outConnect) {
	outPar << "Connectivity matrix - every " << sim.outIntConn << " year";
	if (sim.outIntConn > 1) outPar << "s";
	if (sim.outStartConn > 0) outPar << " starting year " << sim.outStartConn;
	outPar << endl;
}
	if (sim.outPaths) {
		outPar << "SMS paths - every " << sim.outIntPaths << " year";
		if (sim.outIntPaths > 1) outPar << "s";
		if (sim.outStartPaths > 0) outPar << " starting year " << sim.outStartPaths;
		outPar << endl;
	}
outPar << "SAVE MAPS: ";
if (sim.saveMaps) {
	outPar << "yes - every " << sim.mapInt << " year";
	if (sim.mapInt > 1) outPar << "s";
	outPar << endl;
}
else outPar << "no" << endl;
outPar << "SAVE TRAITS MAPS: ";
if (sim.saveTraitMaps) {
	outPar << "yes - every " << sim.traitInt << " year";
	if (sim.traitInt > 1) outPar << "s";
	outPar << endl;
}
else outPar << "no" << endl;
if (trfr.moveModel && trfr.moveType == 1) {
	outPar << "SMS HEAT MAPS: ";
	if (sim.saveVisits) outPar << "yes" << endl;
	else outPar << "no" << endl;
}

outPar.close(); outPar.clear();

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
