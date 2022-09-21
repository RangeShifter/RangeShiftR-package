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

#include "Model.h"

ofstream outPar;

//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma package(smart_init)
#endif
//---------------------------------------------------------------------------
#if RS_ABC
int RunModel(Landscape *pLandscape,int seqsim,ABCmaster* pABCmaster)
#else
#if RS_RCPP && !R_CMD
#if RS_THREADSAFE
Rcpp::List RunModel(Landscape *pLandscape,int seqsim,Rcpp::S4 ParMaster)
#else
Rcpp::List RunModel(Landscape *pLandscape,int seqsim)
#endif // RS_THREADSAFE
#else
int RunModel(Landscape *pLandscape,int seqsim)
#endif // RS_RCPP && !R_CMD
#endif // RS_ABC
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
#if VIRTUALECOLOGIST
virtParams virt = paramsSim->getVirtParams();
#endif

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
#if SEASONAL
		pLandscape->allocatePatches(pSpecies,dem.nSeasons);
#else
		pLandscape->allocatePatches(pSpecies);
#endif // SEASONAL 
	}
	pComm = new Community(pLandscape); // set up community
	// set up a sub-community associated with each patch (incl. the matrix)
	pLandscape->updateCarryingCapacity(pSpecies,0,0);
//	if (ppLand.rasterType <= 2 && ppLand.dmgLoaded) 
//		pLandscape->updateDamageIndices();
#if SPATIALDEMOG
	if (ppLand.rasterType == 2 && ppLand.spatialdemog)
		pLandscape->updateDemoScalings(0); // TODO
#endif // SPATIALDEMOG 

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
#if VIRTUALECOLOGIST
	if (sim.virtualEcologist) {
		pVirt = new VirtualEcologist(paramsSim); // set up virtual ecologist
		if (virt.patchMethod != 2) pVirt->samplePatches(pLandscape,pComm);
	}
#endif
}

#if VCL
// NOTE: GRADIENT AND TRAIT MAPS ARE AVAILABLE FOR GUI VERSION ONLY
SetupVisualOutput();
#endif

#if RS_ABC
int nABCyears = (int)pABCmaster->NYears();
int yearABC = -1;
int ixABC;
bool abcYear;
#endif

#if RS_RCPP && !R_CMD
Rcpp::List list_outPop;
#endif

// Loop through replicates
for (int rep = 0; rep < sim.reps; rep++) {
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): starting simulation=" << sim.simulation << " rep=" << rep << endl;
#endif
#if RS_RCPP && !R_CMD
		Rcpp::Rcout << endl << "starting replicate " << rep << endl;
#else
#if BATCH
		cout << endl << "starting replicate " << rep << endl;
#endif
#endif

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
#if LINUX_CLUSTER
			pComm->addSubComm(ppp.pPatch,ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
#else
			SubCommunity *pSubComm = pComm->addSubComm(ppp.pPatch,ppp.patchNum); // SET UP ALL SUB-COMMUNITIES
#endif
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
#if VIRTUALECOLOGIST
		if (sim.virtualEcologist) {
			pVirt = new VirtualEcologist(paramsSim); // set up virtual ecologist
			if (virt.patchMethod != 2) pVirt->samplePatches(pLandscape,pComm);
		}
#endif
		MemoLine("...completed...");
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): finished generating populations" << endl;
#endif
	}
#if VCL
	else {
		if (sim.batchMode && (sim.saveMaps || sim.saveVisits)) {
			pLandscape->setLandMap();
			pLandscape->drawLandscape(rep,0,0);
		}
	}
#endif
	if (init.seedType == 0 && init.freeType < 2 && init.initFrzYr > 0) {
		// restrict available landscape to initialised region
		pLandscape->setLandLimits(init.minSeedX,init.minSeedY,
			init.maxSeedX,init.maxSeedY);
	}
	else {
		pLandscape->resetLandLimits();
	}

	filesOK = true;
#if RS_THREADSAFE
	if(init.seedType==2 && init.indsFile=="NULL"){ // initialisation from InitInds list of dataframes
		if(rep > 0){
			int error_init = 0;
			Rcpp::S4 InitParamsR("InitialisationParams");
			InitParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("init"));
			Rcpp::List InitIndsList = Rcpp::as<Rcpp::List>(InitParamsR.slot("InitIndsList"));
			error_init = ReadInitIndsFileR(0, pLandscape, Rcpp::as<Rcpp::DataFrame>(InitIndsList[rep]));
			if(error_init>0) {
				MemoLine("UNABLE TO GET NEXT INITIALISATION DATAFRAME");
				filesOK = false;
			}
		}
	}
#endif
	if (rep == 0) {
		// open output files
		if (sim.outRange) { // open Range file
			if (!pComm->outRangeHeaders(pSpecies,ppLand.landNum)) {
				MemoLine("UNABLE TO OPEN RANGE FILE");
				filesOK = false;
			}
#if RS_CONTAIN
			if (ppLand.dmgLoaded) {
				if (!pLandscape->outSummDmgHeaders(ppLand.landNum)) {
					MemoLine("UNABLE TO OPEN SUMMARY DAMAGE FILE");
					filesOK = false;
				}				
			}
#endif // RS_CONTAIN
		}
		if (sim.outOccup && sim.reps > 1)
			if (!pComm->outOccupancyHeaders(0)) {
				MemoLine("UNABLE TO OPEN OCCUPANCY FILE(S)");
				filesOK = false;
			}
#if RS_RCPP
		if (sim.outPop && sim.CreatePopFile) {
#else
		if (sim.outPop) {
#endif
			// open Population file
			if (!pComm->outPopHeaders(pSpecies,ppLand.landNum)) {
				MemoLine("UNABLE TO OPEN POPULATION FILE");
				filesOK = false;
			}
#if RS_CONTAIN
			if (pCull->isCullApplied()) {
				// also open Cull file
				if (!pComm->outCullHeaders(pSpecies,ppLand.landNum)) {
					MemoLine("UNABLE TO OPEN MANAGEMENT CULL FILE");
					filesOK = false;
				}
			}
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
		if (ppLand.dmgLoaded && sim.outDamage) {
			if (!pLandscape->outDamageHeaders(ppLand.landNum)) {
				MemoLine("UNABLE TO OPEN DAMAGE INDICES FILE");
				filesOK = false;
			}				
		}
#endif // RS_CONTAIN
#if VIRTUALECOLOGIST
		if (sim.virtualEcologist) {
			if (!pVirt->outLandGenHeaders(ppLand.landNum,ppLand.patchModel)) {
#if RSDEBUG
DEBUGLOG << "RunModel(): UNABLE TO OPEN LANDSCAPE GENETICS FILE" << endl;
#endif
				MemoLine("UNABLE TO OPEN LANDSCAPE GENETICS FILE");
				filesOK = false;
			}
			if (virt.outGenomes) {
				if (!pVirt->outGenSamplesHeaders(ppLand.landNum,ppLand.patchModel)) {
#if RSDEBUG
DEBUGLOG << "RunModel(): UNABLE TO OPEN GENETIC SAMPLES FILE" << endl;
#endif
					MemoLine("UNABLE TO OPEN GENETIC SAMPLES FILE");
					filesOK = false;
				}
			}
		}
#endif // VIRTUALECOLOGIST 
#if PEDIGREE
		if (!pComm->outGroupHeaders(0)) {
			MemoLine("UNABLE TO OPEN GROUPS FILE");
			filesOK = false;		
		}
#endif
	} // rep==0
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
#if RS_CONTAIN
			if (ppLand.dmgLoaded) pLandscape->outSummDmgHeaders(-999);				
#endif // RS_CONTAIN 
		}
		if (sim.outOccup && sim.reps > 1)
			pComm->outOccupancyHeaders(-999);
		if (sim.outPop) {
			pComm->outPopHeaders(pSpecies,-999);
#if RS_CONTAIN
			if (pCull->isCullApplied()) {
				pComm->outCullHeaders(pSpecies,-999);				
			}
#endif // RS_CONTAIN
		}
		if (sim.outTraitsCells)
			pComm->outTraitsHeaders(pSpecies,-999);
		if (sim.outTraitsRows)
			pComm->outTraitsRowsHeaders(pSpecies,-999);
		if (sim.outConnect && ppLand.patchModel)
			pLandscape->outConnectHeaders(-999);
#if RS_CONTAIN
		if (ppLand.dmgLoaded && sim.outDamage) pLandscape->outDamageHeaders(-999);				
#endif // RS_CONTAIN 
#if VIRTUALECOLOGIST
		if (sim.virtualEcologist) {
			pVirt->outLandGenHeaders(-999,false);
			if (virt.outGenomes) pVirt->outGenSamplesHeaders(-999,false);
		}
#endif
#if PEDIGREE
		pComm->outGroupHeaders(-999); // close groups file
#endif
#if RS_RCPP && !R_CMD
		return Rcpp::List::create(Rcpp::Named("Errors") = 666);
#else
		return 666;
#endif
	}

#if VCL
	ResetVisualOutput();
#endif

	if (env.stoch && !env.local) {
#if BUTTERFLYDISP
		if (env.fromFile) {
			pLandscape->readGlobalStoch(sim.years+1,envstochfilename);
		}
		else {
			// create time series in case of global environmental stochasticity
			pLandscape->setGlobalStoch(sim.years+1);
		}
#else
		// create time series in case of global environmental stochasticity
		pLandscape->setGlobalStoch(sim.years+1);
#endif
#if VCL
		if (v.viewGraph) pLandscape->drawGlobalStoch(sim.years+1);
#endif
	}

	if (grad.gradient) { // set up environmental gradient
#if SEASONAL
		pLandscape->setEnvGradient(pSpecies,dem.nSeasons,true); 
#else
		pLandscape->setEnvGradient(pSpecies,true);
#endif // SEASONAL 
	}

#if RS_ABC
	if (ppLand.patchModel && (sim.outConnect || pABCmaster->obsConnectivity()))
#else
	if (sim.outConnect && ppLand.patchModel)
#endif
		pLandscape->createConnectMatrix();

	// variables to control dynamic landscape
	landChange landChg; landChg.chgnum = 0; landChg.chgyear = 999999;
	if (!ppLand.generated && ppLand.dynamic) {
		landChg = pLandscape->getLandChange(0); // get first change year
	}

#if PEDIGREE
	// create relatedness matrix ready for initial population
	Pedigree *pPed = new Pedigree(sim.relMatSize);
#endif
	
	// set up populations in the community
	pLandscape->updateCarryingCapacity(pSpecies,0,0);
#if RS_CONTAIN
//	pLandscape->resetPrevDamage();
	if (ppLand.rasterType <= 2 && ppLand.dmgLoaded)
		pLandscape->updateDamageIndices();
#endif // RS_CONTAIN
#if SPATIALDEMOG
	if (ppLand.rasterType == 2 && ppLand.spatialdemog)
		pLandscape->updateDemoScalings(0);
#endif // SPATIALDEMOG 
#if RSDEBUG
DEBUGLOG << "RunModel(): completed updating carrying capacity" << endl;
#endif
//	if (init.seedType != 2) {
#if PEDIGREE
		pComm->initialise(pSpecies,pPed,-1);
#else
		pComm->initialise(pSpecies,-1);
#endif
//	}
	bool updateland = false;
	int landIx = 0; // landscape change index

#if RSDEBUG
DEBUGLOG << "RunModel(): completed initialisation, rep=" << rep
	<< " pSpecies=" << pSpecies << endl;
#endif
#if BATCH
#if RS_RCPP && !R_CMD
	Rcpp::Rcout << "RunModel(): completed initialisation " << endl;
#else
	cout << "RunModel(): completed initialisation " << endl;
#endif
#endif

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
#if RS_RCPP
		// open a new movement paths file for each replicate
		if (sim.outPaths)
			pLandscape->outPathsHeaders(rep,0);
#endif

#if RS_ABC
	ixABC = 0;
	if (ixABC < nABCyears) yearABC = pABCmaster->getYear(ixABC++);
#endif

#if SEASONAL && PARTMIGRN
	// determine if specific extreme events will occur
	bool extremeEvents = false;
//	int prevYear,prevSeason;
//	prevYear = prevSeason = 99999999;
	extEvent eEvent;
	int nEvents = pLandscape->numExtEvents();
	int ixevent = 0;
	if (nEvents) {
		extremeEvents = true;
		// read first event
		eEvent = pLandscape->getExtEvent(ixevent++); 
	}
#endif  

	// years loop
	MemoLine("...running...");
	for (yr = 0; yr < sim.years; yr++) {
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): starting simulation=" << sim.simulation 
	<< " rep=" << rep << " yr=" << yr << endl;
#endif
#if RS_RCPP && !R_CMD
		Rcpp::checkUserInterrupt();
#endif
		bool updateCC = false;
		if (yr < 4
		|| (yr < 31 && yr%10 == 0)
		|| (yr < 301 && yr%100 == 0)
		|| (yr < 3001 && yr%1000 == 0)
		|| (yr < 30001 && yr%10000 == 0)
		|| (yr < 300001 && yr%100000 == 0)
		|| (yr < 3000001 && yr%1000000 == 0)
		) {
#if RS_RCPP && !R_CMD
			Rcpp::Rcout << "starting year " << yr << "..." << endl;
#else
			cout << "starting year " << yr << endl;
#endif
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
#if SEASONAL
						commStats s = pComm->getStats(0);
#else
						commStats s = pComm->getStats();
#endif // SEASONAL 
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
#if SEASONAL
					commStats s = pComm->getStats(0);
#else
					commStats s = pComm->getStats();
#endif // SEASONAL 
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
#if SPATIALMORT
		int period = 0;
		if (sim.mortMapLoaded && yr >= paramsSim->getMortChgYear()) period = 1;
#endif
		// environmental gradient, stochasticity & local extinction
		// or dynamic landscape
		updateland = false;
		if (env.stoch || grad.gradient || ppLand.dynamic) {
			if (grad.shifting && yr > grad.shift_begin && yr < grad.shift_stop) {
				paramsGrad->incrOptY();
#if SEASONAL
				pLandscape->setEnvGradient(pSpecies,dem.nSeasons,false);
#else
				pLandscape->setEnvGradient(pSpecies,false);
#endif // SEASONAL 
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
#if VCL
			if (grad.gradient && v.viewGrad) {
				pLandscape->drawGradient();
			}
#endif
		} // end of environmental gradient, etc.

		if (updateCC) {
			pLandscape->updateCarryingCapacity(pSpecies,yr,landIx);
#if RS_CONTAIN
			if (ppLand.rasterType <= 2 && ppLand.dmgLoaded)
				pLandscape->updateDamageIndices();
#endif // RS_CONTAIN 
#if SPATIALDEMOG
			if (ppLand.rasterType == 2 && ppLand.spatialdemog)
				pLandscape->updateDemoScalings((short)landIx);
#endif // SPATIALDEMOG 
		}
#if RS_CONTAIN
		pLandscape->resetDamageLocns();
#endif // RS_CONTAIN 

#if TEMPMORT
//		trfrRules trfr = pSpecies->getTrfr();
		if (trfr.moveModel)
			if (trfr.moveType == 1 && trfr.smType == 2) { 
				// SMS with temporally variable per-step mortality
				pSpecies->updateMortality(yr);
			}
#endif // TEMPMORT 
		
#if GROUPDISP
#if PEDIGREE
#if RSDEBUG
int mdim = (int)pPed->getRelMatDim();
if (mdim > 40) mdim = 40;
Individual *pInd;
int indID;
DEBUGLOG  << endl << "RunModel(): relatedness matrix, mdim=" << mdim 
	<< " yr=" << yr << " nInds=" << pPed->getNInds() << endl;
DEBUGLOG  << "ID";
for (int i = 0; i < mdim; i++) {
	pInd = pPed->getInd(i);
	if (pInd == 0) indID = -9; else indID = pInd->getId();
	DEBUGLOG << "\t" << indID;
}
DEBUGLOG  << endl;
DEBUGLOG  << "Ma P";
for (int i = 0; i < mdim; i++) {
	DEBUGLOG << "\t" << pPed->getParentPosn(i,0);
}
DEBUGLOG  << endl;
DEBUGLOG  << "Pa P";
for (int i = 0; i < mdim; i++) {
	DEBUGLOG << "\t" << pPed->getParentPosn(i,1);
}
DEBUGLOG  << endl;
for (int i = 0; i < mdim; i++) {
	pInd = pPed->getInd(i);
	if (pInd == 0) indID = -9; else indID = pInd->getId();
	DEBUGLOG << indID;
	for (int j = 0; j < mdim; j++) {
		DEBUGLOG << "\t" << pPed->getRelMat(i,j);		
	}
	DEBUGLOG  << endl;
}
DEBUGLOG  << endl;
#endif
#endif
#endif

		
#if RS_ABC
		if (ppLand.patchModel && (sim.outConnect || pABCmaster->obsConnectivity()))
#else
		if (sim.outConnect && ppLand.patchModel)
#endif
			pLandscape->resetConnectMatrix();

		if (ppLand.dynamic && updateland) {
//			trfrRules trfr = pSpecies->getTrfr();
			if (trfr.moveModel && trfr.moveType == 1) { // SMS
				if (!trfr.costMap) pLandscape->resetCosts(); // in case habitats have changed
			}
			// apply effects of landscape change to species present in changed patches
			pComm->patchChanges();
#if RS_ABC
			pComm->dispersal(landIx,false);
#else
#if PEDIGREE
			pComm->dispersal(pPed,rep,yr,0,landIx);
#else
#if SEASONAL
			pComm->dispersal(landIx,0);
#else
#if RS_RCPP
			pComm->dispersal(landIx,yr);
#else
			pComm->dispersal(landIx);
#endif // RS_RCPP
#endif // SEASONAL 
#endif // PEDIGREE
#endif // RS_ABC
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
#if PEDIGREE
			pComm->initialise(pSpecies,pPed,yr);
#else
			pComm->initialise(pSpecies,yr);
#endif
		}

#if RS_CONTAIN
		// update landscape change index for all sub-communities
		pComm->setHabIndex(pSpecies,landIx);
#endif // RS_CONTAIN 
		
#if VIRTUALECOLOGIST
		if (sim.virtualEcologist) {
			if ((yr >= virt.outStart && yr%virt.outInt == 0)) {
				if (virt.landscapeGenetics && virt.patchMethod == 2) { // dynamic sampling
					pVirt->samplePatches(pLandscape,pComm);
				}
			}
		}
#endif

#if SEASONAL
		for(int gen = 0; gen < dem.nSeasons; gen++) // seasonal loop
#else
		for(int gen = 0; gen < dem.repSeasons; gen++) // generation loop
#endif
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

#if VCL
			if (gen == 0 && v.viewGraph) {
				// popn graphics for first generation only
				// NOTE: CURRENTLY SHOWING TOTAL OF ALL SPECIES COMBINED
				DrawPopnGraph(pComm,yr);
			}
#endif

#if RS_ABC
			if (yr == yearABC && gen == 0) abcYear = true;
			else abcYear = false;
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
#if RS_ABC
			if (!dem.stageStruct && (sim.outRange || sim.outPop || abcYear))
				RangePopOutput(pComm,rep,yr,gen,pABCmaster,abcYear);
#else
			if (!dem.stageStruct && (sim.outRange || sim.outPop))
				RangePopOutput(pComm,rep,yr,gen);
#endif
#if RS_RCPP && !R_CMD
			if ( sim.ReturnPopRaster && sim.outPop && yr >= sim.outStartPop && yr%sim.outIntPop == 0) {
				list_outPop.push_back(pComm->addYearToPopList(rep,yr), "rep" + std::to_string(rep) + "_year" + std::to_string(yr));
			}
#endif
#if VIRTUALECOLOGIST
			if (!dem.stageStruct && sim.virtualEcologist) {
				// likewise
				if ((yr >= virt.outStart && yr%virt.outInt == 0 && gen == 0)) {
//					if (virt.landscapeGenetics && virt.patchMethod == 2) { // dynamic sampling
//						pVirt->samplePatches(pComm);
//					}
					pVirt->outLandGen(pLandscape,rep,yr,gen,ppLand.patchModel,false);
				}
			}
#endif

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

#if SEASONAL
			if (sim.outConnect && ppLand.patchModel)
				pLandscape->resetConnectMatrix();
#endif

#if RS_CONTAIN
			damageparams d = pDamageParams->getDamageParams();
			if (d.timing == 0) pComm->updateDamage(pSpecies,pCull);
#endif // RS_CONTAIN 
			
			// reproduction
#if SEASONAL
			pComm->reproduction(yr,gen);
#else
#if GROUPDISP
			pComm->reproduction(pSpecies,yr);
#else
#if BUTTERFLYDISP
			if (dem.stageStruct) {
				pComm->reproduction(pSpecies,yr,0); // classical reproduction
			}
			else {
				if (dem.dispersal == 0) // dispersal during reproduction
					pComm->reproduction(pSpecies,yr,1); // mating only before dispersal
				else
					pComm->reproduction(pSpecies,yr,0); // classical reproduction
			}
#else
			pComm->reproduction(yr);
#endif // BUTTERFLYDISP
#endif // GROUPDISP
#endif // SEASONAL

			if (dem.stageStruct) {
				if (sstruct.survival == 0) { // at reproduction
#if SEASONAL
#if PARTMIGRN
					pComm->survival(gen,0,1,1); // survival of ALL stages
					if (extremeEvents) {
						Patch *pEventPatch;
						Population *pPopn;
						intptr pspecies = (intptr)pSpecies;
						while (eEvent.year == yr && eEvent.season == gen) {
							pEventPatch = pLandscape->findPatch(eEvent.patchID);
							intptr ppopn = 0;
							if (pEventPatch != 0) { // specified patch exists
								ppopn = pEventPatch->getPopn(pspecies);
								if (ppopn != 0) { // population exists
									pPopn = (Population*)ppopn;
									pPopn->extremeEvent(eEvent.probMort);
								}							
							}
#if RSDEBUG
DEBUGLOG << "RunModel(): EXTREME EVENT year=" << eEvent.year << " season=" << eEvent.season 
	<< " patch=" << eEvent.patchID << " probMort=" << eEvent.probMort 
	<< " pEventPatch=" << pEventPatch << " ppopn=" << ppopn 
	<< endl;
#endif
							if (ixevent < nEvents) {
								// read next event
								eEvent = pLandscape->getExtEvent(ixevent++);								
							}
							else {
								eEvent.year = eEvent.season = 99999999;   
							}
						}
					}
#else
					pComm->survival(gen,0,1,2); // survival of all non-juvenile stages
#endif // PARTMIGRN 
#else
#if SPATIALMORT
					pComm->survival(0,period,2,1); // survival of all non-juvenile stages
#else
#if PEDIGREE
					pComm->survival(pPed,0,2,1); // survival of all non-juvenile stages
#else
					pComm->survival(0,2,1); // survival of all non-juvenile stages
#endif // PEDIGREE
#endif // SPATIALMORT 
#endif // SEASONAL 
				}
			}

#if RS_CONTAIN
			if (d.timing == 1) pComm->updateDamage(pSpecies,pCull);
#endif // RS_CONTAIN 
			
			// Output and pop. visualisation AFTER reproduction
#if RS_ABC
			if (dem.stageStruct && (sim.outRange || sim.outPop || abcYear))
				RangePopOutput(pComm,rep,yr,gen,pABCmaster,abcYear);
#else
			if (dem.stageStruct && (sim.outRange || sim.outPop))
				RangePopOutput(pComm,rep,yr,gen);
#endif
#if VIRTUALECOLOGIST
			if (dem.stageStruct && sim.virtualEcologist) {
				// likewise
				if ((yr >= virt.outStart && yr%virt.outInt == 0 && gen == 0)) {
//					if (virt.landscapeGenetics && virt.patchMethod == 2) { // dynamic sampling
//						pVirt->samplePatches(pComm);
//					}
					pVirt->outLandGen(pLandscape,rep,yr,gen,ppLand.patchModel,false);
				}
			}
#endif

#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed reproduction" << endl;
#endif

#if RS_CONTAIN
			if (pCull->isCullApplied()) {
				if (pCull->cullTiming() == 0) {
#if RSDEBUG
#if RS_CONTAIN 
DEBUGLOG << endl << "RunModel(): PERFORMING MANAGEMENT CULL before dispersal" << endl;
#else
DEBUGLOG << endl << "RunModel(): PERFORMING MANAGEMENT CULL after reproduction" << endl;
#endif // RS_CONTAIN  
#endif
					if (dem.stageStruct) ManagementCull(pLandscape,yr,sstruct.nStages);					
					else ManagementCull(pLandscape,yr,1);								
				}
			}
#endif // RS_CONTAIN 

			// Dispersal

#if RS_DISEASE
			pComm->emigration(pSpecies,gen);
#else
#if SEASONAL
			pComm->emigration(gen);
#else
			pComm->emigration();
#endif // SEASONAL 
#endif // RS_DISEASE  
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed emigration" << endl;
#endif
#if RS_ABC
			pComm->dispersal(landIx,pABCmaster->obsConnectivity());
#else
#if PEDIGREE
			pComm->dispersal(pPed,rep,yr,gen,landIx);
#else
#if SEASONAL
			if ((gen+1) < dem.nSeasons) pComm->dispersal(landIx,(gen+1));
			else pComm->dispersal(landIx,0);
#else
#if RS_RCPP
			pComm->dispersal(landIx,yr);
#else
			pComm->dispersal(landIx);
#endif // RS_RCPP
#endif // SEASONAL 
#endif // PEDIGREE
#endif // RS_ABC
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed dispersal" << endl;
#endif
#if VCL
			if (stopRun) break;
#endif

#if RS_CONTAIN 
			if (pCull->isCullApplied()) {
				if (pCull->cullTiming() == 1) {
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): PERFORMING MANAGEMENT CULL after dispersal" << endl;
#endif				
					if (dem.stageStruct) ManagementCull(pLandscape,yr,sstruct.nStages);					
					else ManagementCull(pLandscape,yr,1);								
				}
			}
#endif // RS_CONTAIN 

#if BUTTERFLYDISP
			if (!dem.stageStruct && dem.dispersal == 0) { // dispersal prior to parturition
				// individual-level output must occur before parturition and survival
				// in order that successful dispersal of adults is recorded, as all
				// adults die at fledging
				// output Individuals
				if (sim.outInds && yr >= sim.outStartInd && yr%sim.outIntInd == 0)
					pComm->outInds(rep,yr,gen,-1);
				// output Genetics
				if (sim.outGenetics && yr >= sim.outStartGenetic && yr%sim.outIntGenetic == 0)
					pComm->outGenetics(rep,yr,gen,-1);
				// parturition ...
				pComm->reproduction(pSpecies,yr,2);
				// ... and fledging (which must be deferred until after parturition has
				// occurred in ALL patches, since in a genetic model, a juvenile's father
				// may be in a different patch, and must not be deleted until inheritance
        // is completed
				pComm->fledge();
			}
#endif

#if RS_ABC
			if (abcYear && gen == 0)
				pComm->outABCpreds(rep,yr,pABCmaster,ppLand.resol);
#endif

			// survival part 0
#if SPATIALMORT
			if (dem.stageStruct) {
				if (sstruct.survival == 0) { // at reproduction
					pComm->survival(0,period,0,1); // survival of juveniles only
				}
				if (sstruct.survival == 1) { // between reproduction events
					pComm->survival(0,period,1,1); // survival of all stages
				}
				if (sstruct.survival == 2) { // annually
					pComm->survival(0,period,1,0); // development only of all stages
				}
			}
			else { // non-structured population
				pComm->survival(0,period,1,1);
			}
#else
			if (dem.stageStruct) {
				if (sstruct.survival == 0) { // at reproduction
#if SEASONAL
#if PARTMIGRN
					// no action, as survival was applied to ALL stages prior to dispersal
#else
					pComm->survival(gen,0,0,1); // survival of juveniles only
#endif // PARTMIGRN 
#else
#if PEDIGREE
					pComm->survival(pPed,0,0,1); // survival of juveniles only
#else
					pComm->survival(0,0,1); // survival of juveniles only
#endif // PEDIGREE
#endif // SEASONAL 
				}
				if (sstruct.survival == 1) { // between reproduction events
#if SEASONAL
					pComm->survival(gen,0,1,1); // survival of all stages
#else
#if PEDIGREE
					pComm->survival(pPed,0,1,1); // survival of all stages
#else
					pComm->survival(0,1,1); // survival of all stages
#endif // PEDIGREE
#endif // SEASONAL 
				}
				if (sstruct.survival == 2) { // annually
#if SEASONAL
					pComm->survival(gen,0,1,0); // development only of all stages
#else
#if PEDIGREE
					pComm->survival(pPed,0,1,0); // development only of all stages
#else
					pComm->survival(0,1,0); // development only of all stages
#endif // PEDIGREE
#endif // SEASONAL 
//					pComm->survival(0,1,0); // development only of all stages
				}
			}
			else { // non-structured population
#if SEASONAL
				pComm->survival(gen,0,1,1);
#else
#if PEDIGREE
				pComm->survival(pPed,0,1,1);
#else
				pComm->survival(0,1,1);
#endif // PEDIGREE
#endif // SEASONAL 
			}
#endif // SPATIALMORT
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed survival part 0" << endl;
#endif

#if PEDIGREE
			pPed->updateRelMat();
#endif


#if BUTTERFLYDISP
			if (dem.stageStruct || dem.dispersal == 1) {
#endif
				// output Individuals
				if (sim.outInds && yr >= sim.outStartInd && yr%sim.outIntInd == 0)
					pComm->outInds(rep,yr,gen,-1);
				// output Genetics
				if (sim.outGenetics && yr >= sim.outStartGenetic && yr%sim.outIntGenetic == 0)
					pComm->outGenetics(rep,yr,gen,-1);
#if BUTTERFLYDISP
			}
#endif
#if VIRTUALECOLOGIST
			if (sim.virtualEcologist) {
				if ((yr >= virt.outStart && yr%virt.outInt == 0 && gen == 0)) {
//					if (virt.landscapeGenetics && virt.patchMethod == 2) { // dynamic sampling
//						pVirt->samplePatches(pComm);
//					}
					pVirt->outLandGen(pLandscape,rep,yr,gen,ppLand.patchModel,true);
				}
			}
#endif
#if RS_ABC
			if (yr == yearABC && gen == 0) {
				// get next ABC year
				if (ixABC < nABCyears) yearABC = pABCmaster->getYear(ixABC++);
			}
#endif

			// survival part 1
#if SPATIALMORT
			if (dem.stageStruct) {
//				if (sstruct.survival != 2) { // at reproduction or between reproduction events
					pComm->survival(1,period,0,1);
//				}
			}
			else { // non-structured population
				pComm->survival(1,period,0,1);
			}
#else
			if (dem.stageStruct) {
//				if (sstruct.survival != 2) { // at reproduction or between reproduction events
#if SEASONAL
					pComm->survival(gen,1,0,1);
#else
#if PEDIGREE
					pComm->survival(pPed,1,0,1);
#else
					pComm->survival(1,0,1);
#endif // PEDIGREE
#endif // SEASONAL 
//				}
			}
			else { // non-structured population
#if SEASONAL
				pComm->survival(gen,1,0,1);
#else
#if PEDIGREE
				pComm->survival(pPed,1,0,1);
#else
				pComm->survival(1,0,1);
#endif // PEDIGREE
#endif // SEASONAL 
			}
#endif // SPATIALMORT
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " gen=" << gen << " completed survival part 1" << endl;
#endif

#if RS_CONTAIN 
			// Write data to Cull file to match data written to Pop file
			if (sim.outPop && yr >= sim.outStartPop && yr%sim.outIntPop == 0)
				pComm->outCull(rep,yr,gen);
#endif // RS_CONTAIN  

#if SEASONAL
			// Connectivity Matrix
			if (sim.outConnect && ppLand.patchModel
			&& yr >= sim.outStartConn && yr%sim.outIntConn == 0)
				pLandscape->outConnect(rep,yr,gen);
#endif // SEASONAL 

#if RS_CONTAIN
			if (ppLand.dmgLoaded) {
//				pLandscape->resetPrevDamage();
				if (sim.outRange && yr%sim.outIntRange == 0) 
					pLandscape->outSummDmg(rep,yr,trfr.moveModel && trfr.moveType == 1,v.viewDamage);
				if (sim.outDamage && yr >= sim.outStartDamage && yr%sim.outIntDamage == 0) 
					pLandscape->outDamage(rep,yr,trfr.moveModel && trfr.moveType == 1);
			}	
			pComm->resetCull();			
#endif // RS_CONTAIN 
#if VCL
			if (v.viewCosts && v.viewPaths) {
				RefreshVisualCost();
			}
			Application->ProcessMessages();
			if (stopRun) break;
#endif
		} // end of the generation loop
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed generation loop" << endl;
#endif

#if !RS_ABC
		totalInds =  pComm->totalInds();
		if (totalInds <= 0) { yr++; break; }
#endif
#if VCL
		if (stopRun) break;
#endif

#if !SEASONAL
		// Connectivity Matrix
		if (sim.outConnect && ppLand.patchModel
		&& yr >= sim.outStartConn && yr%sim.outIntConn == 0)
			pLandscape->outConnect(rep,yr);
#endif // !SEASONAL 
#if RS_ABC
		if (ppLand.patchModel && abcYear) {
			// process any connectivity predictions for ABC
			obsdata obs;
			int nobs = (int)pABCmaster->NObs();
			for (int i = 0; i < nobs; i++) {
				obs = pABCmaster->getObsData(i);
#if RSDEBUG
//DEBUGLOG << "RunModel(): i=" << i << " yr=" << yr
//	<< " obs.year=" << obs.year << " obs.type=" << obs.type << " obs.name=" << obs.name
//	<< endl;
#endif
				if (obs.year == yr && obs.type == 4) {
#if RSDEBUG
DEBUGLOG << "RunModel(): i=" << i << " PROCESS Type 4 Connectivity"
	<< " obs.id=" << obs.id << " obs.name=" << obs.name
	<< " obs.x=" << obs.x << " obs.y=" << obs.y
	<< " obs.value=" << obs.value
	<< " obs.weight=" << obs.weight
	<< endl;
#endif
					if (obs.name == "NInds") {
						int npred = pLandscape->outABCconnect(obs.x,obs.y);
						pABCmaster->AddNewPred(sim.simulation,obs.id,rep,obs.value,
							npred,obs.weight);
					}
				}
			}
		} // if (ppLand.patchModel && abcYear)
#endif

		if (dem.stageStruct && sstruct.survival == 2) {  // annual survival - all stages
#if SPATIALMORT
			pComm->survival(0,period,1,2);
			pComm->survival(1,period,0,1);
#else
#if PEDIGREE
			pComm->survival(pPed,0,1,2);
			pComm->survival(pPed,1,0,1);
#else
#if !SEASONAL
			pComm->survival(0,1,2);
			pComm->survival(1,0,1);
#endif // !SEASONAL
#endif // PEDIGREE
#endif // SPATIALMORT
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed annual survival" << endl;
#endif
		}

		if (dem.stageStruct) {
			pComm->ageIncrement(); // increment age of all individuals
			if (sim.outInds && yr >= sim.outStartInd && yr%sim.outIntInd == 0)
				pComm->outInds(rep,yr,-1,-1); // list any individuals dying having reached maximum age
#if SPATIALMORT
			pComm->survival(1,period,0,1);		// delete any such individuals
#else
#if PEDIGREE
			pComm->survival(pPed,1,0,1);						// delete any such individuals
#else
#if SEASONAL
			pComm->survival(888,1,0,1);						// delete any such individuals
#else
			pComm->survival(1,0,1);						// delete any such individuals
#endif // SEASONAL
#endif // PEDIGREE
#endif // SPATIALMORT
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed Age_increment and final survival" << endl;
#endif
#if !RS_ABC
		totalInds =  pComm->totalInds();
		if (totalInds <= 0) { yr++; break; }
#endif
		}

	} // end of the years loop

#if RS_ABC
	if (yr == yearABC) abcYear = true;
	else abcYear = false;
#endif

	// Final output and popn. visualisation
#if VCL
	if (!stopRun) {
		if (v.viewPop || (sim.saveMaps && yr%sim.mapInt == 0)) {
			if (updateland) {
				pLandscape->drawLandscape(rep,landIx,ppLand.landNum);
			}
			pComm->draw(rep,yr,0,ppLand.landNum);
		}
		if (v.viewGraph) {
			DrawPopnGraph(pComm,yr);
		}
	}
#endif
#if BATCH
	if (sim.saveMaps && yr%sim.mapInt == 0) {
		if (updateland) {
			pLandscape->drawLandscape(rep,landIx,ppLand.landNum);
		}
		pComm->draw(rep,yr,0,ppLand.landNum);
	}
#endif
#if VCL
	if (!stopRun) {
#endif
		// produce final summary output
		if (v.viewPop || v.viewTraits || sim.outOccup
		|| 	sim.outTraitsCells || sim.outTraitsRows || sim.saveMaps)
			PreReproductionOutput(pLandscape,pComm,rep,yr,0);
#if RS_ABC
		if (sim.outRange || sim.outPop || abcYear)
			RangePopOutput(pComm,rep,yr,0,pABCmaster,abcYear);
#else
		if (sim.outRange || sim.outPop)
			RangePopOutput(pComm,rep,yr,0);
#endif
#if RSDEBUG
DEBUGLOG << "RunModel(): yr=" << yr << " completed final summary output" << endl;
#endif
#if VCL
	}
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

#if RS_ABC
	if (ppLand.patchModel && (sim.outConnect || pABCmaster->obsConnectivity()))
#else
	if (sim.outConnect && ppLand.patchModel)
#endif
		pLandscape->resetConnectMatrix(); // set connectivity matrix to zeroes

#if VCL
	// redraw maps for the next replicate
	if (!stopRun && (v.viewLand || sim.saveMaps) && rep != sim.reps-1)
		pLandscape->drawLandscape(0,0,0);
#endif

	if (sim.outInds) // close Individuals output file
		pComm->outInds(rep,0,0,-999);
	if (sim.outGenetics) // close Genetics output file
		pComm->outGenetics(rep,0,0,-999);

	if (sim.saveVisits) {
#if VCL		
		if (!sim.batchMode) pLandscape->saveVisits(rep,ppLand.landNum);
#endif
		pLandscape->outVisits(rep,ppLand.landNum);
		pLandscape->resetVisits();
	}

#if RS_RCPP
	if (sim.outPaths)
		pLandscape->outPathsHeaders(rep,-999);
#endif
#if VCL
	if (stopRun) break;
#endif
#if RSDEBUG
DEBUGLOG << endl << "RunModel(): finished rep=" << rep << endl;
#endif

#if GROUPDISP
#if PEDIGREE
	if (pPed != 0) delete pPed;
//	pSpecies->deleteRelMat();
#endif
#endif

} // end of the replicates loop

#if RS_ABC

// Average predictions across replicates
#if RSDEBUG
//DEBUGLOG << "RunModel(): PREDICTIONS" << endl;
#endif

pABCmaster->sortPreds();

preddata pred,predcomp;
predcomp.sample = predcomp.id = predcomp.rep = 0;
predcomp.obsvalue = predcomp.predvalue = 0.0;
int nreps = 0;
double meanpred = 0.0;
int npreds = (int)pABCmaster->NPreds();
//cout << "RunModel(): npreds=" << npreds << endl;
for (int i = 0; i < npreds; i++) {
	pred = pABCmaster->getPredData(i);
#if RSDEBUG
//DEBUGLOG << "RunModel(): sample=" << pred.sample
//	<< " id=" << pred.id << " rep=" << pred.rep
//	<< " obsvalue=" << pred.obsvalue << " predvalue=" << pred.predvalue
//	<< " weight=" << pred.weight
//	<< endl;
#endif
	if (pred.id != predcomp.id) { // save comparison for previous id
		if (nreps > 0) meanpred /= (double)nreps; else meanpred = 0.0;
		pABCmaster->AddNewPredComp(predcomp.sample,predcomp.id,nreps,
			predcomp.obsvalue,meanpred,predcomp.weight);
		meanpred = 0.0; nreps = 0;
	}
	predcomp = pred;
	meanpred += pred.predvalue; nreps++;
}
// save comparison for last id

//cout << endl << endl << "***** npreds=" << npreds
//	<< " nreps=" << nreps
//	<< " meanpred=" << meanpred
//	<< endl << endl;

if (nreps > 0) meanpred /= (double)nreps; else meanpred = 0.0;
pABCmaster->AddNewPredComp(predcomp.sample,predcomp.id,nreps,
	predcomp.obsvalue,meanpred,predcomp.weight);

pABCmaster->deletePreds();

#endif // RS_ABC

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

#if RS_CONTAIN
#if VCL
if (v.viewDamage) {
	pLandscape->outTotDamage(true);
	pLandscape->deleteTotDamage(sim.reps);
}
#endif
#endif // RS_CONTAIN 

if (sim.outRange) {
	pComm->outRangeHeaders(pSpecies,-999); // close Range file
#if RS_CONTAIN	
	if (ppLand.dmgLoaded) pLandscape->outSummDmgHeaders(-999);
#endif // RS_CONTAIN 
}
#if RS_RCPP
	if (sim.outPop && sim.CreatePopFile) {
#else
	if (sim.outPop) {
#endif
	pComm->outPopHeaders(pSpecies,-999); // close Population file
#if RS_CONTAIN
	if (pCull->isCullApplied()) {
		pComm->outCullHeaders(pSpecies,-999); // close Population file		
	}
#endif // RS_CONTAIN 
}
if (sim.outTraitsCells)
	pComm->outTraitsHeaders(pSpecies,-999); // close Traits file
if (sim.outTraitsRows)
	pComm->outTraitsRowsHeaders(pSpecies,-999); // close Traits rows file
// close Individuals & Genetics output files if open
// they can still be open if the simulation was stopped by the user
if (sim.outInds) pComm->outInds(0,0,0,-999);
if (sim.outGenetics) pComm->outGenetics(0,0,0,-999);
#if RS_CONTAIN	
if (ppLand.dmgLoaded && sim.outDamage) pLandscape->outDamageHeaders(-999);
#endif // RS_CONTAIN 
#if VIRTUALECOLOGIST
if (sim.virtualEcologist) {
	pVirt->outLandGenHeaders(-999,false);
	if (virt.outGenomes) pVirt->outGenSamplesHeaders(-999,false);
}
#endif
#if PEDIGREE
pComm->outGroupHeaders(-999); // close groups file
#endif

MemoLine("Deleting community...");
delete pComm; pComm = 0;
#if VIRTUALECOLOGIST
if (sim.virtualEcologist && pVirt != 0) {
	delete pVirt; // delete virtual ecologist
}
#endif
MemoLine("...finished");

#if VCL
Application->ProcessMessages();
#endif

// Write performance data
//t1 = time(0);
//RSlog << "Simulation," << sim.simulation << "," << sim.reps << "," << sim.years
//	<< "," << t1-t0 << endl;

#if RS_RCPP && !R_CMD
	return list_outPop;
#else
	return 0;
#endif

}

#if RS_EMBARCADERO || LINUX_CLUSTER || RS_RCPP 
// Check whether a specified directory path exists
bool is_directory(const char *pathname) {
struct stat info;
if (stat(pathname, &info) != 0) return false; // path does not exist
if (S_ISDIR(info.st_mode)) return true;
return false;
}
#endif

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
#if SEASONAL
	pComm->updateOccupancy(yr/sim.outIntOcc,rep,0);
#else
	pComm->updateOccupancy(yr/sim.outIntOcc,rep);
#endif // SEASONAL 

#if RSDEBUG
//DEBUGLOG << "PreReproductionOutput(): 88888 " << endl;
#endif

// Remaining graphical output actions are performed for GUI only
#if VCL
if (v.viewTraits) Outputs_Visuals_B(rep,yr,gen,ppLand.landNum);
#endif

#if RSDEBUG
//DEBUGLOG << "PreReproductionOutput(): finished " << endl;
#endif
}

//For outputs and population visualisations pre-reproduction
#if RS_ABC
void RangePopOutput(Community *pComm,int rep,int yr,int gen,
	ABCmaster *pABCmaster,bool abcYear)
#else
void RangePopOutput(Community *pComm,int rep,int yr,int gen)
#endif
{
simParams sim = paramsSim->getSim();

#if RS_ABC
if ((sim.outRange && (yr%sim.outIntRange == 0 || pComm->totalInds() <= 0)) || abcYear)
	pComm->outRange(pSpecies,rep,yr,gen,pABCmaster,abcYear);
#else
if (sim.outRange && (yr%sim.outIntRange == 0 || pComm->totalInds() <= 0))
	pComm->outRange(pSpecies,rep,yr,gen);
#endif

#if RS_ABC
bool popOutputYear = false;
if (sim.outPop && yr >= sim.outStartPop && yr%sim.outIntPop == 0)
	popOutputYear = true;
if (popOutputYear || abcYear)
	pComm->outPop(pSpecies,rep,yr,gen,pABCmaster,abcYear,popOutputYear);
#else
#if RS_RCPP
if (sim.outPop && sim.CreatePopFile && yr >= sim.outStartPop && yr%sim.outIntPop == 0)
#else
if (sim.outPop && yr >= sim.outStartPop && yr%sim.outIntPop == 0)
#endif
	pComm->outPop(rep,yr,gen);
#endif

}

#if RS_CONTAIN

void ManagementCull(Landscape *pLandscape,int year,int nstages)
{
culldata c = pCull->getCullData();
int nTargetPatches = pComm->findCullTargets(pCull,year,nstages);
#if RSDEBUG
DEBUGLOG << "ManagementCull(): maxNpatches=" << c.maxNpatches 
	<< " nTargetPatches=" << nTargetPatches << endl;
#endif
if (nTargetPatches <= 0) return;
if (nTargetPatches <= c.maxNpatches) {
	// cull all target patches
	pComm->cullAllTargets(pCull);
}
else {
	// cull the maximum number of patches
	// NB IF THE CULL IS TO BE SPATIALLY CORRELATED,
	//    IT WILL HAVE TO BECOME A FUNCTION OF THE Landscape
	switch (c.method) {	
	case 0: // random
	case 1: // random recently colonised
	case 2: // random closest to damage
	case 3: // random closest to damage X popn size
		pComm->cullRandomTargets(pCull,year);
//		pComm->cullClosestTargets(pCull,year);
		break;
	case 4: // closest to damage
	case 5: // closest to damage X popn size
		pComm->cullTargets(pCull);
		break;
	case 6: // random reactive to previous damage
		pComm->cullRandomTargets(pCull,year);
		break;
	case 7: // reactive to previous damage
		pComm->cullTargets(pCull);
		break;
	}
}
pComm->resetCullTargets();
}

#endif // RS_CONTAIN 

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
#if SEASONAL
//outPar << " PARTIAL MIGRATION MODEL";
outPar << " SEASONAL MODEL ";
#endif // SEASONAL 
#if GOBYMODEL
outPar << " GOBY MODEL ";
#endif // GOBYMODEL
#if SOCIALMODEL
outPar << " SOCIAL PHENOTYPE MODEL ";
#endif // SOCIALMODEL
#if GROUPDISP
outPar << " GROUP DISPERSAL MODEL ";
#endif // GROUPDISP 
#if RS_ABC
outPar << " APPROXIMATE BAYESIAN COMPUTATION ";
#endif // RS_ABC
#if RS_CONTAIN
outPar << " ADAPTIVE MANAGEMENT ";
#endif // RS_CONTAIN 

#if !RS_RCPP
#if RSWIN64
outPar << " - 64 bit implementation";
#else
outPar << " - 32 bit implementation";
#endif
#endif
outPar << endl;

outPar << "================ ";
#if SEASONAL
//outPar << " =======================";
outPar << " ============== ";
#endif // SEASONAL 
#if GOBYMODEL
outPar << " ========== ";
#endif // GOBYMODEL
#if SOCIALMODEL
outPar << " ====================== ";
#endif // SOCIALMODEL
#if GROUPDISP
outPar << " ===================== ";
#endif // GROUPDISP 
#if RS_ABC
outPar << " ================================ ";
#endif // RS_ABC
#if RS_CONTAIN
outPar << " =================== ";
#endif // RS_CONTAIN 

outPar << "   =====================";
outPar << endl << endl;

outPar << "BATCH MODE \t";
if (sim.batchMode) outPar << "yes" << endl; else outPar << "no" << endl;
#if RS_RCPP
outPar << "SEED \t" << RS_random_seed << endl;
#endif
outPar << "REPLICATES \t" << sim.reps << endl;
outPar << "YEARS \t" << sim.years << endl;
#if SEASONAL
outPar << "NO. SEASONS / YEAR\t" << dem.nSeasons << endl;
#else
outPar << "REPRODUCTIVE SEASONS / YEAR\t" << dem.repSeasons << endl;
#endif
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
#if RS_RCPP
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
#else
	if (sim.batchMode) outPar << " (see batch file) " << landFile << endl;
	else {
		outPar << habmapname << endl;
		if (ppLand.rasterType == 1) { // habitat % cover - list additional layers
			for (int i = 0; i < ppLand.nHab-1; i++) {
				outPar  << "           "<< hfnames[i] << endl;
			}
		}
		if (ppLand.patchModel) {
			outPar << "PATCH FILE: " << patchmapname << endl;
		}
	}
#endif
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
#if RS_CONTAIN
if (ppLand.dmgLoaded) {
	outPar << endl << "ECONOMIC / ENVIRONMENTAL DAMAGE: " << endl;
	outPar << "DAMAGE MAP:             ";     
	if (sim.batchMode) outPar << "(see batch file) " << landFile << endl;
	else outPar << dmgmapname << endl;
	outPar << "DISTANCE DECAY (alpha): " << pLandscape->getAlpha() << endl;
}
#endif // RS_CONTAIN 
#if SEASONAL
#if PARTMIGRN
if (pLandscape->numExtEvents() > 0) {
	outPar << endl << "SPECIFIC EXTREME EVENTS: " << endl;
	outPar << "No. of events: " << pLandscape->numExtEvents() << endl;	
}
#endif // PARTMIGRN 
#endif // SEASONAL 

outPar << endl << "SPECIES DISTRIBUTION LOADED: \t";
//if (sim.initDistLoaded)
if (ppLand.spDist)
{
	outPar << "yes" << endl;
	outPar << "RESOLUTION (m)\t" << ppLand.spResol << endl;
	outPar << "FILE NAME: ";
#if !RS_RCPP
	if (sim.batchMode) outPar << " (see batch file) " << landFile << endl;
	else {
		outPar << distnmapname << endl;
	}
#else
	outPar << name_sp_dist << endl;
#endif
}
else outPar << "no" << endl;

#if SPATIALMORT
outPar << endl << "SPATIAL MORTALITY MAPS UPLOADED: \t";
if (sim.mortMapLoaded) {
	outPar << "yes" << endl;
	if (sim.batchMode) {
		outPar << "FILE NAME: (see batch file) " << landFile << endl;
	}
	else {
		outPar << "FILE NAME MAP 1: " << mortmapname[0] << endl;
		outPar << "FILE NAME MAP 2: " << mortmapname[1] << endl;
	}
	outPar << "MAP 2 APPLIED FROM YEAR: " << paramsSim->getMortChgYear() << endl;
}
else outPar << "no" << endl;
#endif

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
#if BUTTERFLYDISP
	if (!env.local && env.fromFile) {
		outPar << "FROM FILE: " << envstochfilename << endl;
	}
	else {
		outPar << "TEMPORAL AUTOCORRELATION (ac)\t" << env.ac << endl;
		outPar << "AMPLITUDE (std)\t" << env.std << endl;
	}
#else
	outPar << "TEMPORAL AUTOCORRELATION (ac)\t" << env.ac << endl;
	outPar << "AMPLITUDE (std)\t" << env.std << endl;
#endif
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
#if GROUPDISP
		if (dem.paternity == 1) outPar << " with random paternity";
		if (dem.paternity == 2) outPar << " with pollen kernel";
//		if (dem.selfing) outPar << " with self-fertilisation";
		outPar << endl;
		if (dem.paternity == 2) {
			outPar << "POLLEN KERNEL: Proportion local: " << dem.propLocal
				<< " Proportion neighbourhood: " << dem.propNghbr << endl;
		}
#else
		outPar << endl;
#endif
		outPar << "PROP. of MALES\t" << dem.propMales << endl;
		break;
	case 2:
		outPar << "Sexual model (explicit mating system)" << endl;
		outPar << "PROP. of MALES\t" << dem.propMales << endl;
		outPar << "MAX. HAREM SIZE (h)\t" << dem.harem << endl;
	break;
#if GROUPDISP
	case 3:
		outPar << "Hermaphrodite";
		if (dem.paternity == 1) outPar << " with random paternity";
		if (dem.paternity == 2) outPar << " with pollen kernel";
		if (dem.selfing) outPar << " with self-fertilisation";
		outPar << endl;
		if (dem.paternity == 2) {
			outPar << "POLLEN KERNEL: Proportion local: " << dem.propLocal
				<< " Proportion neighbourhood: " << dem.propNghbr << endl;
		}
		break;
#endif
}
#if RS_CONTAIN
int nhabitats = ppLand.nHab; 
if (nhabitats > NHABITATS) nhabitats = NHABITATS;
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
#if SEASONAL
		if (dem.habDepDem) {
			for (int h = 0; h < nhabitats; h++) {
				for (int j = 0; j < dem.nSeasons; j++) {
					outPar << endl << "HABITAT: " << pLandscape->getHabCode(h);
					outPar << " SEASON: " << j;
//					if (pSpecies->getBreeding(j)) outPar << " Breeding";
//					else outPar << " Non-breeding";
					outPar << endl;
					outPar << "FECUNDITIES:" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						outPar << "stage\t" << i << ":\t" << pSpecies->getFec(h,j,i,0) << endl;
					}
					outPar << "DEVELOPMENT PROB.:" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						outPar << "stage\t" << i << ":\t" << pSpecies->getDev(h,j,i,0) << endl;
					}
					outPar << "SURVIVAL PROB.:" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(h,j,i,0) << endl;
					}
				}
			}
			outPar << endl;
		}
		else {
			for (int j = 0; j < dem.nSeasons; j++) {
				outPar << endl << "SEASON: " << j;
				if (pSpecies->getBreeding(j)) outPar << " Breeding";
				else outPar << " Non-breeding";
				outPar << endl;
				outPar << "FECUNDITIES:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getFec(0,j,i,0) << endl;
				}
				outPar << "DEVELOPMENT PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getDev(0,j,i,0) << endl;
				}
				outPar << "SURVIVAL PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(0,j,i,0) << endl;
				}
			}
		}
#else
		if (dem.habDepDem) {
			for (int h = 0; h < nhabitats; h++) {
				outPar << endl << "HABITAT: " << pLandscape->getHabCode(h);
				outPar << endl;
				outPar << "FECUNDITIES:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getFec(h,i,0) << endl;
				}
				outPar << "DEVELOPMENT PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getDev(h,i,0) << endl;
				}
				outPar << "SURVIVAL PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(h,i,0) << endl;
				}
			}
			outPar << endl;
		}
		else {
			outPar << "FECUNDITIES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getFec(0,i,0) << endl;
			}
			outPar << "DEVELOPMENT PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getDev(0,i,0) << endl;
			}
			outPar << "SURVIVAL PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(0,i,0) << endl;
			}
		}
#endif // SEASONAL 
#else
#if SEASONAL
		extrmevent e;
		for (int j = 0; j < dem.nSeasons; j++) {
			outPar << endl << "SEASON: " << j;
			if (pSpecies->getBreeding(j)) outPar << " Breeding";
			else outPar << " Non-breeding";
			outPar << endl;
			outPar << "FECUNDITIES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getFec(j,i,0) << endl;
			}
			outPar << "DEVELOPMENT PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getDev(j,i,0) << endl;
			}
			outPar << "SURVIVAL PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(j,i,0) << endl;
			}
			e = pSpecies->getExtreme(j);
			outPar << "RANDOM EXTREME EVENTS: Probability " << e.prob << " Mortality " << e.mort << endl;
		}
		outPar << endl;
#else
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
#endif // SEASONAL 
#endif // RS_CONTAIN 
	}
	// sex-specific demographic parameters
	else {
		outPar << "MIN. AGES:" << endl;
		for (int i = 0; i < sstruct.nStages; i++) {
			outPar << "males " << i << ":\t" << pSpecies->getMinAge(i,1) << " years;\t";
			outPar << "females " << i << ":\t" << pSpecies->getMinAge(i,0) << " years" << endl;
		}
#if RS_CONTAIN
#if SEASONAL
		for (int h = 0; h < nhabitats; h++) {
			for (int j = 0; j < dem.nSeasons; j++) {
				outPar << "HABITAT:" << h;
				outPar << " SEASON:" << j << endl;
				outPar << "FECUNDITIES:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "males   " << i << ":\t" << pSpecies->getFec(h,j,i,1) << endl;
					outPar << "females " << i << ":\t" << pSpecies->getFec(h,j,i,0) << endl;
				}
				outPar << "DEVELOPMENT PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "males   " << i << ":\t" << pSpecies->getDev(h,j,i,1) << endl;
					outPar << "females " << i << ":\t" << pSpecies->getDev(h,j,i,0) << endl;
				}
				outPar << "SURVIVAL PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "males   " << i << ":\t" << pSpecies->getSurv(h,j,i,1) << endl;
					outPar << "females " << i << ":\t" << pSpecies->getSurv(h,j,i,0) << endl;
				}
			}
		}
#else
		for (int h = 0; h < nhabitats; h++) {
			outPar << "HABITAT:" << h << endl;
			outPar << "FECUNDITIES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getFec(h,i,1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getFec(h,i,0) << endl;
			}
			outPar << "DEVELOPMENT PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getDev(h,i,1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getDev(h,i,0) << endl;
			}
			outPar << "SURVIVAL PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getSurv(h,i,1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getSurv(h,i,0) << endl;
			}
		}
#endif // SEASONAL 
#else
#if SEASONAL
		for (int j = 0; j < dem.nSeasons; j++) {
			outPar << "SEASON:" << j << endl;
			outPar << "FECUNDITIES:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getFec(j,i,1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getFec(j,i,0) << endl;
			}
			outPar << "DEVELOPMENT PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getDev(j,i,1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getDev(j,i,0) << endl;
			}
			outPar << "SURVIVAL PROB.:" << endl;
			for (int i = 0; i < sstruct.nStages; i++) {
				outPar << "males   " << i << ":\t" << pSpecies->getSurv(j,i,1) << endl;
				outPar << "females " << i << ":\t" << pSpecies->getSurv(j,i,0) << endl;
			}
		}
#else
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
#endif // SEASONAL 
#endif // RS_CONTAIN 
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

#if PARTMIGRN
outPar << endl << "PROBABILITIES OF DISPERSAL/MIGRATION STRATEGIES:" << endl;
outPar << "1 = Philopatric resident                             " << pSpecies->getPropDispMigrn(1) << endl;
outPar << "2 = Philopatric migrant, non-breeding site fixed     " << pSpecies->getPropDispMigrn(2) << endl;
outPar << "3 = Philopatric migrant, non-breeding site not fixed " << pSpecies->getPropDispMigrn(3) << endl;
outPar << "4 = Dispersed resident                               " << pSpecies->getPropDispMigrn(4) << endl;
outPar << "5 = Dispersed migrant, non-breeding site fixed       " << pSpecies->getPropDispMigrn(5) << endl;
outPar << "6 = Dispersed migrant, non-breeding site not fixed   " << pSpecies->getPropDispMigrn(6) << endl;
#endif // PARTMIGRN 

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
#if SEASONAL
for (int j = 0; j < dem.nSeasons; j++) {
	for (int i = 0; i < nhab; i++) {
		k = pSpecies->getHabK(i,j) * (10000.0/(float)(ppLand.resol*ppLand.resol));
		if (!ppLand.generated && ppLand.rasterType == 0) { // imported & habitat codes
			outPar << "Season " << j << " Habitat " << pLandscape->getHabCode(i) << ": \t";
		}
		else {
			outPar << "Season " << j << " Habitat " << i << ": ";
		}
		if (dem.stageStruct) outPar << "1/b ";
		else outPar << "K ";
		outPar << k << endl;		
	}
}
#else
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
#endif // SEASONAL 

#if GOBYMODEL
outPar << endl << "SOCIALITY PHENOTYPE:" << endl;
socialParams s = pSpecies->getSocialParams();
outPar << "Phenotype initialisation:  mean " << s.socMean << "  s.d. " << s.socSD
	<< "  scaling factor " << s.socScale << endl;
outPar << "F asocial:      " << sstruct.asocF << endl;
outPar << "D asocial:      " << emig.asocD << endl;
outPar << "AlphaS asocial: " << sett.alphaSasoc << endl;
outPar << "BetaS asocial:  " << sett.betaSasoc << endl;
#endif

emigTraits ep0,ep1;
emigParams eparams0,eparams1;
string sexdept 		= "SEX-DEPENDENT:   ";
string stgdept  = "STAGE-DEPENDENT: ";
string indvar = "INDIVIDUAL VARIABILITY: ";
string emigstage = "EMIGRATION STAGE: ";

#if BUTTERFLYDISP
outPar << endl << "SCHEDULING OF DISPERSAL: ";
switch (dem.dispersal) {
case 0:
	outPar << "During reproduction" << endl;
	break;
case 1:
	outPar << "After reproduction" << endl;
	break;
}
#endif

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
#if GROUPDISP
if (emig.groupdisp) {
	outPar << "GROUP DISPERSAL: \tmean group size: " << emig.groupmean;
	if (emig.grouptype == 0) outPar << " at population level";
	else outPar << " at sibling level";
	outPar << endl;
}
#endif

// Transfer

outPar << endl << "DISPERSAL - TRANSFER: \t";

if (trfr.moveModel) {
	bool straigtenPath;
	if (trfr.moveType == 1) { // SMS
		trfrSMSTraits move = pSpecies->getSMSTraits();
		straigtenPath = move.straigtenPath;
		if (trfr.costMap) {
			outPar << "SMS\tcosts from imported cost map" << endl;
#if !RS_RCPP
			outPar << "FILE NAME: " << costmapname << endl;
#endif
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
#if TEMPMORT
	trfrCRWTraits move = pSpecies->getCRWTraits();      
	switch ((trfr.smType)) {	
		case 0:
			outPar << "constant " << move.stepMort << endl;
			break;
		case 1:
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
			break;
		case 2:
			outPar << "temporally variable:\t" << endl;
			outPar << "FILE NAME: " << mortfilename << endl;
			break;
	} 
#else
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
#endif // TEMPMORT 
} // end of movement process
else { // kernel
	string meandist = "MEAN DISTANCE";
	string probkern = "PROB. KERNEL I";
	trfrKernTraits kern0,kern1;
	trfrKernParams k0,k1;
	outPar << "dispersal kernel" << endl << "TYPE: \t";
#if RS_CONTAIN
	if (trfr.kernType == 1) outPar << "double ";
	if (trfr.kernType <= 1) outPar << "negative exponential" << endl;
	if (trfr.kernType == 2) outPar << "2Dt" << endl;
	if (trfr.kernType == 3) outPar << "WALD" << endl;
#else
	if (trfr.twinKern) outPar << "double ";
	outPar << "negative exponential" << endl;
#endif // RS_CONTAIN 

#if RS_CONTAIN
	if (trfr.kernType >= 2) {
		if (trfr.kernType == 2) { // 2DT 
			trfr2Dt t2 = pSpecies->getTrfr2Dt(); 
			outPar << "Kernel 1   U0: " << t2.u0Kernel1 << " P0: " << t2.p0Kernel1 << endl;
			outPar << "Kernel 2   U0: " << t2.u0Kernel2 << " P0: " << t2.p0Kernel2 << endl;
			outPar << "Prop kernel 1: " << t2.propKernel1 << endl;
		}
		else { // WALD 
			trfrWald w = pSpecies->getTrfrWald();
			outPar << "Mean U: " << w.meanU << " sigma_w: " << w.sigma_w << " hc: " << w.hc << " vt: " << w.vt << " kappa: " << w.kappa << endl;
			outPar << "hr:";
			if (dem.stageStruct) {
				for (int i = 1; i < sstruct.nStages; i++) {
					outPar << "  stage " << i << ": " << pSpecies->getTrfrHr(i);					
				}
			}
			outPar << endl;
			outPar << "Mean wind direction: " << w.meanDirn << " s.d.: " << w.sdDirn << endl;
		}		
	}
	else {
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
				if (trfr.kernType == 1)
#else
				if (trfr.twinKern)
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
				if (trfr.kernType == 1)
#else
				if (trfr.twinKern)
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
				if (trfr.kernType == 1)
#else
				if (trfr.twinKern)
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
				if (trfr.kernType == 1)
#else
				if (trfr.twinKern)
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
				if (trfr.kernType == 1) 
#else
				if (trfr.twinKern) 
#endif // RS_CONTAIN 
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
#if RS_CONTAIN
				if (trfr.kernType == 1) 
#else
				if (trfr.twinKern) 
#endif // RS_CONTAIN 
				{
					outPar << meandist << " II: \t" << kern0.meanDist2 <<endl;
					outPar << probkern << ": \t" << kern0.probKern1 <<endl;
				}
			}
		}
	}
#if RS_CONTAIN
	}
#endif // RS_CONTAIN 

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
#if GOBYMODEL
if (true)
#else
#if SOCIALMODEL
if (true)
#else
if (emig.indVar || trfr.indVar || sett.indVar || d.neutralMarkers)
#endif
#endif
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

#if VIRTUALECOLOGIST
// Landscape genetics

virtParams virt = paramsSim->getVirtParams();
if (virt.landscapeGenetics) {
	outPar << endl << "LANDSCAPE GENETICS:" << endl;
	outPar << "Patch sampling method:             " << virt.patchMethod;
	if (virt.patchMethod == 1) {
		sampleLimits lim = paramsSim->getLimits();
		outPar << " \tMin. X: " << lim.minX << " Min. Y: " << lim.minY ;
		outPar << " \tMax. X: " << lim.maxX << " Max. Y: " << lim.maxY ;
	}
	outPar << endl;
	if (virt.patchMethod == 3) {
		outPar << "FILE NAME: " << patchfilename << endl;;		
	}
	else {
		outPar   << "Max. no. of sampled patches:       " << virt.maxNPatches << endl;
		if (virt.maxPatchNum > 0) 
			outPar   << "Maximum patch number:              " << virt.maxPatchNum << endl;
	}
	if (virt.patchMethod == 2) { // dynamic
		outPar << "No. of rows from range front:      " << virt.rowsFront << endl;
		outPar << "Min. no. of individuals per patch: " << virt.minIndsPatch << endl;
	}
	outPar << "Max. no. of individuals per patch: " << virt.maxIndsPatch << endl;
	outPar << "Stage sampling method:             " << virt.stgMethod << endl;
	if (pSpecies->sampleAllLoci()) {
		outPar << "All loci sampled" << endl;
	}
	else {
		outPar << "Sample file:                       " << locfilename << endl;
		traitAllele sample;
		outPar << "Sampled loci (chrom:locus):       ";
		int nloci = pSpecies->nSampleLoci();
		for (int i = 0; i < nloci; i++) {
			sample = pSpecies->getSampleLocus(i);
			outPar << " " << sample.chromo << ":" << sample.locus;
		}
		outPar << endl;
	}
}

#endif // VIRTUALECOLOGIST

#if SOCIALMODEL

// ADDITIONAL PARAMETERS FOR PROBIS SOCIAL POLYMORPHISM MODEL
socialParams s = pSpecies->getSocialParams();
outPar << endl << "SOCIAL POLYMORPHISM PARAMETERS:" << endl;
outPar << "Phenotype initialisation:  mean " << s.socMean << "  s.d. " << s.socSD
	<< "  scaling factor " << s.socScale << endl;
outPar << "Asocial morph ratio for K:    " << s.asocK << endl;
outPar << "Asocial morph ratio for Rmax: " << s.asocRmax << endl;
outPar << "Asocial morph ratio for bc:   " << s.asocBc << endl;
//outPar << "Growth rates (Fogarty et al., 2011): ra " << soc.ra << " rs " << soc.rs << endl;
outPar << "Allee effect thresholds: Ta " << s.Ta << " Ts " << s.Ts << endl;
//outPar << "Strength of Allee effect below T: Ca " << soc.Ca << " Cs " << soc.Cs << endl;
outPar << "Parameters for Allee effect below T: ca " << s.ca << " ba " << s.ba
	<< " cs " << s.cs << " bs " << s.bs << endl;
outPar << "Fitness-independent dispersal rate (dK):       " << s.dK << endl;
outPar << "Change in dispersal rate with fitness (alpha): " << s.alpha << endl;

#endif // SOCIALMODEL

#if RS_CONTAIN

outPar << endl << "MANAGEMENT CULL:" << endl;

if (pCull->isCullApplied()) {
	culldata c;
	cullstagedata cstage;
	c = pCull->getCullData();
	outPar << "Method:                    ";
	switch (c.method) {	
	case 0:
		outPar << "random" << endl;
		break;
	case 1:
		outPar << "random weighted recently colonised" << endl;
		break;
	case 2:
		outPar << "random weighted closest to damage" << endl;
		break;
	case 3:
		outPar << "random weighted closest to damage X population size" << endl;
		break;
	case 4:
		outPar << "weighted closest to damage" << endl;
		break;
	case 5:
		outPar << "weighted closest to damage X population size" << endl;
		break;
	case 6:
		outPar << "random reactive to incurred damage" << endl;
		break;
	case 7:
		outPar << "reactive to incurred damage" << endl;
		break;
	}
	outPar << "Timing of cull:            ";
	switch (c.timing) {	
	case 0:
		outPar << "before dispersal" << endl;
		break;
	case 1:
		outPar << "after dispersal" << endl;
		break;
	}
	outPar << "Max. no. of cells/patches: " << c.maxNpatches << endl;
//	outPar << "Threshold population:      " << c.popnThreshold << endl;
	outPar << "Threshold density:         " << c.densThreshold << endl;
	outPar << "Count c.v. (%):            " << c.countCV << endl;
	if (dem.stageStruct) {
		outPar << "Stages to be culled:       ";
		bool firststage = true;
		for (int i = 0; i < sstruct.nStages; i++) {
			if (pCull->getCullStage(i)) {
				if (!firststage) outPar << ", ";
				outPar << i;
				firststage = false;			
			}
		}
		outPar << endl;	
	}
	switch (c.cullRate) {		
	case 0: // constant
		outPar << "Cull rate function:        constant " << c.cullMaxRate << endl;
		break;
	case 1: // logistic function of density of stages to be culled
		outPar << "Cull rate function:        logistic function of density" << endl;
		outPar << "Max. cull rate:            " << c.cullMaxRate << endl;
		outPar << "Slope:                     " << c.cullAlpha << endl;
		outPar << "Inflection point:          " << c.cullBeta << endl;
		break;
	case 2: // stage-specific logistic functions of density
		outPar << "Cull rate function:        stage-specific logistic functions of density" << endl;
		for (int i = 0; i < NSTAGES; i++) {
			if (pCull->getCullStage(i)) {
				cstage = pCull->getCullStageData(i);
				outPar << "Stage " << i << ":                   " << "max. cull rate: " << cstage.cullMaxRate
					<< "  slope: " << cstage.cullAlpha << "  inflection point: " << cstage.cullBeta << endl;
			}			
		}
		break;
	}
//	outPar << "Delay (years):             " << c.delay << endl;
//	outPar << "Bias towards range edge:   ";
//	if (c.edgeBias) outPar << "yes"; else outPar << "no"; outPar << endl;
}
else {
	outPar << "Not applied" << endl;
}

if (ppLand.dmgLoaded) {
	damageparams d = pDamageParams->getDamageParams();
	outPar << "Occupancy damage assessment: ";
	switch (d.timing) {		
	case 0:
		outPar << "before"; break;
	case 1:
		outPar << "after"; break;
	}
	outPar << " reproduction" << endl;	
	outPar << "Damage type:                 ";
	switch (d.type) {		
	case 0:
		outPar << "catastrophic"; break;
	case 1:
		outPar << "asymptotic"; break;
	case 2:
		outPar << "logistic"; break;
	}
	outPar << endl;	
	if (d.type > 0) {
		outPar << "Damage-occupancy function:   ";
		switch (d.occOption) {		
		case 0:
			outPar << "total population size, "; break;
		case 1:
			outPar << "population density, "; break;
		case 2:
			outPar << "density of culled stages, "; break;
		case 3:
			outPar << "stage-specific density, stage "; 
			outPar << d.stage << ", ";
			break;
		}
		if (d.type == 2) outPar << "alpha " << d.alphaOccupancy << " ";
		outPar << "beta " << d.betaOccupancy << endl;
		if (trfr.moveModel && trfr.moveType == 1) { // SMS
			outPar << "Damage-traversal function:   ";
			if (d.type == 2) outPar << "alpha " << d.alphaTraversal << " ";
			outPar << "beta " << d.betaTraversal << endl;
		}		
	}
}

#endif // RS_CONTAIN 

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
#if VIRTUALECOLOGIST
if (virt.landscapeGenetics) {
	outPar << "Landscape genetics - every " << virt.outInt << " year";
	if (virt.outInt > 1) outPar << "s";
	if (virt.outStart > 0) outPar << " starting year " << virt.outStart;
	outPar << endl;
	outPar << "Output sampled genomes: ";
	if (virt.outGenomes) outPar << "yes"; else outPar << "no";
	outPar << endl;
}
#endif

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
#if RS_RCPP
	if (sim.outPaths) {
		outPar << "SMS paths - every " << sim.outIntPaths << " year";
		if (sim.outIntPaths > 1) outPar << "s";
		if (sim.outStartPaths > 0) outPar << " starting year " << sim.outStartPaths;
		outPar << endl;
	}
#endif
#if RS_CONTAIN
if (sim.outDamage) {
	outPar << "Damage indices - every " << sim.outIntDamage << " year";
	if (sim.outIntDamage > 1) outPar << "s";
	if (sim.outStartDamage > 0) outPar << " starting year " << sim.outStartDamage;
	outPar << endl;
}
#endif // RS_CONTAIN 
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
