/*----------------------------------------------------------------------------
 *	
 *	Copyright (C) 2020 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
 *	
 *	This file is part of RangeShiftR.
 *	
 *	RangeShiftR is free software: you can redistribute it and/or modify
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
 *	along with RangeShiftR. If not, see <https://www.gnu.org/licenses/>.
 *	
 --------------------------------------------------------------------------*/
 
 
/*------------------------------------------------------------------------------

RangeShifter v2.0 Main

Entry level function for the R-package RangeshiftR.

For compilation with GCC g++

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe'er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species' responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Author: Anne-Kathleen Malchow, Humboldt University Berlin
        large parts modified from 'Main.cpp' and 'BatchMode.cpp' created by
        Steve Palmer, University of Aberdeen

------------------------------------------------------------------------------*/

#include "Rinterface.h"

class msghdrs1;
string habmapname, patchmapname, distnmapname; // req'd for compilation, but not used
string costmapname, genfilename;               // ditto
vector<string> hfnames;                        // ditto

paramGrad* paramsGrad;   // pointer to environmental gradient parameters
paramStoch* paramsStoch; // pointer to environmental stochasticity parameters
paramInit* paramsInit;   // pointer to initialisation parameters
paramSim* paramsSim;     // pointer to simulation parameters

Species* pSpecies; // pointer to species
Community* pComm;  // pointer to community
RSrandom* pRandom; // pointer to random number routines

#if RSDEBUG
ofstream DEBUGLOG;
ofstream MUTNLOG;
#endif
// ofstream batchlog;
ofstream rsLog; // performance log for recording simulation times, etc.


// global variables passed between parsing functions...
int batchnum;
int patchmodel, resolution, landtype, maxNhab, speciesdist, distresolution;
int reproductn;
int repseasons;
int stagestruct, stages, transfer;
int sexesDem;  // no. of explicit sexes for demographic model
int sexesDisp; // no. of explicit sexes for dispersal model
int firstsimul;
int fileNtraits; // no. of traits defined in genetic architecture file
rasterdata landraster,patchraster,spdistraster,costsraster;
// rasterdata landraster;
// ...including names of the input files
// string parameterFile;
// string landFile;
string name_landscape, name_patch, name_sp_dist, name_costfile;
//string name_dynland;
//#if RS_CONTAIN
// string name_damagefile;
// string name_managefile;
//#endif // RS_CONTAIN
//#if SPATIALMORT
// string name_mortfile[2];
//#endif // SPATIALMORT
//#if SEASONAL
// string seasonFile;
//#endif // SEASONAL
// string stageStructFile,transMatrix;
// string emigrationFile,transferFile,settleFile,geneticsFile,initialFile;
//#if VIRTUALECOLOGIST
// string virtEcolFile;
//#endif // VIRTUALECOLOGIST
//#if RS_ABC
// string abcParamsFile,abcObsFile;
// int nABCsamples;
//#endif // RS_ABC
// string prevInitialIndsFile = " ";
//
string msgnlines = "No. of lines for final Simulation ";
string msgshldbe = " should be ";
string msgresol0 = "*** Resolution of ";
string msgresol1 = " does not match set Resolution ";
string msgresol2 = " does not match set Distribution Resolution ";
string msghdrs0 = "*** Headers of ";
string msghdrs1 = " do not match headers of LandscapeFile";
string msgpatch = " is required for patch-based model";
string msgmatch = " must match the specification exactly";
string msgcase = " case-sensitive parameter names";
//
// float **matrix = NULL;	// temporary matrix used in batch mode
// int matrixsize = 0; 		// size of temporary matrix

//----------------------------

//-------- R interface globals
// static Rcpp::Environment global = Rcpp::Environment::global_env();
// static Rcpp::S4 ParMaster;
// static Landscape *LandParamsR;
static bool anyIndVar;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


/* # // [[Rcpp::export(name = ".run_from_R")]]
// # // [[Rcpp::export]]
Rcpp::List run_from_R(Rcpp::S4 ParMaster, Rcpp::String dirpath)
{
	return BatchMainR(dirpath, ParMaster);
}*/

/* this function is for calling the file batch version from R;
 * it is excluded atm since it is not supposed to be used before the GUI can not generate the batch parameter files from the settings made;
 * to include, all file parsing and reading functions have to be included from the BatchMode.h/.cpp (like in the standalone version)
#if !RS_RCPP
Rcpp::List BatchMainFile(string dirpath, Rcpp::S4 ParMaster)
{
	// int i,t0,t1,Nruns;
	int i, t0, t1;
	int nSimuls, nLandscapes; // no. of simulations and landscapes in batch

	t0 = time(0);

	// set up parameter objects
	paramsGrad = new paramGrad;
	paramsStoch = new paramStoch;
	paramsInit = new paramInit;
	paramsSim = new paramSim;

	// set up working directory and control file name

	paramsSim->setDir(dirpath); // full path name of directory passed as a parameter
	// NOT IMPLEMENTED: // control file number also passed as a parameter:
	// int i = atoi(argv[2]);
	// cname = paramsSim->getDir(0) + "Inputs/CONTROL" + Int2Str(i) + ".txt";
	cname = paramsSim->getDir(0) + "Inputs/CONTROL.txt";

#if RSDEBUG
	Rcpp::Rcout << endl << "Working directory: " << paramsSim->getDir(0) << endl;
	Rcpp::Rcout << endl << "Control file:      " << cname << endl << endl;
#endif

	bool errorfolder = CheckDirectory();
	if(errorfolder) {
		Rcpp::Rcout << endl << "***** Invalid working directory: " << paramsSim->getDir(0) << endl << endl;
		Rcpp::Rcout << "***** Working directory must contain Inputs, Outputs and Output_Maps folders" << endl << endl;
		Rcpp::Rcout << "*****" << endl;
		Rcpp::Rcout << "***** Simulation ABORTED " << endl;
		Rcpp::Rcout << "*****" << endl;
		return Rcpp::List::create(666);
	}

#if RSDEBUG
	// set up debugging log file
	string name = paramsSim->getDir(2) + "DebugLog.txt";
	DEBUGLOG.open(name.c_str());
	name = paramsSim->getDir(2) + "MutnLog.txt";
	MUTNLOG.open(name.c_str());
	// DEBUGLOG << "Main(): paramsSim = " << paramsSim << endl;
	if(DEBUGLOG.is_open())
		Rcpp::Rcout << endl << "Main(): DEBUGLOG is open" << endl << endl;
	else
		Rcpp::Rcout << endl << "Main(): DEBUGLOG is NOT open" << endl << endl;
#endif

	// set up species
	// FOR MULTI-SPECIES MODEL, THERE WILL BE AN ARRAY OF SPECIES POINTERS
	// OR A COMMUNITY CLASS TO HOLD THE SPECIES
	pSpecies = new Species;
	demogrParams dem = pSpecies->getDemogr();
	stageParams sstruct = pSpecies->getStage();
	trfrRules trfr = pSpecies->getTrfr();

	batchfiles b;
	string indir = paramsSim->getDir(1);
	string outdir = paramsSim->getDir(2);
	b = ParseControlFile(cname, indir, outdir);
	if(b.ok) {
		nSimuls = b.nSimuls;
		nLandscapes = b.nLandscapes;
		dem.repType = b.reproductn;
		dem.repSeasons = b.repseasons;
		if(b.stagestruct == 0)
			dem.stageStruct = false;
		else
			dem.stageStruct = true;
		sstruct.nStages = b.stages;
		if(b.transfer == 0)
			trfr.moveModel = false;
		else {
			trfr.moveModel = true;
			trfr.moveType = b.transfer;
		}
		Rcpp::Rcout << endl << "Batch input files OK" << endl;
		pSpecies->setDemogr(dem);
		pSpecies->setStage(sstruct);
		pSpecies->setTrfr(trfr);
		simParams sim = paramsSim->getSim();
		sim.batchMode = true;
		sim.batchNum = b.batchNum;
		paramsSim->setSim(sim);
	} else {
		Rcpp::Rcout << endl << "Error in parsing batch input files - " << endl;
	}
#if RSDEBUG
	DEBUGLOG << "Main(): dem.repType = " << dem.repType << endl;
#endif

	clear_outPop();

// set up random number stream
#if RSDEBUG
	pRandom = new RSrandom(666);
#else
	pRandom = new RSrandom(0);
#endif

#if RANDOMCHECK
	randomCheck();
#else
	if(b.ok) {
		RunBatch(nSimuls, nLandscapes);
	}
#endif

	delete pRandom;

#if RSDEBUG
	if(DEBUGLOG.is_open()) {
		DEBUGLOG.close();
		DEBUGLOG.clear();
	}
	if(MUTNLOG.is_open()) {
		MUTNLOG.close();
		MUTNLOG.clear();
	}
#endif

	delete paramsGrad;
	delete paramsStoch;
	delete paramsInit;
	delete paramsSim;
	delete pSpecies;

	t1 = time(0);
	Rcpp::Rcout << endl << "***** Elapsed time " << t1 - t0 << " seconds" << endl << endl;

	Rcpp::Rcout << "*****" << endl;
	Rcpp::Rcout << "***** Simulation completed " << endl;
	Rcpp::Rcout << "***** Outputs folder: " << outdir << endl;
	Rcpp::Rcout << "*****" << endl;


	return Rcpp::List::create(0);
}
#endif // !RS_RCPP
*/


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//------ R interface --------------------------------------------------------
//---------------------------------------------------------------------------


//[[Rcpp::export(name = "run_from_R")]]
Rcpp::List BatchMainR(std::string dirpath, Rcpp::S4 ParMaster)
{
	// int i,t0,t1,Nruns;
	int t0, t1;
	int nSimuls, nLandscapes; // no. of simulations and landscapes in batch

	t0 = time(0);
	//clear_outPop();

	// set up parameter objects
	paramsGrad = new paramGrad;
	paramsStoch = new paramStoch;
	paramsInit = new paramInit;
	paramsSim = new paramSim;

	// set up working directory and log file names
	paramsSim->setDir(dirpath);

#if RSDEBUG
	Rcpp::Rcout << endl << "Working directory: " << paramsSim->getDir(0) << endl;
#endif

	bool errorfolder = CheckDirectory();
	if(errorfolder) {
		Rcpp::Rcout << endl << "***** Invalid working directory: " << paramsSim->getDir(0) << endl << endl;
		Rcpp::Rcout << "***** Working directory must contain Inputs, Outputs and Output_Maps folders" << endl << endl;
		Rcpp::Rcout << "*****" << endl;
		Rcpp::Rcout << "***** Simulation ABORTED " << endl;
		Rcpp::Rcout << "*****" << endl;
		return Rcpp::List::create(Rcpp::Named("Errors") = 666);
	}

#if RSDEBUG
	// set up debugging log file
	string logname = paramsSim->getDir(2) + "DebugLog.txt";
	DEBUGLOG.open(logname.c_str());
	logname = paramsSim->getDir(2) + "MutnLog.txt";
	MUTNLOG.open(logname.c_str());
	if(DEBUGLOG.is_open())
		Rcpp::Rcout << endl << "Main(): DEBUGLOG is open" << endl << endl;
	else
		Rcpp::Rcout << endl << "Main(): DEBUGLOG is NOT open" << endl << endl;
#else
#endif


	// set global variables
	Rcpp::S4 control("ControlParams");
	// ParMaster = global.get(ParMaster_name);
	control = Rcpp::as<Rcpp::S4>(ParMaster.slot("control"));
	setglobalvarsR(control);

	nSimuls = 1;
	nLandscapes = 1;

	//////////////////////// batchfiles ParseControlR(string dirpath) ///////////////////

	//string indir  = dirpath + "Inputs/";
	string outdir = dirpath + "Outputs/";

	// Doublecheck global variables
	Rcpp::Rcout << "Checking Control parameters " << endl;

	int errors = 0;
	string filetype = "Control parameters";

	if(batchnum < 0) {
		BatchErrorR(filetype, -999, 19, "BatchNum");
		errors++;
	}

	if(patchmodel < 0 || patchmodel > 1) {
		BatchErrorR(filetype, -999, 1, "PatchModel");
		errors++;
	}

	if(resolution < 1) {
		BatchErrorR(filetype, -999, 11, "Resolution");
		errors++;
	}

	if(landtype != 0 && landtype != 2 && landtype != 9) {
		BatchErrorR(filetype, -999, 0, "LandType");
		Rcpp::Rcout << "LandType must be 0, 2 or 9" << endl;
		errors++;
	} else {
		if(landtype == 9 && patchmodel) {
			BatchErrorR(filetype, -999, 0, "LandType");
			Rcpp::Rcout << "LandType may not be 9 for a patch-based model" << endl;
			errors++;
		}
	}

	if(landtype == 0) { // raster with unique habitat codes
		if(maxNhab < 2) {
			BatchErrorR(filetype, -999, 12, "MaxHabitats");
			errors++;
		}
	} else { // raster with habitat quality OR artificial landscape
		if(maxNhab != 1) {
			BatchErrorR(filetype, -999, 0, " ");
			errors++;
			Rcpp::Rcout << "MaxHabitats must be 1 for LandType = " << landtype << endl;
		}
		/*else
		{
		if (landtype == 9) // artificial landscape
		// although the user enters 1, the actual number of habitats is 2
		b.maxNhab = 2; // only changes maxNhab in batchfile b, the global variable remains maxNhab=1
		}*/
	}

	if(speciesdist < 0 || speciesdist > 1) {
		BatchErrorR(filetype, -999, 1, "SpeciesDist");
		errors++;
	} else {
		if(speciesdist != 0 && landtype == 9) {
			BatchErrorR(filetype, -999, 0, "SpeciesDist");
			Rcpp::Rcout << "SpeciesDist must be 0 for an artificial landscape" << endl;
			errors++;
		}
	}

	if(speciesdist == 1) { // distribution resolution is required
		if(distresolution < resolution) {
			BatchErrorR(filetype, -999, 0, "DistResolution");
			Rcpp::Rcout << "DistResolution may not be less than Resolution" << endl;
			errors++;
		} else {
			if(distresolution % resolution) {
				BatchErrorR(filetype, -999, 0, "DistResolution");
				Rcpp::Rcout << "DistResolution must be an integer multiple of Resolution" << endl;
				errors++;
			}
		}
	}

	if(reproductn < 0 || reproductn > 2) {
		BatchErrorR(filetype, -999, 2, "Reproduction");
		errors++;
	}
	else {
		switch(reproductn) {
		case 0: {
			sexesDem = 1;
			sexesDisp = 1;
			break;
		}
		case 1: {
			sexesDem = 1;
			sexesDisp = 2;
			break;
		}
		case 2: {
			sexesDem = 2;
			sexesDisp = 2;
			break;
		}
		}
	}

	if(repseasons < 1) {
		BatchErrorR(filetype, -999, 11, "RepSeasons");
		errors++;
	}

	if(stagestruct < 0 || stagestruct > 1) {
		BatchErrorR(filetype, -999, 1, "StageStruct");
		errors++;
	}

	if(stagestruct) {
		if(stages < 2 || stages > NSTAGES) {
			BatchErrorR(filetype, -999, 0, " ");
			errors++;
			Rcpp::Rcout << "Stages must be between 2 and " << NSTAGES << endl;
		}
	} else { // non-stage-structured model must have 2 stages
		stages = 2;
	}

	if(transfer < 0 || transfer > 2) {
		BatchErrorR(filetype, -999, 2, "Transfer");
		errors++;
	}

	if(errors > 0) { // terminate batch error checking

		Rcpp::Rcout << endl << "*** Control parameters must be corrected before proceeding" << endl;
		return Rcpp::List::create(Rcpp::Named("Errors") = errors);
	}

	//////////////////////////////////////


	// set up species
	// FOR MULTI-SPECIES MODEL, THERE WILL BE AN ARRAY OF SPECIES POINTERS
	// OR A COMMUNITY CLASS TO HOLD THE SPECIES
	pSpecies = new Species;
	demogrParams dem = pSpecies->getDemogr();
	stageParams sstruct = pSpecies->getStage();
	trfrRules trfr = pSpecies->getTrfr();

	if(errors == 0) {
		//   nSimuls = b.nSimuls;
		//   nLandscapes = b.nLandscapes;
		dem.repType = reproductn;
		dem.repSeasons = repseasons;
		if(stagestruct == 0)
			dem.stageStruct = false;
		else
			dem.stageStruct = true;
		sstruct.nStages = stages;
		if(transfer == 0)
			trfr.moveModel = false;
		else {
			trfr.moveModel = true;
			trfr.moveType = transfer;
		}
		pSpecies->setDemogr(dem);
		pSpecies->setStage(sstruct);
		pSpecies->setTrfr(trfr);

		simParams sim = paramsSim->getSim();
		sim.batchMode = true;
		sim.batchNum = batchnum;
		paramsSim->setSim(sim);
		Rcpp::Rcout << endl << "Control Parameters checked" << endl;
	} else {
		Rcpp::Rcout << endl << "Error in parsing parameters - " << endl;
	}

#if RSDEBUG
	DEBUGLOG << "BatchMainR(): dem.repType = " << dem.repType << endl;
#endif

	// set up random number stream
	std::int64_t seed = Rcpp::as<std::int64_t>(control.slot("seed"));

	if(seed == 0) { // don't set seed from R
#if RSDEBUG
		pRandom = new RSrandom(666);  // fixed debug seed
#else
		pRandom = new RSrandom(-1);  // random seed
#endif
	}
	else pRandom = new RSrandom(seed);

	//Rcpp::RNGScope rngScope;

	Rcpp::List list_outPop;
	if(errors == 0) {
		Rcpp::Rcout << endl << "Run Simulation(s)";
		if(seed > 0) {
			Rcpp::Rcout << " with seed " << seed;
		}else{
			Rcpp::Rcout << " with random seed";
		}
		Rcpp::Rcout << " ..." << endl;

		list_outPop = RunBatchR(nSimuls, nLandscapes, ParMaster);
	}

#if RSDEBUG
	if(DEBUGLOG.is_open()) {
		DEBUGLOG.close();
		DEBUGLOG.clear();
	}
	if(MUTNLOG.is_open()) {
		MUTNLOG.close();
		MUTNLOG.clear();
	}
#endif

	simParams sim = paramsSim->getSim();

	delete paramsGrad;
	delete paramsStoch;
	delete paramsInit;
	delete paramsSim;
	delete pSpecies;

	delete pRandom;

	t1 = time(0);
	Rcpp::Rcout << endl << "***** Elapsed time: " << t1 - t0 << " seconds" << endl << endl;

	Rcpp::Rcout << "*****" << endl;
	Rcpp::Rcout << "***** Simulation completed " << endl;
	Rcpp::Rcout << "***** Outputs folder: " << outdir << endl;
	Rcpp::Rcout << "*****" << endl;

	if(sim.ReturnPopRaster && sim.outIntPop > 0) {
		// return Rcpp::List::create(Rcpp::Named("runs") = errors);
		return list_outPop;
	} else {
		return Rcpp::List::create(Rcpp::Named("Errors") = errors);
	}
}

//---------------------------------------------------------------------------

bool ReadLandParamsR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{
	landParams ppLand = pLandscape->getLandParams();
	genLandParams ppGenLand = pLandscape->getGenLandParams();

	Rcpp::S4 LandParamsR("LandParams");
	LandParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("land"));

	int errors = 0;

	if(landtype == 9) { // artificial landscape
#if RSDEBUG
		DEBUGLOG << "ReadLandParamsR(): artificial: " << endl;
#endif
		ppLand.nHab = 2;
		ppLand.rasterType = 9;
		ppLand.landNum = Rcpp::as<int>(LandParamsR.slot("LandNum"));
		ppGenLand.fractal = Rcpp::as<bool>(LandParamsR.slot("fractal"));
		ppGenLand.continuous = Rcpp::as<bool>(LandParamsR.slot("continuous"));
		ppLand.dimX = Rcpp::as<int>(LandParamsR.slot("dimX"));
		ppLand.dimY = Rcpp::as<int>(LandParamsR.slot("dimY"));
		ppGenLand.minPct = Rcpp::as<float>(LandParamsR.slot("minPct"));
		ppGenLand.maxPct = Rcpp::as<float>(LandParamsR.slot("maxPct"));
		ppGenLand.propSuit = Rcpp::as<float>(LandParamsR.slot("propSuit"));
		ppGenLand.hurst = Rcpp::as<float>(LandParamsR.slot("hurst"));
		ppLand.maxX = ppLand.dimX - 1;
		ppLand.maxY = ppLand.dimY - 1;

		// check dimensions  -->  these shouldn't occur, as already checked for on R-level
		if(ppGenLand.fractal && ppLand.maxX > ppLand.maxY) {
			//return -901;
			// fix it by swapping X and Y:
			ppLand.maxX = ppLand.dimY - 1;
			ppLand.maxY = ppLand.dimX - 1;
			ppLand.dimX = ppLand.maxX + 1;
			ppLand.dimY = ppLand.maxY + 1;
		}
		if(ppGenLand.fractal) {
			if((ppLand.dimX < 3 || ppLand.dimX % 2 != 1) // why only check for uneven number, not ((power of 2)-1) ?
			   || (ppLand.dimY < 3 || ppLand.dimY % 2 != 1)) {
				//return -902;  // -> no action, should be covered on R-level
			}
		}
		// SCFP 26/9/13 - min and max habitat percentages need to be set for all types of
		// fractal landscape (including discrete), as they are passed to the fractal generator
		// NOTE that will not have been checked for a discrete landscape
		if(ppGenLand.fractal && !ppGenLand.continuous) {
			ppGenLand.minPct = 1;
			ppGenLand.maxPct = 100;
		}
	} else { // imported raster map
		ppLand.landNum = Rcpp::as<int>(LandParamsR.slot("LandNum"));
		ppLand.nHab = Rcpp::as<int>(LandParamsR.slot("Nhabitats")); // no longer necessary to read no. of habitats from landFile

		Rcpp::IntegerVector dynland_years;
		Rcpp::StringVector habitatmaps, patchmaps, costmaps;
		dynland_years = Rcpp::as<Rcpp::IntegerVector>(LandParamsR.slot("DynamicLandYears"));
		if(dynland_years.size() == 1 && dynland_years[0] == 0 ) ppLand.dynamic = false;
		else ppLand.dynamic = true;
		if(ppLand.dynamic) {
			habitatmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("LandscapeFile"));
			patchmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("PatchFile"));
			costmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("CostsFile"));
			name_landscape = habitatmaps(0);
			name_patch = patchmaps(0);
			name_costfile = costmaps(0);
		} else {
			name_landscape = Rcpp::as<string>(LandParamsR.slot("LandscapeFile"));
			name_patch = Rcpp::as<string>(LandParamsR.slot("PatchFile"));
			name_costfile = Rcpp::as<string>(LandParamsR.slot("CostsFile"));
		}
		if(!patchmodel && name_patch != "NULL") Rcpp::Rcout << "PatchFile must be NULL in a cell-based model!" << endl;

		name_sp_dist = Rcpp::as<string>(LandParamsR.slot("SpDistFile"));

		if(landtype == 2)
			ppLand.nHab = 1; // habitat quality landscape has one habitat class

		// CHECK IMPORTED RASTER FILES
		string indir = paramsSim->getDir(1);
		string ftype, fname;
		string filetype = "LandFile";
		//rasterdata patchraster, spdistraster; //declared globally

		// check landscape filename
		ftype = "LandscapeFile";
		fname = indir + name_landscape;
		landraster = ParseRasterHead(fname);
		if(landraster.ok) {
			if(landraster.cellsize == resolution)
				Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
			else {
				errors++;
				Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol1 << endl;
			}
		} else {
			errors++;
			if(landraster.errors == -111)
				OpenErrorR(ftype, fname);
			else
				FormatErrorR(fname, landraster.errors);
		}

		// check patch map filename
		ftype = "PatchFile";
		if(name_patch == "NULL") {
			if(patchmodel) {
				BatchErrorR(filetype, -999, 0, " ");
				errors++;
				Rcpp::Rcout << ftype << msgpatch << endl;
			}
		} else {
			if(patchmodel) {
				fname = indir + name_patch;
				patchraster = ParseRasterHead(fname);
				if(patchraster.ok) {
					// check resolutions match
					if(patchraster.cellsize == resolution) {
						if(!errors) {
							if(patchraster.cellsize == landraster.cellsize) {
								// check that extent matches landscape extent
								if(patchraster.ncols == landraster.ncols && patchraster.nrows == landraster.nrows) {
									// check origins match
									if((int)patchraster.xllcorner == (int)landraster.xllcorner &&
									   (int)patchraster.yllcorner == (int)landraster.yllcorner) {
										Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
									} else {
										Rcpp::Rcout << "*** Origin co-ordinates of " << ftype << msghdrs1 << endl;
										errors++;
									}
								} else {
									Rcpp::Rcout << "*** Extent of " << ftype << " " << fname << msghdrs1 << endl;
									errors++;
								}
							} else {
								Rcpp::Rcout << msgresol0 << ftype << " " << fname << msghdrs1 << endl;
								errors++;
							}
						}
					} else {
						Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol1 << endl;
						errors++;
					}
				} else {
					errors++;
					if(patchraster.errors == -111)
						OpenErrorR(ftype, fname);
					else
						FormatErrorR(fname, patchraster.errors);
				}
			}
		}

		// check cost map filename
		ftype = "CostMapFile";
		if (name_costfile == "NULL") {
			if (transfer == 1) { // SMS
				if (landtype == 2) { // habitat quality
					BatchErrorR(filetype, -999, 0, " ");
					errors++;
					Rcpp::Rcout << ftype << " is required for a habitat quality landscape" << endl;
				}
			}
		}
		else {
			if (transfer == 1) { // SMS
				fname = indir + name_costfile;
				costsraster = ParseRasterHead(fname);
				if(costsraster.ok) {
					// check resolutions match
					if(costsraster.cellsize == resolution) {
						if(!errors) {
							if(costsraster.cellsize == landraster.cellsize) {
								// check that extent matches landscape extent
								if(costsraster.ncols == landraster.ncols && costsraster.nrows == landraster.nrows) {
									// check origins match
									if((int)costsraster.xllcorner == (int)landraster.xllcorner &&
									   (int)costsraster.yllcorner == (int)landraster.yllcorner) {
										Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
									} else {
										Rcpp::Rcout << "*** Origin co-ordinates of " << ftype << msghdrs1 << endl;
										errors++;
									}
								} else {
									Rcpp::Rcout << "*** Extent of " << ftype << " " << fname << msghdrs1 << endl;
									errors++;
								}
							} else {
								Rcpp::Rcout << msgresol0 << ftype << " " << fname << msghdrs1 << endl;
								errors++;
							}
						}
					} else {
						Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol1 << endl;
						errors++;
					}
				} else {
					errors++;
					if(costsraster.errors == -111)
						OpenErrorR(ftype, fname);
					else
						FormatErrorR(fname, costsraster.errors);
				}
			}
			else {
				BatchErrorR(filetype, -999, 0, " ");
				errors++;
				Rcpp::Rcout << ftype << " must be NULL if transfer model is not SMS" << endl;
			}
		}

		// check dynamic landscape filename
		// ...most checks are done at reading time in ReadDynLandR()
		ftype = "Dynamic landscape";
		if(ppLand.dynamic) {
			// check valid years
			if(dynland_years[0]!=0) {
				errors++;
				Rcpp::Rcout << "First year in dynamic landscape must be 0." << endl;
			} else {
				for(int i=1; i<dynland_years.size(); i++ ) {
					if(dynland_years[i-1] >= dynland_years[i]) {
						errors++;
						Rcpp::Rcout << "Year in dynamic landscape must strictly increase." << endl;
					}
				}
			}
			if(dynland_years.size() != habitatmaps.size()) {
				errors++;
				Rcpp::Rcout << "Dynamic landscape: Years must have as many elements as habitat maps." << endl;
			}
			if(patchmodel) {
				if( dynland_years.size() != patchmaps.size() ||
				    habitatmaps.size()   != patchmaps.size() ) {
					errors++;
					Rcpp::Rcout << "Dynamic landscape: Patchmaps must have as many elements as Years and habitat maps." << endl;
				}
			}
			if (name_costfile != "NULL") {
				if( dynland_years.size() != costmaps.size() ||
				    habitatmaps.size()   != costmaps.size() ) {
					errors++;
					Rcpp::Rcout << "Dynamic landscape: Costmaps must have as many elements as Years and habitat maps." << endl;
				}
			}
			if(errors==0) {
				// store land changes
				string landchangefile,patchchangefile;
				landChange chg;
				for(int i=1; i<dynland_years.size(); i++ ) {
					chg.chgnum = i;
					chg.chgyear = dynland_years[i];
					chg.habfile = indir + habitatmaps(i);
					if(patchmodel) chg.pchfile = indir + patchmaps(i);
					else chg.pchfile = "NULL";
					if (name_costfile == "NULL") chg.costfile = "none";
					else chg.costfile = indir + costmaps(i);
					pLandscape->addLandChange(chg);
				}
			}
		}

		// check initial distribution map filename
		ftype = "Species Distribution map";
		if(name_sp_dist == "NULL") {
			if(speciesdist) {
				BatchErrorR(filetype, -999, 0, " ");
				errors++;
				Rcpp::Rcout << ftype << " is required as SpeciesDist is 1 in Control" << endl;
			}
		} else {
			if(speciesdist) {
				fname = indir + name_sp_dist;
				spdistraster = ParseRasterHead(fname);
				if(spdistraster.ok) {
					if(spdistraster.cellsize == distresolution) {
						if(!errors) {
							// check origins match
							if((int)spdistraster.xllcorner == (int)landraster.xllcorner && (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
								// check extents match
								if(spdistraster.cellsize == landraster.cellsize) {
									// same resolution
									if(spdistraster.ncols == landraster.ncols && spdistraster.nrows == landraster.nrows) {
										Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
									} else {
										Rcpp::Rcout << "*** Extents of " << ftype << msghdrs1 << endl;
										errors++;
									}
								} else { // different resolution
									if((spdistraster.cellsize % landraster.cellsize)==0) {
										double coarse = spdistraster.cellsize/landraster.cellsize;
										double rightx = ceil(landraster.ncols/coarse);
										double righty = ceil(landraster.nrows/coarse);
										if( spdistraster.ncols == int(rightx) && spdistraster.nrows == int(righty) ) {
											Rcpp::Rcout << ftype << " headers OK: " << fname << endl;
										} else {
											Rcpp::Rcout << "*** Extents of " << ftype << msghdrs1 << endl;
											errors++;
										}
									} else {
										BatchErrorR(filetype, -999, 0, "DistResolution");
										Rcpp::Rcout << "Resolution of initial distribution must be an integer multiple of Landscape resolution" << endl;
										errors++;
									}
								}
							} else {
								Rcpp::Rcout << "*** Origin co-ordinates of " << ftype << msghdrs1 << endl;
								errors++;
							}
						}
					} else {
						Rcpp::Rcout << msgresol0 << ftype << " " << fname << msgresol2 << endl;
						errors++;
					}
				} else {
					errors++;
					if(spdistraster.errors == -111)
						OpenErrorR(ftype, fname);
					else
						FormatErrorR(fname, spdistraster.errors);
				}
			}
		}

	}

	pLandscape->setLandParams(ppLand, true);
	pLandscape->setGenLandParams(ppGenLand);

#if RSDEBUG
// DEBUGLOG << "ReadLandParamsR(): NHab=" << ppLand.nHab << endl;
	DEBUGLOG << "ReadLandParamsR(): ppLand.landNum=" << ppLand.landNum << endl;
#endif

	if(errors) return false; //=landOK
	else return true;
}

//---------------------------------------------------------------------------

int ReadDynLandR(Landscape *pLandscape, Rcpp::S4 LandParamsR)
{

#if RSDEBUG
	DEBUGLOG << "ReadDynLandR(): pLandscape=" << pLandscape << endl;
#endif

	Rcpp::StringVector habitatmaps, patchmaps, costmaps;
	habitatmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("LandscapeFile"));
	if (patchmodel) {
		patchmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("PatchFile"));
	}
	bool costs = false;
	costmaps = Rcpp::as<Rcpp::StringVector>(LandParamsR.slot("CostsFile"));
	if (costmaps(0) != "NULL") costs = true;


	//------------ int ParseDynamicFile(string indir) {

	string indir = paramsSim->getDir(1);
	string fname,ftype;
	wstring header;
	wifstream hfile,pfile,cfile;
	int errors,ncols,nrows,cellsize,habnodata,pchnodata,costnodata;
	double xllcorner,yllcorner;

	errors = ncols = nrows = cellsize = habnodata = pchnodata = costnodata = 0;
	xllcorner = yllcorner = 0.0;

	if (patchmodel) {
		pLandscape->createPatchChgMatrix();
	}
	if (costs) {
		pLandscape->createCostsChgMatrix();
	}

	for(int i=1; i < habitatmaps.size(); i++ ) {

		// Habitat change file
		fname = indir + habitatmaps(i);

		// open file
#if RSWIN64
		hfile.open(fname.c_str());
#else
		hfile.open(fname, std::ios::binary);
#endif
		if(!hfile.is_open()) {
			OpenErrorR("Dynamic landscape habitat map ",  fname);
#if RSDEBUG
			DEBUGLOG << "Dynamic landscape habitat map failed to open: " << fname << std::endl;
#endif
			hfile.clear();
			return -212;
		} else {
#if RSDEBUG
			DEBUGLOG << "Dynamic landscape habitat map #" << i << " open to read" << std::endl;
#endif
#if !RSWIN64
			// check BOM for UTF-16
			if(check_bom(fname) == "utf16")
				// apply BOM-sensitive UTF-16 facet
				hfile.imbue(std::locale(hfile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif
			// ASCII header
			hfile >> header;
			if (!hfile.good()) {
#if RSDEBUG
				DEBUGLOG << "ReadDynLandR(): failed to read landscape habitat map #" << i << ": " << fname << std::endl;
#endif
				errors = -1112;
				hfile.close();
				hfile.clear();
				return errors;
			}
			if (header != L"ncols" && header != L"NCOLS") errors++;
			hfile >> ncols;

			hfile >> header >> nrows;
			if (header != L"nrows" && header != L"NROWS") errors++;

			hfile >> header >> xllcorner;
			if (header != L"xllcorner" && header != L"XLLCORNER") errors++;

			hfile >> header >> yllcorner;
			if (header != L"yllcorner" && header != L"YLLCORNER") errors++;

			hfile >> header >> cellsize;
			if (header != L"cellsize" && header != L"CELLSIZE") errors++;

			hfile >> header >> habnodata;
			if (header != L"NODATA_value" && header != L"NODATA_VALUE") errors++;

			if (errors > 0)  {
				FormatErrorR(fname,errors);
#if RSDEBUG
				DEBUGLOG << "ReadDynLandR(): failed to read Raster header of landscape habitat map #" << i << ": " << fname << std::endl;
#endif
				hfile.close();
				hfile.clear();
			} else {
				// check resolution match
				if (cellsize == resolution) {
					// check that extent matches landscape extent
					if(ncols == landraster.ncols && nrows == landraster.nrows) {
						// check origins match
						if((int)xllcorner == (int)landraster.xllcorner &&
						   (int)yllcorner == (int)landraster.yllcorner) {
							Rcpp::Rcout << "Dynamic landscape habitat map #" << i << ", headers OK: " << fname << endl;
						} else {
							Rcpp::Rcout << "*** Origin co-ordinates of " << fname << msghdrs1 << endl;
							errors++;
						}
					} else {
						Rcpp::Rcout << "*** Extent of " << fname << msghdrs1 << endl;
						errors++;
					}
				} else {
					Rcpp::Rcout << "*** Resolution of " << fname << msghdrs1 << endl;
					errors++;
				}
				if (errors > 0)  {
					FormatErrorR(fname,errors);
#if RSDEBUG
					DEBUGLOG << "ReadDynLandR(): Errors in raster header of landscape habitat map #" << i << ": " << fname << std::endl;
#endif
					hfile.close();
					hfile.clear();
				}
			} // end of reading ASCII header
		}

		// Do the same for corresponding patch map, if applicable
		if (patchmodel) {
			// Patch change file
			fname = indir + patchmaps(i);

			// open file
#if RSWIN64
			pfile.open(fname.c_str());
#else
			pfile.open(fname, std::ios::binary);
#endif
			if(!pfile.is_open()) {
				OpenErrorR("Dynamic landscape patch map ",  fname);
#if RSDEBUG
				DEBUGLOG << "Dynamic landscape patch map failed to open: " << fname << std::endl;
#endif
				pfile.clear();
				return -213;
			} else {
#if RSDEBUG
				DEBUGLOG << "Dynamic landscape patch map #" << i << " open to read" << std::endl;
#endif
#if !RSWIN64
				// check BOM for UTF-16
				if(check_bom(fname) == "utf16")
					// apply BOM-sensitive UTF-16 facet
					pfile.imbue(std::locale(pfile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif
				// ASCII header
				pfile >> header;
				if (!pfile.good()) {
#if RSDEBUG
					DEBUGLOG << "ReadDynLandR(): failed to read landscape patch map #" << i << ": " << fname << std::endl;
#endif
					errors = -1113;
					pfile.close();
					pfile.clear();
					return errors;
				}
				if (header != L"ncols" && header != L"NCOLS") errors++;
				pfile >> ncols;

				pfile >> header >> nrows;
				if (header != L"nrows" && header != L"NROWS") errors++;

				pfile >> header >> xllcorner;
				if (header != L"xllcorner" && header != L"XLLCORNER") errors++;

				pfile >> header >> yllcorner;
				if (header != L"yllcorner" && header != L"YLLCORNER") errors++;

				pfile >> header >> cellsize;
				if (header != L"cellsize" && header != L"CELLSIZE") errors++;

				pfile >> header >> pchnodata;
				if (header != L"NODATA_value" && header != L"NODATA_VALUE") errors++;

				if (errors > 0)  {
					FormatErrorR(fname,errors);
#if RSDEBUG
					DEBUGLOG << "ReadDynLandR(): failed to read Raster header of landscape patch map #" << i << ": " << fname << std::endl;
#endif
					pfile.close();
					pfile.clear();
				} else {
					// check resolution match
					if (cellsize == resolution) {
						// check that extent matches landscape extent
						if(ncols == landraster.ncols && nrows == landraster.nrows) {
							// check origins match
							if((int)xllcorner == (int)landraster.xllcorner &&
							   (int)yllcorner == (int)landraster.yllcorner) {
								Rcpp::Rcout << "Dynamic landscape patch map #" << i << ", headers OK: " << fname << endl;
							} else {
								Rcpp::Rcout << "*** Origin co-ordinates of " << fname << msghdrs1 << endl;
								errors++;
							}
						} else {
							Rcpp::Rcout << "*** Extent of " << fname << msghdrs1 << endl;
							errors++;
						}
					} else {
						Rcpp::Rcout << "*** Resolution of " << fname << msghdrs1 << endl;
						errors++;
					}
					if (errors > 0)  {
						FormatErrorR(fname,errors);
#if RSDEBUG
						DEBUGLOG << "ReadDynLandR(): Errors in raster header of landscape patch map #" << i << ": " << fname << std::endl;
#endif
						pfile.close();
						pfile.clear();
					}
				} // end of reading ASCII header
			}
		} // end of if(patchmodel)
		else {
			pfile.clear();
		}

		// Do the same for corresponding cost map, if applicable
		if (costs) {
			// Cost change file
			fname = indir + costmaps(i);

			// open file
#if RSWIN64
			cfile.open(fname.c_str());
#else
			cfile.open(fname, std::ios::binary);
#endif
			if(!cfile.is_open()) {
				OpenErrorR("Dynamic SMS cost map ",  fname);
#if RSDEBUG
				DEBUGLOG << "Dynamic SMS cost map failed to open: " << fname << std::endl;
#endif
				cfile.clear();
				return -214;
			} else {
#if RSDEBUG
				DEBUGLOG << "Dynamic SMS cost map #" << i << " open to read" << std::endl;
#endif
#if !RSWIN64
				// check BOM for UTF-16
				if(check_bom(fname) == "utf16")
					// apply BOM-sensitive UTF-16 facet
					cfile.imbue(std::locale(cfile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif
				// ASCII header
				cfile >> header;
				if (!cfile.good()) {
#if RSDEBUG
					DEBUGLOG << "ReadDynLandR(): failed to read SMS cost map #" << i << ": " << fname << std::endl;
#endif
					errors = -1113;
					cfile.close();
					cfile.clear();
					return errors;
				}
				if (header != L"ncols" && header != L"NCOLS") errors++;
				cfile >> ncols;

				cfile >> header >> nrows;
				if (header != L"nrows" && header != L"NROWS") errors++;

				cfile >> header >> xllcorner;
				if (header != L"xllcorner" && header != L"XLLCORNER") errors++;

				cfile >> header >> yllcorner;
				if (header != L"yllcorner" && header != L"YLLCORNER") errors++;

				cfile >> header >> cellsize;
				if (header != L"cellsize" && header != L"CELLSIZE") errors++;

				cfile >> header >> pchnodata;
				if (header != L"NODATA_value" && header != L"NODATA_VALUE") errors++;

				if (errors > 0)  {
					FormatErrorR(fname,errors);
#if RSDEBUG
					DEBUGLOG << "ReadDynLandR(): failed to read Raster header of SMS cost map #" << i << ": " << fname << std::endl;
#endif
					cfile.close();
					cfile.clear();
				} else {
					// check resolution match
					if (cellsize == resolution) {
						// check that extent matches landscape extent
						if(ncols == landraster.ncols && nrows == landraster.nrows) {
							// check origins match
							if((int)xllcorner == (int)landraster.xllcorner &&
							   (int)yllcorner == (int)landraster.yllcorner) {
								Rcpp::Rcout << "Dynamic SMS cost map #" << i << ", headers OK: " << fname << endl;
							} else {
								Rcpp::Rcout << "*** Origin co-ordinates of " << fname << msghdrs1 << endl;
								errors++;
							}
						} else {
							Rcpp::Rcout << "*** Extent of " << fname << msghdrs1 << endl;
							errors++;
						}
					} else {
						Rcpp::Rcout << "*** Resolution of " << fname << msghdrs1 << endl;
						errors++;
					}
					if (errors > 0)  {
						FormatErrorR(fname,errors);
#if RSDEBUG
						DEBUGLOG << "ReadDynLandR(): Errors in raster header of SMS cost map #" << i << ": " << fname << std::endl;
#endif
						cfile.close();
						cfile.clear();
					}
				} // end of reading ASCII header
			}
		} // end of if(costs)
		else {
			cfile.clear();
		}

		// Now read raster data of Habitat and, if applicable, Patch and/or Cost maps:
		int imported = 0;
		if (errors == 0)  {
			imported = pLandscape->readLandChange(i-1, costs, hfile, pfile, cfile, habnodata, pchnodata, costnodata);
			if (imported != 0) {
				if(hfile.is_open()) hfile.close();
				hfile.clear();
				if(patchmodel) {
					if(pfile.is_open()) pfile.close();
					pfile.clear();
				}
				if (costs) {
					if(cfile.is_open()) cfile.close();
					cfile.clear();
				}
				return imported;
			}
			if (patchmodel) {
				pLandscape->recordPatchChanges(i);
			}
			if (costs) {
				pLandscape->recordCostChanges(i);
			}
		}

		// Close files
		if(hfile.is_open()) hfile.close();
		hfile.clear();
		if(patchmodel) {
			if(pfile.is_open()) pfile.close();
			pfile.clear();
		}
		if (costs) {
			if(cfile.is_open()) cfile.close();
			cfile.clear();
		}
	} // end of loop over landscape changes i

	if(patchmodel) {
		// record changes back to original landscape for multiple replicates
		pLandscape->recordPatchChanges(0);
		pLandscape->deletePatchChgMatrix();
	}
	if (costs) {
		pLandscape->recordCostChanges(0);
		pLandscape->deleteCostsChgMatrix();
	}
#if RSDEBUG
	DEBUGLOG << "ReadDynLandR(): finished" << endl;
#endif
	return 0;
}


//---------------------------------------------------------------------------

int ReadParametersR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{
	Rcpp::S4 ParamParamsR("SimulationParams");
	ParamParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("simul"));

	Rcpp::S4 DemogParamsR("DemogParams");
	DemogParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("demog"));

	Rcpp::S4 LandParamsR("LandParams");
	LandParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("land"));

	int error = 0;
	landParams paramsLand = pLandscape->getLandParams();
	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogr();
	simParams sim = paramsSim->getSim();

	sim.simulation = Rcpp::as<int>(ParamParamsR.slot("Simulation"));
	sim.reps = Rcpp::as<int>(ParamParamsR.slot("Replicates"));
	sim.years = Rcpp::as<int>(ParamParamsR.slot("Years"));
#if RSDEBUG
	DEBUGLOG << "ReadParametersR(): paramsSim = " << paramsSim << endl;
	DEBUGLOG << "ReadParametersR(): simulation = " << sim.simulation << " reps = " << sim.reps
	         << " years = " << sim.years << endl;
#endif

	int iiii, gradType, shift_begin, shift_stop;
	float grad_inc, opt_y, f, optEXT, shift_rate;
	bool shifting;

	// Landscape boundary and any 'no-data' areas are absorbing?:
	sim.absorbing = Rcpp::as<bool>(ParamParamsR.slot("Absorbing"));

	// Environmental gradient: 0 = none, 1  = carrying capacity (or 1/b), 2 = growth rate (or fecundity), 3 = local
	// extinction probability. Must be 0 for patch-based models. NOTE: change in code numbers from v1.0
	gradType = Rcpp::as<int>(ParamParamsR.slot("Gradient"));
	// following not required for Gradient = 0
	grad_inc = Rcpp::as<float>(ParamParamsR.slot("GradSteep"));
	opt_y = Rcpp::as<float>(ParamParamsR.slot("Optimum"));
	f = Rcpp::as<float>(ParamParamsR.slot("f"));
	optEXT = Rcpp::as<float>(ParamParamsR.slot("ExtinctOptim"));
	paramsGrad->setGradient(gradType, grad_inc, opt_y, f, optEXT);

	shifting = Rcpp::as<bool>(ParamParamsR.slot("Shifting"));
	shift_rate = Rcpp::as<float>(ParamParamsR.slot("ShiftRate"));
	shift_begin = Rcpp::as<int>(ParamParamsR.slot("ShiftStart"));
	shift_stop = Rcpp::as<int>(ParamParamsR.slot("ShiftEnd"));
	if(shifting)
		paramsGrad->setShifting(shift_rate, shift_begin, shift_stop);
	else
		paramsGrad->noShifting();

	// Environmental stochasticity: 0 = none, 1 = global, 2 = local:
	iiii = Rcpp::as<int>(ParamParamsR.slot("EnvStoch"));
	if(iiii == 0)
		env.stoch = false;
	else {
		env.stoch = true;
		if(iiii == 2)
			env.local = true;
		else
			env.local = false;
	}
	// For a patch-based model, EnvStoch = 2 is NOT allowed:
	if(paramsLand.patchModel && env.local)
		error = 101;

	// Environmental stochasticity type: 0 = in growth rate, 1 = in carrying capacity
	// Not required for EnvStoch = 0; stochasticity in carrying capacity is allowed for an artificial landscape only
	env.inK = Rcpp::as<bool>(ParamParamsR.slot("EnvStochType"));

	// as from v1.1, there is just one pair of min & max values,
	// which are attributes of the species
	// ULTIMATELY, THE PARAMETER FILE SHOULD HAVE ONLY TWO COLUMNS ...
	env.ac = Rcpp::as<float>(ParamParamsR.slot("ac")); // Temporal autocorrelation coefficient
	env.std = Rcpp::as<float>(ParamParamsR.slot(
	                              "std")); // Amplitude of stochastic fluctuations: standard deviation of a normal distribution (having mean 0)
	float minR, maxR, minK, maxK;
	minR = Rcpp::as<float>(ParamParamsR.slot("minR")); // not required for EnvStoch = 0 or EnvStochType = 1
	maxR = Rcpp::as<float>(ParamParamsR.slot("maxR")); // -"-
	minK = Rcpp::as<float>(ParamParamsR.slot("minK")); // not required for EnvStoch = 0 or EnvStochType = 0
	maxK = Rcpp::as<float>(ParamParamsR.slot("maxK")); // -"-
	if(env.inK) {
		float minKK, maxKK;
		minKK = minK * (((float)(pow(paramsLand.resol, double(2)))) / 10000.0);
		maxKK = maxK * (((float)(pow(paramsLand.resol, double(2)))) / 10000.0);
		pSpecies->setMinMax(minKK, maxKK);
	} else
		pSpecies->setMinMax(minR, maxR);

	env.localExt = Rcpp::as<bool>(ParamParamsR.slot("LocalExt"));
	// For a patch-based model no LocalExt allowed:
	if(paramsLand.patchModel && env.localExt)
		error = 102;
	env.locExtProb = Rcpp::as<float>(ParamParamsR.slot("LocalExtProb"));
	paramsStoch->setStoch(env);

	dem.propMales = Rcpp::as<float>(DemogParamsR.slot("PropMales"));
	dem.harem = Rcpp::as<float>(DemogParamsR.slot("Harem"));
	dem.bc = Rcpp::as<float>(DemogParamsR.slot("bc"));
	dem.lambda = Rcpp::as<float>(DemogParamsR.slot("Rmax"));

	pSpecies->setDemogr(dem);

	float k;
	if(landtype == 9) { // artificial landscape
		// only one value of K is read, but the first 'habitat' is the matrix where K = 0
		pSpecies->createHabK(2);
		k = Rcpp::as<float>(LandParamsR.slot("K_or_DensDep"));
		k *= ((double)(pow(paramsLand.resol, double(2)))) / 10000.0;
		pSpecies->setHabK(0,0);
		pSpecies->setHabK(1,k);
	} else {
		pSpecies->createHabK(paramsLand.nHabMax);
		Rcpp::NumericVector k_vec;
		k_vec = Rcpp::as<Rcpp::NumericVector>(LandParamsR.slot("K_or_DensDep"));
		for (int i = 0; i < paramsLand.nHabMax; i++) {
			k = k_vec[i] * ((double)(pow(paramsLand.resol, double(2)))) / 10000.0;
			pSpecies->setHabK(i,k);
		}
	}

#if RSDEBUG
	DEBUGLOG << "ReadParametersR(): dem.lambda = " << dem.lambda
	         << " habK[0] = " << pSpecies->getHabK(0)
	         << " nHabMax = " << paramsLand.nHabMax << endl;
#endif

	// Output start years
	sim.outStartPop = Rcpp::as<int>(ParamParamsR.slot("OutStartPop"));
	sim.outStartInd = Rcpp::as<int>(ParamParamsR.slot("OutStartInd"));
	sim.outStartGenetic = Rcpp::as<int>(ParamParamsR.slot("OutStartGenetic"));
	sim.outStartTraitCell = Rcpp::as<int>(ParamParamsR.slot("OutStartTraitCell"));
	sim.outStartTraitRow = Rcpp::as<int>(ParamParamsR.slot("OutStartTraitRow"));
	sim.outStartConn = Rcpp::as<int>(ParamParamsR.slot("OutStartConn"));
	sim.outStartPaths = Rcpp::as<int>(ParamParamsR.slot("OutStartPaths"));
	// Output intervals
	sim.outIntRange = Rcpp::as<int>(ParamParamsR.slot("OutIntRange"));
	sim.outIntOcc = Rcpp::as<int>(ParamParamsR.slot("OutIntOcc"));
	sim.outIntPop = Rcpp::as<int>(ParamParamsR.slot("OutIntPop"));
	sim.outIntInd = Rcpp::as<int>(ParamParamsR.slot("OutIntInd"));
	sim.outIntGenetic = Rcpp::as<int>(ParamParamsR.slot("OutIntGenetic"));

	if(sim.outIntRange > 0)
		sim.outRange = true;
	else
		sim.outRange = false;
	if(sim.outIntOcc > 0)
		sim.outOccup = true;
	else
		sim.outOccup = false;
	if(sim.outIntPop > 0)
		sim.outPop = true;
	else
		sim.outPop = false;
	if(sim.outIntInd > 0)
		sim.outInds = true;
	else
		sim.outInds = false;
	if(sim.outIntGenetic > 0)
		sim.outGenetics = true;
	else
		sim.outGenetics = false;

	// Output genetics for: 0 = juveniles only, 1 =  all individuals, 2 = adults only:
	sim.outGenType = Rcpp::as<int>(ParamParamsR.slot("OutGenType"));

	// Output genetics as a cross table?
	sim.outGenXtab = Rcpp::as<bool>(ParamParamsR.slot("OutGenCrossTab"));

	sim.outIntTraitCell = Rcpp::as<int>(ParamParamsR.slot("OutIntTraitCell"));
	sim.outIntTraitRow = Rcpp::as<int>(ParamParamsR.slot("OutIntTraitRow"));
	sim.outIntConn = Rcpp::as<int>(ParamParamsR.slot("OutIntConn"));
	sim.outIntPaths = Rcpp::as<int>(ParamParamsR.slot("OutIntPaths"));
	if(sim.outIntTraitCell > 0)
		sim.outTraitsCells = true;
	else
		sim.outTraitsCells = false;
	if(sim.outIntTraitRow > 0)
		sim.outTraitsRows = true;
	else
		sim.outTraitsRows = false;
	if(sim.outIntConn > 0)
		sim.outConnect = true;
	else
		sim.outConnect = false;
	if(sim.outIntPaths > 0)
		sim.outPaths = true;
	else
		sim.outPaths = false;

	if(sim.outOccup && sim.reps < 2)
		error = 103;
	if(paramsLand.patchModel) {
		if(sim.outTraitsRows)
			error = 104;
	} else {
		if(sim.outConnect)
			error = 105;
	}
#if RSDEBUG
	DEBUGLOG << "ReadParametersR(): outRange = " << sim.outRange << " outInt = " << sim.outIntRange << endl;
#endif

	sim.saveMaps = Rcpp::as<bool>(ParamParamsR.slot("SaveMaps"));
	sim.mapInt = Rcpp::as<int>(ParamParamsR.slot("MapsInterval"));
	sim.saveVisits = Rcpp::as<bool>(ParamParamsR.slot("SMSHeatMap"));
	sim.drawLoaded = Rcpp::as<bool>(ParamParamsR.slot("DrawLoadedSp"));
// sim.saveInitMap = false;
	sim.ReturnPopRaster = Rcpp::as<bool>(ParamParamsR.slot("ReturnPopRaster"));
	sim.CreatePopFile = Rcpp::as<bool>(ParamParamsR.slot("CreatePopFile"));

	paramsSim->setSim(sim);

	return error;
}

//---------------------------------------------------------------------------


int ReadStageStructureR(Rcpp::S4 ParMaster)
{
	Rcpp::S4 DemogParamsR("DemogParams");
	DemogParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("demog"));
	Rcpp::S4 StagesParamsR("StagesParams");
	StagesParamsR = Rcpp::as<Rcpp::S4>(DemogParamsR.slot("StageStruct"));

	demogrParams dem = pSpecies->getDemogr();
	stageParams sstruct = pSpecies->getStage();
	int matrixsize, i, j, stg;
	float ss, dd, devCoeff, survCoeff;
	Rcpp::NumericMatrix trmatrix, wtsmatrix;
	Rcpp::IntegerVector minAge;

#if RSDEBUG
	DEBUGLOG << "ReadStageStructureR(): sstruct.nStages = " << sstruct.nStages << endl;
#endif
	// int simulation;
	// simulation = Rcpp::as<int>(StagesParamsR.slot("Simulation")); //Must match simulation numbers in ParamParams
	sstruct.disperseOnLoss = Rcpp::as<bool>(StagesParamsR.slot("PostDestructn"));
	sstruct.probRep = Rcpp::as<float>(StagesParamsR.slot("PRep"));
	sstruct.repInterval = Rcpp::as<int>(StagesParamsR.slot("RepInterval"));
	sstruct.maxAge = Rcpp::as<int>(StagesParamsR.slot("MaxAge"));
	trmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("TransMatrix"));
	// parsing is done on R-level
	minAge = Rcpp::as<Rcpp::IntegerVector>(StagesParamsR.slot("MinAge"));

	// Store Transition matrix:
	if(dem.repType != 2) { // asexual or implicit sexual model
		matrixsize = sstruct.nStages;
		for(i = 0; i < matrixsize; i++) {
			// set minimum ages
			pSpecies->setMinAge(i, 0, minAge(i));
			// set fecundities
			pSpecies->setFec(i, 0, (float)trmatrix(0, i));
			// set survival and development probalities
			ss = (float)trmatrix(i, i); // survival prob
			if((i + 1) != matrixsize) {
				dd = (float)trmatrix(i + 1, i);
			} // development prob
			else {
				dd = 0.0;
			}
			pSpecies->setSurv(i, 0, ss + dd);
			pSpecies->setDev(i, 0, dd / (ss + dd));
		}
	} else {                              // complex sexual model
		matrixsize = sstruct.nStages * 2; // juv (i=0) sex-independent as ending stage -> columns = rows + 1
		stg = 1;
		for(i = 1; i < (matrixsize - 1); i++) { // loop over rows
			// set minAge and fecundities
			if(i % 2) {                         // odd columns  -> males (1)
				pSpecies->setMinAge(stg, 1, minAge(i));
				pSpecies->setFec(stg, 1, (float)trmatrix(0, i + 1));
			} else { // even columns -> females (0)
				pSpecies->setMinAge(stg, 0, minAge(i));
				pSpecies->setFec(stg, 0, (float)trmatrix(0, i + 1));
				stg++;
			}
		}
		stg = 0;
		for(i = 0; i < matrixsize; i++) { // loop over columns
			// survival and development
			if(i != 0)
				ss = (float)trmatrix(i - 1, i); // survival prob
			else
				ss = (float)trmatrix(i, i);
			if((i + 2) != matrixsize)
				dd = (float)trmatrix(i + 1, i); // development prob
			else
				dd = 0.0;
			if(i % 2 == 0) { // even rows -> males (1)
				pSpecies->setSurv(stg, 1, ss + dd);
				pSpecies->setDev(stg, 1, dd / (ss + dd));
			} else { // odd rows -> females (0)
				pSpecies->setSurv(stg, 0, ss + dd);
				pSpecies->setDev(stg, 0, dd / (ss + dd));
				stg++;
			}
		}
	}
#if RSDEBUG
	DEBUGLOG << "Read_Transition Matrix: matrix = " << trmatrix << endl;
#endif

	// Survival schedule: 0 = At reproduction; 1 = Between reproductive events; 2 = Annually
	sstruct.survival = Rcpp::as<int>(StagesParamsR.slot("SurvSched"));

	if(dem.repType != 2)
		matrixsize = sstruct.nStages;
	else
		matrixsize = sstruct.nStages * NSEXES;
	// Fecundity
	sstruct.fecDens = Rcpp::as<bool>(StagesParamsR.slot("FecDensDep")); // Density-dependence in reproduction
	sstruct.fecStageDens =
	    Rcpp::as<bool>(StagesParamsR.slot("FecStageWts")); // stage-specific density dependence of fecundity
	if(sstruct.fecStageDens) {
		wtsmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("FecStageWtsMatrix"));
		pSpecies->createDDwtFec(matrixsize);
		for(i = 0; i < matrixsize; i++) {
			for(j = 0; j < matrixsize; j++) {
				pSpecies->setDDwtFec(i, j, wtsmatrix(i, j));
			}
		}
#if RSDEBUG
		DEBUGLOG << "Read_StageWeights(): completed reading fecundity weights matrix " << endl;
#endif
	}

	// Development
	sstruct.devDens = Rcpp::as<bool>(StagesParamsR.slot("DevDensDep")); // Density-dependence in development
	devCoeff = Rcpp::as<float>(StagesParamsR.slot("DevDensCoeff"));
	sstruct.devStageDens =
	    Rcpp::as<bool>(StagesParamsR.slot("DevStageWts")); // stage-specific density dependence of development
	if(sstruct.devStageDens) {
		wtsmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("DevStageWtsMatrix"));
		pSpecies->createDDwtDev(matrixsize);
		for(i = 0; i < matrixsize; i++) {
			for(j = 0; j < matrixsize; j++) {
				pSpecies->setDDwtDev(i, j, wtsmatrix(i, j));
			}
		}
#if RSDEBUG
		DEBUGLOG << "Read_StageWeights(): completed reading development weights matrix " << endl;
#endif
	}

	// Survival
	sstruct.survDens = Rcpp::as<bool>(StagesParamsR.slot("SurvDensDep")); // Density-dependence in survival
	survCoeff = Rcpp::as<float>(StagesParamsR.slot("SurvDensCoeff"));
	sstruct.survStageDens =
	    Rcpp::as<bool>(StagesParamsR.slot("SurvStageWts")); // stage-specific density dependence of survival
	if(sstruct.survStageDens) {
		wtsmatrix = Rcpp::as<Rcpp::NumericMatrix>(StagesParamsR.slot("SurvStageWtsMatrix"));
		pSpecies->createDDwtSurv(matrixsize);
		for(i = 0; i < matrixsize; i++) {
			for(j = 0; j < matrixsize; j++) {
				pSpecies->setDDwtSurv(i, j, wtsmatrix(i, j));
			}
		}
#if RSDEBUG
		DEBUGLOG << "Read_StageWeights(): completed reading survival weights matrix " << endl;
#endif
	}

	pSpecies->setStage(sstruct);
	if(sstruct.devDens || sstruct.survDens) {
		pSpecies->setDensDep(devCoeff, survCoeff);
	}

	return 0;
}

//---------------------------------------------------------------------------

int ReadEmigrationR(Rcpp::S4 ParMaster)
{

	Rcpp::S4 DispParamsR("DispersalParams");
	DispParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("dispersal"));
	Rcpp::S4 EmigParamsR("EmigrationParams");
	EmigParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Emigration"));

	int error = 0;
	int Nlines, emigstage, stage, sex, offset;
	Rcpp::NumericMatrix EmigMatrix;
	Rcpp::NumericVector EmigScalesVec;
	demogrParams dem = pSpecies->getDemogr();
	stageParams sstruct = pSpecies->getStage();
	emigRules emig = pSpecies->getEmig();
	emigTraits etraits;
	emigParams eparams;
	emigScales scale = pSpecies->getEmigScales();

	// int simulation;
	// simulation = Rcpp::as<int>(EmigParamsR.slot("Simulation")); // REMOVED in R-interface // Must match simulation
	// numbers in ParamParams
	emig.densDep = Rcpp::as<bool>(EmigParamsR.slot("DensDep"));
	pSpecies->setFullKernel(Rcpp::as<bool>(EmigParamsR.slot("UseFullKern"))); // Emigration rate derived from kernel. Only for kernel-based transfer and DensDep = 0
	emig.stgDep = Rcpp::as<bool>(EmigParamsR.slot("StageDep")); // Stage-dependent emigration. Must be 0 if IndVar is 1
	emig.sexDep = Rcpp::as<bool>(EmigParamsR.slot("SexDep"));   // Sex-dependent emigration.
	emig.indVar =
	    Rcpp::as<bool>(EmigParamsR.slot("IndVar")); // Inter-individual variability. Must be 0 if StageDep is 1
	if(stagestruct && emig.indVar){
		emigstage = Rcpp::as<int>(EmigParamsR.slot("EmigStage")); // Stage which emigrates. Required for stage-strucutred population having IndVar = 1
		if(emigstage >= 0 && emigstage < sstruct.nStages)
			emig.emigStage = emigstage;
		else
			emig.emigStage = 0;
	} else {
		emig.emigStage = 0;
	}

	pSpecies->setEmig(emig);

	if(!dem.repType && emig.sexDep)
		error = 301;
	if(!dem.stageStruct && emig.stgDep)
		error = 303;

	// no.of lines according to known stage- and sex-dependency and corresponding column offset
	if(emig.stgDep) {
		if(emig.sexDep) {
			Nlines = sstruct.nStages * sexesDisp;
			offset = 2;
		} else {
			Nlines = sstruct.nStages;
			offset = 1;
		}
	} else {
		if(emig.sexDep) {
			Nlines = sexesDisp;
			offset = 1;
		} else {
			Nlines = 1;
			offset = 0;
		}
	}

#if RSDEBUG
	DEBUGLOG << "ReadEmigrationR(): Nlines = " << Nlines << " emig.densDep = " << emig.densDep
	         << " emig.indVar = " << emig.indVar << " sexesDisp = " << sexesDisp << endl;
#endif

	EmigMatrix = Rcpp::as<Rcpp::NumericMatrix>(EmigParamsR.slot("EmigProb"));

	for(int line = 0; line < Nlines; line++) {

		if(emig.stgDep) {
			if(emig.sexDep) {
				stage = (int)EmigMatrix(line, 0);
				sex = (int)EmigMatrix(line, 1);
			} else {
				stage = (int)EmigMatrix(line, 0);
				sex = 0;
			}
		} else {
			if(emig.sexDep) {
				stage = 0;
				sex = (int)EmigMatrix(line, 0);
			} else {
				stage = 0;
				sex = 0;
			}
		}

		if(emig.densDep) {
			if(emig.indVar) {
				eparams.d0Mean = (float)EmigMatrix(line, offset + 0);
				eparams.d0SD = (float)EmigMatrix(line, offset + 1);
				eparams.alphaMean = (float)EmigMatrix(line, offset + 2);
				eparams.alphaSD = (float)EmigMatrix(line, offset + 3);
				eparams.betaMean = (float)EmigMatrix(line, offset + 4);
				eparams.betaSD = (float)EmigMatrix(line, offset + 5);

				pSpecies->setEmigParams(stage, sex, eparams);
			} else {
				etraits.d0 = (float)EmigMatrix(line, offset + 0);
				etraits.alpha = (float)EmigMatrix(line, offset + 1);
				etraits.beta = (float)EmigMatrix(line, offset + 2);

				pSpecies->setEmigTraits(stage, sex, etraits);
			}
		} else { // !emig.densDep
			if(emig.indVar) {
				eparams.d0Mean = (float)EmigMatrix(line, offset + 0);
				eparams.d0SD = (float)EmigMatrix(line, offset + 1);
				eparams.alphaMean = eparams.betaMean = 0.0;
				eparams.alphaSD = eparams.betaSD = 0.00000001;

				pSpecies->setEmigParams(stage, sex, eparams);
			} else {
				etraits.d0 = (float)EmigMatrix(line, offset + 0);
				etraits.alpha = etraits.beta = 0.0;

				pSpecies->setEmigTraits(stage, sex, etraits);
			}
		}
	} // end of Nlines for loop

	EmigScalesVec = Rcpp::as<Rcpp::NumericVector>(EmigParamsR.slot("TraitScaleFactor"));

	if(emig.indVar) {
		if(emig.densDep) {
			scale.d0Scale = (float)EmigScalesVec(0);
			scale.alphaScale = (float)EmigScalesVec(1);
			scale.betaScale = (float)EmigScalesVec(2);
		} else {
			scale.d0Scale = (float)EmigScalesVec(0);
			scale.alphaScale = scale.betaScale = 0.00000001;
		}
		pSpecies->setEmigScales(scale);
	}

	if(emig.indVar) anyIndVar = true;

	return error;
}

//---------------------------------------------------------------------------

int ReadTransferR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{
	Rcpp::S4 DispParamsR("DispersalParams");
	DispParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("dispersal"));

	int Nlines, offset, stage, sex;
	int error = 0;

	landParams paramsLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogr();
	stageParams sstruct = pSpecies->getStage();
	trfrRules trfr = pSpecies->getTrfr();

#if RSDEBUG
	DEBUGLOG << "ReadTransferR(): TransferType=" << trfr.moveModel
	         << " paramsLand.generated=" << paramsLand.generated
	         << " paramsLand.rasterType=" << paramsLand.rasterType
	         << " trfr.moveModel=" << trfr.moveModel
	         << " trfr.twinKern=" << trfr.twinKern
	         << endl;
#endif

	// Create Costs vector of species
	if(trfr.moveModel) {
#if RSDEBUG
		DEBUGLOG << "ReadTransferR(): creating cost/mortality matrix, dimension=";
		if(paramsLand.generated)
			DEBUGLOG << paramsLand.nHab;
		else
			DEBUGLOG << paramsLand.nHabMax;
		DEBUGLOG << endl;
#endif
		if(paramsLand.generated) {
			pSpecies->createHabCostMort(paramsLand.nHab);
		} else {
			pSpecies->createHabCostMort(paramsLand.nHabMax);
		}
	}

	int TransferType; // new local variable to replace former global variable
	if(trfr.moveModel)
		TransferType = trfr.moveType;
	else
		TransferType = 0;

	trfrKernTraits k;
	trfrMovtTraits movt;
	trfrKernParams kparams;
	trfrMortParams mort;
	trfrScales scale;
	string CostsFile;
	trfrSMSParams smsparams;
	trfrCRWParams mparams;
	Rcpp::NumericVector ReadVec;

	switch(TransferType) {
	case 0: { // dispersal kernel

		Rcpp::S4 TransParamsR("DispersalKernel");
		TransParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Transfer"));

		Rcpp::NumericMatrix DispMatrix;
		Rcpp::NumericVector DispScalesVec;

		// simulation = Rcpp::as<int>(TransParamsR.slot("Simulation")); // REMOVED in R-interface //Must match
		// simulation numbers in ParamParams
		trfr.stgDep = Rcpp::as<bool>(TransParamsR.slot("StageDep")); // Stage-dependent transfer. Must be 0 if IndVar is 1
		trfr.sexDep = Rcpp::as<bool>(TransParamsR.slot("SexDep")); // Sex-dependent transfer.
		trfr.twinKern = Rcpp::as<bool>(TransParamsR.slot("DoubleKernel")); // 0 = negative exponential; 1 = double negative exponential
		trfr.distMort = Rcpp::as<bool>(TransParamsR.slot("DistMort")); // Distance-dependent mortality
		trfr.indVar = Rcpp::as<bool>(TransParamsR.slot("IndVar"));

		// some error checks
		if(dem.repType == 0) {
			if(trfr.sexDep) {
				error = 401;
			}
		}
		if(dem.stageStruct) {
		} // if (trfr.indVar) error = 402;
		else {
			if(trfr.stgDep) {
				error = 403;
			}
		}

		// set no. of lines according to known stage- and sex-dependency and corresponding column offset
		if(trfr.stgDep) {
			if(trfr.sexDep) {
				Nlines = sstruct.nStages * sexesDisp;
				offset = 2;
			} else {
				Nlines = sstruct.nStages;
				offset = 1;
			}
		} else {
			if(trfr.sexDep) {
				Nlines = sexesDisp;
				offset = 1;
			} else {
				Nlines = 1;
				offset = 0;
			}
		}

		DispMatrix = Rcpp::as<Rcpp::NumericMatrix>(TransParamsR.slot("Distances"));

		for(int line = 0; line < Nlines; line++) {

			if(trfr.stgDep) {
				if(trfr.sexDep) {
					stage = (int)DispMatrix(line, 0);
					sex = (int)DispMatrix(line, 1);
				} else {
					stage = (int)DispMatrix(line, 0);
					sex = 0;
				}
			} else {
				if(trfr.sexDep) {
					stage = 0;
					sex = (int)DispMatrix(line, 0);
				} else {
					stage = 0;
					sex = 0;
				}
			}

			if(trfr.twinKern) {
				if(trfr.indVar) {
					kparams.dist1Mean = (float)DispMatrix(line, offset + 0);
					kparams.dist1SD = (float)DispMatrix(line, offset + 1);
					kparams.dist2Mean = (float)DispMatrix(line, offset + 2);
					kparams.dist2SD = (float)DispMatrix(line, offset + 3);
					kparams.PKern1Mean = (float)DispMatrix(line, offset + 4);
					kparams.PKern1SD = (float)DispMatrix(line, offset + 5);
					// MAXDist = kparams.maxDist1;
					pSpecies->setKernParams(stage, sex, kparams, paramsLand.resol);
				} else { // const kernel parameters
					k.meanDist1 = (float)DispMatrix(line, offset + 0);
					k.meanDist2 = (float)DispMatrix(line, offset + 1);
					k.probKern1 = (float)DispMatrix(line, offset + 2);

					pSpecies->setKernTraits(stage, sex, k, paramsLand.resol);
				}
			} else { // single kernel
				if(trfr.indVar) {
					kparams.dist1Mean = (float)DispMatrix(line, offset + 0);
					kparams.dist1SD = (float)DispMatrix(line, offset + 1);
					kparams.dist2Mean = kparams.dist1Mean;
					kparams.dist2SD = kparams.dist1SD;
					kparams.PKern1Mean = 0.999;
					kparams.PKern1SD = 0.001;

					pSpecies->setKernParams(stage, sex, kparams, paramsLand.resol);
				} else { // const kernel parameters
					k.meanDist1 = (float)DispMatrix(line, offset + 0);
					k.meanDist2 = k.meanDist1;
					k.probKern1 = 1.0;

					pSpecies->setKernTraits(stage, sex, k, paramsLand.resol);
				}
			}
		} // end of Nlines for-loop

		// Mutation scales
		scale = pSpecies->getTrfrScales();
		DispScalesVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("TraitScaleFactor"));
		if(trfr.indVar) {
			//	if (!trfr.indVar) error = 411;
			//	if (dem.stageStruct) error = 412;
			if(trfr.twinKern) {
				scale.dist1Scale = (float)DispScalesVec(0);
				scale.dist2Scale = (float)DispScalesVec(1);
				scale.PKern1Scale = (float)DispScalesVec(2);
			} else {
				scale.dist1Scale = (float)DispScalesVec(0);
				scale.dist2Scale = scale.dist1Scale;
				scale.PKern1Scale = 0.00000001;
			}
			pSpecies->setTrfrScales(scale);
		}

		// mortality
		mort.fixedMort = Rcpp::as<float>(TransParamsR.slot("MortProb"));
		mort.mortAlpha = Rcpp::as<float>(TransParamsR.slot("Slope"));
		mort.mortBeta  = Rcpp::as<float>(TransParamsR.slot("InflPoint"));
		pSpecies->setMortParams(mort);

		pSpecies->setTrfr(trfr);

	}
	break; // end of dispersal kernel

	case 1: { // SMS

		Rcpp::S4 TransParamsR("StochMove");
		TransParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Transfer"));

		// simulation = Rcpp::as<int>(TransParamsR.slot("Simulation")); // REMOVED in R-interface // Must match simulation numbers in ParamParams
		movt.pr = Rcpp::as<short>(TransParamsR.slot("PR")); // Perceptual range (cells)
		movt.prMethod = Rcpp::as<short>(TransParamsR.slot("PRMethod")); // Perceptual range method: 1 = arithmetic mean; 2 = harmonic mean; 3 = weighted arithmtic mean
		movt.memSize = Rcpp::as<short>(TransParamsR.slot("MemSize"));    // No. of previous steps over which to calculate current direction to apply DP [1-14]
		movt.goalType = Rcpp::as<short>(TransParamsR.slot("GoalType"));   // Goal type: 0 (none) or 2 (dispersal bias)
		trfr.indVar = Rcpp::as<bool>(TransParamsR.slot("IndVar"));
		if(trfr.indVar) {
			smsparams = pSpecies->getSMSParams(0,0);
			scale = pSpecies->getTrfrScales();
		}
		ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("DP")); // Directional persistence, Must be >= 1.0
		if(trfr.indVar) {
			if(ReadVec.size() == 3) {
				smsparams.dpMean = (float)ReadVec[0]; //  Directional persistence initial mean (m); Required for IndVar = 1
				smsparams.dpSD   = (float)ReadVec[1]; //  Directional persistence initial SD (m); Required for IndVar = 1
				scale.dpScale  = (float)ReadVec[2]; //  Directional persistence scaling factor; Required for IndVar = 1
			} else {
				error = 435;
			}
		} else {
			if(ReadVec.size() == 1) {
				movt.dp = (float)ReadVec[0];
			} else {
				error = 436;
			}
		}
		ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("GoalBias")); // Goal bias strength, Must be >= 1.0
		if(trfr.indVar) {
			if(ReadVec.size() == 3) {
				smsparams.gbMean = (float)ReadVec[0]; //  Goal bias strength initial mean (m); Required for IndVar = 1
				smsparams.gbSD   = (float)ReadVec[1]; // Goal bias strength initial SD (m); Required for IndVar = 1
				scale.gbScale  = (float)ReadVec[2]; //  Goal bias strength scaling factor; Required for IndVar = 1
			} else {
				error = 435;
			}
		} else {
			if(ReadVec.size() == 1) {
				movt.gb = (float)ReadVec[0];
			} else {
				error = 436;
			}
		}
		if(movt.goalType == 2) {	// dispersal bias
			ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("AlphaDB")); // Dispersal bias decay rate (> 0)
			if(trfr.indVar) {
				if(ReadVec.size() == 3) {
					smsparams.alphaDBMean = (float)ReadVec[0]; // decay rate initial mean (m); Required for IndVar = 1
					smsparams.alphaDBSD   = (float)ReadVec[1]; // decay rate initial SD (m); Required for IndVar = 1
					scale.alphaDBScale    = (float)ReadVec[2]; // decay rate scaling factor; Required for IndVar = 1
				} else {
					error = 435;
				}
			} else {
				if(ReadVec.size() == 1) {
					movt.alphaDB = (float)ReadVec[0];
				} else {
					error = 436;
				}
			}
			ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("BetaDB")); // Dispersal bias decay inflection point (no. of steps) (> 0)
			if(trfr.indVar) {
				if(ReadVec.size() == 3) {
					smsparams.alphaDBMean = (float)ReadVec[0]; // decay inflection point initial mean (m); Required for IndVar = 1
					smsparams.alphaDBSD   = (float)ReadVec[1]; // decay inflection point initial SD (m); Required for IndVar = 1
					scale.alphaDBScale    = (float)ReadVec[2]; // decay inflection point scaling factor; Required for IndVar = 1
				} else {
					error = 435;
				}
			} else {
				if(ReadVec.size() == 1) {
					movt.betaDB = (float)ReadVec[0];
				} else {
					error = 436;
				}
			}
		}
#if RSDEBUG
		DEBUGLOG << "ReadTransferR(): SMS" << endl
		         << " indVar=" << trfr.indVar << " PR=" << movt.pr << " PRmethod=" << movt.prMethod << endl;
		DEBUGLOG << "ReadTransferR(): dp=" << movt.dp << " MemSize=" << movt.memSize << " gb=" << movt.gb
		         << " goaltype=" << movt.goalType << endl;
#endif

		//smsparams.dpMean = Rcpp::as<float>(TransParamsR.slot("DPMean"));
		//smsparams.dpSD = Rcpp::as<float>(TransParamsR.slot("DPSD"));
		//smsparams.gbMean = Rcpp::as<float>(TransParamsR.slot("GBMean"));
		//smsparams.gbSD = Rcpp::as<float>(TransParamsR.slot("GBSD"));
		//smsparams.alphaDBMean = Rcpp::as<float>(TransParamsR.slot("AlphaDBMean"));
		//smsparams.alphaDBSD = Rcpp::as<float>(TransParamsR.slot("AlphaDBSD"));
		//smsparams.betaDBMean = Rcpp::as<float>(TransParamsR.slot("BetaDBMean"));
		//smsparams.betaDBSD = Rcpp::as<float>(TransParamsR.slot("BetaDBSD"));

		//scale.dpScale = Rcpp::as<float>(TransParamsR.slot("DPScale"));
		//scale.gbScale = Rcpp::as<float>(TransParamsR.slot("GBScale"));
		//scale.alphaDBScale = Rcpp::as<float>(TransParamsR.slot("AlphaDBScale"));
		//scale.betaDBScale = Rcpp::as<float>(TransParamsR.slot("BetaDBScale"));

		if (trfr.indVar) {
			pSpecies->setSMSParams(0,0,smsparams);
			pSpecies->setTrfrScales(scale);
		}

		movt.straigtenPath = Rcpp::as<bool>(TransParamsR.slot("StraightenPath")); // Straighten path after decision not to settle?

		// Mortality
		Rcpp::NumericVector HabMortVec;
		HabMortVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("StepMort"));
		if(HabMortVec.size() == 1) {
			trfr.habMort = false;                 // Per-step mortality type: 0 = constant
			movt.stepMort = (float)HabMortVec[0]; // Constant per-step mortality probability
		} else {
			trfr.habMort = true; // Per-step mortality type: 1 = habitat-dependent
			movt.stepMort = -9;
			if(paramsLand.generated) {
				// values are for habitat (hab=1) then for matrix (hab=0)
				pSpecies->setHabMort(1, (double)HabMortVec[1]);
				pSpecies->setHabMort(0, (double)HabMortVec[0]);
#if RSDEBUG
				DEBUGLOG << "ReadTransferR(): Generated Landscpae with MortHabitat=" << pSpecies->getHabMort(1)
				         << " MortMatrix=" << pSpecies->getHabMort(0) << endl;
#endif
			}
			if(paramsLand.rasterType == 0) {
#if RSDEBUG
				DEBUGLOG << "ReadTransferR(): nHabMax = " << paramsLand.nHabMax << endl;
#endif
				for(int i = 0; i < paramsLand.nHabMax; i++) {
					pSpecies->setHabMort(i, (double)HabMortVec[i]);
#if RSDEBUG
					DEBUGLOG << "ReadTransferR(): Habitat #" << i << ": mortality = " << pSpecies->getHabMort(i)
					         << endl;
#endif
				}
			}
		}


#if RSDEBUG
		DEBUGLOG << "ReadTransferR(): SMtype=" << trfr.habMort << " SMconst=" << movt.stepMort << endl;
#endif

		// Costs
		trfr.costMap = Rcpp::as<bool>(TransParamsR.slot("CostMap"));

		// read habitat costs for land types
		Rcpp::NumericVector HabCostVec;
		if(!trfr.costMap) {
			HabCostVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("Costs"));
#if RSDEBUG
			DEBUGLOG << "ReadTransferR(): Read Habitat Cost vector of length " << HabCostVec.size() << endl;
#endif
		}

		if(!paramsLand.generated) {          // real landscape
			if(paramsLand.rasterType == 0) { // habitat codes
				if(!trfr.costMap) {          // habitat costs
					for(int i = 0; i < paramsLand.nHabMax; i++) {
						pSpecies->setHabCost(i, (int)HabCostVec[i]);
#if RSDEBUG
						DEBUGLOG << "ReadTransferR(): Habitat #" << i << ": cost = " << pSpecies->getHabCost(i) << endl;
#endif
					}
				}
			} else { // habitat quality
				// should have trfr.costMap = 1
			}
		} else {               // artificial landscape
			if(trfr.costMap) { // should not occur
				// should have trfr.costMap = 0
			} else { // habitat costs
				// costs are for habitat (hab=1) then for matrix (hab=0)
				pSpecies->setHabCost(1, (int)HabCostVec[1]);
				pSpecies->setHabCost(0, (int)HabCostVec[0]);
			}
		}
		pSpecies->setTrfr(trfr);
		pSpecies->setMovtTraits(movt);
	}
	break; // end of SMS

	case 2: { // CRW

		Rcpp::S4 TransParamsR("CorrRW");
		TransParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Transfer"));

		// simulation = Rcpp::as<int>(TransParamsR.slot("Simulation"));// REMOVED in R-interface  //Must match
		// simulation numbers in ParamParams
		trfr.indVar = Rcpp::as<bool>(TransParamsR.slot("IndVar"));

		ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("StepLength")); // Step length params
		if(trfr.indVar) {
			if(ReadVec.size() == 3) {
				mparams.stepLgthMean = (float)ReadVec[0]; // Step length initial mean (m); Required for IndVar = 1
				mparams.stepLgthSD = (float)ReadVec[1];   // Step length initial SD (m); Required for IndVar = 1
				scale.stepLScale = (float)ReadVec[2];     // Step length scaling factor; Required for IndVar = 1
			} else {
				error = 435;
			}
		} else {
			if(ReadVec.size() == 1) {
				movt.stepLength = (float)ReadVec[0]; // Step length (m); Required for IndVar = 0; must be > 0
			} else {
				error = 436;
			}
		}

		ReadVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("Rho")); // Step correlation coefficient params
		if(trfr.indVar) {
			if(ReadVec.size() == 3) {
				mparams.rhoMean = (float)ReadVec[0]; // Step correlation coefficient initial mean (m); Required for IndVar = 1
				mparams.rhoSD   = (float)ReadVec[1]; // Step correlation coefficient initial SD (m); Required for IndVar = 1
				scale.rhoScale  = (float)ReadVec[2]; // Step correlation coefficient scaling factor; Required for IndVar = 1
			} else {
				error = 435;
			}
		} else {
			if(ReadVec.size() == 1) {
				movt.rho =
				    (float)ReadVec[0]; // Step correlation coefficient; Required for IndVar = 0; must be > 0.0 and < 1.0
			} else {
				error = 436;
			}
		}

		movt.straigtenPath = Rcpp::as<bool>(TransParamsR.slot("StraightenPath")); // Straighten path after decision not to settle?
		pSpecies->setTrfrScales(scale);

#if RSDEBUG
		DEBUGLOG << "ReadTransferR():"
		         << " paramsLand.rasterType=" << paramsLand.rasterType << " trfr.indVar=" << trfr.indVar
		         << " move.stepLength=" << movt.stepLength << " move.rho=" << movt.rho
		         << " mparams.stepLgthMean=" << mparams.stepLgthMean << " mparams.rhoMean=" << mparams.rhoMean
		         << " move.straigtenPath=" << movt.straigtenPath << endl;
#endif

		// Mortality
		Rcpp::NumericVector HabMortVec;
		HabMortVec = Rcpp::as<Rcpp::NumericVector>(TransParamsR.slot("StepMort"));

		if(HabMortVec.size() == 1) {
			trfr.habMort = 0;                     // Per-step mortality type: 0 = constant
			movt.stepMort = (float)HabMortVec[0]; // Constant per-step mortality probability
		} else {
			trfr.habMort = 1; // Per-step mortality type: 1 = habitat-dependent
			movt.stepMort = -9;
			if(paramsLand.generated) {
				// values are for habitat (hab=1) then for matrix (hab=0)
				pSpecies->setHabMort(1, (double)HabMortVec[1]);
				pSpecies->setHabMort(0, (double)HabMortVec[0]);
#if RSDEBUG
				DEBUGLOG << "ReadTransferR(): Generated Landscpae with MortHabitat=" << pSpecies->getHabMort(1)
				         << " MortMatrix=" << pSpecies->getHabMort(0) << endl;
#endif
			}
			if(paramsLand.rasterType == 0) {
#if RSDEBUG
				DEBUGLOG << "ReadTransferR(): nHabMax = " << paramsLand.nHabMax << endl;
#endif
				for(int i = 0; i < paramsLand.nHabMax; i++) {
					pSpecies->setHabMort(i, (double)HabMortVec[i]);
#if RSDEBUG
					DEBUGLOG << "ReadTransferR(): Habitat #" << i << ": mortality = " << pSpecies->getHabMort(i)
					         << endl;
#endif
				}
			}
		}
		if(trfr.habMort && paramsLand.rasterType != 0)
			error = 434; // habitat percentage landscape cant have habitat-dependent mortality

#if RSDEBUG
		DEBUGLOG << "ReadTransferR(): SMtype=" << trfr.habMort << " SMconst=" << movt.stepMort << endl;
#endif

		pSpecies->setTrfr(trfr);
		pSpecies->setMovtTraits(movt);
		pSpecies->setCRWParams(0, 0, mparams);

	}
	break; // end of CRW

	default:
		error = 440;
	} // end of switch (TransferType)

	if(trfr.indVar) anyIndVar = true;

	return error;
}

//---------------------------------------------------------------------------
// NOTE that stage- and sex-dependent settlement parameters are set for
// ALL stage/sex combinations, even if the species has stage- and/or
// sex-independent settlement rules
int ReadSettlementR(Rcpp::S4 ParMaster)
{

	Rcpp::S4 DispParamsR("DispersalParams");
	DispParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("dispersal"));
	Rcpp::S4 SettleParamsR("SettlementParams");
	SettleParamsR = Rcpp::as<Rcpp::S4>(DispParamsR.slot("Settlement"));

	int Nlines, stage, sex, offset, sexSettle, settType;
	int error = 0;
	bool findmate, densdep;

	demogrParams dem = pSpecies->getDemogr();
	stageParams sstruct = pSpecies->getStage();
	trfrRules trfr = pSpecies->getTrfr();
	settleType sett = pSpecies->getSettle();
	settleRules srules;
	settleSteps ssteps = pSpecies->getSteps(0,0);
	settleTraits settleDD = pSpecies->getSettTraits(0,0);
	settParams sparams = pSpecies->getSettParams(0,0);

	// int simulation;
	// simulation = Rcpp::as<int>(SettleParamsR.slot("Simulation")); // REMOVED in R-interface  // Must match simulation
	// numbers in ParamParams
	sett.stgDep = Rcpp::as<bool>(SettleParamsR.slot("StageDep")); // Stage-dependent settlement.
	sett.sexDep = Rcpp::as<bool>(SettleParamsR.slot("SexDep"));   // Sex-dependent settlement.
	if (transfer == 0) {
		// dispersal kernel                                         // dispersal kernel
		sett.indVar = false;
		if(dem.repType == 0) {
			if(sett.sexDep)
				error = 501; // sex-dependent settlement is not possible with asexual models
		}
		if(!dem.stageStruct) {
			if(sett.stgDep)
				error = 502; // stage-dependent settlement is not possible without stage structure
		}
	} else { // movement process
		sett.indVar = Rcpp::as<bool>(SettleParamsR.slot("IndVar"));
		densdep = Rcpp::as<bool>(SettleParamsR.slot("DensDep"));
		if(dem.repType == 0) {
			if(sett.sexDep)
				error = 508; // sex-dependent settlement is not possible with asexual models
		}
		if(!dem.stageStruct) {
			if(sett.stgDep)
				error = 509; // stage-dependent settlement is not possible without stage structure
		}
	}
	pSpecies->setSettle(sett);

	// no.of lines according to known stage- and sex-dependency and corresponding column offset
	if(sett.stgDep) {
		if(sett.sexDep) {
			Nlines = sstruct.nStages * sexesDisp;
			offset = 2;
		} else {
			Nlines = sstruct.nStages;
			offset = 1;
		}
	} else {
		if(sett.sexDep) {
			Nlines = sexesDisp;
			offset = 1;
		} else {
			Nlines = 1;
			offset = 0;
		}
	}
#if RSDEBUG
	DEBUGLOG << "ReadSettlementR(): sett.stgDep = " << sett.stgDep << ", sett.sexDep = " << sett.sexDep
	         << ", sett.indVar = " << sett.indVar << endl;
#endif

	Rcpp::NumericMatrix SettleCondMatrix = Rcpp::as<Rcpp::NumericMatrix>(SettleParamsR.slot("Settle"));
	Rcpp::LogicalVector FindMate = Rcpp::as<Rcpp::LogicalVector>(SettleParamsR.slot("FindMate"));
	bool constFindMate = false;
	if(FindMate.length() == 1) {
		constFindMate = true;
	}
	Rcpp::NumericVector MutationCoeffs = Rcpp::as<Rcpp::NumericVector>(SettleParamsR.slot("TraitScaleFactor"));
	Rcpp::IntegerVector MinSteps = Rcpp::as<Rcpp::IntegerVector>(SettleParamsR.slot("MinSteps"));
	bool constMinSteps = false;
	if(MinSteps.length() == 1) {
		constMinSteps = true;
	}
	Rcpp::IntegerVector MaxSteps = Rcpp::as<Rcpp::IntegerVector>(SettleParamsR.slot("MaxSteps"));
	bool constMaxSteps = false;
	if(MaxSteps.length() == 1) {
		constMaxSteps = true;
	}
	Rcpp::IntegerVector MaxStepsYr = Rcpp::as<Rcpp::IntegerVector>(SettleParamsR.slot("MaxStepsYear"));
	bool constMaxStepsYr = false;
	if(MaxStepsYr.length() == 1) {
		constMaxStepsYr = true;
	}

	sexSettle = 2 * sett.stgDep + sett.sexDep;

	for(int line = 0; line < Nlines; line++) {

		// determine stage and sex of this line
		if(sett.stgDep) {
			if(sett.sexDep) {
				stage = (int)SettleCondMatrix(line, 0);
				sex = (int)SettleCondMatrix(line, 1);
			} else {
				stage = (int)SettleCondMatrix(line, 0);
				sex = 0;
			}
		} else {
			if(sett.sexDep) {
				stage = 0;
				sex = (int)SettleCondMatrix(line, 0);
			} else {
				stage = 0;
				sex = 0;
			}
		}

		// read settlement conditions for...
		if(trfr.moveModel) { // ...movement process            //  /!\ different to ReadSettlement() : all cases are
			// covered here (in a way that parameters for IIV (i.e. 'settParams sparams;') are only set
			// for stage #0,
			//                                      however on R-level (IndVar && StageDep) is not
			//                                      admissible

			// densdep = (bool)SettleCondMatrix(line, offset+0);
			if(constFindMate) {
				findmate = (bool)FindMate(0);
			} else {
				findmate = (bool)FindMate(line);
			}
			if(findmate && dem.repType == 0)
				error = 504;

			if(constMinSteps) {
				ssteps.minSteps = (int)MinSteps(0);
			} else {
				ssteps.minSteps = (int)MinSteps(line);
			}
			if(constMaxSteps) {
				ssteps.maxSteps = (int)MaxSteps(0);
			} else {
				ssteps.maxSteps = (int)MaxSteps(line);
			}
			if(constMaxStepsYr) {
				ssteps.maxStepsYr = (int)MaxStepsYr(0);
			} else {
				ssteps.maxStepsYr = (int)MaxStepsYr(line);
			}
			if(densdep) {
				if(sett.indVar) {
					sparams.s0Mean = (float)SettleCondMatrix(
					                     line, offset + 0); // Required for DensDep = 1 and IndVar = 1. 0.0 < S0 <= 1.0
					sparams.s0SD = (float)SettleCondMatrix(
					                   line, offset + 1); // Required for DensDep = 1 and IndVar = 1. 0.0 < S0 <= 1.0
					sparams.alphaSMean = (float)SettleCondMatrix(line, offset + 2);
					sparams.alphaSSD = (float)SettleCondMatrix(line, offset + 3);
					sparams.betaSMean = (float)SettleCondMatrix(line, offset + 4);
					sparams.betaSSD = (float)SettleCondMatrix(line, offset + 5);
					if(stage == 0) {
						sparams.s0Scale = (float)MutationCoeffs(0);
						sparams.alphaSScale = (float)MutationCoeffs(1);
						sparams.betaSScale = (float)MutationCoeffs(2);
					}
				} else {
					settleDD.s0 = (float)SettleCondMatrix(
					                  line, offset + 0); // Max. settlement probability for density reaction norm. Required for
					// DensDep = 1 and IndVar = 0; 0.0 < S0 <= 1.0
					settleDD.alpha =
					    (float)SettleCondMatrix(line, offset + 1); // Required for DensDep = 1 and IndVar = 0
					settleDD.beta =
					    (float)SettleCondMatrix(line, offset + 2); // Required for DensDep = 1 and IndVar = 0
				}
			}

			switch(sexSettle) {

			case 0: { // no sex- / stage-dependence     // why do the parameters for stage=0, sex=1 not get set if
				// (dem.stageStruct) ??? (remove the else?? ) and why the different rules for the sexes regarding
				// setSettTraits() for stages>0
				srules = pSpecies->getSettRules(0, 0);
				srules.densDep = densdep;
				srules.findMate = findmate;
				pSpecies->setSettRules(0, 0, srules);
				pSpecies->setSteps(0, 0, ssteps);
				if(srules.densDep) {
					if(sett.indVar)
						pSpecies->setSettParams(0, 0, sparams);
					else
						pSpecies->setSettTraits(0, 0, settleDD);
				}
				if(dem.stageStruct) { // model is structured - also set parameters for all stages
					for(int i = 1; i < sstruct.nStages; i++) {
						pSpecies->setSettRules(i, 0, srules);
						pSpecies->setSteps(i, 0, ssteps);
						if(srules.densDep && !sett.indVar)
							pSpecies->setSettTraits(i, 0, settleDD); //  /!\ different to ReadSettlement()
						if(dem.repType > 0) {                        // model is sexual - also set parameters for males
							pSpecies->setSettRules(i, 1, srules);
							pSpecies->setSteps(i, 1, ssteps);
							if(srules.densDep && !sett.indVar)
								pSpecies->setSettTraits(i, 1, settleDD);
						}
					}
				} else {                  // see comment above (at case label)
					if(dem.repType > 0) { // model is sexual - also set parameters for males
						pSpecies->setSettRules(0, 1, srules);
						pSpecies->setSteps(0, 1, ssteps);
						if(srules.densDep) {
							if(sett.indVar)
								pSpecies->setSettParams(0, 1, sparams);
							else
								pSpecies->setSettTraits(0, 1, settleDD);
						}
					}
				}
			}
			break;

			case 1: { // sex-dependent
				srules = pSpecies->getSettRules(0, sex);
				srules.densDep = densdep;
				srules.findMate = findmate;
				pSpecies->setSettRules(0, sex, srules);
				pSpecies->setSteps(0, sex, ssteps);
#if RSDEBUG
				DEBUGLOG << "ReadSettlementR(): stage=" << stage << " sex=" << sex
				         << " ssteps.maxStepsYr =" << ssteps.maxStepsYr << endl;
#endif
				if(srules.densDep) {
					if(sett.indVar)
						pSpecies->setSettParams(0, sex, sparams);
					else
						pSpecies->setSettTraits(0, sex, settleDD);
				}
				if(dem.stageStruct) { // model is structured - also set parameters for all stages
					for(int i = 1; i < sstruct.nStages; i++) {
						pSpecies->setSettRules(i, sex, srules);
						pSpecies->setSteps(i, sex, ssteps);
						if(srules.densDep && !sett.indVar)
							pSpecies->setSettTraits(i, sex, settleDD);
					}
				}
			}
			break;

			case 2: { // stage-dependent
				srules = pSpecies->getSettRules(stage, 0);
				srules.densDep = densdep;
				srules.findMate = findmate;
				pSpecies->setSettRules(stage, 0, srules);
				pSpecies->setSteps(stage, 0, ssteps);
				if(srules.densDep) {
					if(sett.indVar) {
						if(stage == 0)
							pSpecies->setSettParams(0, 0, sparams);
					} else
						pSpecies->setSettTraits(stage, 0, settleDD);
				}
				if(dem.repType > 0) { // model is sexual - also set parameters for males
					pSpecies->setSettRules(stage, 1, srules);
					pSpecies->setSteps(stage, 1, ssteps);
					if(srules.densDep) {
						if(sett.indVar) {
							if(stage == 0)
								pSpecies->setSettParams(0, 1, sparams);
						} else
							pSpecies->setSettTraits(stage, 1, settleDD);
					}
				}
			}
			break;

			case 3: { // sex- & stage-dependent
				srules = pSpecies->getSettRules(stage, sex);
				srules.densDep = densdep;
				srules.findMate = findmate;
				pSpecies->setSettRules(stage, sex, srules);
				pSpecies->setSteps(stage, sex, ssteps);
				if(srules.densDep) {
					if(sett.indVar) {
						if(stage == 0)
							pSpecies->setSettParams(0, sex, sparams);
					} else
						pSpecies->setSettTraits(stage, sex, settleDD);
				}
			}
			break;
			} // end sexSettle

		} // end of movement model

		// read settlement conditions for...
		else { // ...dispersal kernel

			settType = (int)SettleCondMatrix(line,
			                                 offset); // Settlement rule if the arrival cell/patch is unsuitable: 0 = die, 1 = wait, 2 = randomly
			// choose a suitable cell/patch or die, 3 = randomly choose a suitable cell/patch or wait.
			// Options 1 and 3 may be chosen for a stage-structured population only
			if(constFindMate) {
				findmate = (bool)FindMate(0);
			} // Mating requirements to settle, required for a sexual population only
			else {
				findmate = (bool)FindMate(line);
			}
			if(findmate && dem.repType == 0)
				error = 504;

			switch(sexSettle) {
			case 0: { // no sex / stage dependence
				if((settType == 1 || settType == 3) && !dem.stageStruct)
					error = 503;
				if(findmate && dem.repType == 0)
					error = 504;
				srules = pSpecies->getSettRules(0, 0);
				switch(settType) {
				case 0:
					srules.wait = false;
					srules.go2nbrLocn = false;
					break;
				case 1:
					srules.wait = true;
					srules.go2nbrLocn = false;
					break;
				case 2:
					srules.wait = false;
					srules.go2nbrLocn = true;
					break;
				case 3:
					srules.wait = true;
					srules.go2nbrLocn = true;
					break;
				}
				srules.findMate = findmate;
				if(dem.stageStruct) { // model is structured - also set parameters for all stages
					for(int i = 0; i < sstruct.nStages; i++) {
						pSpecies->setSettRules(i, 0, srules);
						if(dem.repType > 0) { // model is sexual - also set parameters for males
							pSpecies->setSettRules(i, 1, srules);
						}
					}
				} else {
					pSpecies->setSettRules(0, 0, srules);
					if(dem.repType > 0) { // model is sexual - also set parameters for males
						pSpecies->setSettRules(0, 1, srules);
					}
				}
			}
			break;

			case 1: { // sex dependent
				if((settType == 1 || settType == 3) && dem.stageStruct == false)
					error = 505;
				srules = pSpecies->getSettRules(0, sex);
				switch(settType) {
				case 0:
					srules.wait = false;
					srules.go2nbrLocn = false;
					break;
				case 1:
					srules.wait = true;
					srules.go2nbrLocn = false;
					break;
				case 2:
					srules.wait = false;
					srules.go2nbrLocn = true;
					break;
				case 3:
					srules.wait = true;
					srules.go2nbrLocn = true;
					break;
				}
				srules.findMate = findmate;
				pSpecies->setSettRules(0, sex, srules);
				if(dem.stageStruct) { // model is structured - also set parameters for all stages
					for(int i = 1; i < sstruct.nStages; i++) {
						pSpecies->setSettRules(i, sex, srules);
					}
				}
			}
			break;

			case 2: { // stage dependent
				if(findmate && dem.repType == 0)
					error = 507;
				srules = pSpecies->getSettRules(stage, 0);
				switch(settType) {
				case 0:
					srules.wait = false;
					srules.go2nbrLocn = false;
					break;
				case 1:
					srules.wait = true;
					srules.go2nbrLocn = false;
					break;
				case 2:
					srules.wait = false;
					srules.go2nbrLocn = true;
					break;
				case 3:
					srules.wait = true;
					srules.go2nbrLocn = true;
					break;
				}
				srules.findMate = findmate;
				pSpecies->setSettRules(stage, 0, srules);
				if(dem.repType > 0) { // model is sexual - also set parameters for males
					pSpecies->setSettRules(stage, 1, srules);
				}
			}
			break;

			case 3: { // sex & stage dependent
				srules = pSpecies->getSettRules(stage, sex);
				switch(settType) {
				case 0:
					srules.wait = false;
					srules.go2nbrLocn = false;
					break;
				case 1:
					srules.wait = true;
					srules.go2nbrLocn = false;
					break;
				case 2:
					srules.wait = false;
					srules.go2nbrLocn = true;
					break;
				case 3:
					srules.wait = true;
					srules.go2nbrLocn = true;
					break;
				}
				srules.findMate = findmate;
				pSpecies->setSettRules(stage, sex, srules);
			}
			break;

			} // end of switch (sexSettle)

		} // end of dispersal kernel

	} // end of for line loop

	if(sett.indVar) anyIndVar = true;

	return error;
}

//---------------------------------------------------------------------------

int ReadInitialisationR(Landscape* pLandscape, Rcpp::S4 ParMaster)
{

	Rcpp::S4 InitParamsR("InitialisationParams");
	InitParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("init"));

	Rcpp::NumericVector PropStages;

	landParams paramsLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogr();
	stageParams sstruct = pSpecies->getStage();
	initParams init = paramsInit->getInit();
	string Inputs = paramsSim->getDir(1);

	int maxcells; //,simulation
	float p, check;
	int error = 0;

	// simulation = Rcpp::as<int>(InitParamsR.slot("Simulation"));// REMOVED in R-interface
	init.seedType = Rcpp::as<int>(
	                    InitParamsR.slot("InitType")); // Initialisation type: 0 = free initialisation, 1 = from species distribution, 2
	// = from initial individuals file. Must be 0 for artificial landscapes.
	init.freeType = Rcpp::as<int>(InitParamsR.slot(
	                                  "FreeType")); // If SeedType = 0: Type of free initialisation (i.e. not from specified distribution): 0 = Random
	// (given no. of cells/patches), 1 = All suitable cells/patches.
	init.spDistType = Rcpp::as<int>(
	                      InitParamsR.slot("SpType")); // If SeedType = 1: Type of initialisation from species distribution: 0 = All
	// suitable cells/patches within the distribution cells, 1 = Some randomly chosen
	// suitable cells/patches within distribution cells
	if(init.seedType == 1 && paramsLand.spDist == false)
		error = 601;

	if(paramsLand.patchModel) {
		init.initDens = Rcpp::as<int>(InitParamsR.slot(
		                                  "InitDens")); // Initial density of individuals in each patch: 0 = at carrying capacity; 1
		// = at half carrying capacity; 2 = specified density. Patch carrying capacity
		// is determined by the carrying capacities of the cells forming the patch.
		init.indsHa =
		    Rcpp::as<float>(InitParamsR.slot("IndsHaCell")); // If initDens = 2: Density of individuals to initialise
		// per hectare; required for initDens = 2.
	} else {
		init.initDens = Rcpp::as<int>(
		                    InitParamsR.slot("InitDens")); // Initial no. of individuals in each cell: 0 = at carrying capacity; 1 = at
		// half carrying capacity; 2 = specified no. of individuals
		init.indsCell =
		    Rcpp::as<int>(InitParamsR.slot("IndsHaCell")); // If initDens = 2: No. of individuals to initialise in each
		// cell; required for initDens = 2, otherwise set to -9
	}

	init.minSeedX =
	    Rcpp::as<int>(InitParamsR.slot("minX")); // Required for SeedType = 0 only; minX and minY may not be less than
	// 0, maxX and maxY may not be less than corresponding minimum
	init.maxSeedX = Rcpp::as<int>(InitParamsR.slot("maxX")); // -"-
	init.minSeedY = Rcpp::as<int>(InitParamsR.slot("minY")); // -"-
	init.maxSeedY = Rcpp::as<int>(InitParamsR.slot("maxY")); // -"-
	init.nSeedPatches = Rcpp::as<int>(
	                        InitParamsR.slot("NrCells")); // If SeedType = 0 and FreeType = 0: No. of cells/patches to initialise
	init.nSpDistPatches = Rcpp::as<int>(InitParamsR.slot(
	                                        "NrCells")); // If  SeedType = 1 and SpType = 1: No. of species dist cells/patches to initialise
	init.initFrzYr =
	    Rcpp::as<int>(InitParamsR.slot("InitFreezeYear")); // Year until which species is confined to initial range
	// limits. Must be >= 0 for SeedType = 0.
	init.restrictRows =
	    Rcpp::as<int>(InitParamsR.slot("RestrictRows")); // No. of rows at northern front to restrict range. Must be >0
	// if applied for SeedType = 0, otherwise 0.
	init.restrictFreq =
	    Rcpp::as<int>(InitParamsR.slot("RestrictFreq")); // Frequency (years) at which range is restricted to northern
	// front. Must be > 0 if RestrictRows is > 0.
	init.finalFrzYr =
	    Rcpp::as<int>(InitParamsR.slot("FinalFreezeYear")); // Year after which species is confined to current range
	// limits. Must be > InitFreezeYear if applied, otherwise 0.
	init.indsFile = Rcpp::as<string>(InitParamsR.slot("InitIndsFile")); // Name of the initial individuals file (*.txt). Required for SeedType = 2, otherwise NULL.
#if RSDEBUG
	DEBUGLOG << "ReadInitialisationR():"
	         //<< " simulation=" << simulation
	         << " seedType=" << init.seedType
	         << " freeType=" << init.freeType << " spDistType=" << init.spDistType << " maxSeedX=" << init.maxSeedX
	         << " maxSeedY=" << init.maxSeedY << " initFrzYr=" << init.initFrzYr
	         << " restrictRows=" << init.restrictRows << " restrictFreq=" << init.restrictFreq
	         << " finalFrzYr=" << init.finalFrzYr << " indsFile=" << init.indsFile << endl;
#endif
	init.restrictRange = false;
	if(init.seedType == 0 && init.restrictRows > 0)
		init.restrictRange = true;

	if(dem.stageStruct) {
		init.initAge = Rcpp::as<int>(InitParamsR.slot("InitAge")); // Initial age distribution within each stage: 0 = lowest possible age, 1 = randomised, 2 =
		// quasi-equilibrium. Required for StageStruct = 1 only - otherwise OMIT COLUMNS.
		PropStages = Rcpp::as<Rcpp::NumericVector>(InitParamsR.slot("PropStages")); // Proportion of the initial individuals in stage class i>0 (No juveniles
		// are initialized). Required for StageStruct = 1 only (number of columns
		// is one less than number of stages) - otherwise OMIT COLUMNS
		if(init.seedType != 2) {
			check = 0.0;
			for(int i = 1; i < sstruct.nStages; i++) {
				p = (float)PropStages(i);
				check += p;
				paramsInit->setProp(i, p);
			}
			if(check < 1.0 || check > 1.0) {
				// this condition should not occur - WHAT COULD BE DONE?? ABORT WITH ERROR CODE ...
#if RSDEBUG
				DEBUGLOG << "ReadInitialisation(): check = " << check << endl;
#endif
			}
		}
	}
	paramsInit->setInit(init);
	switch(init.seedType) {
	case 0: { // free initialisation
		if(init.minSeedX < 0)
			init.minSeedX = 0;
		if(init.minSeedY < 0)
			init.minSeedY = 0;
		if(init.maxSeedX < 0 || init.maxSeedX > paramsLand.maxX)
			init.maxSeedX = paramsLand.maxX; // added upper boundary checks
		if(init.maxSeedY < 0 || init.maxSeedY > paramsLand.maxY)
			init.maxSeedY = paramsLand.maxY; // added upper boundary checks
		if(init.minSeedY > init.maxSeedY || init.minSeedX > init.maxSeedX) {
#if RSDEBUG
			DEBUGLOG << "ReadInitialisationR(): maxSeedX=" << init.maxSeedX << " paramsLand.maxX=" << paramsLand.maxX
			         << " maxSeedY=" << init.maxSeedY << " paramsLand.maxY=" << paramsLand.maxY << endl;
#endif
			error = 603;
		}
		maxcells = (init.maxSeedY - init.minSeedY) * (init.maxSeedX - init.minSeedX);
		if(init.freeType == 0 && init.nSeedPatches > maxcells)
			error = 602;
	}
	break;
	case 1: // from species distribution
		break;
	case 2: { // from initial individuals file
		// if (init.indsFile != prevInitialIndsFile) {
		// read and store the list of individuals to be initialised
		error = ReadInitIndsFileR(0, pLandscape); //open, parse, read header and lines, store in vector "initinds"
		// prevInitialIndsFile = init.indsFile;
		//}
	}
	break;
	default:
		;
	}
	return error;
}

//---------------------------------------------------------------------------

int ReadGeneticsR(Rcpp::S4 GeneParamsR)
{
	emigRules emig = pSpecies->getEmig();
	trfrRules trfr = pSpecies->getTrfr();
	settleType sett = pSpecies->getSettle();
	demogrParams dem = pSpecies->getDemogr();

	int arch;
	string archfile;
	int error = 0;
	genomeData g;

	arch = Rcpp::as<int>(GeneParamsR.slot("Architecture")); // 0 = One chromosome per trait, 1 = Read from file
	g.nLoci = Rcpp::as<int>(GeneParamsR.slot("NLoci")); // No. of loci per chromosome, Required for Architecture=0; > 0
	archfile = Rcpp::as<string>(GeneParamsR.slot("ArchFile")); // Name of the genetic architecture file (*.txt),Required for Architecture=1, otherwise NULL
	g.probMutn = Rcpp::as<float>(GeneParamsR.slot("ProbMutn")); // Probability of mutation of an individual allele at meiosis, 0 <= prob <= 1
	g.probCrossover = Rcpp::as<float>(GeneParamsR.slot("ProbCross")); // Probability of crossover at an individual locus at meiosis, 0 <= prob <= 1
	g.alleleSD = Rcpp::as<float>(GeneParamsR.slot("AlleleSD")); // S.d. of initial allelic values around phenotypic value, > 0
	g.mutationSD = Rcpp::as<float>(GeneParamsR.slot("MutationSD")); // S.d. of mutation magnitude, > 0

	if (dem.repType == 0) g.diploid = false;
	else g.diploid = true;

#if RSDEBUG
	DEBUGLOG << "ReadGeneticsR(): arch=" << arch
	         << " g.nLoci=" << g.nLoci << " archfile=" << archfile
	         << " g.probMutn=" << g.probMutn << " g.probCrossover=" << g.probCrossover
	         << " g.alleleSD=" << g.alleleSD << " g.mutationSD=" << g.mutationSD
	         << endl;
#endif

	g.neutralMarkers = false;
	if (arch == 0) { // no architecture file
		g.trait1Chromosome = true;
		pSpecies->set1ChromPerTrait(g.nLoci);
	} else { // architecture file
		g.trait1Chromosome = false;
		g.nLoci = 0;
		if (!(emig.indVar || trfr.indVar || sett.indVar)) {
			g.neutralMarkers = true;
		}
		//check architecture file
		Rcpp::Rcout << "Checking Archiecture file " << archfile << endl;
		string fname = paramsSim->getDir(1) + archfile;
		wifstream archFile;
		//open file
#if RSWIN64
		archFile.open(fname.c_str());
#else
		archFile.open(fname, std::ios::binary);
#endif
		if(!archFile.is_open()) {
			OpenErrorR("Architecture file", fname);
#if RSDEBUG
			DEBUGLOG << "Architecture file failed to open: " << fname << std::endl;
#endif
			return -217;
		} else {
#if RSDEBUG
			DEBUGLOG << "Architecture file open to read" << std::endl;
#endif
#if !RSWIN64
			// check BOM for UTF-16
			if(check_bom(fname) == "utf16")
				// apply BOM-sensitive UTF-16 facet
				archFile.imbue(std::locale(archFile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif

			error = ReadArchFileR(archFile);

			if (archFile.is_open()) archFile.close();
			archFile.clear();
		}
	}

	pSpecies->setGenomeData(g);

	return error;
}

//---------------------------------------------------------------------------

int ReadArchFileR(wifstream& archFile)
{
	//int ParseArchFile(void){

	wstring paramname;
	int traitnum,prevtrait,nchromosomes,chrom,locus,nloci;
	int errors = 0;
	int fileNtraits = 0;
	bool formatError = false;
	int *chromsize = 0;
	string filetype = "ArchFile";

// check no. of chromosomes, and terminate if in error
	archFile >> paramname >> nchromosomes;
	if (paramname == L"NChromosomes") {
		if (nchromosomes < 1) {
			BatchErrorR(filetype,-999,11,"NChromosomes");
			errors++;
			return -111;
		}
	} else {
		ArchFormatErrorR();
		return -111;
	}
	chromsize = new int[nchromosomes];
	for (int i = 0; i < nchromosomes; i++) chromsize[i] = 0;

// check no. of loci on each chromosome, and terminate if in error
	archFile >> paramname;
	if (paramname != L"NLoci") formatError = true;
	int locerrors = 0;
	for (int i = 0; i < nchromosomes; i++) {
		nloci = -999;
		archFile >> nloci;
		if (archFile.eof()){ // iOS problem: if at end of file: problems clearing the status flags later; therefore: clear it now (eof-flag will also be cleared)
		    archFile.clear();
		    }
		if (nloci < 1) locerrors++;
		else chromsize[i] = nloci;
	}
	if (locerrors) {
		BatchErrorR(filetype,-999,11,"NLoci");
		return -111;
	}
//if (formatError) batchlog << "formatError is TRUE" << endl;
//else batchlog << "formatError is FALSE" << endl;

// check unspecified no. of traits
	fileNtraits = 0;
	traitnum = prevtrait = -1;
	bool traitError = false;
	bool lociError = false;
	bool chromError = false;
	bool locusError = false;
	paramname = L"XXXyyyZZZ";
//batchlog << "paramname=" << paramname << endl;
	if (!archFile.eof()) archFile >> paramname; // only read another parameter if not end of file (if clear() was called before, it will try reading in another character, but won't succeed -> not going into while-loop
//batchlog << "paramname=" << paramname << endl;
	while (paramname != L"XXXyyyZZZ") {
		archFile >> traitnum;
//	batchlog << "traitnum=" << traitnum << endl;
		if (paramname != L"Trait") formatError = true;
		if (traitnum == (prevtrait+1)) prevtrait = traitnum;
		else traitError = true;
		archFile >> paramname >> nloci;
//	batchlog << "paramname=" << paramname << " nloci=" << nloci << endl;
		if (paramname != L"NLoci") formatError = true;
		if (nloci < 1) lociError = true;
		for (int i = 0; i < nloci; i++) {
			chrom = locus = -999999;
			archFile >> chrom >> locus;
//		batchlog << "chrom=" << chrom << " locus=" << locus << endl;
			if (chrom == -999999 || locus == -999999) {
				BatchErrorR(filetype,-999,0," ");
				errors++;
				Rcpp::Rcout << "Too few loci listed for trait " << traitnum << endl;
			} else {
				if (chrom >= 0 && chrom < nchromosomes) {
//				batchlog << "chromsize[" << chrom << "]=" << chromsize[chrom] << endl;
					if (locus < 0 || locus >= chromsize[chrom]) locusError = true;
				} else chromError = true;
			}
		}
		fileNtraits++;
		paramname = L"XXXyyyZZZ";
		if (archFile.eof()){ // iOS problem: if at end of file: problems clearing the status flags later; therefore: clear it now (eof-flag will also be cleared)
		    archFile.clear();
		}
		if (!archFile.eof()) archFile >> paramname; // only read another parameter if not end of file (if clear() was called before, it will try reading in another character, but won't succeed go out of while loop
//	batchlog << "paramname=" << paramname << " (end of loop)" << endl;
	}
//batchlog << "paramname=" << paramname << " (after loop)" << endl;

	if (traitError) {
		BatchErrorR(filetype,-999,0," ");
		errors++;
		Rcpp::Rcout << "Traits must be sequentially numbered starting at 0 " << endl;
	}
	if (lociError) {
		BatchErrorR(filetype,-999,11,"Trait NLoci");
		errors++;
	}
	if (chromError) {
		BatchErrorR(filetype,-999,0," ");
		errors++;
		Rcpp::Rcout << "Chromosome no. must be from 0 to " << (nchromosomes-1) << endl;
	}
	if (locusError) {
		BatchErrorR(filetype,-999,0," ");
		errors++;
		Rcpp::Rcout << "Locus no. must not exceed no. of loci on specified chromosome " << endl;
	}

	if (chromsize != 0) delete[] chromsize;

// check if parsing was successful before starting to read
	if (formatError || errors > 0) { // terminate batch error checking
		if (formatError) ArchFormatErrorR();
		return -111;
	} else {
		// final read should have hit EOF
		if (!archFile.eof()) {
			EOFerrorR(filetype);
		}
		archFile.clear();
		archFile.seekg(0);
		archFile.sync();
		if(!archFile.good()) {
			Rcpp::Rcout << "Error re-opening Architecture file with state " << archFile.rdstate() << endl;
			return -331;
		}
		else Rcpp::Rcout << "Architecture file OK" << endl;
	}

	// READING
	//int ReadArchFile(string archfile){

	emigRules emig = pSpecies->getEmig();
	trfrRules trfr = pSpecies->getTrfr();
	settleType sett = pSpecies->getSettle();

// set no. of chromosomes
	archFile >> paramname >> nchromosomes;
	pSpecies->setNChromosomes(nchromosomes);
	int nchromset = pSpecies->getNChromosomes();

	if (nchromset <= 0) errors = 1;
	if (emig.indVar || trfr.indVar || sett.indVar) {
		pSpecies->setTraitData(fileNtraits);
	} else { // neutral markers only
		pSpecies->setTraitData(0);
	}
// set no. of loci for each chromosome
	archFile >> paramname;
	for (int i = 0; i < nchromosomes; i++) {
		archFile >> nloci;
		pSpecies->setNLoci(i,nloci);
	}
	if (emig.indVar || trfr.indVar || sett.indVar) {
		// set trait maps
		paramname = L"XXXyyyZZZ";
		archFile >> paramname;
		while (paramname != L"XXXyyyZZZ") {
			archFile >> traitnum >> paramname >> nloci;
			pSpecies->setTraitMap(traitnum,nloci);
			for (int allele = 0; allele < nloci; allele++) {
				chrom = locus = -999999;
				archFile >> chrom >> locus;
				pSpecies->setTraitAllele(traitnum,allele,chrom,locus);
			}
			paramname = L"XXXyyyZZZ";
			archFile >> paramname;
		};
	}

// any loci not contributing to a trait are recorded as neutral
	if (emig.indVar || trfr.indVar || sett.indVar) {
		pSpecies->setNeutralLoci(false);
	} else { // model has neutral markers only
		pSpecies->setNeutralLoci(true);
	}

	return errors;
}

//---------------------------------------------------------------------------

Rcpp::List RunBatchR(int nSimuls, int nLandscapes, Rcpp::S4 ParMaster)
{
	int land_nr;
	int t0, t1, t00, t01;
	int read_error;
	bool params_ok;
	simParams sim = paramsSim->getSim();

	Rcpp::List list_outPop;
	Landscape* pLandscape = NULL; // pointer to landscape

#if RSDEBUG
	DEBUGLOG << endl;
	DEBUGLOG << "RunBatchR(): nSimuls=" << nSimuls << " nLandscapes=" << nLandscapes << endl;
	DEBUGLOG << "RunBatchR(): landtype=" << landtype << " maxNhab=" << maxNhab << endl;
#endif

	t0 = time(0);

	// int batch_line = 0;

	string name = paramsSim->getDir(2) + "Batch" + Int2Str(sim.batchNum) + "_RS_log.csv";
	if(rsLog.is_open()) {
		rsLog.close();
		rsLog.clear();
	}
	rsLog.open(name.c_str());
	if(!rsLog.is_open()) {
		Rcpp::Rcout << endl
		            << "Error - unable to open Batch" << sim.batchNum << "_RS_log.csv file - aborting batch run"
		            << endl;
		return Rcpp::List::create(Rcpp::Named("Errors") = -678);
	}
	rsLog << "Event,Number,Reps,Years,Time" << endl;
#if RSDEBUG
	rsLog << "WARNING,***** RSDEBUG mode is active *****,,," << endl;
#endif
	rsLog << "RNG SEED,,,," << RS_random_seed << endl;

	// loop over landscpaes

	for(int j = 0; j < nLandscapes; j++) {

#if RSDEBUG
		DEBUGLOG << endl;
#endif
		// create new landscape
		if(pLandscape != NULL)
			delete pLandscape;
		pLandscape = new Landscape;
		bool landOK = true;

		t00 = time(0);

		landOK = ReadLandParamsR(pLandscape, ParMaster);
		//land_nr = ReadLandParamsR(pLandscape, ParMaster);
		land_nr = j; // TODO: ReadLandParamsR() is supposed to return land_nr; this is a temporary replacement

		if(!landOK) {
			rsLog << "Error reading landscape ASCII haeders - aborting" << endl;
			Rcpp::Rcout << "Error reading landscape ASCII haeders - aborting" << endl;
		} else {

			MemoLine(("Starting landscape " + Int2Str(land_nr) + "...").c_str());

#if RSDEBUG
			DEBUGLOG << endl << "RunBatchR(): j=" << j << " land_nr=" << land_nr << " landtype=" << landtype;
			if(landtype != 9)
				DEBUGLOG << " name_landscape=" << name_landscape
				         << " name_patch=" << name_patch
				         << " name_costfile=" << name_costfile
				         << " name_sp_dist=" << name_sp_dist;
			DEBUGLOG << endl;
#endif
			landParams paramsLand = pLandscape->getLandParams();
			paramsLand.patchModel = patchmodel;
			paramsLand.resol = resolution;
			paramsLand.rasterType = landtype;
			if(landtype == 9) {
				paramsLand.generated = true;
				paramsLand.nHab = 2;
			} else {
				paramsLand.generated = false;
				/*if(name_dynland == "NULL")
					paramsLand.dynamic = false;
				else
					paramsLand.dynamic = true;*/
			}
			paramsLand.nHabMax = maxNhab;
			paramsLand.spDist = speciesdist;
			paramsLand.spResol = distresolution;
			pLandscape->setLandParams(paramsLand, sim.batchMode);

			if(landtype != 9) { // imported landscape
				string hname = paramsSim->getDir(1) + name_landscape;
				int landcode;
				string cname;
				if (name_costfile == "NULL" || name_costfile == "none") cname = "NULL";
				else cname = paramsSim->getDir(1) + name_costfile;
				if(paramsLand.patchModel) {
					string pname = paramsSim->getDir(1) + name_patch;
					landcode = pLandscape->readLandscape(0, hname, pname, cname);
				} else
					landcode = pLandscape->readLandscape(0, hname, " ", cname);
				if(landcode != 0) {
					rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
					Rcpp::Rcout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
					landOK = false;
				}
				if(paramsLand.dynamic) {
					Rcpp::S4 LandParamsR("LandParams");
					LandParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("land"));
					landcode = ReadDynLandR(pLandscape, LandParamsR);
					if(landcode != 0) {
						rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
						Rcpp::Rcout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
						landOK = false;
					}
				}
				if(landtype == 0) {
					pLandscape->updateHabitatIndices();
				}
#if RSDEBUG
				landParams tempLand = pLandscape->getLandParams();
				DEBUGLOG << "RunBatchR(): j=" << j << " land_nr=" << land_nr << " landcode=" << landcode
				         << " nHab=" << tempLand.nHab << endl;
#endif

				// species distribution

				if(paramsLand.spDist) { // read initial species distribution
					// WILL NEED TO BE CHANGED FOR MULTIPLE SPECIES ...
					string distname = paramsSim->getDir(1) + name_sp_dist;
					landcode = pLandscape->newDistribution(pSpecies, distname);
					if(landcode == 0) {
					} else {
						rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
						Rcpp::Rcout << endl
						            << "Error reading initial distribution for landscape " << land_nr << " - aborting"
						            << endl;
						landOK = false;
					}
				}
				paramsSim->setSim(sim);
#if RSDEBUG
				DEBUGLOG << "RunBatchR(): j=" << j << " spDist=" << paramsLand.spDist << endl;
#endif

				if(landOK) {
					t01 = time(0);
					rsLog << "Landscape," << land_nr << ",,," << t01 - t00 << endl;

				} // end of landOK condition

			} // end of imported landscape
		}
		if(landOK) {

			// Open all other batch files and read header records

			// nSimuls is the total number of lines (simulations) in
			// the batch and is set in the control function
			string msgsim = "Simulation,";
			string msgerr = ",ERROR CODE,";
			string msgabt = ",simulation aborted";
			for(int i = 0; i < nSimuls; i++) { // this loop is useless at the moment since nSimuls is set to one in R entry function BatchMainR()
				t00 = time(0);
				params_ok = true;
				anyIndVar = false;
				read_error = ReadParametersR(pLandscape, ParMaster);
				simParams sim = paramsSim->getSim();
				if(read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				if(stagestruct) {
					ReadStageStructureR(ParMaster);
				}
				read_error = ReadEmigrationR(ParMaster);
				if(read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				read_error = ReadTransferR(pLandscape, ParMaster);
				if(read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				read_error = ReadSettlementR(ParMaster);
				if(read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				if(params_ok) {
#if RSDEBUG
					DebugGUI("RunBatchR(): simulation i=" + Int2Str(i));
#endif
					pSpecies->setNChromosomes(0);
					pSpecies->setTraits();
				}
				Rcpp::S4 GeneParamsR("GeneticsParams");
				GeneParamsR = Rcpp::as<Rcpp::S4>(ParMaster.slot("gene"));
				if (anyIndVar || Rcpp::as<int>(GeneParamsR.slot("Architecture")) == 1) {
					read_error = ReadGeneticsR(GeneParamsR);
					if(read_error) {
						rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
						params_ok = false;
					}
				} else {
					// use default genetics parameters
					// (by setting illegal values except for diploid)
					genomeData g;
					g.nLoci = -1;
					g.probMutn = g.probCrossover = g.alleleSD = g.mutationSD = -1.0;
					if(reproductn == 0)
						g.diploid = false;
					else
						g.diploid = true;
					g.neutralMarkers = g.trait1Chromosome = false;

					pSpecies->setGenomeData(g);
				}
				read_error = ReadInitialisationR(pLandscape, ParMaster);
				if(read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}

				if(params_ok) {
					simParams sim = paramsSim->getSim();

#if RSDEBUG
					DEBUGLOG << endl
					         << "RunBatchR(): i=" << i << " simulation=" << sim.simulation
							 //<< " landFile=" << landFile
					         << " outRange=" << sim.outRange << " outIntRange=" << sim.outIntRange << endl;
#endif

					Rcpp::Rcout << endl
					            << "Running simulation nr. " << Int2Str(sim.simulation)
					            //<< " on landscape no. " << Int2Str(land_nr)
					            << endl;

					MemoLine(("Starting simulation " + Int2Str(sim.simulation) + "...").c_str());

					// for batch processing, include landscape number in parameter file name
					OutParameters(pLandscape);

					// run the model
					list_outPop = RunModel(pLandscape, i);
#if RSDEBUG
					// DEBUGLOG << endl << "RunBatchR(): real landscape, i = " << i
					//	<< " simulation = " << sim.simulation << " landFile = " << landFile
					//	<< endl;
#endif

					t01 = time(0);
					rsLog << msgsim << sim.simulation << "," << sim.reps << "," << sim.years << "," << t01 - t00
					      << endl;
				} // end of if (params_ok)
				else {
					Rcpp::Rcout << endl << "Error in reading parameter file(s)... see RS log." << endl;
				}
			} // end of nSimuls for loop

			// close input files

//		if (landtype != 9) {
			if (pLandscape != NULL) {
				delete pLandscape;
				pLandscape = NULL;
			}

		} // end of landOK condition

	} // end of nLandscapes loop

// Write performance data to log file
	t1 = time(0);
	rsLog << endl << "Batch,,,," << t1 - t0 << endl;

	if(rsLog.is_open()) {
		rsLog.close();
		rsLog.clear();
	}

	return list_outPop;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void setglobalvarsR(Rcpp::S4 control)
{

	// nSimuls = 2;                                   //read from parameterfile
	// nLandscapes = 1;                               //read from landfile, but not yet used in R-class RSparams!

	batchnum = Rcpp::as<int>(control.slot("batchnum"));
	patchmodel = Rcpp::as<int>(control.slot("patchmodel"));
	resolution = Rcpp::as<int>(control.slot("resolution"));
	landtype = Rcpp::as<int>(control.slot("landtype"));
	maxNhab = Rcpp::as<int>(control.slot("maxNhab"));
	speciesdist = Rcpp::as<int>(control.slot("speciesdist"));
	distresolution = Rcpp::as<int>(control.slot("distresolution"));
	reproductn = Rcpp::as<int>(control.slot("reproductn"));
	repseasons = Rcpp::as<int>(control.slot("repseasons"));
	stagestruct = Rcpp::as<int>(control.slot("stagestruct"));
	stages = Rcpp::as<int>(control.slot("stages"));
	transfer = Rcpp::as<int>(control.slot("transfer"));

#if RSDEBUG
	/*
	Function slotNames("slotNames");
	CharacterVector snames = slotNames(obj);
	for (int i = 0; i < snames.size(); i++) {
	SEXP slot = obj.slot(Rcpp::as<std::string>(snames[i]));
	// do something with slot
	}
	*/
#endif
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#if !RSWIN64
// check for UTF-16 encoding
string check_bom(string file)
{
	/*
	const char *UTF_16_BE_BOM = "\xFE\xFF";
	const char *UTF_16_LE_BOM = "\xFF\xFE";
	const char *UTF_8_BOM     = "\xEF\xBB\xBF";
	const char *UTF_32_BE_BOM = "\x00\x00\xFE\xFF";
	const char *UTF_32_LE_BOM = "\xFF\xFE\x00\x00";
	*/
	string enc = "undef";

	ifstream infile(file, ios::in | ios::binary);
	if(!infile) {
		enc = "error";
		infile.clear();
	} else {
		char buffer[20];
		infile.read(buffer, 4);
		if(buffer[0] == '\xFE' && buffer[1] == '\xFF') { // UTF-16 LE
			enc = "utf16";
		}
		if(buffer[0] == '\xFF' && buffer[1] == '\xFE') { // UTF-16 BE
			enc = "utf16";
		}
		if(buffer[0] == '\xEF' && buffer[1] == '\xBB' && buffer[2] == '\xBF') { // UTF-8
			enc = "utf8";
		}
		if(buffer[0] == '\x00' && buffer[1] == '\x00' && buffer[2] == '\xFE' && buffer[3] == '\xFF') { // UTF-32 BE
			enc = "utf32";
		}
		if(buffer[0] == '\xFF' && buffer[1] == '\xFE' && buffer[2] == '\x00' && buffer[3] == '\x00') { // UTF-32 BE
			enc = "utf32";
		}
		/*
		if (size >= 3) {
		    if (memcmp(data, UTF_8_BOM, 3) == 0)
		        return "UTF-8";
		}
		if (size >= 4) {
		    if (memcmp(data, UTF_32_LE_BOM, 4) == 0)
		        return "UTF-32-LE";
		    if (memcmp(data, UTF_32_BE_BOM, 4) == 0)
		        return "UTF-32-BE";
		}
		if (size >= 2) {
		    if (memcmp(data, UTF_16_LE_BOM, 2) == 0)
		        return "UTF-16-LE";
		    if (memcmp(data, UTF_16_BE_BOM, 2) == 0)
		        return "UTF-16-BE";
		}*/
		infile.close();
	}
	return enc;
}
#endif

//---------------------------------------------------------------------------

rasterdata ParseRasterHead(string file)
{
	wifstream infile;
	rasterdata r;
	wstring header;
	int inint;

	r.ok = true;
	r.errors = r.ncols = r.nrows = r.cellsize = 0;
	r.xllcorner = r.yllcorner = 0.0;
	r.utf = false;

	// open file
#if RSWIN64
	infile.open(file.c_str());
#else
	infile.open(file, std::ios::binary);
#endif
	if (infile.is_open()) {
#if RSDEBUG
		DEBUGLOG << "Parsing raster file " << file << std::endl;
#endif
#if !RSWIN64
		// check BOM for UTF-16
		if(check_bom(file) == "utf16") {
			// apply BOM-sensitive UTF-16 facet
			infile.imbue(std::locale(infile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
			r.utf = true;
		} else {
			r.utf = false;
		}
#endif
		// parse ASCII header
		infile >> header;
		if (!infile.good()) {
#if RSDEBUG
			DEBUGLOG << "ParseRasterHead(): failed to read " << file << std::endl;
#endif
			r.ok = false;
			r.errors = -211;
			infile.close();
			infile.clear();
			return r;
		}
		if (header != L"ncols" && header != L"NCOLS") r.errors++;
		infile >> r.ncols;

		infile >> header >> r.nrows;
		if (header != L"nrows" && header != L"NROWS") r.errors++;

		infile >> header >> r.xllcorner;
		if (header != L"xllcorner" && header != L"XLLCORNER") r.errors++;

		infile >> header >> r.yllcorner;
		if (header != L"yllcorner" && header != L"YLLCORNER") r.errors++;

		infile >> header >> r.cellsize;
		if (header != L"cellsize" && header != L"CELLSIZE") r.errors++;

		infile >> header >> inint;
		if (header != L"NODATA_value" && header != L"NODATA_VALUE") r.errors++;

		if (r.errors > 0) r.ok = false;

	} else {
		r.ok = false;
		r.errors = -111;
		//OpenErrorR("Raster file ", file);
#if RSDEBUG
		DEBUGLOG << "Raster file failed to open: " << file << std::endl;
#endif
	}
	infile.close();
	infile.clear();
	return r;
}

//----------------------------------------------------------------------------------------------

int ReadInitIndsFileR(int option, Landscape* pLandscape)
{
	landParams paramsLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogr();
	//stageParams sstruct = pSpecies->getStage();
	initParams init = paramsInit->getInit();

	string indsfile = paramsSim->getDir(1) + init.indsFile;
	wifstream initIndsFile;

	wstring header;
	string filetype = "InitIndsFile";
	string filename, ftype2, fname;

	int line, prevyear;
	//int year, sex, species, patchID, x, y, ninds, age, stage;

	int errors = 0;

	if(option == 0) { // open file, parse and read header and lines
		// open file
#if RSWIN64
		initIndsFile.open(indsfile.c_str());
#else
		initIndsFile.open(indsfile, std::ios::binary);
#endif
		if(!initIndsFile.is_open()) {
			OpenErrorR("Initial individuals file", indsfile);
#if RSDEBUG
			DEBUGLOG << "Initial individuals file failed to open: " << indsfile << std::endl;
#endif
			return -21;
		} else {
#if RSDEBUG
			DEBUGLOG << "Initial individuals file open to read" << std::endl;
#endif
#if !RSWIN64
			// check BOM for UTF-16
			if(check_bom(indsfile) == "utf16")
				// apply BOM-sensitive UTF-16 facet
				initIndsFile.imbue(std::locale(initIndsFile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
#endif

			// Check right headers format
			initIndsFile >> header;
			if(!initIndsFile.good()) {
				Rcpp::Rcout << "Initial individuals failed to read" << std::endl;
#if RSDEBUG
				DEBUGLOG << "Initial individuals failed to read" << std::endl;
#endif
				return -211;
			}
			if(header != L"Year")
				errors++;
			initIndsFile >> header;
			if(header != L"Species")
				errors++;
			if(patchmodel) {
				initIndsFile >> header;
				if(header != L"PatchID")
					errors++;
			} else {
				initIndsFile >> header;
				if(header != L"X")
					errors++;
				initIndsFile >> header;
				if(header != L"Y")
					errors++;
			}
			initIndsFile >> header;
			if(header != L"Ninds")
				errors++;
			if(reproductn > 0) {
				initIndsFile >> header;
				if(header != L"Sex")
					errors++;
			}
			if(stagestruct) {
				initIndsFile >> header;
				if(header != L"Age")
					errors++;
				initIndsFile >> header;
				if(header != L"Stage")
					errors++;
			}
			// Report any errors in headers, and if so, terminate validation
			if(errors > 0) {
				FormatErrorR(filetype, errors);
				return -111;
			}

			paramsInit->resetInitInds();

			// Read data lines
			initInd iind;
			int ninds;
			int totinds = 0;

			line = 1;
			iind.year = prevyear = -98765;
			initIndsFile >> iind.year;

			while(iind.year != -98765) {
				// Year
				if(iind.year < 0) {
					BatchErrorR(filetype, line, 19, "Year");
					errors++;
				} else {
					if(iind.year < prevyear) {
						BatchErrorR(filetype, line, 2, "Year", "previous Year");
						errors++;
					}
				}
				prevyear = iind.year;

				// Species
				initIndsFile >> iind.species;
				if(iind.species != 0) {
					BatchErrorR(filetype, line, 0, " ");
					errors++;
					Rcpp::Rcout << "Species must be 0" << endl;
				}

				// Patch | Coordinates
				if(paramsLand.patchModel) {
					initIndsFile >> iind.patchID;
					if(iind.patchID < 1) {
						BatchErrorR(filetype, line, 11, "PatchID");
						errors++;
						iind.x = iind.y = 0;
					}
				} else {
					initIndsFile >> iind.x >> iind.y;
					if(iind.x < 0 || iind.y < 0) {
						BatchErrorR(filetype, line, 19, "X and Y");
						errors++;
						iind.patchID = 0;
					}
				}

				// No of individuals
				initIndsFile >> ninds;
				if(ninds < 1) {
					BatchErrorR(filetype, line, 11, "Ninds");
					errors++;
				}

				// Sex
				if(dem.repType > 0){
					initIndsFile >> iind.sex;
					if(iind.sex < 0 || iind.sex > 1) {
						BatchErrorR(filetype, line, 1, "Sex");
						errors++;
					}
				}
				else iind.sex = 0;

				// Stage
				if(dem.stageStruct) {
					initIndsFile >> iind.age >> iind.stage;
					if(iind.age < 1) {
						BatchErrorR(filetype, line, 11, "Age");
						errors++;
					}
					if(iind.stage < 1) {
						BatchErrorR(filetype, line, 11, "Stage");
						errors++;
					}
					if(iind.stage >= stages) {
						BatchErrorR(filetype, line, 4, "Stage", "no. of stages");
						errors++;
					}
				} else {
					iind.age = iind.stage = 0;
				}

				for(int i = 0; i < ninds; i++) {
					totinds++;
					paramsInit->addInitInd(iind);
				}

				line++;
				iind.year = -98765;				// finished current line
				if(!errors){					// check for format errors
					if(!initIndsFile.eof())		// check for EOF to end loop
					initIndsFile >> iind.year;		// read next year
				}
			} // end of while loop over lines
			if(!initIndsFile.eof()) {
				EOFerrorR(filetype);
				errors++;
			}

		} // end of "file is open"

		if(initIndsFile.is_open())
			initIndsFile.close();
		initIndsFile.clear();

		Rcpp::Rcout << "Initial individuals file OK:" << indsfile << std::endl;

		return errors; //totinds;

	} // end of option 0

	if(option == 9) { // close file
		if(initIndsFile.is_open()) {
			initIndsFile.close();
		}
		initIndsFile.clear();
		return 0;
	}
	return -1;
}

//---------------------------------------------------------------------------

void BatchErrorR(string filename, int line, int option, string fieldname)
{
	if(line == -999) { // message does not cite line number
		Rcpp::Rcout << "*** Error in " << filename << ": ";
	} else {
		Rcpp::Rcout << "*** Error in " << filename << " at line " << line << ": ";
	}
	switch(option) {
	case 0:
		break;
	case 1:
		Rcpp::Rcout << fieldname << " must be 0 or 1";
		break;
	case 2:
		Rcpp::Rcout << fieldname << " must be 0, 1 or 2";
		break;
	case 3:
		Rcpp::Rcout << fieldname << " must be 0, 1, 2 or 3";
		break;
	case 4:
		Rcpp::Rcout << fieldname << " must be from 0 to 4";
		break;
	case 5:
		Rcpp::Rcout << fieldname << " must be from 0 to 5";
		break;
	case 6:
		Rcpp::Rcout << fieldname << " must be from 0 to 6";
		break;
	case 7:
		Rcpp::Rcout << fieldname << " must be from 0 to 7";
		break;
	case 10:
		Rcpp::Rcout << fieldname << " must be greater than zero";
		break;
	case 11:
		Rcpp::Rcout << fieldname << " must be 1 or more";
		break;
	case 12:
		Rcpp::Rcout << fieldname << " must be 2 or more";
		break;
	case 13:
		Rcpp::Rcout << fieldname << " must be 3 or more";
		break;
	case 18:
		Rcpp::Rcout << fieldname << " must be greater than 1.0";
		break;
	case 19:
		Rcpp::Rcout << fieldname << " must be 0 or more";
		break;
	case 20:
		Rcpp::Rcout << fieldname << " must be between 0 and 1";
		break;
	case 21:
		Rcpp::Rcout << fieldname << " must be greater than 1";
		break;
	case 33:
		Rcpp::Rcout << fieldname << " must be 1, 2 or 3";
		break;
	case 44:
		Rcpp::Rcout << fieldname << " must be from 1 to 4";
		break;
	case 55:
		Rcpp::Rcout << fieldname << " must be from 1 to 5";
		break;
	case 66:
		Rcpp::Rcout << fieldname << " must be from 1 to 6";
		break;
	case 100:
		Rcpp::Rcout << fieldname << " must be between 0 and 100";
		break;
	case 111:
		Rcpp::Rcout << fieldname << " must match the first Simulation in ParameterFile";
		break;
	case 222:
		Rcpp::Rcout << "Simulation numbers must be sequential integers";
		break;
	case 333:
		Rcpp::Rcout << "No. of " << fieldname << " columns must equal max. no. of habitats (" << maxNhab
		            << ") and be sequentially numbered starting from 1";
		break;
	case 444:
		Rcpp::Rcout << "No. of " << fieldname << " columns must be one fewer than no. of stages, i.e. " << stages - 1
		            << ", and be sequentially numbered starting from 1";
		break;
	case 555:
		Rcpp::Rcout << "No. of " << fieldname << " columns must equal no. of stages, i.e. " << stages
		            << ", and be sequentially numbered starting from 0";
		break;
	case 666:
		Rcpp::Rcout << fieldname << " must be a unique positive integer";
		break;
	default:
		Rcpp::Rcout << "*** Unspecified error regarding parameter " << fieldname;
	}
	if(option != 0)
		Rcpp::Rcout << endl;
}

void BatchErrorR(string filename,int line,int option,string fieldname,string fieldname2)
{
	if (line == -999) { // message does not cite line number
		Rcpp::Rcout << "*** Error in " << filename << ": ";
	} else {
		Rcpp::Rcout << "*** Error in " << filename << " at line " << line <<": ";
	}
	switch (option) {
	case 0:
		break;
	case 1:
		Rcpp::Rcout << fieldname << " must be greater than " << fieldname2;
		break;
	case 2:
		Rcpp::Rcout << fieldname << " must be greater than or equal to " << fieldname2;
		break;
	case 3:
		Rcpp::Rcout << fieldname << " must be less than or equal to " << fieldname2;
		break;
	case 4:
		Rcpp::Rcout << fieldname << " must be less than " << fieldname2;
		break;
	default:
		Rcpp::Rcout << "*** Unspecified error regarding parameters " << fieldname
		            << " and " << fieldname2;
		;
	}
	if (option != 0) Rcpp::Rcout << endl;
}

void ArchFormatErrorR(void)
{
//batchlog << "*** Format error in ArchFile: case-sensitive parameter names "
//	<< "must match the specification exactly" << endl;
	Rcpp::Rcout << "*** Format error in ArchFile:" << msgcase << msgmatch << endl;
}

void FormatErrorR(string filename, int errors)
{
	Rcpp::Rcout << "*** Format error in header line of ";
	if(errors == 0) {
		Rcpp::Rcout << filename << endl;
	} else {
		Rcpp::Rcout << filename << ": " << errors << " error";
		if(errors > 1)
			Rcpp::Rcout << "s";
		Rcpp::Rcout << " detected" << endl;
	}
}

void OpenErrorR(string ftype, string fname)
{
	Rcpp::Rcout << "*** Unable to open " << ftype << " " << fname << std::endl;
}

void EOFerrorR(string filename)
{
	Rcpp::Rcout << "*** Did not read to EOF in " << filename << std::endl;
}

void StreamErrorR(string filename)
{
	Rcpp::Rcout << "*** Corrupted file stream in " << filename << std::endl << "Too few entries? Unsuppoerted file encoding? (You might try to use a different one, like UTF-8.)" << std::endl;
#if RSDEBUG
	DEBUGLOG << "Corrupted file stream in " << filename << std::endl;
#endif
}

//---------------------------------------------------------------------------

// Dummy functions corresponding to those used in GUI version

/* Batch mode of v2.0 currently has no facility to save maps (unless initiated from GUI).
 */

const string Int2Str(const int x)
{
	ostringstream o;
	if(!(o << x))
		return "ERROR";
	return o.str();
}
const string Int2Str(const int x, unsigned int width)
{
	ostringstream o;
	if(!(o << std::setfill('0') << std::setw(width) << x))
		return "ERROR";
	return o.str();
}
const string Float2Str(const float x)
{
	ostringstream o;
	if(!(o << x))
		return "ERROR";
	return o.str();
}
const string Double2Str(const double x)
{
	ostringstream o;
	if(!(o << x))
		return "ERROR";
	return o.str();
}

void MemoLine(string msg)
{
	// dummy function for batch version
}

#if RSDEBUG
void DebugGUI(string msg)
{
	// dummy function for batch version
}
#endif

traitCanvas SetupTraitCanvas(void)
{
	traitCanvas tcanv;
	for(int i = 0; i < NTRAITS; i++) {
		tcanv.pcanvas[i] = 0;
	}
	return tcanv;
}

void Landscape::setLandMap(void)
{
}
void Landscape::drawLandscape(int rep, int yr, int landnum)
{
}
void Community::viewOccSuit(int year, double mn, double se)
{
}
void Community::draw(int rep, int yr, int gen, int landNum)
{
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
