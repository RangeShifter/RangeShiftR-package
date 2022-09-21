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


/*------------------------------------------------------------------------------

RangeShifter v2.0 Species

Implements the Species class

There is ONE instance of a Species for each species within the Community
AND THIS IS CURRENTLY LIMITED TO A SINGLE SPECIES.
The class holds all the demographic and dispersal parameters of the species.

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe�er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species� responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 28 July 2021 by Greta Bocedi

------------------------------------------------------------------------------*/

#ifndef SpeciesH
#define SpeciesH

#if VIRTUALECOLOGIST
#include <vector>
#endif

#include "Version.h"
#include "Parameters.h"

#if SEASONAL
// structure for seasonal extreme event parameters
struct extrmevent {
	float prob; float mort;
};
#endif // SEASONAL

// structures for demographic parameters

struct demogrParams {
	short repType;
#if SEASONAL
	short nSeasons;
#else
	short repSeasons;
#endif
#if RS_CONTAIN
	bool habDepDem;
#endif // RS_CONTAIN
	float propMales; float harem; float bc; float lambda;
	bool stageStruct;
#if BUTTERFLYDISP
	short dispersal;
#endif
#if GROUPDISP
	bool selfing; short paternity; float propLocal; float propNghbr;
#endif
};
struct stageParams {
	short nStages; short repInterval; short maxAge; short survival;
	float probRep;
#if GOBYMODEL
	float asocF;
#endif
	bool fecDens;  bool fecStageDens; bool devDens; bool devStageDens;
	bool survDens; bool survStageDens; bool disperseOnLoss;
};
struct densDepParams {
	float devCoeff; float survCoeff;
};
#if GOBYMODEL
struct socialParams {
	double socMean; double socSD; double socScale;
};
#endif
#if SOCIALMODEL
// ADDITIONAL STRUCTURES FOR PROBIS SOCIAL POLYMORPHISM MODEL
// for convenience, all social parameters are handled by a single structure
struct socialParams {
	double socMean; double socSD; double socScale;
	float asocK; float asocRmax; float asocBc;
//	float rs; float ra; float Ts; float Ta; float dK; float alpha;
//	float Ts; float Ta; float Cs; float Ca; float dK; float alpha;
	float Ts; float Ta; float cs; float ca; float bs; float ba; float dK; float alpha;
};
#endif

// structures for genetics

struct genomeData {
	int nLoci;
	bool diploid; bool neutralMarkers; bool pleiotropic; bool trait1Chromosome;
	double probMutn,probCrossover,alleleSD,mutationSD;
} ;

struct traitAllele {
	short chromo; short locus;
} ;

struct traitMap {
	short nAlleles;
	traitAllele **traitalleles;
} ;

struct traitData {
	short nTraitMaps;
	traitMap **traitmaps;
	traitMap *neutralloci;
} ;

// structures for emigration parameters

struct emigRules {
	bool densDep; bool stgDep; bool sexDep; bool indVar;
	short emigStage;
	short emigTrait[2];
#if GOBYMODEL
	float asocD;
#endif
#if GROUPDISP
	float groupmean; bool groupdisp; short grouptype;
#endif
};
struct emigTraits {
	float d0; float alpha; float beta;
};
struct emigParams {
	double d0Mean;    double d0SD;    double d0Scale;
	double alphaMean; double alphaSD; double alphaScale;
	double betaMean;  double betaSD;  double betaScale;
};
struct emigScales {
	double d0Scale; double alphaScale; double betaScale;
};

// structures for transfer parameters

struct trfrRules {
	bool moveModel; bool stgDep; bool sexDep;
	bool distMort; bool indVar;
#if RS_CONTAIN
	short kernType;
#else
	bool twinKern;
#endif // RS_CONTAIN
#if TEMPMORT
	short smType;
#else
	bool habMort;
#endif // TEMPMORT
	short moveType; bool costMap;
	short movtTrait[2];
};
struct trfrKernTraits {
	float	meanDist1; float	meanDist2; float	probKern1;
};
struct trfrMortParams {
	float fixedMort; float mortAlpha; float mortBeta;
};
struct trfrMovtTraits {
	short	pr; short	prMethod; short	memSize; short goalType;
	float	dp; float	gb; float alphaDB; int betaDB;
	float	stepMort; float	stepLength; float	rho;
	bool straigtenPath;
};
struct trfrCRWTraits {
	float	stepMort; float	stepLength; float	rho; bool straigtenPath;
};
struct trfrSMSTraits {
	short	pr; short	prMethod; short	memSize; short goalType;
	float	dp; float	gb; float alphaDB; int betaDB; float	stepMort;
	bool straigtenPath;
};
struct trfrKernParams {
	double dist1Mean; double dist1SD; double dist1Scale;
	double dist2Mean; double dist2SD; double dist2Scale;
	double PKern1Mean; double PKern1SD; double PKern1Scale;
};
struct trfrSMSParams {
	double dpMean; double dpSD; double gbMean; double gbSD;
	double alphaDBMean; double alphaDBSD; double betaDBMean; double betaDBSD;
	double dpScale; double gbScale; double alphaDBScale; double betaDBScale;
};
#if TEMPMORT
struct mortChange {
	int chgyear; double gradient;
};
#endif
struct trfrCRWParams {
	double stepLgthMean; double stepLgthSD; double stepLScale;
	double rhoMean; double rhoSD; double rhoScale;
};
struct trfrScales {
	float dist1Scale; float dist2Scale; float	PKern1Scale;
	float dpScale; float gbScale; float alphaDBScale; float betaDBScale;
	float stepLScale; float rhoScale;
};
#if RS_CONTAIN
struct trfr2Dt {
	float u0Kernel1,p0Kernel1,u0Kernel2,p0Kernel2,propKernel1;
};
struct trfrWald {
//	float mu,gamma;
	float meanU,sigma_w,hc,vt,kappa;
	float meanDirn,sdDirn;
};
#endif // RS_CONTAIN

// structures for settlement parameters

struct settleType {
	bool stgDep; bool sexDep; bool indVar;
	short settTrait[2];
#if GOBYMODEL
	float alphaSasoc; float betaSasoc;
#endif
};
struct settleRules {
	 bool densDep; bool wait; bool go2nbrLocn; bool findMate;
};
struct settleSteps {
	int minSteps; int maxSteps; int maxStepsYr;
};
struct settleTraits {
	float s0; float alpha; float beta;
};
struct settParams {
	double s0Mean;     double s0SD;     double s0Scale;
	double alphaSMean; double alphaSSD; double alphaSScale;
	double betaSMean;  double betaSSD;  double betaSScale;
};
struct settScales {
	double s0Scale; double alphaSScale; double betaSScale;
};


//---------------------------------------------------------------------------

class Species {

public:
	Species(void);
	~Species(void);
	short getSpNum(void);

#if SEASONAL
	void setBreeding(short,bool);	// set seasonal breeding indicator
	bool getBreeding(short);			// get seasonal breeding indicator
	void setExtreme(short,extrmevent);		// set seasonal extreme event
	extrmevent getExtreme(short);					// get seasonal extreme event
#if PARTMIGRN
	void setPropDispMigrn(				// set probability of dispersal/migration strategy
		short,		// migration status
		float			// probability
	);
	float getPropDispMigrn(				// get probability of dispersal/migration strategy
		short			// migration status
	);
	void setResetMigrn(bool);			// set reset dispersal/migration strategy indicator
	bool getResetMigrn(void);			// get reset dispersal/migration strategy indicator
#endif // PARTMIGRN
#endif // SEASONAL

	// demographic parameter functions

#if SEASONAL
	void createHabK( // Create habitat carrying capacity table
		short,	// no. of habitats
		short		// no. of seasons
	);
	void setHabK(
		short,	// habitat index no. (NB may differ from habitat no. supplied by user)
		short,	// season
		float		// carrying capacity (inds/cell)
	);
	float getHabK(
		short,	// habitat index no. (NB may differ from habitat no. supplied by user)
		short		// season
	);
	float getMaxK( // return highest carrying capacity over all habitats
		short		// no. of seasons
	);
#else
	void createHabK( // Create habitat carrying capacity table
		short	// no. of habitats
	);
	void setHabK(
		short,	// habitat index no. (NB may differ from habitat no. supplied by user)
		float		// carrying capacity (inds/cell)
	);
	float getHabK(
		short		// habitat index no. (NB may differ from habitat no. supplied by user)
	);
	float getMaxK(void); // return highest carrying capacity over all habitats
#endif // SEASONAL
	void deleteHabK(void); // Delete habitat carrying capacity table
	void setStage( // Set stage structure parameters
		const stageParams	// structure holding stage structure parameters
	);
	stageParams getStage(void); // Get stage structure parameters
	void setDemogr( // Set general demographic parameters
		const demogrParams	// structure holding general demographic parameters
	);
	demogrParams getDemogr(void); // Get general demographic parameters
	short getRepType(void);
	bool stageStructured(void);
	void setDensDep( // Set demographic density dependence coefficients
		float,	// development coefficient
		float		// survival coefficient
	);
	densDepParams getDensDep(void); // Get development and survival coefficients

#if RS_CONTAIN

	void resetDem( // Reset demographic rates to zero for a specified habitat
		short		// habitat index (if < 0 then all habitats are reset)
	);
#if SEASONAL
	void setFec( // Set fecundity
		short,	// habitat index (must be >= 0)
		short,	// season (must be >= 0)
		short,	// stage (must be > 0)
		short,	// sex
		float		// fecundity
	);
	float getFec( // Get fecundity
		short,	// habitat index (must be >= 0)
		short,	// season (must be >= 0)
		short,	// stage
		short		// sex
	);
	void setDev( // Set development probability
		short,	// habitat index (must be >= 0)
		short,	// season (must be >= 0)
		short,	// stage
		short,	// sex
		float		// development probability
	);
	float getDev( // Get development probability
		short,	// habitat index (must be >= 0)
		short,	// season (must be >= 0)
		short,	// stage
		short		// sex
	);
	void setSurv( // Set survival probability
		short,	// habitat index (must be >= 0)
		short,	// season (must be >= 0)
		short,	// stage
		short,	// sex
		float		// survival probability
	);
	float getSurv( // Get survival probability
		short,	// habitat index (must be >= 0)
		short,	// season (must be >= 0)
		short,	// stage
		short		// sex
	);
#else
	void setFec( // Set fecundity
		short,	// habitat index (must be >= 0)
		short,	// stage (must be > 0)
		short,	// sex
		float		// fecundity
	);
	float getFec( // Get fecundity
		short,	// habitat index (must be >= 0)
		short,	// stage
		short		// sex
	);
	void setDev( // Set development probability
		short,	// habitat index (must be >= 0)
		short,	// stage
		short,	// sex
		float		// development probability
	);
	float getDev( // Get development probability
		short,	// habitat index (must be >= 0)
		short,	// stage
		short		// sex
	);
	void setSurv( // Set survival probability
		short,	// habitat index (must be >= 0)
		short,	// stage
		short,	// sex
		float		// survival probability
	);
	float getSurv( // Get survival probability
		short,	// habitat index (must be >= 0)
		short,	// stage
		short		// sex
	);
#endif // SEASONAL

#else

#if SEASONAL
	void setFec( // Set fecundity
		short,	// season (must be >= 0)
		short,	// stage (must be > 0)
		short,	// sex
		float		// fecundity
	);
	float getFec( // Get fecundity
		short,	// season (must be >= 0)
		short,	// stage
		short		// sex
	);
	void setDev( // Set development probability
		short,	// season (must be >= 0)
		short,	// stage
		short,	// sex
		float		// development probability
	);
	float getDev( // Get development probability
		short,	// season (must be >= 0)
		short,	// stage
		short		// sex
	);
	void setSurv( // Set survival probability
		short,	// season (must be >= 0)
		short,	// stage
		short,	// sex
		float		// survival probability
	);
	float getSurv( // Get survival probability
		short,	// season (must be >= 0)
		short,	// stage
		short		// sex
	);

#else

	void setFec( // Set fecundity
		short,	// stage (must be > 0)
		short,	// sex
		float		// fecundity
	);
	float getFec( // Get fecundity
		short,	// stage
		short		// sex
	);
	void setDev( // Set development probability
		short,	// stage
		short,	// sex
		float		// development probability
	);
	float getDev( // Get development probability
		short,	// stage
		short		// sex
	);
	void setSurv( // Set survival probability
		short,	// stage
		short,	// sex
		float		// survival probability
	);
	float getSurv( // Get survival probability
		short,	// stage
		short		// sex
	);
#endif // SEASONAL
#endif // RS_CONTAIN

#if SPATIALDEMOG
	void setFecSpatial(bool);
	bool getFecSpatial(void){return fecSpatial;};
	void setDevSpatial(bool);
	bool getDevSpatial(void){return devSpatial;};
	void setSurvSpatial(bool);
	bool getSurvSpatial(void){return survSpatial;};
	void setFecLayer( // set the layer of the spatial demographic scaling used for fecundity
		short, // stage
		short, // sex
		short // layer
	);

	short getFecLayer( // get the layer of the spatial demographic scaling used for fecundity
		short, // stage
		short // sex
	);

	void setDevLayer( // set the layer of the spatial demographic scaling used for development
		short, // stage
		short, // sex
		short // layer
	);

	short getDevLayer( // get the layer of the spatial demographic scaling used for development
		short, // stage
		short // sex
	);

	void setSurvLayer( // set the layer of the spatial demographic scaling used for survival
		short, // stage
		short, // sex
		short // layer
	);

	short getSurvLayer( // get the layer of the spatial demographic scaling used for survival
		short, // stage
		short // sex
	);
#endif // SPATIALDEMOG

	float getMaxFec(void); // Get highest fecundity of any stage
	void setMinAge( // Set minimum age
		short,	// stage
		short,	// sex
		int			// minimum age (years) (must be zero for stages 0 and 1)
	);
	short getMinAge( // Get minimum age
		short,	// stage
		short		// sex
	);
	void createDDwtFec( // Create fecundity weights matrix
		short		// matrix dimension - no. of stages * no. of sexes
	);
	void setDDwtFec( // Set fecundity weights matrix element
		short,	// row
		short,	// column
		float		// weight
	);
	float getDDwtFec( // Get fecundity weights matrix element
		short,	// row
		short		// column
	);
	void deleteDDwtFec(void); // Delete fecundity weights matrix
	void createDDwtDev( // Create development weights matrix
		short		// matrix dimension - no. of stages * no. of sexes
	);
	void setDDwtDev( // Set development weights matrix element
		short,	// row
		short,	// column
		float		// weight
	);
	float getDDwtDev( // Get development weights matrix element
		short,	// row
		short		// column
	);
	void deleteDDwtDev(void); // Delete development weights matrix
	void createDDwtSurv( // Create survival weights matrix
		short		// matrix dimension - no. of stages * no. of sexes
	);
	void setDDwtSurv( // Set survival weights matrix element
		short,	// row
		short,	// column
		float		// weight
	);
	float getDDwtSurv( // Get survival weights matrix element
		short,	// row
		short		// column
	);
	void deleteDDwtSurv(void); // Delete survival weights matrix
	// Functions to handle min/max R or K (under environmental stochasticity)
	void setMinMax( // Set min and max values
		float,	// min
		float		// max
	);
	float getMinMax( // Get min/max value
		short		// option: 0 = return minimum, otherwise = return maximum
	);
#if GOBYMODEL
	void setSocialParams( // Set social phenotype initialisation parameters
//		const short,			// stage (NB implemented for stage 0 only)
//		const short,			// sex
		const socialParams	// structure holding parameters
	);
	socialParams getSocialParams( // Get social phenotype initialisation parameters
//		short,	// stage (NB implemented for stage 0 only)
//		short		// sex
		void
	);
#endif
#if SOCIALMODEL
	// ADDITIONAL FUNCTIONS FOR PROBIS SOCIAL POLYMORPHISM MODEL
	void setSocialParams( // Set social phenotype initialisation parameters
//		const short,			// stage (NB implemented for stage 0 only)
//		const short,			// sex
		const socialParams	// structure holding parameters
	);
	socialParams getSocialParams( // Get social phenotype initialisation parameters
//		short,	// stage (NB implemented for stage 0 only)
//		short		// sex
		void
	);
#endif


	// genome functions

	void setGenomeData(genomeData);
	genomeData getGenomeData(void);
	bool isDiploid(void);
	void setNChromosomes( // Set no. of chromosomes
		int // no. of chromosomes
	);
	int getNChromosomes(void);
	void setNLoci(
		const short,	// chromosome no.
		const short		// locus no.
	);
	int getNLoci(
		const short		// chromosome no.
	);
	void deleteLoci(void);
	void set1ChromPerTrait( // Set 1:1 mapping of trait to chromosome
		const int	// no. of loci on each chromosome
	);
	bool has1ChromPerTrait(void);
	void setTraits(void); // Set trait attributes for the species
	void setTraitNames(void);
	void deleteTraitNames(void);
	string getTraitName(
		const int	// trait no.
	);
	int getNTraits(void);
	void setTraitData(
		const int	// no. of traits
	);
	void deleteTraitData(void);
	int getNTraitMaps(void);
	void setTraitMap(
		const short,	// trait no.
		const short		// no. of alleles
	);
	int getNTraitAlleles(
		const int		// trait no.
	);
	void setTraitAllele(
		const short,	// trait no.
		const short,	// trait allele no.
		const short,	// chromosome no.
		const short		// chromosome locus no.
	);
	traitAllele getTraitAllele(
		const short,	// trait no.
		const short		// trait allele no.
	);
	void setNeutralLoci(bool);
	void deleteNeutralLoci(void);
	int getNNeutralLoci(void);
	traitAllele getNeutralAllele(
		const short		// allele no.
	);
#if VIRTUALECOLOGIST
	// functions for sampling genome
	bool sampleAllLoci(void);
	void setSampleAllLoci(bool);
	void resetSampleLoci(void);
	void addSampleLocus(
		const short,	// chromosome no.
		const short		// locus no.
	);
	int nSampleLoci(void);
	traitAllele getSampleLocus(
		const short	// sample no.
	);
#endif

	// emigration parameter functions

	void setEmig( // Set emigration rules
		const emigRules	// structure holding emigration rules
	);
	emigRules getEmig(void); // Get emigration rules
	void setEmigTraits( // Set emigration trait parameters
		const short,			// stage
		const short,			// sex
		const emigTraits	// structure holding emigration trait parameters
	);
	emigTraits getEmigTraits( // Get emigration trait parameters
		short,	// stage
		short		// sex
	);
	float getEmigD0( // Get (maximum) emigration probability
		short,	// stage
		short		// sex
	);
	void setEmigParams( // Set emigration initialisation parameters
		const short,			// stage (NB implemented for stage 0 only)
		const short,			// sex
		const emigParams	// structure holding parameters
	);
	emigParams getEmigParams( // Get emigration initialisation parameters
		short,	// stage (NB implemented for stage 0 only)
		short		// sex
	);
	void setEmigScales( // Set emigration mutation parameters
		const emigScales	// structure holding emigration mutation parameters
	);
	emigScales getEmigScales(void); // Get emigration mutation parameters

	// transfer parameter functions

	void setTrfr( // Set transfer rules
		const trfrRules	// structure holding transfer rules
	);
	trfrRules getTrfr(void); // Get transfer rules
	void setFullKernel( // Set fullKernel condition
		bool						// fullKernel value
	);
	bool useFullKernel(void);
	void setKernTraits( // Set transfer by kernel parameters
		const short,					// stage
		const short,					// sex
		const trfrKernTraits,	// structure holding transfer by kernel parameters
		const int							// Landscape resolution
	);
	trfrKernTraits getKernTraits( // Get transfer by kernel parameters
		short,	// stage
		short		// sex
	);
	void setMortParams( // Set transfer mortality parameters
		const trfrMortParams	// structure holding transfer mortality parameters
	);
	trfrMortParams getMortParams(void); // Get transfer mortality parameters
	void setMovtTraits( // Set transfer movement model parameters
		const trfrMovtTraits	// structure holding transfer movement model parameters
	);
	trfrMovtTraits getMovtTraits(void); // Get transfer movement model traits
	trfrCRWTraits getCRWTraits(void);		// Get CRW traits
	trfrSMSTraits getSMSTraits(void);		// Get SMS traits
	void setKernParams( // Set initial transfer by kernel parameter limits
		const short,					// stage (NB implemented for stage 0 only)
		const short,					// sex
		const trfrKernParams, // structure holding min and max values
		const double						// Landscape resolution
	);
	trfrKernParams getKernParams( // Get initial transfer by kernel parameter limits
		short,	// stage (NB implemented for stage 0 only)
		short		// sex
	);
	void setSMSParams( // Set initial transfer by SMS parameter limits
		const short,					// stage (NB implemented for stage 0 only)
		const short,					// sex   (NB implemented for sex 0 only)
		const trfrSMSParams		// structure holding min and max values
	);
	trfrSMSParams getSMSParams( // Get initial transfer by SMS parameter limits
		short,	// stage (NB implemented for stage 0 only)
		short		// sex   (NB implemented for sex 0 only)
	);
#if TEMPMORT
	void clearMortalities(void);
	void addMortChange(
		int,		// year of change
		double	// mortality change gradient
	);
	void updateMortality(int);
	double getMortality(void);
#endif // TEMPMORT
	void setCRWParams( // Set initial transfer by CRW parameter limits
		const short,					// stage (NB implemented for stage 0 only)
		const short,					// sex   (NB implemented for sex 0 only)
		const trfrCRWParams		// structure holding min and max values
	);
	trfrCRWParams getCRWParams( // Get initial transfer by CRW parameter limits
		short,	// stage (NB implemented for stage 0 only)
		short		// sex   (NB implemented for sex 0 only)
	);
	void setTrfrScales( // Set transfer mutation parameters
		const trfrScales	// structure holding transfer mutation parameters
	);
	trfrScales getTrfrScales(void); // Get transfer mutation parameters
	// Return dimension of habitat-dependent step mortality and costs matrices
	short getMovtHabDim(void);
	void createHabCostMort( // Create habitat-dependent costs and mortality matrices
		short	// no. of habitats
	);
	void setHabCost( // Set habitat-dependent cost
		short,	// habitat index no.
		int			// cost value
	);
	void setHabMort( // Set habitat-dependent per-step mortality
		short,	// habitat index no.
		double	// mortality rate
	);
	int getHabCost( // Get habitat-dependent cost
		short		// habitat index no.
	);
	double getHabMort( // Get habitat-dependent per-step mortality
		short		// habitat index no.
	);
	void deleteHabCostMort(void); // Delete habitat-dependent costs and mortality matrices
#if RS_CONTAIN
	void setTrfr2Dt(trfr2Dt);
	trfr2Dt getTrfr2Dt(void);
	void setTrfrWald(trfrWald);
	trfrWald getTrfrWald(void);
	void setTrfrHr(float,unsigned short);
	float getTrfrHr(unsigned short);
#endif // RS_CONTAIN

	// settlement parameter functions

	void setSettle( // Set settlement type
		const settleType	// structure holding settlement type (stage- and/or sex-dependent)
	);
	settleType getSettle(void); // Get settlement type
	void setSettRules( // Set settlement rules
		const short,			// stage
		const short,			// sex
		const settleRules	// structure holding settlement rules
	);
	settleRules getSettRules( // Get settlement rules
		short,	// stage
		short		// sex
	);
	void setSteps( // Set path step limit parameters
		const short,			// stage
		const short,			// sex
		const settleSteps	// structure holding path step limit parameters
	);
	settleSteps getSteps( // Set path step limit parameters
		short,	// stage
		short		// sex
	);
	void setSettTraits( // Set settlement density dependence traits
		const short,					// stage
		const short,					// sex
		const settleTraits	// structure holding density dependence traits
	);
	settleTraits getSettTraits( // Get settlement density dependence traits
		short,	// stage
		short		// sex
	);
	void setSettParams( // Set settlement initialisation parameters
		const short,			// stage (NB implemented for stage 0 only)
		const short,			// sex
		const settParams	// structure holding parameters
	);
	settParams getSettParams( // Get settlement initialisation parameters
		short,	// stage (NB implemented for stage 0 only)
		short		// sex
	);
	void setSettScales( // Set settlement mutation parameters
		const settScales	// structure holding settlement mutation parameters
	);
	settScales getSettScales(void); // Get settlement mutation parameters

private:

	// NOTE: SEQUENCE OF PARAMETER VARIABLES MAY NEED TO BE ALTERED FOR EFFICIENTCY ...
	// ... but that is of low importance, as there will only be one (or a few) instance(s)

	// demographic parameters

	short repType;			// 0 = asexual, 1 = simple two sex, 2 = complex two sex
#if GROUPDISP
											// 3 = hermaphrodite
#endif
#if BUTTERFLYDISP
	short dispersal; 		// dispersal timing: 0 = during reprodn (i.e. prior to parturition),
											// 1 = after reprodn
#endif
	short nStages;      // no. of stages (incl. juveniles) in structured population
	float propMales;    // proportion of males at birth in sexual model
	float harem;        // max harem size in complex sexual model
	float bc;						// competition coefficient for non-structured population
	float lambda;       // max intrinsic growth rate for non-structured population
	float probRep; 			// probability of reproducing in subsequent seasons
#if GOBYMODEL
	float asocF;				// fecundity adjustment parameter for asocial phenotype
#endif
#if SEASONAL
	short nSeasons;			// no. of seasons per year
#else
	short repSeasons;		// no. of reproductive seasons per year
#endif
	short repInterval;	// no. of reproductive seasons between subsequent reproductions
	short maxAge;       // max age in structured population
	short survival;			// survival timing: 0 = at reprodn, 1 = between reprodns, 2 = anually
	bool stageStruct;
	bool fecDens;
	bool fecStageDens;
	bool devDens;
	bool devStageDens;
	bool survDens;
	bool survStageDens;
	bool disperseOnLoss;	// individuals disperse on complete loss of patch
												// (otherwise they die)
#if SPATIALDEMOG
	bool fecSpatial, devSpatial, survSpatial;
#endif //SPATIALDEMOG
#if GROUPDISP
	bool selfing;       // self-fertilisation possible
	short paternity;		// 0 = fixed (one father per breeding attempt), 1 = assigned at random
											// 2 = pollen ratio - local : neighbour : global
	float propLocal;    // proportion of local pollen (from same population)
	float propNghbr;    // proportion of pollen from neighbouring populations
#endif
	short habDimK;			// dimension of carrying capacities matrix
#if SEASONAL
	float **habK;				// seasonal habitat-specific carrying capacities (inds/cell)
#else
	float *habK;				// habitat-specific carrying capacities (inds/cell)
#endif // SEASONAL
	float devCoeff; 		// density-dependent development coefficient
	float survCoeff; 		// density-dependent survival coefficient
	float **ddwtFec;    // density-dependent weights matrix for fecundity
	float **ddwtDev;    // density-dependent weights matrix for development
	float **ddwtSurv;   // density-dependent weights matrix for survival
	// NB for the following arrays, sex 0 is females, sex 1 is males
#if RS_CONTAIN
#if SEASONAL
	float fec[NHABITATS][NSEASONS][NSTAGES][NSEXES];		// fecundities
	float dev[NHABITATS][NSEASONS][NSTAGES][NSEXES];		// development probabilities
	float surv[NHABITATS][NSEASONS][NSTAGES][NSEXES];		// survival probabilities
#else
	float fec[NHABITATS][NSTAGES][NSEXES];		// fecundities
	float dev[NHABITATS][NSTAGES][NSEXES];		// development probabilities
	float surv[NHABITATS][NSTAGES][NSEXES];		// survival probabilities
#endif // SEASONAL
	bool habDepDem;														// habitat-dependent demography
#else
#if SEASONAL
	float fec[NSEASONS][NSTAGES][NSEXES];			// fecundities
	float dev[NSEASONS][NSTAGES][NSEXES];			// development probabilities
	float surv[NSEASONS][NSTAGES][NSEXES];		// survival probabilities
#else
	float fec[NSTAGES][NSEXES];			// fecundities
	float dev[NSTAGES][NSEXES];			// development probabilities
	float surv[NSTAGES][NSEXES];		// survival probabilities
#endif // SEASONAL
#endif // RS_CONTAIN
#if SPATIALDEMOG
	int fecLayer[NSTAGES][NSEXES]; // layer for spatial varying fecundity
	int devLayer[NSTAGES][NSEXES]; // layer for spatial varying development
	int survLayer[NSTAGES][NSEXES]; // layer for spatial varying survival
#endif // SPATIALDEMOG
	short minAge[NSTAGES][NSEXES];	// minimum age to enter stage
	// NOTE - IN THEORY, NEXT 3 VARIABLES COULD BE COMMON, BUT WE WOULD NEED TO ENSURE THAT
	// ALL MATRICES ARE DELETED IF THERE IS A CHANGE IN NO. OF STAGES OR REPRODUCTION TYPE
	// ***** TO BE RECONSIDERED LATER *****
	short ddwtFecDim;		// dimension of density-dependent weights matrix for fecundity
	short ddwtDevDim;		// dimension of density-dependent weights matrix for fecundity
	short ddwtSurvDim;	// dimension of density-dependent weights matrix for fecundity
	float minRK; 				// minimum ) growth rate OR carrying capacity
	float maxRK; 				// maximum ) (under environmental stochasticity)
#if GOBYMODEL
	double socMean;
	double socSD;
	double socScale;
#endif
#if SOCIALMODEL
	double socMean;
	double socSD;
	double socScale;
#endif

	// genome parameters

	short nTraits;						// no. of inheritable traits
	short emigTrait[2];				// to record first and no. of emigration traits
	short movtTrait[2];				// to record first and no. of transfer traits
	short settTrait[2];				// to record first and no. of settlement traits
	bool diploid;
	bool neutralMarkers;			// neutral markers in absence of any adaptive traits
	bool pleiotropic;
	bool trait1Chromosome;		// 1:1 mapping of chromosome to trait
	short nChromosomes;				// no. of chromosomes
	double probMutn;					// allelic mutation probability
	double probCrossover; 		// crossover probability at meiosis
	double alleleSD;					// s.d. of initial allelic values around phenotypic value
	double mutationSD;				// s.d. of mutation magnitude
	short nNLoci;							// no. of nLoci set
	short *nLoci;							// no. of loci per chromosome
	short nTraitNames;				// no. of trait names set
	traitData *traitdata;			// for mapping of chromosome loci to traits
	string *traitnames;				// trait names for parameter output
#if VIRTUALECOLOGIST
	// parameters for sampling genome
	bool sampleAll;
	std::vector <traitAllele> samples;
#endif

	// emigration parameters

	bool	densDepEmig;	// density-dependent emigration
	bool	stgDepEmig;   // stage-dependent emigration
	bool	sexDepEmig;   // sex-dependent emigration
	bool	indVarEmig;   // individual variation in emigration
	short emigStage;		// stage which emigrates (used for stage-strucutred population
											// having individual variability in emigration probability)
	// NB for the following arrays, sex 0 is females, sex 1 is males
	float	d0[NSTAGES][NSEXES];				 // maximum emigration probability
	float	alphaEmig[NSTAGES][NSEXES];	 // slope of density-dependent reaction norm
	float	betaEmig[NSTAGES][NSEXES];	 // inflection point of reaction norm (in terms of N/K)
#if GOBYMODEL
	float asocD;				// emigration adjustment parameter for asocial phenotype
#endif // GOBYMODEL
#if GROUPDISP
	float groupmean;    // mean group size (> 1)
	bool groupdisp;     // group dispersal occurs
	short grouptype;    // 0 = population-level grouping, 1 = sibling-level grouping
											// NB OPTION 1 IS NOT CURRENTLY IMPLEMENTED
#endif // GROUPDISP
	// NB Initialisation parameters are made double to avoid conversion errors (reason unclear)
	// on traits maps using FloatToStr()
	// As evolving traits are not stage-dependent, no. of rows can be 1
	// Indeed, they could be 1-D arrays
	double d0Mean[1][NSEXES];
	double d0SD[1][NSEXES];
	double alphaMean[1][NSEXES];
	double alphaSD[1][NSEXES];
	double betaMean[1][NSEXES];
	double betaSD[1][NSEXES];
	double d0Scale;											// scaling factor for d0
	double alphaScale;									// scaling factor for alpha
	double betaScale;										// scaling factor for beta

	// transfer parameters

	bool moveModel;
	bool stgDepTrfr;
	bool sexDepTrfr;
	bool distMort;
	bool indVarTrfr;
#if RS_CONTAIN
	short kernType;	// 0 = single negative exponential, 1 = double negative exponential
									// 2 = 2Dt, 3 = WALD (inverse Gaussian)
#else
	bool twinKern;
#endif // RS_CONTAIN
#if TEMPMORT
	short smType;		// per-step mortality type: 0 = constant, 1 = habitat-dependent
									// 2 = temporally variable
#else
	bool habMort;		// habitat-dependent mortality
#endif // TEMPMORT
	float	meanDist1[NSTAGES][NSEXES];	// mean of 1st dispersal kernel (m)
	float	meanDist2[NSTAGES][NSEXES]; // mean of 2nd dispersal kernel (m)
	float	probKern1[NSTAGES][NSEXES]; // probability of dispersing with the 1st kernel
	// NB INITIAL limits are made double to avoid conversion errors (reason unclear)
	// on traits maps using FloatToStr()
	// As evolving traits are are not stage-dependent, no. of rows can be 1
	// Indeed, as they are INITIAL limits, which may subsequently be exceeded, they could be
	// 1-D arrays
	double dist1Mean[1][NSEXES];	// mean of initial mean of the 1st dispersal kernel (m)
	double dist1SD[1][NSEXES]; 		// s.d. of initial mean of the 1st dispersal kernel (m)
	double dist2Mean[1][NSEXES]; 	// mean of initial mean of the 2nd dispersal kernel (m)
	double dist2SD[1][NSEXES]; 		// s.d. of initial mean of the 2nd dispersal kernel (m)
	double PKern1Mean[1][NSEXES];	// mean of initial prob. of dispersing with 1st kernel
	double PKern1SD[1][NSEXES];	  // s.d. of initial prob. of dispersing with 1st kernel
	float dist1Scale;		// scaling factor for mean of 1st dispersal kernel (m)
	float dist2Scale;		// scaling factor for mean of 2nd dispersal kernel (m)
	float PKern1Scale;	// scaling factor for prob. of dispersing with 1st kernel
	float fixedMort;		// constant mortality probability
	float mortAlpha;		// slope for mortality distance dependence function
	float mortBeta;			// inflection point for mortality distance dependence function
#if RS_CONTAIN
	// parameters for 2Dt kernel
	float u0Kernel1;		// scaling parameter for 1st kernel
	float p0Kernel1;		// shape parameter for 1st kernel
	float u0Kernel2;		// scaling parameter for 2nd kernel
	float p0Kernel2;		// shape parameter for 2nd kernel
	float propKernel1;	// probability of dispersing with the 1st kernel
	// parameters for WALD (inverse Gaussian) kernel
//	float mu;						// scaling parameter
//	float gamma;				// shape parameter
	float meanU;				// mean horizontal wind speed (m/s)
	float sigma_w;			// s.d. of vertical wind speed (m/s)
	float hc;						// canopy height (m)
	float hr[NSTAGES];	// stage-dependent seed release height (m)
	float vt;						// seed terminal velocity (m/s)
	float kappa;				// turbulence coefficient
	float meanDirn;			// mean wind direction (degrees)
	float sdDirn;				// s.d. of wind direction (degrees)
#endif // RS_CONTAIN
	short moveType; 		// 1 = SMS, 2 = CRW
	short pr;						// SMS perceptual range (cells)
	short prMethod;			// SMS perceptual range evaluation method:
											// 1 = arith. mean, 2 = harmonic mean, 3 = inverse weighted arith. mean
	short memSize;			// SMS memory size (1-14 steps)
	short goalType;			// SMS goal bias type: 0 = none, 1 = towards goal, 2 = dispersal bias
	float dp;						// SMS directional persistence
	float gb;						// SMS goal bias
	float alphaDB; 			// SMS dispersal bias decay rate
	int betaDB; 				// SMS dispersal bias decay inflection point (no. of steps)
	float stepMort;			// constant per-step mortality probability for movement models
	double *habStepMort;	// habitat-dependent per-step mortality probability
	float stepLength;		// CRW step length (m)
	float rho;					// CRW correlation coefficient
	double dpMean[1][NSEXES];				// mean of initial SMS directional persistence
	double dpSD[1][NSEXES];	 				// s.d. of initial SMS directional persistence
	double gbMean[1][NSEXES];				// mean of initial SMS goal bias
	double gbSD[1][NSEXES];	 				// s.d. of initial SMS goal bias
	double alphaDBMean[1][NSEXES];	// mean of initial SMS dispersal bias decay rate
	double alphaDBSD[1][NSEXES];	 	// s.d. of initial SMS dispersal bias decay rate
	double betaDBMean[1][NSEXES];		// mean of initial SMS dispersal bias decay infl. pt.
	double betaDBSD[1][NSEXES];	 		// s.d. of initial SMS dispersal bias decay infl. pt.
	float dpScale;									// scaling factor for SMS directional persistence
	float gbScale;									// scaling factor for SMS goal bias
	float alphaDBScale;							// scaling factor for SMS dispersal bias decay rate
	float betaDBScale;							// scaling factor for SMS dispersal bias decay infl. pt.
#if TEMPMORT
	// list of dynamic per-step mortality changes
	std::vector <mortChange> mortchanges;
	int nextChange,nextYear;
	double currentMortality,currentGradient,nextGradient;
#endif // TEMPMORT
	double stepLgthMean[1][NSEXES];	// mean of initial step length (m)
	double stepLgthSD[1][NSEXES];		// s.d. of initial step length (m)
	double rhoMean[1][NSEXES];			// mean of initial correlation coefficient
	double rhoSD[1][NSEXES];	 			// s.d. of initial correlation coefficient
	float stepLScale;								// scaling factor for step length (m)
	float rhoScale;									// scaling factor for correlation coefficient
	short habDimTrfr;		// dimension of habitat-dependent step mortality and costs matrices
	int *habCost;				// habitat costs
	bool costMap;				// import cost map from file?
	bool straigtenPath;	// straighten path after decision not to settle
	bool fullKernel;		// used to indicate special case when density-independent emigration
											// is 1.0, and kernel-based movement within the natal cell is used
											// to determine philopatry

	// settlement parameters

	bool stgDepSett;
	bool sexDepSett;
	bool indVarSett;   								// individual variation in settlement
	bool densDepSett[NSTAGES][NSEXES];
	bool wait[NSTAGES][NSEXES];				// wait to continue moving next season (stage-structured model only)
	bool go2nbrLocn[NSTAGES][NSEXES];	// settle in neighbouring cell/patch if available (ditto)
	bool findMate[NSTAGES][NSEXES];
	int minSteps;     								// minimum no. of steps
	int maxSteps;											// maximum total no. of steps
	int maxStepsYr[NSTAGES][NSEXES]; 	// maximum no. of steps in any one dispersal period
	float	s0[NSTAGES][NSEXES];				// maximum settlement probability
	float alphaS[NSTAGES][NSEXES];		// slope of the settlement reaction norm to density
	float betaS[NSTAGES][NSEXES];			// inflection point of the settlement reaction norm to density
#if GOBYMODEL
	float alphaSasoc;									// alphaS adjustment parameter for asocial phenotype
	float betaSasoc;									// betaS adjustment parameter for asocial phenotype
#endif
	double s0Mean[1][NSEXES];					// mean of initial maximum settlement probability
	double s0SD[1][NSEXES]; 					// s.d. of initial maximum settlement probability
	double alphaSMean[1][NSEXES];			// mean of initial settlement reaction norm slope
	double alphaSSD[1][NSEXES]; 	 		// s.d. of initial settlement reaction norm slope
	double betaSMean[1][NSEXES]; 			// mean of initial settlement reaction norm inflection point
	double betaSSD[1][NSEXES];	  		// s.d. of initial settlement reaction norm inflection point
	float s0Scale;										// scaling factor for maximum settlement probability
	float alphaSScale;								// scaling factor for settlement reaction norm slope
	float betaSScale;									// scaling factor for settlement reaction norm inflection point

	// other attributes

	int spNum;

#if SEASONAL
	std::vector <bool> breeding;			// seasonal breeding
	std::vector <extrmevent> extreme;	// seasonal random extreme events
#if PARTMIGRN
	float propDispMigrn[7];						// probabilities of dispersal/migration strategies
	bool resetMigrn;									// reset dispersal/migration strategies every year
#endif // PARTMIGRN
#endif // SEASONAL

#if SOCIALMODEL
	// ADDITIONAL PARAMETERS FOR PROBIS SOCIAL POLYMORPHISM MODEL
	// see also Fogarty et al. (2011), Am. Nat., 177, 273-287

	float asocK;			// ratio of K for asocial morph to K for social morph
	float asocRmax;		// ratio of Rmax for asocial morph to Rmax for social morph
	float asocBc;			// ratio of bc for asocial morph to bc for social morph
//	float rs,ra; 			// Fogarty's r for social and asocial morphs
	float Ts,Ta; 			// Allee effect thresholds (relative to K) for social and asocial morphs
	float cs,ca; 			// Parameters of Allee effect below Ts,Ta
	float bs,ba; 			// Parameters of Allee effect below Ts,Ta
	float dK;					// dispersal rate independent of reproduction expectation
	float alpha;			// change in dispersal rate with reproduction expectation
#endif

};

/* IMPORTANT NOTE:
At the time of writing (24/4/14) the stage- and sex-dependent parameters for emigration
and dispersal (e.g. d0[NSTAGES][NSEXES]) are set according to the appropriate stage- and
sex-dependent settings; thus a species could be structured and sexual, but parameter values
are set for elements [0][0] only if emigration/transfer is stage- and sex-independent.
However, the parameters for settlement are set for ALL stages and sexes appropriate to the
species, regardless of settlement dependency. The rationale for this is that settlement
parameters need to be accessed many more times for a movement model (at each step) than
emigration or transfer parameters, and therefore there will be a performance gain in
avoiding nested if statements in Individual::moveStep() which would otherwise be required
to get the correct parameter values from the settlement arrays. Whether that particular
rationale is justified remains to be tested!
*/

//---------------------------------------------------------------------------

#if RSDEBUG
//extern ofstream DEBUGLOG;
extern void DebugGUI(string);
#endif

//---------------------------------------------------------------------------
#endif
