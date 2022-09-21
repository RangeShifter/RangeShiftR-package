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

#include "Genome.h"
//---------------------------------------------------------------------------

#if RS_EMBARCADERO 
#pragma package(smart_init)
#endif

ofstream outGenetic;

//---------------------------------------------------------------------------

#ifdef RSDEBUG
//std::ofstream chromdebug("chromdebug.txt");
//std::ofstream chromdebug1("chromdebug1.txt");
#endif

//---------------------------------------------------------------------------

Chromosome::Chromosome(int nloc)
{
#if RSDEBUG
//DEBUGLOG << "Chromosome::Chromosome(): this=" << this
//	<< " nloc=" << nloc
////	<< " maternalLoci.size()=" << maternalLoci.size()
////	<< " paternalLoci.size()=" << paternalLoci.size()
//	<< endl;
//DebugGUI("Chromosome::Chromosome(): this=" + Int2Str((int)this)
//	+ " nloc=" + Int2Str(nloc)
////	+ " maternalLoci.size()=" + Int2Str((int)maternalLoci.size())
////	+ " paternalLoci.size()=" + Int2Str((int)paternalLoci.size())
//	);
#endif
if (nloc > 0) nloci = nloc; else nloci = 1;
pLoci = new locus[nloci];
for (int i = 0; i < nloci; i++) {
	pLoci[i].allele[0] = pLoci[i].allele[1] = 0;
}
}

Chromosome::~Chromosome() {
#if RSDEBUG
//DEBUGLOG << "Chromosome::~Chromosome(): this=" << this << endl;
#endif
#if RSDEBUG
//DEBUGLOG << "Chromosome::Chromosome(): deleting this=" << this << endl;
//DebugGUI("Chromosome::Chromosome(): deleting this=" + Int2Str((int)this));
#endif
if (pLoci != 0) {
//	for (int i = 0; i < nloci; i++) {
//		delete pLoci[i]; pLoci[i] = NULL;
//	}
	delete[] pLoci; pLoci = NULL;
}
}

short Chromosome::nLoci(void) { return nloci; }

locus Chromosome::alleles(const int loc) { // return allele values at a specified locus
locus l; l.allele[0] = l.allele[1] = 0;
if (loc >= 0 && loc < nloci) {
	l.allele[0] = pLoci[loc].allele[0]; l.allele[1] = pLoci[loc].allele[1];
}
return l;
}

double Chromosome::additive(const bool diploid) {
int sum = 0;
for (int i = 0; i < nloci; i++) {
	sum += pLoci[i].allele[0];
	if (diploid) sum += pLoci[i].allele[1];
}
#if RSDEBUG
//DEBUGLOG << "Chromosome::additive(): this=" << this
//	<< " sum=" << sum
//	<< endl;
#endif
return (double)sum / INTBASE;
}

double Chromosome::meanvalue(const bool diploid) {
int sum = 0;
double mean;
for (int i = 0; i < nloci; i++) {
	sum += pLoci[i].allele[0];
	if (diploid) sum += pLoci[i].allele[1];
}
mean = (double)sum / (double)nloci;
if (diploid) mean /= 2.0;
mean /= INTBASE;
#if RSDEBUG
//DEBUGLOG << "Chromosome::meanvalue(): this=" << this
//	<< " sum=" << sum
//	<< " mean=" << mean
//	<< endl;
#endif
return  mean;
}

double Chromosome::additive(const short loc,const bool diploid) {
int sum = 0;
sum += pLoci[loc].allele[0];
if (diploid) sum += pLoci[loc].allele[1];
#if RSDEBUG
//DEBUGLOG << "Chromosome::additive(): this=" << this
//	<< " sum=" << sum
//	<< endl;
#endif
return (double)sum / INTBASE;
}

/*
double Chromosome::probval(const bool diploid) {
double genval,phenval;
int sum = 0;
for (int i = 0; i < nloci; i++) {
	sum += pLoci[i].allele[0];
	if (diploid) sum += pLoci[i].allele[1];
}
genval = (double)sum / INTBASE;
phenval = 1.0 / (1.0 + exp(-genval));
#if RSDEBUG
//DEBUGLOG << "Chromosome::probval(): this=" << this
//	<< " sum=" << sum
//	<< " genval=" << genval
//	<< " phenval=" << phenval
//	<< endl;
#endif
return phenval;
}
*/

// Set up chromosome at simulation initialisation
//void Chromosome::initialise(const float mean,const float sd,
//	const bool diploid) {
void Chromosome::initialise(const double mean, const double sd, 
		const bool diploid) {
double avalue;
double intbase = INTBASE;
//// adjust mean and s.d. allowing for number of alleles determining phenotype
//double adjmean,adjsd,factor;
//if (diploid) factor = (double)(nloci * 2); else factor = (double)(nloci);
//adjmean = mean / factor;
//adjsd = sqrt(sd * sd / factor);
for (int i = 0; i < nloci; i++) {
//	avalue = pRandom->Normal(adjmean,adjsd);
	avalue = pRandom->Normal(mean,sd);
	if (avalue > 0.0)
		pLoci[i].allele[0] = (int)(avalue * intbase + 0.5);
	else
		pLoci[i].allele[0] = (int)(avalue * intbase - 0.5);
#if RSDEBUG
//DEBUGLOG << "Chromosome::initialise(): this=" << this
//	<< " mean=" << mean << " sd=" << sd
//	<< " i=" << i << " avalue=" << avalue
//	<< " allele[0]=" << pLoci[i].allele[0];
//DEBUGLOG << endl;
#endif
	if (diploid) {
//		avalue = pRandom->Normal(adjmean,adjsd);
		avalue = pRandom->Normal(mean,sd);
		if (avalue > 0.0)
			pLoci[i].allele[1] = (int)(avalue * intbase + 0.5);
		else
			pLoci[i].allele[1] = (int)(avalue * intbase - 0.5);
#if RSDEBUG
//DEBUGLOG << "Chromosome::initialise(): this=" << this
//	<< " mean=" << mean << " sd=" << sd
//	<< " i=" << i << " avalue=" << avalue
//	<< " allele[1]=" << pLoci[i].allele[1];
//DEBUGLOG << endl;
#endif
	}
}

}

// Set up specified locus at simulation initialisation
void Chromosome::initialise(const short locus,const short posn,const int aval)
{
// note that initialising value is ADDED to current value to allow for pleiotropy
pLoci[locus].allele[posn] += aval;
}

// Inherit from specified parent
void Chromosome::inherit(const Chromosome *parentChr,const short posn,const short nloc,
	const double probmutn,const double probcross,const double mutnSD,const bool diploid)
{

// NOTE: At present for diploid genome, presence of crossover is determined at each
// locus (except first). However, Roslyn has shown that it is more efficient to sample
// crossover locations from geometric distribution if number of loci is large.
// HOW LARGE IS 'LARGE' IN THIS CASE?...

//// adjust mutation variance for number of loci
//double mutnsd;
//if (diploid)
//	mutnsd = sqrt(mutnSD * mutnSD / (double)(2*nloc));
//else
//	mutnsd = sqrt(mutnSD * mutnSD / (double)nloc);
int ix = 0; // indexes maternal and paternal strands
if (diploid) ix = pRandom->Bernoulli(0.5); // start index at random
for (int i = 0; i < nloc; i++) {
	if (diploid) {
		pLoci[i].allele[posn] = parentChr->pLoci[i].allele[ix];
		if (pRandom->Bernoulli(probcross)) { // crossover occurs
			if (ix == 0) ix = 1; else ix = 0;
		}
	}
	else
		pLoci[i].allele[posn] = parentChr->pLoci[i].allele[0];
#if RSDEBUG
//DEBUGLOG << "Chromosome::inherit(): this=" << this
//	<< " posn=" << posn << " nloc=" << nloc << " pmutn=" << pmutn
//	<< " i=" << i << " allele=" << pLoci[i].allele[posn]
//	<< endl;
#endif
	if (pRandom->Bernoulli(probmutn)) { // mutation occurs
		double intbase = INTBASE;
#if RSDEBUG
		int oldval = pLoci[i].allele[posn];
#endif
//		double mutnvalue = pRandom->Normal(0,mutnsd);
		double mutnvalue = pRandom->Normal(0,mutnSD);
		if (mutnvalue > 0.0)
			pLoci[i].allele[posn] += (int)(intbase * mutnvalue + 0.5);
		else
			pLoci[i].allele[posn] += (int)(intbase * mutnvalue - 0.5);
#if RSDEBUG
//DEBUGLOG << "Chromosome::inherit(): this=" << this
//	<< " probmutn=" << probmutn << " nloc=" << nloc
////	<< " mutnsd=" << mutnsd
//	<< " mutnSD=" << mutnSD
//	<< " posn=" << posn << " locus=" << i << " MUTATED "
//	<< " old=" << oldval << " new=" << pLoci[i].allele[posn]
//	<< endl;
#endif
#if RSDEBUG
//MUTNLOG << 1 << endl;
MUTNLOG << mutnvalue << " " << oldval << " " << pLoci[i].allele[posn] << " " << endl;
#endif
	}
#if RSDEBUG
//	else {
//MUTNLOG << 0 << endl;
//	}
#endif
}
}


//---------------------------------------------------------------------------

// NB THIS FUNCTION IS CURRENTLY NOT BEING CALLED TO CONSTRUCT AN INSTANCE OF Genome
// Genome(int) IS USED INSTEAD

Genome::Genome(){
pChromosome = NULL;
nChromosomes = 0;
}

// Set up new genome at initialisation for 1 chromosome per trait
Genome::Genome(int nchromosomes,int nloci,bool d) {
#if RSDEBUG
//DEBUGLOG << "Genome::Genome(): this=" << this
//	<< " nchromosomes=" << nchromosomes << " nLoci=" << nLoci
//	<< endl;
#endif

diploid = d;
if (nchromosomes > 0) nChromosomes = nchromosomes; else nChromosomes = 1;
pChromosome = new Chromosome *[nChromosomes];
for (int i = 0; i < nChromosomes; i++) {
	pChromosome[i] = new Chromosome(nloci);
//	pChromosome[i]->initialise(alleleMean,alleleSD);
}

}

// Set up new genome at initialisation for trait mapping
Genome::Genome(Species *pSpecies) {
int nloci;
nChromosomes = pSpecies->getNChromosomes();
diploid = pSpecies->isDiploid();
#if RSDEBUG
//DEBUGLOG << "Genome::Genome(): this=" << this
//	<< " nChromosomes=" << nChromosomes
//	<< endl;
#endif
pChromosome = new Chromosome *[nChromosomes];
for (int i = 0; i < nChromosomes; i++) {
	nloci = pSpecies->getNLoci(i);
#if RSDEBUG
//DEBUGLOG << "Genome::Genome(): this=" << this
//	<< " i=" << i << " nloci=" << nloci << endl;
#endif
	pChromosome[i] = new Chromosome(nloci);
}
}

// Inherit genome from parent(s)
Genome::Genome(Species *pSpecies,Genome *mother,Genome *father)
{
genomeData gen = pSpecies->getGenomeData();

nChromosomes = mother->nChromosomes;
diploid = mother->diploid;
pChromosome = new Chromosome *[nChromosomes];

for (int i = 0; i < nChromosomes; i++) {
	pChromosome[i] = new Chromosome(mother->pChromosome[i]->nLoci());
	inherit(mother,0,i,gen.probMutn,gen.probCrossover,gen.mutationSD);
	if (diploid) {
		if (father == 0) { // species is hermaphrodite - inherit again from mother
			inherit(mother,1,i,gen.probMutn,gen.probCrossover,gen.mutationSD);
		}
		else inherit(father,1,i,gen.probMutn,gen.probCrossover,gen.mutationSD);
	}
}

}

Genome::~Genome(){

if (pChromosome == NULL) return;

for (int i = 0; i < nChromosomes; i++) {
	delete pChromosome[i];
}
delete[] pChromosome;

}

/*
Genome::~Genome(void)
{
if (genome != 0) {
	delete genome; genome = NULL;
}
}
*/

//---------------------------------------------------------------------------

void Genome::setDiploid(bool dip) { diploid = dip; }
bool Genome::isDiploid(void) { return diploid; }
short Genome::getNChromosomes(void) { return nChromosomes; }
#if VIRTUALECOLOGIST
int Genome::getChromosomeNloci(short i) {
if (i >= 0 && i < nChromosomes) {
	return pChromosome[i]->nLoci();
}
else {
	return 0;
}
}
#endif

//---------------------------------------------------------------------------

// Inherit from specified parent
void Genome::inherit(const Genome *parent,const short posn,const short chr,
	const double probmutn,const double probcross,const double mutnSD)
{
// adjust mutation variance for number of loci
//double mutnSD = sqrt(mutationSD * mutationSD / (double)nLoci);
//for (int i = 0; i < nChromosomes; i++) {
//	pChromosome[i]->inherit(parent->pChromosome[i],pChromosome[i]->nLoci());
//}
pChromosome[chr]->inherit(parent->pChromosome[chr],posn,parent->pChromosome[chr]->nLoci(),
	probmutn,probcross,mutnSD,diploid);

}

#if GROUPDISP || ROBFITT
void Genome::outGenHeaders(const int rep,const int landNr,const bool patchmodel,
	const bool xtab)
#else
void Genome::outGenHeaders(const int rep,const int landNr,const bool xtab)
#endif
{

if (landNr == -999) { // close file
	if (outGenetic.is_open()) {
		outGenetic.close(); outGenetic.clear();
	}
	return;
}

string name;
simParams sim = paramsSim->getSim();

if (sim.batchMode) {
	name = paramsSim->getDir(2)
		+ "Batch" + Int2Str(sim.batchNum) + "_"
		+ "Sim" + Int2Str(sim.simulation)
		+ "_Land" + Int2Str(landNr) + "_Rep" + Int2Str(rep) + "_Genetics.txt";
}
else {
	name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation)
		+ "_Rep" + Int2Str(rep) +"_Genetics.txt";
}
outGenetic.open(name.c_str());

outGenetic << "Rep\tYear\tSpecies\tIndID";
if (xtab) {
#if GROUPDISP || ROBFITT
	if (patchmodel) outGenetic << "\tPatchID";
	else outGenetic << "\tX\tY";
#endif
	for (int i = 0; i < nChromosomes; i++) {
		int nloci = pChromosome[i]->nLoci();
		for (int j = 0; j < nloci; j++) {
			outGenetic << "\tChr" << i << "Loc" << j << "Allele0";
			if (diploid) outGenetic << "\tChr" << i << "Loc" << j << "Allele1";
		}
	}
	outGenetic << endl;
}
else {
	outGenetic << "\tChromosome\tLocus\tAllele0";
	if (diploid) outGenetic << "\tAllele1";
	outGenetic << endl;
}

}

#if GROUPDISP || ROBFITT
void Genome::outGenetics(const int rep,const int year,const int spnum,
	const int indID,const int X,const int Y,const bool patchmodel,const bool xtab)
#else
void Genome::outGenetics(const int rep,const int year,const int spnum,
	const int indID,const bool xtab)
#endif
	{
locus l;
if (xtab) {
	outGenetic << rep << "\t" << year << "\t" << spnum << "\t" << indID;
#if GROUPDISP || ROBFITT
	if (patchmodel) outGenetic << "\t" << X;
	else outGenetic << "\t" << X << "\t" << Y;
#endif
	for (int i = 0; i < nChromosomes; i++) {
		int nloci = pChromosome[i]->nLoci();
		for (int j = 0; j < nloci; j++) {
			l = pChromosome[i]->alleles(j);
			outGenetic << "\t" << l.allele[0];
			if (diploid) outGenetic << "\t" << l.allele[1];
		}
	}
	outGenetic << endl;
}
else {
	for (int i = 0; i < nChromosomes; i++) {
		int nloci = pChromosome[i]->nLoci();
		for (int j = 0; j < nloci; j++) {
			outGenetic << rep << "\t" << year << "\t" << spnum << "\t"
				<< indID << "\t" << i << "\t" << j;
			l = pChromosome[i]->alleles(j);
			outGenetic << "\t" << l.allele[0];
			if (diploid) outGenetic << "\t" << l.allele[1];
			outGenetic << endl;
		}
	}
}
}

//---------------------------------------------------------------------------

// Set up new gene at initialisation for 1 chromosome per trait
//void Genome::setGene(const short chr,const short exp,
//	const float traitval,const float alleleSD)
void Genome::setGene(const short chr, const short exp,
		const double traitval, const double alleleSD) 
// NB PARAMETER exp FOR EXPRESSION TYPE IS NOT CURRENTLY USED...
{
#if RSDEBUG
//DEBUGLOG << "Genome::setGene(): this=" << this
//	<< " chr=" << chr
//	<< " exp=" << exp << " traitval=" << traitval
//	<< " alleleSD=" << alleleSD
//	<< endl;
#endif
if (chr >= 0 && chr < nChromosomes) {
	pChromosome[chr]->initialise(traitval,alleleSD,diploid);
}
}

// Set up trait at initialisation for trait mapping
//void Genome::setTrait(Species *pSpecies,const int trait,
//	const float traitval,const float alleleSD)
void Genome::setTrait(Species* pSpecies, const int trait,
		const double traitval, const double alleleSD)	
{
traitAllele allele;
int nalleles = pSpecies->getNTraitAlleles(trait);
int ntraitmaps = pSpecies->getNTraitMaps();
#if RSDEBUG
//DEBUGLOG << "Genome::setTrait(): this=" << this << " ntraitmaps=" << ntraitmaps
//	<< " trait=" << trait << " traitval=" << traitval << " alleleSD=" << alleleSD
//	<< " nalleles=" << nalleles
//	<< endl;
#endif

int avalue;
double intbase = INTBASE;
if (trait < ntraitmaps) {
//	// adjust mean and s.d. allowing for number of alleles determining phenotype
//	double adjmean,adjsd,factor;
//	if (diploid) factor = (double)(nalleles * 2); else factor = (double)(nalleles);
//	adjmean = traitval / factor;
//	adjsd = sqrt(alleleSD * alleleSD / factor);

	for (int i = 0; i < nalleles; i++) {
		allele = pSpecies->getTraitAllele(trait,i);
//		avalue = (int)(pRandom->Normal(adjmean,adjsd) * intbase);
		avalue = (int)(pRandom->Normal(traitval,alleleSD) * intbase);
#if RSDEBUG
//DEBUGLOG << "Genome::setTrait(): this=" << this
//	<< " i=" << i << " chromo=" << allele.chromo << " locus=" << allele.locus
//	<< " posn=" << 0 << " avalue=" << avalue
//	<< endl;
#endif
		pChromosome[allele.chromo]->initialise(allele.locus,0,avalue);
		if (diploid) {
//			avalue = (int)(pRandom->Normal(adjmean,adjsd) * intbase);
			avalue = (int)(pRandom->Normal(traitval,alleleSD) * intbase);
			pChromosome[allele.chromo]->initialise(allele.locus,1,avalue);
		}
	}
}
else { // insufficient traits were defined
	// alleles cannot be initialised - all individuals have mean phenotype
#if RSDEBUG
//DEBUGLOG << "Genome::setTrait(): this=" << this
//	<< " *** unable to initialise undefined trait ***" << endl;
#endif
}

}

// Set up trait at initialisation for trait mapping
void Genome::setNeutralLoci(Species *pSpecies,const double alleleSD)
{
traitAllele allele;
int nneutral = pSpecies->getNNeutralLoci();
#if RSDEBUG
//DEBUGLOG << "Genome::setNeutralLoci(): this=" << this
//	<< " nneutral=" << nneutral << " alleleSD=" << alleleSD
//	<< endl;
#endif

//int avalue;
double avalue;
double intbase = INTBASE;
//// adjust allele s.d. for diploid genotype
//double adjsd,factor;
//if (diploid) factor = 2.0; else factor = 1.0;
//adjsd = sqrt(alleleSD * alleleSD / factor);

for (int i = 0; i < nneutral; i++) {
	allele = pSpecies->getNeutralAllele(i);
//	avalue = (int)(pRandom->Normal(0.0,adjsd) * intbase);
//	pChromosome[allele.chromo]->initialise(allele.locus,0,avalue);
	avalue = pRandom->Normal(0.0,alleleSD);
	if (avalue > 0.0)
		pChromosome[allele.chromo]->initialise(allele.locus,0,(int)(avalue * intbase + 0.5));
	else
		pChromosome[allele.chromo]->initialise(allele.locus,0,(int)(avalue * intbase - 0.5));
#if RSDEBUG
//DEBUGLOG << "Genome::setNeutralLoci(): this=" << this
//	<< " i=" << i << " chromo=" << allele.chromo << " locus=" << allele.locus
//	<< " avalue=" << avalue
//	<< endl;
#endif
	if (diploid) {
//		avalue = (int)(pRandom->Normal(0.0,adjsd) * intbase);
//		pChromosome[allele.chromo]->initialise(allele.locus,1,avalue);
		avalue = pRandom->Normal(0.0,alleleSD);
		if (avalue > 0.0)
			pChromosome[allele.chromo]->initialise(allele.locus,1,(int)(avalue * intbase + 0.5));
		else
			pChromosome[allele.chromo]->initialise(allele.locus,1,(int)(avalue * intbase - 0.5));
	}
}
}

// Return a single allele, having applied mutation
// Either allele is returned with equal probability, unless the gene is sex-linked,
// in which case a male returns the allele inherited from his father and a female
// returns the allele inherited from her mother
/*
double Genome::copy(int i,int sexx)
{
double genevalue = 0.0;
return genevalue;
}
*/

// Return the expressed value of a gene when species has one chromosome per trait
double Genome::express(short chr,short expr,short indsex)
{
double genevalue = 0.0;
//if (expr == 1) {
//	genevalue = pChromosome[chr]->probval(diploid);
//}
//else
//genevalue = pChromosome[chr]->additive(diploid);
genevalue = pChromosome[chr]->meanvalue(diploid);
#if RSDEBUG
//DEBUGLOG << "Genome::express(): this=" << this
//	<< " chr=" << chr
//	<< " genevalue=" << genevalue
//	<< endl;
#endif
return genevalue;
}
/*
// Return the expressed value of an allele
double Genome::express(short chr,short loc)
{
double genevalue = 0.0;
//if (expr == 1) {
//	genevalue = pChromosome[chr]->probval(diploid);
//}
//else
genevalue = pChromosome[chr]->additive(loc,diploid);
#if RSDEBUG
//DEBUGLOG << "Genome::express(): this=" << this
//	<< " expr=" << expr
//	<< " genevalue=" << genevalue
//	<< endl;
#endif
return genevalue;
}
*/

// Return the expressed value of a trait when genetic architecture is defined
double Genome::express(Species *pSpecies,short traitnum)
{
double genevalue = 0.0;

traitAllele allele;
int nalleles = pSpecies->getNTraitAlleles(traitnum);
if (nalleles > 0) {
	for (int i = 0; i < nalleles; i++) {
		allele = pSpecies->getTraitAllele(traitnum,i);
		genevalue += pChromosome[allele.chromo]->additive(allele.locus,diploid);
	}
	genevalue /= (double)nalleles;
	if (diploid) genevalue /= 2.0;
}

//	genevalue = pChromosome[chr]->additive(diploid);
#if RSDEBUG
//DEBUGLOG << "Genome::express(): this=" << this
//	<< " traitnum=" << traitnum << " nalleles=" << nalleles
//	<< " genevalue=" << genevalue
//	<< endl;
#endif
return genevalue;
}


locusOK Genome::getAlleles(short chr,short loc) {
locusOK l;
l.allele[0] = l.allele[1] = 0; l.ok = false;
//double intbase = INTBASE;
if (chr >= 0 && chr < nChromosomes) {
	if (pChromosome[chr] != 0) {
		if (loc >= 0 && loc < pChromosome[chr]->nLoci()) {
			locus a = pChromosome[chr]->alleles(loc);
			l.allele[0] = a.allele[0]; l.allele[1] = a.allele[1]; l.ok = true;
		}
	}
}
//	l.allele[0] = (int)(intbase * pChromosome[0]->additive(diploid));

return l;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

// OLD CODE WHICH MIGHT BE USEFUL FOR IMPLEMENTING ALTERNATIVE FORMS OF EXPRESSION

/*

// Set up new genome at initialisation
Genome::Genome(int g) {
if (g < 1) return;

genomesize = g;
genome = new gene *[genomesize];
for (int i = 0; i < genomesize; i++) {
	genome[i] = new gene;
	genome[i]->expression = 0; genome[i]->mutntype = 0;
	genome[i]->Pmutn = 0.0; genome[i]->mutnSize = 0.0;
	genome[i]->allele[0] = genome[i]->allele[1] = 0.0;
}

}

// Inherit genome from parent(s)
Genome::Genome(Genome *mother,Genome *father)
{
#if RSDEBUG
//locn currloc = currCell->getLocn();
//DEBUGLOG << "Genome::setGenes(): indId=" << indId
//	<< " x=" << currloc.x << " y=" << currloc.y
//	<< " pSpecies=" << pSpecies
//	<< " mother=" << mother
//	<< " father=" << father
//	<< endl;
#endif

genomesize = mother->genomesize;
genome = new gene *[genomesize];
for (int i = 0; i < genomesize; i++) {
	genome[i] = new gene;
	genome[i]->expression = mother->genome[i]->expression;
	genome[i]->mutntype   = mother->genome[i]->mutntype;
	genome[i]->Pmutn      = mother->genome[i]->Pmutn;
	genome[i]->mutnSize   = mother->genome[i]->mutnSize;
	genome[i]->allele[0] = mother->copy(i,0);
	if (father == 0) genome[i]->allele[1] = 0.0;
	else genome[i]->allele[1] = father->copy(i,1);
}

}

Genome::~Genome(void)
{
if (genome != 0) {
	for (int i = 0; i < genomesize; i++) {
		delete genome[i]; genome[i] = NULL;
	}
	delete genome; genome = NULL;
}

}

// Set up new gene
void Genome::setGene(int i,short exp,short mtype,float pmut,float msize,double a0,double a1)
{
if (exp >= 0 && exp <= 4) genome[i]->expression = exp;
if (mtype >= 0 && mtype <= 5) genome[i]->mutntype = mtype;
if (pmut >= 0.0 && pmut <= 1.0) genome[i]->Pmutn = pmut;
if (msize > 0.0) genome[i]->mutnSize = msize;
genome[i]->allele[0] = a0;
genome[i]->allele[1] = a1;
}

// Return a single allele, having applied mutation
// Either allele is returned with equal probability, unless the gene is sex-linked,
// in which case a male returns the allele inherited from his father and a female
// returns the allele inherited from her mother
double Genome::copy(int i,int sexx)
{
double genevalue,logP,logitP;

switch (genome[i]->expression) {
case 0: // haploid
	genevalue = genome[i]->allele[0];
	break;
case 1: // sex-linked gene
	if (sexx == 0 || sexx == 1) genevalue = genome[i]->allele[sexx];
	else return -999999.999; // indicates error in sex of parent animal
	break;
default: // otherwise - select one of the alleles at random
	genevalue = genome[i]->allele[pRandom->Bernoulli(0.5)];
	;
}
#if RSDEBUG
//DEBUGLOG << "Gene::copy(): expression=" << expression
//	<< " mutntype=" << mutntype
//	<< " Pmutn=" << Pmutn << " mutnSize=" << mutnSize
//	<< " allele[0]=" << allele[0] << " allele[1]=" << allele[1]
//	<< " genevalue=" << genevalue
//	<< endl;
#endif
// apply mutation to the gene
if (pRandom->Random() < genome[i]->Pmutn) {
	switch (genome[i]->mutntype) {
		case 0: genevalue = pRandom->Random(); break;
		case 1: genevalue += pRandom->Normal(0.0,genome[i]->mutnSize); break;
		case 2:
			logitP = log(genevalue/(1.0-genevalue));
			logitP += pRandom->Normal(0.0,genome[i]->mutnSize);
			genevalue = 1.0 / (1 + exp(-logitP));
			if (genevalue >= 1.0) genevalue = 0.999999999;
			if (genevalue <= 0.0) genevalue = 0.000000001;
			break;
		case 3:
			logP = log(genevalue);
			logP += pRandom->Normal(0.0,genome[i]->mutnSize);
			genevalue = exp(logP);
			break;
		case 4:
			genevalue += pRandom->Random()*(2.0*genome[i]->mutnSize) - genome[i]->mutnSize;
			break;
		case 5:
			genevalue += pRandom->Random()*(2.0*genome[i]->mutnSize) - genome[i]->mutnSize;
			if (genevalue > 1.0) genevalue = 1.0;
			if (genevalue < 0.0) genevalue = 0.0;
			break;
		default: genevalue = -888888.888; // error in mutation type
	}
}
#if RSDEBUG
//DEBUGLOG << "Gene::copy():"
//	<< " allele[0]=" << allele[0] << " allele[1]=" << allele[1]
//	<< " genevalue=" << genevalue
//	<< endl;
#endif
return genevalue;
}

// Return the expressed value of a gene
double Genome::express(int i,int sexx)
{
double genevalue;

switch (genome[i]->expression) {
	case 0: // sex-linked gene
		genevalue = genome[i]->allele[0];
		break;
	case 1: // sex-linked gene
		if (sexx == 0 || sexx == 1) genevalue = genome[i]->allele[sexx];
		else genevalue = -999999.999; // indicates error in sex of parent animal
		break;
	case 2: // dominance of greater value
		genevalue = genome[i]->allele[0];
		if (genome[i]->allele[1] > genevalue) genevalue = genome[i]->allele[1];
		break;
	case 3: // co-dominance (mean value)
		genevalue = (genome[i]->allele[0] + genome[i]->allele[1]) / 2.0;
		break;
	case 4: // sex-expressed co-dominance
// NB does not differ from co-dominance other than when genes interact to determine
// phenotypic level of expression
		genevalue = (genome[i]->allele[0] + genome[i]->allele[1]) / 2.0;
//    if (sexx == 0 || sexx == 1) genevalue = genome[i]->allele[sexx];
//    else genevalue = -999999.999; // indicates error in sex of parent animal
    break;
default: genevalue = -777777.777; // error in expression type
}

return genevalue;

}

locus Genome::getAlleles(int g) {
locus l;
if (g >= 0 && g < genomesize) {
	l.allele[0] = genome[g]->allele[0];
	l.allele[1] = genome[g]->allele[1];
}
else {
	l.allele[0] = l.allele[1] = 0.0;
}
return l;
}

*/

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

