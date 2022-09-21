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

RangeShifter v2.0 Genome

Implements the Genome class

Author: Steve Palmer & Roslyn Henry, University of Aberdeen

Last updated: 28 July 2021 by Greta Bocedi

------------------------------------------------------------------------------*/

#ifndef GenomeH
#define GenomeH

#include <vector>
#include <algorithm>
//#include "maths.h"

#include "Parameters.h"
#include "Species.h"

#define INTBASE 100.0; // to convert integer alleles into continuous traits

struct locus { short allele[2]; };
struct locusOK { short allele[2]; bool ok; };

//---------------------------------------------------------------------------

class Chromosome {

public:
	Chromosome(int);
	~Chromosome();
	short nLoci(void);
	double additive( // Return trait value on normalised genetic scale
		const bool	// diploid
	);
	double meanvalue( // Return trait value on normalised genetic scale
		const bool	// diploid
	);
	double additive( // Return trait value on normalised genetic scale
		const short,	// locus
		const bool		// diploid
	);
//	double probval(const bool);
	locus alleles(	// Return allele values at a specified locus
		const int			// position of locus on chromosome
	);
	void initialise( // Set up chromosome at simulation initialisation
		const double,	// normalised phenotypic trait value			
		const double,	// s.d. of allelic variance (genetic scale)
		const bool		// diploid
	);
	void initialise( // Set up specified locus at simulation initialisation
		const short,	// locus
		const short,	// position: 0 from mother, 1 from father
		const int			// allele value
	);
	void inherit( // Inherit chromosome from specified parent
		const Chromosome*,	// pointer to parent's chromosome
		const short,				// position: 0 from mother, 1 from father
		const short,				// no. of loci
		const double,				// mutation probability
		const double,				// crossover probability
		const double,				// s.d. of mutation magnitude (genetic scale)
		const bool					// diploid
	);

protected:

private:
	short nloci;
	locus *pLoci;

};

//---------------------------------------------------------------------------

class Genome{

public:
//
//	static float delPMutation, delEffectSize, delDominance, delBackMutation ,genomeMeanRecombination;
//	static int delMaxSize, delNMutations;
//	static bool genomeCanRecombine, genomeCompletelyUnlinked;
//

	Genome();
	Genome(int,int,bool);
	Genome(Species*);
	Genome(Species*,Genome*,Genome*);
	~Genome();
	void setGene( // Set up new gene at initialisation for 1 chromosome per trait
		const short,	// chromosome number
		const short,	// expression type (NOT CURRENTLY USED)
		const double,	// normalised trait value		
		const double		// s.d. of allelic variance
	);
	void setTrait( // Set up trait at initialisation for trait mapping
		Species*,			// pointer to Species
		const int,		// trait number			
		const double,	// normalised trait value		
		const double		// s.d. of allelic variance
	);
	void setNeutralLoci( // Set up neutral loci at initialisation
		Species*,			// pointer to Species
		const double		// s.d. of allelic variance
	);
//	double copy(int,int);
	double express(
	// Return the expressed value of a gene when species has one chromosome per trait
		short,	// chromosome number
		short,	// expression type (NOT CURRENTLY USED)
		short		// individual's sex (NOT CURRENTLY USED)
	);
//	double express(
//		short,	// chromosome number
//		short		// locus on chromosome
//	);
	double express(
	// Return the expressed value of a trait when genetic architecture is defined
		Species*,	// pointer to Species
		short			// trait number
//		bool			// true if trait is sex-dependent
	);
	locusOK getAlleles( // Get allele values at a specified locus
		short,	// chromosome number
		short		// locus position on chromosome
	);
	// SCFP NEW DECLARATIONS
	void setDiploid(bool);
	bool isDiploid(void);
	void inherit( // Inherit from specified parent
		const Genome*,	// pointer to parent's genome
		const short,		// position: 0 from mother, 1 from father
		const short,		// chromasome number
		const double,		// mutation probability
		const double,		// crossover probability
		const double			// s.d. of mutation magnitude (genetic scale)
	);
//	void setStaticData(genomeData);
//	genomeData getStaticData(void);
	short getNChromosomes(void);
#if VIRTUALECOLOGIST
	int getChromosomeNloci(short);
#endif
#if GROUPDISP || ROBFITT
	void outGenHeaders(
		const int, 	 // replicate
		const int,	 // landscape number
		const bool,  // patch-based landscape
		const bool	 // output as cross table?
	);
	void outGenetics(
		const int,		// replicate
		const int,		// year
		const int,		// species number
		const int,		// individual ID
		const int,		// X co-ordinate OR patch ID
		const int,		// Y co-ordinate
		const bool,  	// patch-based landscape
		const bool		// output as cross table?
	);
#else
	void outGenHeaders(
		const int,	// replicate
		const int,	// landscape number
		const bool	// output as cross table?
	);
	void outGenetics(
		const int,	// replicate
		const int,	// year
		const int,	// species number
		const int, 	// individual ID
		const bool 	// output as cross table?
	);
#endif


private:
	short nChromosomes;						// no. of chromosomes
	bool diploid;
	Chromosome **pChromosome;

};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

/*

struct gene {
	short expression;	// used to control how gene is expressed:
		// 0 = haploid
		// 1 = sex-linked
		// 2 = dominance of greater value
		// 3 = co-dominance (mean value)
		// 4 = sex-expressed (but not sex-linked)
	short mutntype;		// mutation type:
		// 0 = random in [0,1]
		// 1 = normal
		// 2 = logit scale
		// 3 = log-normal (gene must be positive)
		// 4 = RS method - random in specified range
		// 5 = RS method - random in specified range and constrained to [0,1]
	float Pmutn;			// mutation probability
	float mutnSize;		// mutation size (s.d. of normal distribution)
	double allele[2]; // allele pair: [0] from mother, [1] from father
};
;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class Genome
{
public:
	Genome(int);
	Genome(Genome*,Genome*);
	~Genome(void);
	void setGene(int);
	void setGene(int,short,short,float,float,double,double);  // new gene
	double copy(int,int);
	double express(int,int);
	locus getAlleles( // Get allele values at a specified locus
		int	// locus position within genome
	);

private:
	int genomesize;
	gene **genome;

};

*/

//---------------------------------------------------------------------------

extern paramSim *paramsSim;
extern RSrandom *pRandom;

#if RSDEBUG
extern ofstream DEBUGLOG;
extern ofstream MUTNLOG;
extern void DebugGUI(string);
#endif

//---------------------------------------------------------------------------

#endif
