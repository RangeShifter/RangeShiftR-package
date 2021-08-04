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
 

#include "RSrandom.h"

//--------------- 2.) New version of RSrandom.cpp

//--------------- 3.) R package version of RSrandom.cpp

	#if RSDEBUG
	#include "Parameters.h"
	extern paramSim *paramsSim;
	//ofstream RSRANDOMLOG;
	#endif

	std::uint32_t RS_random_seed = 0;

	// C'tor
	// if parameter seed is negative, a random seed will be generated, else it is used as seed
	RSrandom::RSrandom(std::int64_t seed)
	{
		// get seed
		std::vector<std::uint32_t> random_seed(3);
		random_seed[0] = 1967593562;
		random_seed[1] = 3271254416;
		if (seed < 0) {
			// random seed
			#if RSWIN64
			random_seed[2] = std::time(NULL) + ( seed * (-17) );
			#else
			std::random_device device;
			random_seed[2] = device();
			#endif
			#if BATCH && RSDEBUG
				DEBUGLOG << "RSrandom::RSrandom(): Generated random seed = ";
			#endif
		}
		else{
			// fixed seed
			random_seed[2] = seed;
			#if BATCH && RSDEBUG
				DEBUGLOG << "RSrandom::RSrandom(): Use fixed seed = ";
			#endif
		}

		RS_random_seed = random_seed[2];
		#if BATCH && RSDEBUG
			DEBUGLOG << RS_random_seed << endl;
		#endif

		// set up Mersenne Twister random number generator with seed sequence
		std::seed_seq seq(random_seed.begin(),random_seed.end());
		gen = new mt19937(seq);

		// Set up standard uniform distribution
		pRandom01 = new uniform_real_distribution<double> (0.0,1.0);
		// Set up standard normal distribution
		pNormal = new normal_distribution<double> (0.0,1.0);
	}

	RSrandom::~RSrandom(void) {
		delete gen;
		if (pRandom01 != 0) delete pRandom01;
		if (pNormal != 0) delete pNormal;
	}

	mt19937 RSrandom::getRNG(void) {
		// return random number generator
		return *gen;
	}

	double RSrandom::Random(void) {
		// return random number between 0 and 1
		return pRandom01->operator()(*gen);
	}

	int RSrandom::IRandom(int min,int max) {
		// return random integer in the interval min <= x <= max
		uniform_int_distribution<int> unif(min,max);
		return unif(*gen);
	}

	int RSrandom::Bernoulli(double p) {
		return Random() < p;
	}

	double RSrandom::Normal(double mean,double sd) {
		return mean + sd * pNormal->operator()(*gen);
	}

	int RSrandom::Poisson(double mean) {
		poisson_distribution<int> poiss(mean);
		return poiss(*gen);
	}


	/* ADDITIONAL DISTRIBUTIONS

	// Beta distribution - sample from two gamma distributions
	double RSrandom::Beta(double p0,double p1) {
		double g0,g1,beta;
		if (p0 > 0.0 && p1 > 0.0) { // valid beta parameters
			gamma_distribution<double> gamma0(p0,1.0);
			gamma_distribution<double> gamma1(p1,1.0);
			g0 = gamma0(*gen);
			g1 = gamma1(*gen);
			beta = g0 / (g0 + g1);
		}
		else { // return invalid value
			beta = -666.0;
		}
		return beta;
	}

	// Gamma distribution
	double RSrandom::Gamma(double p0,double p1) {  // using shape (=p0) and scale (=p1)
		double p2,gamma;
		if (p0 > 0.0 && p1 > 0.0) { // valid gamma parameters
			p2 = 1.0 / p1;
			gamma_distribution<double> gamma0(p0,p2);  // using shape/alpha (=p0) and rate/beta (=p2=1/p1)
			gamma = gamma0(*gen);
		}
		else { // return invalid value
			gamma = -666.0;
		}
	return  gamma;
	}

	// Cauchy distribution
	double RSrandom::Cauchy(double loc, double scale) {
		double res;
		if (scale > 0.0) { // valid scale parameter
			cauchy_distribution<double> cauchy(loc,scale);
			res = cauchy(*gen);
		}
		else { // return invalid value
			res = -666.0;
		}
		return res;
	}


	*/

//---------------------------------------------------------------------------
