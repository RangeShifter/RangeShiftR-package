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

#if RS_EMBARCADERO

//--------------- 1.) Former version of RSrandom.cpp

#pragma hdrstop
#pragma package(smart_init)


#if RSDEBUG
#include "Parameters.h"
extern paramSim *paramsSim;
#include <fstream>
//ofstream RSRANDOMLOG;
#endif

#if !RS_ABC

#if !LINUX_CLUSTER
// set up Mersenne Twister random generator
#if RSWIN64
#if RSDEBUG
int RS_random_seed = 666;
tr1::mt19937 gen(RS_random_seed);
#else
int RS_random_seed = time(0);
tr1::mt19937 gen(RS_random_seed);
#endif // RSDEBUG 
#else
// for some unknown reason, 32-bit compile fails if RS_random_seed is passed as parameter (as above)...
#if RSDEBUG
int RS_random_seed = 666;
tr1::mt19937 gen(666);
#else
// ...there is thus the extremely low possibiity that the value in RS_random_seed and
// the value used to seed the random number stream may differ by 1
int RS_random_seed = time(0);
tr1::mt19937 gen(time(0));
#endif // RSDEBUG 
#endif // RSWIN64 
#endif // !LINUX_CLUSTER

#endif // !RS_ABC

#if RS_ABC
int RS_random_seed;
RSrandom::RSrandom(int seed)
#else
RSrandom::RSrandom(void)
#endif
{

#if RS_ABC

#if !LINUX_CLUSTER
// set up Mersenne Twister random generator
// NB see comments above
// presumably 32-bit compile would fail, but ABC version is likely to always be
// compiled as a 64-bit version 
#if RSDEBUG
RS_random_seed = 666;
gen = new tr1::mt19937 (RS_random_seed);
// RS random initialisation log
#else
RS_random_seed = time(0) + 17 * seed;
gen = new tr1::mt19937 (RS_random_seed);
#endif // RSDEBUG
#endif // !LINUX_CLUSTER

#else

#if BATCH && RSDEBUG
DEBUGLOG << "RSrandom::RSrandom(): RS_random_seed=" << RS_random_seed
	<< endl;
#endif // RSDEBUG 

#endif // RS_ABC 

#if RSDEBUG
// RS random initialisation log added by SCFP 25/8/16
//string name = paramsSim->getDir(2) + "RandomLog.txt";
//RSRANDOMLOG.open(name.c_str());
//RSRANDOMLOG << "RSrandom::RSrandom(): creating new random stream" << endl;
#endif
/*
#  if defined(BOOST_HAS_CPP_0X)
cout << "BOOST_HAS_CPP_0X" << endl;
#  else
cout << "***NOT*** BOOST_HAS_CPP_0X" << endl;
#  endif
#  if defined(BOOST_HAS_TR1_RANDOM)
cout << "BOOST_HAS_TR1_RANDOM" << endl;
#     ifdef BOOST_HAS_INCLUDE_NEXT
cout << "BOOST_HAS_INCLUDE_NEXT" << endl;
//#        include_next <random>
#     else
cout << "***NOT*** BOOST_HAS_INCLUDE_NEXT" << endl;
cout << "BOOST_TR1_STD_HEADER" << endl;
//#        include BOOST_TR1_STD_HEADER(random)
#     endif
#  else
cout << "***NOT*** BOOST_HAS_TR1_RANDOM" << endl;
#  endif
*/

normal_x2_valid = 0;
#if !LINUX_CLUSTER
// Set up standard uniform distribution
pRandom01 = new	tr1::uniform_real<> (0.0,1.0);
// Set up standard normal distribution
pNormal = new	tr1::normal_distribution<> (0.0,1.0);
#endif

//cout << endl;
//for (int i = 0; i < 5; i++) {
//	cout << gen() << " ";
//}
//cout << endl << endl;

#if RSDEBUG
//RSRANDOMLOG.close(); RSRANDOMLOG.clear();
#endif
}

RSrandom::~RSrandom(void) {
#if RS_ABC
delete gen;
#endif
#if !LINUX_CLUSTER
if (pRandom01 != 0) delete pRandom01;
if (pNormal != 0) delete pNormal;
#endif
}

double RSrandom::Random(void) {
#if LINUX_CLUSTER
return unif_rand();
#else
#if RS_ABC
return pRandom01->operator()(*gen);
#else
return pRandom01->operator()(gen);
#endif
//double result_Random;
//result_Random = pRandom01->operator()(*gen);
//#if RSDEBUG
//DEBUGLOG << "RSrandom::Random(): result_Random=" << result_Random
//	<< endl;
//#endif
//return result_Random;
#endif
}

int RSrandom::IRandom(int min,int max) {
//#if LINUX_CLUSTER
//return irand(min,max);
//#else
//tr1::uniform_int<> dist(min,max);
//return dist(gen);
//#endif
// output random integer in the interval min <= x <= max (copied from mersenne.cpp)
int r;
r = int((max - min + 1) * Random()) + min; // multiply interval with random and truncate
if (r > max) r = max;
if (max < min) return 0x80000000;
return r;
}

int RSrandom::Bernoulli(double p) {
#if RSDEBUG
//DEBUGLOG << "RSrandom::Bernoulli(): p=" << p
//	<< endl;
#endif
#if LINUX_CLUSTER
return unif_rand() < p;
#else
return Random() < p;
#endif
}

double RSrandom::Normal(double mean,double sd) {
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): mean=" << mean << " sd=" << sd
//	<< endl;
#endif
#if LINUX_CLUSTER
// normal distribution derived from Agner Fog's method
double normal_x1;          // first random coordinate (normal_x2 is member of class)
double w;                  // radius
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): normal_x2_valid=" << normal_x2_valid
//	<< endl;
#endif
if (normal_x2_valid) {     // we have a valid result from last call
	normal_x2_valid = 0;
	return normal_x2 * sd + mean;
}
// make two normally distributed variates by Box-Muller transformation
do {
	normal_x1 = 2. * Random() - 1.;
	normal_x2 = 2. * Random() - 1.;
	w = normal_x1*normal_x1 + normal_x2*normal_x2;
} while (w >= 1. || w < 1E-30);
w = sqrt(log(w)*(-2./w));
normal_x1 *= w;  normal_x2 *= w;     // normal_x1 and normal_x2 are independent normally distributed variates
normal_x2_valid = 1;                 // save normal_x2 for next call
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): normal_x1=" << normal_x1
//	<< endl;
#endif
return normal_x1 * sd + mean;
#else
//double norm = pNormal->operator()(gen);
//#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): norm=" << norm
//	<< endl;
//#endif
//return mean + sd * norm;
#if RS_ABC
return mean + sd * pNormal->operator()(*gen);
#else
return mean + sd * pNormal->operator()(gen);
#endif
#endif
}

int RSrandom::Poisson(double mean) {
#if LINUX_CLUSTER
return rpois(mean);
#else
if (mean > 50.0) {
	return this->Normal(mean,sqrt(mean));
}
else {
	tr1::poisson_distribution<> poiss(mean);
	return poiss(gen);
}
#endif
}

#if RS_ABC

// Beta distribution - sample from two gamma distributions
double RSrandom::Beta(double p0,double p1) {
double g0,g1,beta;
#if RSDEBUG
//DEBUGLOG << "RSrandom::Beta(): p0=" << p0 << " p1=" << p1
//	<< endl;
#endif
if (p0 > 0.0 && p1 > 0.0) { // valid beta parameters
#if LINUX_CLUSTER
	g0 = rgamma(p0,1.0);
	g1 = rgamma(p1,1.0);
#else
	tr1::gamma_distribution<> gamma0(p0,1.0);
	tr1::gamma_distribution<> gamma1(p1,1.0);
	g0 = gamma0(*gen);
	g1 = gamma1(*gen);
#endif // LINUX_CLUSTER
	beta = g0 / (g0 + g1);
#if RSDEBUG
//DEBUGLOG << "RSrandom::Beta(): g0=" << g0 << " g1=" << g1 << " beta=" << beta
//	<< endl;
#endif
}
else { // return invalid value
	beta = -666.0;
}
return  beta;
}

// Gamma distribution
double RSrandom::Gamma(double p0,double p1) {
double p2,gamma;
if (p0 > 0.0 && p1 > 0.0) { // valid gamma parameters
	p2 = 1.0 / p1;
#if LINUX_CLUSTER
	gamma = rgamma(p0,p2);
#else
	tr1::gamma_distribution<> gamma0(p0,p2);
	gamma = gamma0(*gen);
#endif // LINUX_CLUSTER
}
else { // return invalid value
	gamma = -666.0;
}
return  gamma;
}

#endif // RS_ABC

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

#else // not RS_EMBARCADERO

//--------------- 2.) New version of RSrandom.cpp

#if !RS_RCPP 


#if RSDEBUG
#include "Parameters.h"
extern paramSim* paramsSim;
// ofstream RSRANDOMLOG;
#endif

int RS_random_seed = 0;

// C'tor
RSrandom::RSrandom()
{
#if RSDEBUG
    // fixed seed
    RS_random_seed = 666;
#else
    // random seed
#if LINUX_CLUSTER
    std::random_device device;
    RS_random_seed = device(); // old versions of g++ on Windows return a constant value within a given Windows
                                   // session; in this case better use time stamp
#else
    RS_random_seed = std::time(NULL);
#endif
#endif // RSDEBUG

#if BATCH && RSDEBUG
    DEBUGLOG << "RSrandom::RSrandom(): RS_random_seed=" << RS_random_seed << endl;
#endif // RSDEBUG

    // set up Mersenne Twister RNG
    gen = new mt19937(RS_random_seed);

    // Set up standard uniform distribution
    pRandom01 = new uniform_real_distribution<double>(0.0, 1.0);
    // Set up standard normal distribution
    pNormal = new normal_distribution<double>(0.0, 1.0);
}

RSrandom::~RSrandom(void)
{
    delete gen;
    if(pRandom01 != 0)
	delete pRandom01;
    if(pNormal != 0)
	delete pNormal;
}

mt19937 RSrandom::getRNG(void)
{
    return *gen;
}

double RSrandom::Random(void)
{
    // return random number between 0 and 1
    return pRandom01->operator()(*gen);
}

int RSrandom::IRandom(int min, int max)
{
    // return random integer in the interval min <= x <= max
    uniform_int_distribution<int> unif(min, max);
    return unif(*gen);
}

int RSrandom::Bernoulli(double p)
{
    return Random() < p;
}

double RSrandom::Normal(double mean, double sd)
{
    return mean + sd * pNormal->operator()(*gen);
}

int RSrandom::Poisson(double mean)
{
    poisson_distribution<int> poiss(mean);
    return poiss(*gen);
}


//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

#else // if RS_RCPP 

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

#endif // RS_RCPP

#endif // RS_EMBARCADERO

//---------------------------------------------------------------------------
