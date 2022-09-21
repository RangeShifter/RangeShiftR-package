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

RangeShifter v2.0 RSrandom

Implements the RSrandom class

Authors: Steve Palmer, University of Aberdeen
				 Anne-Kathleen Malchow, Potsdam University

Last updated: 12 January 2021 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef RSrandomH
#define RSrandomH

#include <stdlib.h>
#include <fstream>
//#include <iostream>

#include "Version.h"

using namespace std;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif



#if RS_EMBARCADERO

//--------------- 1.) Former version of RSrandom.cpp


	#if LINUX_CLUSTER
	//#include <random>
	//#include <tr1/random>
	#include "maths.h"
	#else
	#if RSWIN64
	#include <dinkumware64/random>
	#else
	#include <dinkumware/random>
	#endif
	#endif

	class RSrandom {

	public:
	#if RS_ABC
		RSrandom(int);
	#else
		RSrandom(void);
	#endif
		~RSrandom(void);
		double Random(void);
		int IRandom(int,int);
		int Bernoulli(double);
		double Normal(double,double);
		int Poisson(double);
	#if RS_ABC
		double Beta(double,double);
		double Gamma(double,double);
	#endif

	private:
		double normal_x2; int normal_x2_valid; // variables used by Normal distribution
	#if !LINUX_CLUSTER
		tr1::uniform_real<> *pRandom01;
		tr1::normal_distribution<> *pNormal;
	#if RS_ABC
		tr1::mt19937 *gen;
	#endif
	#endif
	};

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

#else // not RS_EMBARCADERO



#if !RS_RCPP

//--------------- 2.) New version of RSrandom.cpp


	#include <cmath>
	#include <random>
	#if !LINUX_CLUSTER
	#include <ctime>
	#endif

	class RSrandom
	{

	public:
	#if RS_ABC
		RSrandom(int);
	#else
		RSrandom(void);
	#endif
		~RSrandom(void);
		double Random(void);
		int IRandom(int, int);
		int Bernoulli(double);
		double Normal(double, double);
		int Poisson(double);
		mt19937 getRNG(void);
  #if RS_ABC
		double Beta(double,double);
		double Gamma(double,double);
	#endif

	private:
		mt19937* gen;
		std::uniform_real_distribution<>* pRandom01;
		std::normal_distribution<>* pNormal;
	};


//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------



//--------------- 3.) R package version of RSrandom.cpp


#else // if RS_RCPP 


	#include <cmath>
	#include <random>
	#if RSWIN64
	#include <ctime>
	#endif

	class RSrandom {

	public:
		RSrandom(std::int64_t);       // if int is negative, a random seed will be generated, else it is used as seed
		~RSrandom(void);
		mt19937 getRNG(void);
		double Random(void);
		int IRandom(int,int);
		int Bernoulli(double);
		double Normal(double,double);
		int Poisson(double);
	/* ADDITIONAL DISTRIBUTIONS
		double Beta(double,double);
		double Gamma(double,double); // !! make sure correct definition is used: using shape and scale (as defined here) OR using shape/alpha and rate/beta (=1/scale)
		double Cauchy(double,double);
	*/

	private:
		mt19937 *gen;
		std::uniform_real_distribution<> *pRandom01;
		std::normal_distribution<> *pNormal;
	};



#endif // !RS_RCPP

#endif // RS_EMBARCADERO


//---------------------------------------------------------------------------

#endif // RSrandomH



