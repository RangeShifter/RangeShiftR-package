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

#pragma hdrstop

#include "RandomCheck.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

ifstream inRandom;
ofstream outRandom;
ofstream outBernoulli;
ofstream outNormal;
ofstream outPoisson;
ofstream outIRandom;

void randomCheck(void)
{

	int samplesize,irandMin,irandMax;
	double bernMean,normMean,normSD,poisMean;
	string name,header;
	simParams sim = paramsSim->getSim();

	name = paramsSim->getDir(1) + "RandomCheck.txt";
	inRandom.open(name.c_str());
	if (!inRandom.is_open()) {
		inRandom.clear();
		return;
	}
	for (int i = 0; i < 7; i++) {
		inRandom >> header;
	}
	inRandom >> samplesize >> bernMean >> normMean >> normSD >> poisMean >> irandMin >> irandMax;

	name = paramsSim->getDir(2) + "Random.txt";
	outRandom.open(name.c_str());
	name = paramsSim->getDir(2) + "Bernoulli.txt";
	outBernoulli.open(name.c_str());
	name = paramsSim->getDir(2) + "Normal.txt";
	outNormal.open(name.c_str());
	name = paramsSim->getDir(2) + "Poisson.txt";
	outPoisson.open(name.c_str());
	name = paramsSim->getDir(2) + "IRandom.txt";
	outIRandom.open(name.c_str());

	for (int i = 0; i < samplesize; i++) {
		outRandom << pRandom->Random() << endl;
		outBernoulli << pRandom->Bernoulli(bernMean) << endl;
		outNormal << pRandom->Normal(normMean,normSD) << endl;
		outPoisson << pRandom->Poisson(poisMean) << endl;
		outIRandom << pRandom->IRandom(irandMin,irandMax) << endl;
	}

	inRandom.close();
	inRandom.clear();
	outRandom.close();
	outRandom.clear();
	outBernoulli.close();
	outBernoulli.clear();
	outNormal.close();
	outNormal.clear();
	outPoisson.close();
	outPoisson.clear();
	outIRandom.close();
	outIRandom.clear();

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
