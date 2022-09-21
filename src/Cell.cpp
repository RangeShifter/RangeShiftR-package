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

#include "Cell.h"
#include "algorithm"

//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma package(smart_init)
#endif

//---------------------------------------------------------------------------

// Cell functions

Cell::Cell(int xx,int yy,intptr patch,int hab)
{
x = xx; y = yy;
pPatch = patch;
envVal = 1.0; // default - no effect of any gradient
envDev = eps = 0.0;
habIxx.push_back(hab);
#if RSDEBUG
//DebugGUI(("Cell::Cell(): this=" + Int2Str((int)this)
//	+ " x=" + Int2Str(x) + " y=" + Int2Str(y)
//	+ " habIndex=" + Int2Str(habIndex)
//).c_str());
#endif
visits = 0;
#if RS_CONTAIN
pDamage = 0;
#endif // RS_CONTAIN
smsData = 0;
}

Cell::Cell(int xx,int yy,intptr patch,float hab)
{
x = xx; y = yy;
pPatch = patch;
envVal = 1.0; // default - no effect of any gradient
envDev = eps = 0.0;
habitats.push_back(hab);
visits = 0;
#if RS_CONTAIN
pDamage = 0;
#endif // RS_CONTAIN
smsData = 0;
}

Cell::~Cell() {
#if RSDEBUG
//DEBUGLOG << "Cell::~Cell(): this = " << this << " smsData = " << smsData << endl;
#endif
habIxx.clear();
habitats.clear();
if (smsData != 0) {
	if (smsData->effcosts != 0) delete smsData->effcosts;
	delete smsData;
}
#if SPATIALDEMOG
demoScalings.clear();
#endif
#if RSDEBUG
//DEBUGLOG << "Cell::~Cell(): deleted" << endl;
#endif
}

void Cell::setHabIndex(short hx) {
#if RSDEBUG
//DebugGUI(("Cell::setHabIndex(): this=" + Int2Str((int)this)
//	+ " x=" + Int2Str(x) + " y=" + Int2Str(y)
//	+ " habIx=" + Int2Str(habIx)
//).c_str());
#endif
if (hx < 0) habIxx.push_back(0);
else habIxx.push_back(hx);
}

void Cell::changeHabIndex(short ix,short hx) {
if (ix >= 0 && ix < (short)habIxx.size() && hx >= 0) habIxx[ix] = hx;
else habIxx[ix] = 0;
}

int Cell::getHabIndex(int ix) {
if (ix < 0 || ix >= (int)habIxx.size())
	// nodata cell OR should not occur, but treat as such
	return -1;
else return habIxx[ix];
}
int Cell::nHabitats(void) {
int nh = (int)habIxx.size();
if ((int)habitats.size() > nh) nh = (int)habitats.size();
return nh;
}

void Cell::setHabitat(float q) {
if (q >= 0.0 && q <= 100.0) habitats.push_back(q);
else habitats.push_back(0.0);
}

float Cell::getHabitat(int ix) {
if (ix < 0 || ix >= (int)habitats.size())
	// nodata cell OR should not occur, but treat as such
	return -1.0;
else return habitats[ix];
}

void Cell::setPatch(intptr p) {
pPatch = p;
}
intptr Cell::getPatch(void)
{
#if RSDEBUG
//DebugGUI(("Cell::getPatch(): this=" + Int2Str((int)this)
//	+ " x=" + Int2Str(x) + " y=" + Int2Str(y)
//	+ " habIxx[0]=" + Int2Str(habIxx[0]) + " pPatch=" + Int2Str(pPatch)
//).c_str());
#endif
return pPatch;
}

locn Cell::getLocn(void) { locn q; q.x = x; q.y = y; return q; }

void Cell::setEnvDev(float d) { envDev = d; }

float Cell::getEnvDev(void) { return envDev; }

void Cell::setEnvVal(float e) {
if (e >= 0.0) envVal = e;
}

float Cell::getEnvVal(void) { return envVal; }

void Cell::updateEps(float ac,float randpart) {
eps = eps*ac + randpart;
}

float Cell::getEps(void) { return eps; }

#if SPATIALMORT
// Functions to handle additional spatial mortality

void Cell::setMort(float mort0,float mort1) {
if (mort0 >= 0.0 && mort0 <= 1.0) mort[0] = mort0;
if (mort1 >= 0.0 && mort1 <= 1.0) mort[1] = mort1;
}

float Cell::getMort(short period) {
float m = 0.0;
if (period == 0 || period == 1) m = mort[period];
return m;
}

#endif // SPATIALMORT

// Functions to handle costs for SMS

int Cell::getCost(void) {
int c;
if (smsData == 0) c = 0; // costs not yet set up
else c = smsData->cost;
return c;
}

void Cell::setCost(int c) {
if (smsData == 0) {
	smsData = new smscosts;
	smsData->effcosts = 0;
}
smsData->cost = c;
}

// Reset the cost and the effective cost of the cell
void Cell::resetCost(void) {
if (smsData != 0) { resetEffCosts(); delete smsData; }
smsData = 0;
}

array3x3f Cell::getEffCosts(void) {
array3x3f a;
if (smsData == 0 || smsData->effcosts == 0) { // effective costs have not been calculated
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			a.cell[i][j] = -1.0;
		}
	}
}
else
	a = *smsData->effcosts;
return a;
}

void Cell::setEffCosts(array3x3f a) {
if (smsData->effcosts == 0) smsData->effcosts = new array3x3f;
*smsData->effcosts = a;
}

// Reset the effective cost, but not the cost, of the cell
void Cell::resetEffCosts(void) {
if (smsData != 0) {
	if (smsData->effcosts != 0) {
		delete smsData->effcosts;
		smsData->effcosts = 0;
	}
}
}

void Cell::resetVisits(void) { visits = 0; }
void Cell::incrVisits(void) { visits++; }
unsigned long int Cell::getVisits(void) { return visits; }

#if RS_CONTAIN

//void Cell::setDamage(unsigned int d) { damage = d; }
//unsigned int Cell::getDamage(void) { return damage; }
void Cell::setDamage(DamageLocn *pDmg) { pDamage = pDmg; }
DamageLocn* Cell::getDamage(void) { return pDamage; }

#endif // RS_CONTAIN

#if SPATIALDEMOG

void Cell::addchgDemoScaling(std::vector<float> ds) {
	std::for_each(ds.begin(), ds.end(), [](float& perc){ if(perc < 0.0 || perc > 100.0) perc=100; });
	demoScalings.push_back(ds);
	return;
}

std::vector<float> Cell::getDemoScaling(short chgyear) {
	if (chgyear < 0 || chgyear >= (int)demoScalings.size()) {
		std::vector<float> ret(1, -1);
		return ret;
	}
	else return demoScalings[chgyear];
}

#endif // SPATIALDEMOG


//---------------------------------------------------------------------------

// Initial species distribution cell functions

DistCell::DistCell(int xx,int yy) {
x = xx; y = yy; initialise = false;
}

DistCell::~DistCell() {

}

void DistCell::setCell(bool init) {
initialise = init;
}

bool DistCell::toInitialise(locn loc) {
if (loc.x == x && loc.y == y) return initialise;
else return false;
}

bool DistCell::selected(void) { return initialise; }

locn DistCell::getLocn(void) {
	locn loc; loc.x = x; loc.y = y; return loc;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------




