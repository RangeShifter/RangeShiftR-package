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

#include "Patch.h"
//---------------------------------------------------------------------------
#if RS_EMBARCADERO
#pragma package(smart_init)
#endif

//---------------------------------------------------------------------------

#if SEASONAL
Patch::Patch(int seqnum,int num,short nseasons) 
#else
Patch::Patch(int seqnum,int num) 
#endif // SEASONAL
{
patchSeqNum = seqnum; patchNum = num; nCells = 0;
xMin = yMin = 999999999; xMax = yMax = 0; x = y = 0;
subCommPtr = 0;
#if RS_CONTAIN
damageIndex = 0.0;
//prevDamage = 0.0;		
damageLocns = false;
#endif // RS_CONTAIN 
#if SEASONAL
for (int i = 0; i < nseasons; i++) localK.push_back(0.0);
#else
localK = 0.0;
#endif // SEASONAL 
for (int sex = 0; sex < NSEXES; sex++) {
	nTemp[sex] = 0;
}
#if SPATIALDEMOG
//for (int i = 0; i < nDSlayer; i++) localDemoScaling.push_back(1.0);
localDemoScaling.assign(nDSlayer,1.0);
#endif

changed = false;
}

Patch::~Patch() {
cells.clear();
popns.clear();
#if SEASONAL
localK.clear();
#endif // SEASONAL 
#if SPATIALDEMOG
localDemoScaling.clear();
#endif // SPATIALDEMOG

}

int Patch::getSeqNum(void) { return patchSeqNum; }

int Patch::getPatchNum(void) { return patchNum; }

int Patch::getNCells(void) { return nCells; }

patchLimits Patch::getLimits(void) {
patchLimits p;
p.xMin = xMin; p.xMax = xMax; p.yMin = yMin; p.yMax = yMax;
return p;
}

// Does the patch fall (partially) within a specified rectangle?
bool Patch::withinLimits(patchLimits rect){
locn loc;
if (xMin <= rect.xMax && xMax >= rect.xMin &&  yMin <= rect.yMax && yMax >= rect.yMin) {
			// patch is within the rectangle UNLESS it is irregular in shape and lies at a corner
			// of the rectangle
			if ((xMin >= rect.xMin && xMax <= rect.xMax)
			||  (yMin >= rect.yMin && yMax <= rect.yMax)) {
				// patch lies within or along an edge of the initialistaion rectangle
				return true;
			}
			else {
				// check for any cell of the patch lying within the rectangle
				int ncells = (int)cells.size();
				for (int i = 0; i < ncells; i++) {
					loc = getCellLocn(i);
					if (loc.x >= rect.xMin && loc.x <= rect.xMax
					&&  loc.y >= rect.yMin && loc.y <= rect.yMax) {
						// cell lies within the rectangle
						return true;
					}
				}
			}
		}
return false;
}

// Reset minimum and maximum co-ordinates of the patch if it has been changed
void Patch::resetLimits(void) {
if (changed) {
	// remove any deleted cells
	std::vector <Cell*> newcells; // for all retained and added cells
	int ncells = (int)cells.size();
	for (int i = 0; i < ncells; i++) {
		if (cells[i] != NULL) {
			newcells.push_back(cells[i]);
		}
	}
	cells.clear();
	cells = newcells;
	// reset patch limits
	locn loc;
	xMin = yMin = 999999999; xMax = yMax = 0;
	ncells = (int)cells.size();
	for (int i = 0; i < ncells; i++) {
		loc = getCellLocn(i);
		if (loc.x < xMin) xMin = loc.x;
		if (loc.x > xMax) xMax = loc.x;
		if (loc.y < yMin) yMin = loc.y;
		if (loc.y > yMax) yMax = loc.y;
	}
	changed = false;
}
}

// Add a cell to the patch
void Patch::addCell(Cell* pCell,int x,int y) {
	cells.push_back(pCell);
	nCells++;
	if (x < xMin) xMin = x;
	if (x > xMax) xMax = x;
	if (y < yMin) yMin = y;
	if (y > yMax) yMax = y;
}

// Calculate the total carrying capacity (no. of individuals) and
// centroid co-ordinates of the patch
void Patch::setCarryingCapacity(Species *pSpecies,patchLimits landlimits,
	float epsGlobal,short nHab,short rasterType,short landIx,bool gradK) {
envStochParams env = paramsStoch->getStoch();
//Cell *pCell;
locn loc;
int xsum,ysum;
short hx;
float k,q,envval;

#if SEASONAL
int nseasons = (int)localK.size();
int *nsuitable = new int[nseasons];
for (int s = 0; s < nseasons; s++) {
	localK[s] = 0.0;
	nsuitable[s] = 0;
}
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " nseasons=" << nseasons 
//	<< endl;
#endif
#else
localK = 0.0; // no. of suitable cells (unadjusted K > 0) in the patch
int nsuitable = 0;
#endif // SEASONAL 
double mean;

#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " xMin=" << xMin << " yMin=" << yMin << " xMax=" << xMax << " yMax=" << yMax
//	<< endl;
#endif

if (xMin > landlimits.xMax || xMax < landlimits.xMin
||  yMin > landlimits.yMax || yMax < landlimits.yMin) {
	// patch lies wholely outwith current landscape limits
#if !SEASONAL
	// NB the next statement is unnecessary, as localK has been set to zero above
	//    retained only for consistency in standard variant
	localK = 0.0;
#endif // !SEASONAL 
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " localK=" << localK
//	<< endl;
#endif
	return;
}

int ncells = (int)cells.size();
xsum = ysum = 0;
for (int i = 0; i < ncells; i++) {
	if (gradK) // gradient in carrying capacity
		envval = cells[i]->getEnvVal(); // environmental gradient value
	else envval = 1.0; // no gradient effect
	if (env.stoch && env.inK) { // environmental stochasticity in K
		if (env.local) {
//			pCell = getRandomCell();
//			if (pCell != 0) envval += pCell->getEps();
			envval += cells[i]->getEps();
		}
		else { // global stochasticity
			envval += epsGlobal;
    }
	}
#if SEASONAL
	for (int s = 0; s < nseasons; s++) {
		switch (rasterType) {
		case 0: // habitat codes
			hx = cells[i]->getHabIndex(landIx);
			k = pSpecies->getHabK(hx,s);
			if (k > 0.0) {
				(nsuitable[s])++;
				localK[s] += envval * k;
			}
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " i=" << i << " s=" << s << " hx=" << hx << " k=" << k << " localK[s]=" << localK[s]
//	<< endl;
#endif
			break;
		case 1: // cover %
			k = 0.0;
			for (int j = 0; j < nHab; j++) { // loop through cover layers
				q = cells[i]->getHabitat(j);
				k += q * pSpecies->getHabK(j,s) / 100.0;
			}
			if (k > 0.0) {
				(nsuitable[s])++;
				localK[s] += envval * k;
			}
			break;
		case 2: // habitat quality
			q = cells[i]->getHabitat(landIx);
			if (q > 0.0) {
				(nsuitable[s])++;
				localK[s] += envval * pSpecies->getHabK(0,s) * q / 100.0;
			}
			break;
		}
	}
#else
	switch (rasterType) {
	case 0: // habitat codes
		hx = cells[i]->getHabIndex(landIx);
		k = pSpecies->getHabK(hx);
		if (k > 0.0) {
			nsuitable++;
			localK += envval * k;
		}
		break;
	case 1: // cover %
		k = 0.0;
		for (int j = 0; j < nHab; j++) { // loop through cover layers
			q = cells[i]->getHabitat(j);
			k += q * pSpecies->getHabK(j) / 100.0f;
		}
		if (k > 0.0) {
			nsuitable++;
			localK += envval * k;
		}
		break;
	case 2: // habitat quality
		q = cells[i]->getHabitat(landIx);
		if (q > 0.0) {
			nsuitable++;
			localK += envval * pSpecies->getHabK(0) * q / 100.0f;
		}
		break;
	}
#endif // SEASONAL 
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " i=" << i << " hx=" << hx << " q=" << q << " k=" << k << " localK=" << localK
//	<< endl;
#endif
	loc = cells[i]->getLocn();
	xsum += loc.x; ysum += loc.y;
}
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " epsGlobal=" << epsGlobal << " localK=" << localK
//	<< endl;
#endif
// calculate centroid co-ordinates
if (ncells > 0) {
	mean = (double)xsum / (double)ncells;
	x = (int)(mean + 0.5);
	mean = (double)ysum / (double)ncells;
	y = (int)(mean + 0.5);
}
if (env.stoch && env.inK) { // environmental stochasticity in K
	// apply min and max limits to K over the whole patch
	// NB limits have been stored as N/cell rather than N/ha
	float limit;
#if SEASONAL
	for (int s = 0; s < nseasons; s++) {
		limit = pSpecies->getMinMax(0) * (float)(nsuitable[s]);
		if (localK[s] < limit) localK[s] = limit;
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " limit=" << limit << " localK=" << localK
//	<< endl;
#endif
		limit = pSpecies->getMinMax(1) * (float)(nsuitable[s]);
		if (localK[s] > limit) localK[s] = limit;		
	}
#else
	limit = pSpecies->getMinMax(0) * (float)nsuitable;
	if (localK < limit) localK = limit;
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " limit=" << limit << " localK=" << localK
//	<< endl;
#endif
	limit = pSpecies->getMinMax(1) * (float)nsuitable;
	if (localK > limit) localK = limit;
#endif // SEASONAL 
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " limit=" << limit << " localK=" << localK
//	<< endl;
#endif
}
#if RSDEBUG
//DEBUGLOG << "Patch::setCarryingCapacity(): patchNum=" << patchNum
//	<< " localK=" << localK
//	<< endl;
#endif
}

#if RS_CONTAIN

void Patch::resetDamageIndex(void) { damageIndex = 0.0; }

void Patch::updateDamageIndex(int dmgX,int dmgY,int damage,double alpha) {
// the damage index for the patch is based on Hanski's index, i.e. it is the sum
// of damage values of each cell in the landscape weighted by a negative exponential
// function of distance between the cell and the patch (approx.) centroid

// if the location of the damage is within the patch, then take the distance to
// be zero - otherwise use Euclidean distance from centroid

bool inpatch = false;
int ncells = (int)cells.size();
locn loc;
double dist = 0.0;
for (int i = 0; i < ncells; i++) {
	loc = cells[i]->getLocn();
	if (loc.x == dmgX && loc.y == dmgY) { inpatch = true; i = ncells+1; break; }
}
if (!inpatch) dist = sqrt((double)((dmgX-loc.x)*(dmgX-loc.x) + (dmgY-loc.y)*(dmgY-loc.y)));
damageIndex += (double)damage * exp(-1.0*alpha*dist);

}	

double Patch::getDamageIndex(void) { return damageIndex; }
/*
void Patch::resetDamageLocns(void) {
int ndlocns = (int)dmglocns.size();
for (int i = 0; i < ndlocns; i++)
	if (dmglocns[i] != NULL) dmglocns[i]->resetDamageLocn();
}

void Patch::setDamage(int xx,int yy,int dmg) {
dmglocns.push_back(new DamageLocn(x,y,dmg));
}

double Patch::totalDamage(void) {
double totdmg = 0.0;
int ndlocns = (int)dmglocns.size();
for (int i = 0; i < ndlocns; i++)
	if (dmglocns[i] != NULL) totdmg += dmglocns[i]->getDamageIndex();
return totdmg;
}
*/
void Patch::setDamageLocns(bool d) { damageLocns = d; }
bool Patch::hasDamageLocns(void) { return damageLocns; }

//void Patch::setPrevDamage(double dmg) { if (dmg >= 0.0) prevDamage = dmg; }
//double Patch::getPrevDamage(void) { return prevDamage; }      
double Patch::getPrevDamage(void) {
double damage = 0.0;
DamageLocn *pDamage; 
int ncells = (int)cells.size();
for (int i = 0; i < ncells; i++) {
	pDamage = cells[i]->getDamage();
	if (pDamage != 0) {
		damage += pDamage->getDamageIndex(false);
	}
}
return damage; 
}      
	
#endif // RS_CONTAIN 


#if SEASONAL
float Patch::getK(int season) { 
if (season >= 0 && season < (int)localK.size()) return localK[season];
else return 0.0; 
}
bool Patch::suitableInAllSeasons(void) {
int nseasons = (int)localK.size();
bool ok = true;
for (int i = 0; i < nseasons; i++) {
	if (localK[i] <= 0.0) ok = false;
}
return ok;
}
#else
float Patch::getK(void) { return localK; }
#endif // SEASONAL

#if SPATIALDEMOG
void Patch::setDemoScaling(std::vector <float> ds) {

	std::for_each(ds.begin(), ds.end(), [](float& perc){ if(perc < 0.0 || perc > 1.0) perc=1; });

	localDemoScaling.assign(ds.begin(), ds.end());

	return;
}

std::vector <float> Patch::getDemoScaling(void) { return localDemoScaling; }

void Patch::setPatchDemoScaling(short landIx, patchLimits landlimits) { 

	// if patch wholly outside current landscape boundaries
	if (xMin > landlimits.xMax || xMax < landlimits.xMin
	||  yMin > landlimits.yMax || yMax < landlimits.yMin) {
		localDemoScaling.assign(nDSlayer,0.0); // set all local scales to zero
		return;
	}
	
	// loop through constituent cells of the patch
	int ncells = (int)cells.size();
	std::vector<float> patchDS(nDSlayer, 0.0);
	std::vector<float> cellDS(nDSlayer, 0.0);

	for (int i = 0; i < ncells; i++) {
		cellDS = cells[i]->getDemoScaling(landIx); // is that ok?

		//add cell value to patch value 
		for (int ly = 0; ly < nDSlayer; ly++) {
			patchDS[ly] += cellDS[ly];
		}
	}
	
	// take mean over cells and divide by 100 to scale to range [0,1]
	for (int ly = 0; ly < nDSlayer; ly++) {
		patchDS[ly] = patchDS[ly] / ncells / 100.0f;
	}
	
	// set values
	setDemoScaling(patchDS);

	return;
}
#endif //SPATIALDEMOG

// Return co-ordinates of a specified cell
locn Patch::getCellLocn(int ix) {
locn loc; loc.x = -666; loc.y = -666;
int ncells = (int)cells.size();
if (ix >= 0 && ix < ncells) {
	loc = cells[ix]->getLocn();
}
return loc;
}
// Return pointer to a specified cell
Cell* Patch::getCell(int ix) {			
int ncells = (int)cells.size();
if (ix >= 0 && ix < ncells) return cells[ix];
else return 0;
}
// Return co-ordinates of patch centroid
locn Patch::getCentroid(void) {
locn loc; loc.x = x; loc.y = y;
return loc;
}

// Select a Cell within the Patch at random, and return pointer to it
// For a cell-based model, this will be the only Cell
Cell* Patch::getRandomCell(void) {
Cell *pCell = 0;
int ix;
int ncells = (int)cells.size();
if (ncells > 0) {
	if (ncells == 1) ix = 0;
	else ix = pRandom->IRandom(0,ncells-1);
	pCell = cells[ix];
}
return pCell;
}

// Remove a cell from the patch
void Patch::removeCell(Cell* pCell) {
int ncells = (int)cells.size();
for (int i = 0; i < ncells; i++) {
	if (pCell == cells[i]) {
		cells[i] = NULL; i = ncells;
		nCells--;
		changed = true;
	}
}
}

void Patch::setSubComm(intptr sc)
{ subCommPtr = sc; }

// Get pointer to corresponding Sub-community (cast as an integer)
intptr Patch::getSubComm(void)
{ return subCommPtr; }

void Patch::addPopn(patchPopn pop) {
popns.push_back(pop);
}

// Return pointer (cast as integer) to the Population of the specified Species
intptr Patch::getPopn(intptr sp)
{
int npops = (int)popns.size();
for (int i = 0; i < npops; i++) {
	if (popns[i].pSp == sp) return popns[i].pPop;
}
return 0;
}

void Patch::resetPopn(void) {
popns.clear();
}

void Patch::resetPossSettlers(void) {
for (int sex = 0; sex < NSEXES; sex++) {
	nTemp[sex] = 0;
}
}

// Record the presence of a potential settler within the Patch
void Patch::incrPossSettler(Species *pSpecies,int sex) {
#if RSDEBUG
//DEBUGLOG << "Patch::incrPossSettler(): 5555: patchNum = " << patchNum
//	<< " sex = " << sex << endl;
#endif
// NOTE: THE FOLLOWING OPERATION WILL NEED TO BE MADE SPECIES-SPECIFIC...
if (sex >= 0 && sex < NSEXES) {
	nTemp[sex]++;
}
}

// Get number of a potential settlers within the Patch
int Patch::getPossSettlers(Species *pSpecies,int sex) {
#if RSDEBUG
//DEBUGLOG << "Patch::getPossSettlers(): 5555: patchNum = " << patchNum
//	<< " sex = " << sex << endl;
#endif
// NOTE: THE FOLLOWING OPERATION WILL NEED TO BE MADE SPECIES-SPECIFIC...
if (sex >= 0 && sex < NSEXES) return nTemp[sex];
else return 0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------



