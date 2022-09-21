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

RangeShifter v2.0 Cell

Implements the following classes:

Cell - Landscape cell

DistCell - Initial species distribution cell

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species’ responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 14 January 2021 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef CellH
#define CellH

#if SPATIALDEMOG
#include <algorithm>
#endif
#include <vector>

using namespace std;

#include "Parameters.h"
#if RS_CONTAIN
#include "Control.h"
#endif // RS_CONTAIN 

//---------------------------------------------------------------------------

struct array3x3f { float cell[3][3]; }; 	// neighbourhood cell array (SMS)
struct smscosts { int cost; array3x3f *effcosts; };	// cell costs for SMS

// Landscape cell

class Cell{
public:
	Cell( // Constructor for habitat codes
		int,				// x co-ordinate
		int,				// y co-ordinate
		intptr,			// pointer (cast as integer) to the Patch to which Cell belongs
		int					// habitat index number
	);
	Cell( // Constructor for habitat % cover or habitat quality
		int,				// x co-ordinate
		int,				// y co-ordinate
		intptr,			// pointer (cast as integer) to the Patch to which Cell belongs
		float				// habitat proportion or cell quality score
	);
	~Cell();
	void setHabIndex(
		short	// habitat index number
	);
	void changeHabIndex(
		short,	// landscape change number
		short		// habitat index number
	);
	int getHabIndex(
		int		// landscape change number
	);
	int nHabitats(void);
	void setHabitat(
		float	// habitat proportions or cell quality score
	);
	float getHabitat( // Get habitat proportion / quality score
		int		// habitat index number / landscape change number
	);
	void setPatch(
		intptr		// pointer (cast as integer) to the Patch to which Cell belongs
	);
	intptr getPatch(void);
	locn getLocn(void);
	void setEnvDev(
		float	// local environmental deviation
	);
	float getEnvDev(void);
	void setEnvVal(
		float	// environmental value
	);
	float getEnvVal(void);
	void updateEps( // Update local environmental stochasticity (epsilon)
		float,	// autocorrelation coefficient
		float		// random adjustment
	);
	float getEps(void);
#if SPATIALMORT
	void setMort(
		float,	// mortality time period 0
		float		// mortality time period 1
	);
	float getMort(
		short		// time period (0 or 1)
	);
#endif
	void setCost(
		int		// cost value for SMS
	);
	int getCost(void);
	void resetCost(void);
	array3x3f getEffCosts(void);
	void setEffCosts(
		array3x3f	// 3 x 3 array of effective costs for neighbouring cells
	);
	void resetEffCosts(void); // Reset the effective cost, but not the cost, of the cell
	void resetVisits(void);
	void incrVisits(void);
	unsigned long int getVisits(void);
#if RS_CONTAIN
//	void setDamage(unsigned int);
//	unsigned int getDamage(void);
	void setDamage(DamageLocn*);
	DamageLocn* getDamage(void);
#endif // RS_CONTAIN
#if SPATIALDEMOG
	void addchgDemoScaling(std::vector<float>);
	void setDemoScaling(std::vector<float>, short);
	std::vector<float> getDemoScaling(short);
#endif // SPATIALDEMOG 

private:
	int x,y;			// cell co-ordinates
	intptr pPatch; 	// pointer (cast as integer) to the Patch to which cell belongs
	// NOTE: THE FOLLOWING ENVIRONMENTAL VARIABLES COULD BE COMBINED IN A STRUCTURE
	// AND ACCESSED BY A POINTER ...
	float envVal; // environmental value, representing one of:
								// gradient in K, r or extinction probability
	float envDev;	// local environmental deviation (static, in range -1.0 to +1.0)
	float eps;		// local environmental stochasticity (epsilon) (dynamic, from N(0,std))
	unsigned long int visits; // no. of times square is visited by dispersers
#if RS_CONTAIN
	DamageLocn *pDamage;	// pointer to damage location (if any)
#endif // RS_CONTAIN 
#if SPATIALMORT
	float mort[2];	// additional mortality in two periods
#endif // SPATIALMORT 
	smscosts *smsData;

	vector <short> habIxx; 		// habitat indices (rasterType=0)
		// NB initially, habitat codes are loaded, then converted to index nos.
		//    once landscape is fully loaded
	vector <float> habitats;	// habitat proportions (rasterType=1) or quality (rasterType=2)
#if SPATIALDEMOG
	std::vector<std::vector<float>> demoScalings;	// demographic scaling layers (only if rasterType==2)
#endif
};

//---------------------------------------------------------------------------

// Initial species distribution cell

class DistCell{
public:
	DistCell(
		int,	// x co-ordinate
		int		// y co-ordinate
	);
	~DistCell();
	void setCell(
		bool	// TRUE if cell is to be initialised, FALSE if not
	);
	bool toInitialise(
		locn	// structure holding co-ordinates of cell
	);
	bool selected(void);
	locn getLocn(void);

private:
	int x,y;					// cell co-ordinates
	bool initialise;  // cell is to be initialised

};

#if RSDEBUG
extern void DebugGUI(string);
#endif

#if SPATIALDEMOG
extern short nDSlayer;
#endif // SPATIALDEMOG

//---------------------------------------------------------------------------

#endif
