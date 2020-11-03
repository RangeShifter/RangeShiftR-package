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

#include "FractalGenerator.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

vector<land> patches;

//----- Landscape creation --------------------------------------------------

land::land(): x_coord(0), y_coord(0), value(0.0), avail(0) {}

bool compare(const land& z, const land& zz) //compares only the values of the cells
{
return z.value < zz.value;
}

vector<land>& fractal_landscape(int X,int Y,double Hurst,double prop,
	double maxValue,double minValue)
{

int ii, jj, x, y;
int ix, iy;
//int x0, y0, size, kx, kx2, ky, ky2;
int kx,kx2,ky,ky2;

double range; //range to draw random numbers at each iteration
double nx, ny;
double i, j;
int Nx = X;
int Ny = Y;

double ran[5]; // to store each time the 5 random numbers for the random displacement

int Nno; // number of cells NON suitable as habitat

// exponents used to obtain the landscape dimensions
double pow2x = log((double)(X-1))/log(2.0);
double pow2y = log((double)(Y-1))/log(2.0);

double **arena = new double *[X];
for(ii = 0; ii < X; ii++) {
	arena[ii] = new double[Y];
}

patches.clear();
// initialise all the landscape with zeroes
for (jj = 0; jj < X; jj++) {
	for (ii = 0; ii < Y; ii++) {
		arena[jj][ii]=0;
	}
}

// initialisation of the four corners
arena[0][0]     = 1.0 + pRandom->Random() * (maxValue-1.0);
arena[0][Y-1]   = 1.0 + pRandom->Random() * (maxValue-1.0);
arena[X-1][0]   = 1.0 + pRandom->Random() * (maxValue-1.0);
arena[X-1][Y-1] = 1.0 + pRandom->Random() * (maxValue-1.0);

/////////////MIDPOINT DISPLACEMENT ALGORITHM//////////////////////////////////
kx = (Nx-1) / 2;
kx2 = 2 * kx;
ky = (Ny-1) / 2;
ky2 = 2 * ky;

for (ii = 0; ii < 5; ii++)  //random displacement
{
	ran[ii] = 1.0 + pRandom->Random() * (maxValue-1.0);
}

//The diamond step:
arena[kx][ky] = ((arena[0][0] + arena[0][ky2] + arena[kx2][0] + arena[kx2][ky2])/4) + ran[0];

//The square step:
//left
arena[0][ky] = ((arena[0][0] +arena[0][ky2] + arena[kx][ky]) / 3) + ran[1];
//top
arena[kx][0] = ((arena[0][0] + arena[kx][ky] + arena[kx2][0]) / 3) + ran[2];
//right
arena[kx2][ky] = ((arena[kx2][0] + arena[kx][ky] + arena[kx2][ky2]) / 3) + ran[3];
//bottom
arena[kx][ky2] = ((arena[0][ky2] + arena[kx][ky] +arena[kx2][ky2]) / 3) + ran[4];

range = maxValue*pow(2,-Hurst);

i = pow2x-1;
j = pow2y-1;

while (i > 0) {
	nx = pow(2,i)+1;
	kx = (nx-1) / 2;
	kx2 = 2 * kx;

	ny = pow(2,j)+1;
	ky = (ny-1) / 2;
	ky2 = 2 * ky;

	ix = 0;
	while (ix <= (Nx-nx)) {
		iy = 0;
		while (iy <= (Ny-ny)) {
			for (ii = 0; ii < 5; ii++)  //random displacement
			{
				ran[ii] = (int)(pRandom->Random() * 2.0 * range - range);
			}
			//The diamond step:

			arena[ix+kx][iy+ky] = ((arena[ix][iy] + arena[ix][iy+ky2] + arena[ix+ky2][iy]
						  + arena[ix+kx2][iy+ky2])/ 4) + ran[0];
			if (arena[ix+kx][iy+ky] < 1) 	arena[ix+kx][iy+ky] = 1;

			//The square step:
			//left
			arena[ix][iy+ky] =((arena[ix][iy] +arena[ix][iy+ky2] + arena[ix+kx][iy+ky])/3)
					   + ran[1];
			if (arena[ix][iy+ky] < 1) arena[ix][iy+ky] = 1;
			//top
			arena[ix+kx][iy] =((arena[ix][iy] + arena[ix+kx][iy+ky] + arena[ix+kx2][iy])/3)
					  + ran[2];
			if (arena[ix+kx][iy] < 1) arena[ix+kx][iy] = 1;
			//right
			arena[ix+kx2][iy+ky] = ((arena[ix+kx2][iy] + arena[ix+kx][iy+ky] +
								arena[ix+kx2][iy+ky2]) / 3) + ran[3];
			if (arena[ix+kx2][iy+ky] < 1) arena[ix+kx2][iy+ky] = 1;
			//bottom
			arena[ix+kx][iy+ky2] = ((arena[ix][iy+ky2] + arena[ix+kx][iy+ky] +
								arena[ix+kx2][iy+ky2]) / 3) + ran[4];
			if (arena[ix+kx][iy+ky2] < 1) arena[ix+kx][iy+ky2] = 1;

			iy += (ny-1);
		}
		ix += (nx-1);
	}
	if (i==j) j--;
	i--;

	range = range*pow(2,-Hurst); //reduce the random number range
}

// Now all the cells will be sorted and the Nno cells with the lower carrying
// capacity will be set as matrix, i.e. with K = 0

land *patch;

for (x = 0; x < X; x++) // put all the cells with their values in a vector
{
	for (y = 0; y < Y; y++)
	{
		patch = new land;
		patch->x_coord = x;
		patch->y_coord = y;
		patch->value = arena[x][y];
		patch->avail = 1;

		patches.push_back(*patch);

		delete patch;
	}
}


sort(patches.begin(),patches.end(),compare);  // sorts the vector

Nno = (prop*X*Y);
for (ii = 0; ii < Nno; ii++)
{
	patches[ii].value = 0.0;
	patches[ii].avail = 0;
}

double min = patches[Nno].value;        // variables for the rescaling
double max = patches[X*Y-1].value;

double diff = max - min;
double diffK = maxValue-minValue;
double new_value;

vector<land>::iterator iter = patches.begin();
while (iter != patches.end())
{
	if (iter->value > 0) // rescale to a range of K between Kmin and Kmax
	{
		new_value = maxValue - diffK * (max - iter->value) / diff;

		iter->value = new_value;
	}
	else iter->value = 0;

	iter++;
}

if (arena != NULL) {
#if RSDEBUG
//DebugGUI(("fractal_landscape(): arena=" + Int2Str((int)arena)
//	+ " X=" + Int2Str(X) + " Y=" + Int2Str(Y)
//	).c_str());
#endif
	for(ii = 0; ii < X; ii++) {
#if RSDEBUG
//DebugGUI(("fractal_landscape(): ii=" + Int2Str(ii)
//	+ " arena[ii]=" + Int2Str((int)arena[ii])
//	).c_str());
#endif
		delete[] arena[ii];
	}
	delete[] arena;
}

return patches;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

