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
 
//Last updated: 26 November 2020 by Greta Bocedi

//---------------------------------------------------------------------------

#ifndef VersionH
#define VersionH

#define RSDEBUG 0

#define LINUX_CLUSTER 0
//#define RSWIN64 0

#define RANDOMCHECK 0

#define BATCH 1
#define VCL 0

#define RS_RCPP 1
//#define R_CMD 0

#define RS_EMBARCADERO 0

//---------------------------------------------------------------------------

#define RS_THREADSAFE 1

#define SPATIALDEMOG 1

	#if SPATIALDEMOG
		#define RS_THREADSAFE 1
	#endif // SPATIALDEMOG
	#if RS_THREADSAFE
		#define RS_RCPP 1
		#define RSDEBUG 0
	#endif // RS_THREADSAFE

//---------------------------------------------------------------------------


#endif
