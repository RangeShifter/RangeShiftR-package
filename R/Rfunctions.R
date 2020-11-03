#---------------------------------------------------------------------------
#	
#	Copyright (C) 2020 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
#	
#	This file is part of RangeShiftR.
#	
#	RangeShiftR is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#	
#	RangeShiftR is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#	GNU General Public License for more details.
#	
#	You should have received a copy of the GNU General Public License
#	along with RangeShiftR. If not, see <https://www.gnu.org/licenses/>.
#	
#----------------------------------------------------------------------------
 
 

# -----
#
# R-level functions
#
# -----


#-- define ClassUnions
#----------------------

#for Integer and Numerical so that 'Integer'-slots can also accept 'Numerical' input
setClassUnion("integer_OR_numeric", c("integer", "numeric"))

#for Matrix and Numerical so that 'Matrix'-slots can also accept 'Numerical' input when 1x1-Matrix is expected
setClassUnion("matrix_OR_numeric", c("matrix", "numeric"))

#for cost layer to accept habitat codes or raster map file
setClassUnion("numeric_OR_character", c("numeric", "character"))


#-- define error and warning messages
#----------------------
warn_msg_ignored = " will be ignored "


#-- density dependence function
#----------------------
densdep <- function(x, A0 = 1.0, alpha = 1.0, beta = 0.0) {
    A0/(1+exp(alpha*(beta-x)))
}


#-- validation function
#----------------------

#' Validate a given RS parameter object
#'
#' @export
validateRSparams <- function(x){validObject(x)}
