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
# Plotting parameterised probabilities functions
# (the respective methods are defined together with their class)
#
# -----

#' Plot parameterised probabilities
#'
#' Visualise the dependencies or distributions of probabilities of the various processes that are defined by the different parameter modules,
#' like e.g. fecundity or mortality.
#'
#' @param x a RangeShiftR parameter object
#' @param xmax,ymax upper plotting limits (lower limits are fixed to \eqn{0})
#' @param stage,sex stage(s) and sexe(s) to plot, default: all
#' @param ... various options, depend on the given parameter module type
#' @details
#' Available methods and their options:
#' \itemize{
#'   \item \code{\link[RangeShiftR]{Emigration}}: plot emigration probability
#'   \item \code{\link[RangeShiftR]{DispersalKernel}}:
#'   \itemize{
#'     \item \code{mortality=FALSE} - plot dispersal distance probability density  (default)
#'     \item \code{mortality= TRUE} - plot mortality probability
#'   }
#'   If a mixed kernel was defined (i.e. \code{DoubleKernel=TRUE}), plot the resulting dispersal probability by...
#'   \itemize{
#'     \item \code{combinekernels=FALSE} - ...plotting both kernels separately (default)
#'     \item \code{combinekernels= TRUE} - ...combining both kernels, i.e. \ifelse{html}{ \out{p(d; &delta;<sub>1</sub>,&delta;<sub>2</sub>) = p<sub>I</sub> p(d;&delta;<sub>1</sub>) + (1-p<sub>I</sub>) p(d;&delta;<sub>1</sub>) } }{\deqn{ p(d; δ_1,δ_2) = p_I p(d;δ_1) + (1-p_I) p(d;δ_2)} }
#'   }
#'   \item \code{\link[RangeShiftR]{StageStructure}}: plot fecundity as well as survival and development probabilities
#' }
#' @export
setGeneric("plotProbs", function(x, ...) standardGeneric("plotProbs") )