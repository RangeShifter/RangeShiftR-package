#---------------------------------------------------------------------------
#
#	Copyright (C) 2020-2022 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
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
#'   \item \code{\link[RangeShiftR]{StageStructure}}: plot fecundity as well as survival and development probabilities.
#'
#'   Survival and development are calculated based on the transition matrix. Consider e.g. a 3 stage matrix with a transition matrix of
#'
#'   \tabular{ccccc}{0 \tab | \tab 0 \tab | \tab \ifelse{html}{\out{&phi;<sub>2</sub>}}{\eqn{φ_2}} \cr
#'   \eqn{1.0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1</sub> (1-&gamma;<sub>1-2</sub>)}}{\eqn{σ_1 (1 − γ_(1-2))}} \tab | \tab \eqn{0} \cr
#'   \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1</sub> &gamma;<sub>1-2</sub>}}{\eqn{σ_1 γ_(1-2)}} \tab | \tab \ifelse{html}{\out{&sigma;<sub>2</sub>}}{\eqn{σ_2}} \cr}
#'
#'   The survival probability is calculated as the sum of the probability to stay in the same stage and the probability to reach the next stage.
#'   E.g. for stage 1:  \ifelse{html}{\out{&sigma;<sub>1</sub> = sum( &sigma;<sub>1</sub> (1-&gamma;<sub>1-2</sub>),  &sigma;<sub>1</sub> &gamma;<sub>1-2</sub>}}{\eqn{σ_1 = sum(σ_1 (1 − γ_(1-2)), σ_1 γ_(1-2))}}
#'
#'   The development probability of stage 1 to stage 2 is the ratio of the probability to reach stage 2 and the previously calculated survival probability.
#'   E.g. for stage 1: \ifelse{html}{\out{&gamma;<sub>1-2</sub> = &sigma;<sub>1</sub> &gamma;<sub>1-2</sub> / &sigma;<sub>1</sub>}}{\eqn{γ_(1-2) = σ_1 / (σ_1 γ_(1-2) )}}
#' }
#' @export
setGeneric("plotProbs", function(x, ...) standardGeneric("plotProbs") )

