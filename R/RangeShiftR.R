#---------------------------------------------------------------------------
#
#	Copyright (C) 2020-2021 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell
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


#' @details
#' \code{RangeShiftR} provides flexible and fast simulations of spatial eco-evolutionary dynamics.
#' It implements the individual-based simulation software \emph{RangeShifter},
#' \insertCite{bocedi2014,bocedi2020}{RangeShiftR},
#' making it available for all machines that can build and run \emph{R}.
#'
#' \code{RangeShiftR} models the processes of demography, dispersal and evolution in an
#' inter-dependent way, offering substantial complexity in the corresponding modelling choices.
#' It is entirely open-source and aims to facilitate the application
#' of individual-based and mechanistic modelling to eco-evolutionary questions.
#'
#' The \emph{RangeShifter} project is presented on our website
#' \url{https://rangeshifter.github.io/}.
#'
#' To get started with the package, please find an overview and a range of tutorials here:
#' \url{https://rangeshifter.github.io/RangeshiftR-tutorials/}.
#'
#' @references
#'         \insertAllCited{}
#' @keywords internal
"_PACKAGE"


#' @useDynLib RangeShiftR
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
NULL


.onAttach <- function(libname, pkgname) {
    packageStartupMessage("RangeshiftR version 1.0.2 (11.01.2021)\n",
                          "Copyright (C) 2020-2021 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell\n\n",
                          "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.\n",
                          "You are welcome to redistribute it and/or modify it under certain conditions;\n",
                          "type 'RangeShiftR_license()' for details.\n")
}

.onUnload <- function (libpath) {
    library.dynam.unload("RangeShiftR", libpath)
}


#' Display RangeShiftR license information
#'
#' @export
RangeShiftR_license <- function ()
{
    cat("\nRangeshiftR version 1.0.3 (13.04.2021)\n")
    cat("Copyright (C) 2020-2021 Anne-Kathleen Malchow, Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Damaris Zurell\n\n")
    cat("This program is free software: you can redistribute it and/or modify\n")
    cat("it under the terms of the GNU General Public License as published by\n")
    cat("the Free Software Foundation, either version 3 of the License, or\n")
    cat("(at your option) any later version.\n\n")
    cat("This software is distributed in the hope that it will be useful,\n")
    cat("but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
    cat("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n")
    cat("GNU General Public License for more details.\n\n")
    cat("You should have received a copy of the GNU General Public License\n")
    cat("along with this software. Copies of the license can also be found at\n")
    cat("http://www.gnu.org/licenses/\n")
    cat("\n")
    cat("'Share and Enjoy.'\n\n")
}


#' Display citation
#'
#' @export
RangeShiftR_citation <- function ()
{
    citation(package = "RangeShiftR", lib.loc = NULL, auto = NULL)
}
