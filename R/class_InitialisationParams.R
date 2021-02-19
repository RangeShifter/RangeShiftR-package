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
 
 

#### CLASS INITIALISATIONPARAMS

# from RS 'Initialisation' file

#' Set Initialisation Parameters
#'
#' Initialisation rules define the way in which initial individuals are placed in the landscape at the start of the simulation.
#' There are three initialisation types:\cr
#' Free initialisation according to the habitat map (\code{InitType}\eqn{= 0}),
#' initialisation from a loaded species distribution map (\code{InitType} \eqn{= 1}),
#' and the initialisation from an initial individuals list (\code{InitType} \eqn{= 2}),
#' with additional options for each type, see the Details.\cr
#' Additionally, initial density and, if applicable, initial stage and age distributions can be set.
#'
#' @include plotProbs.R
#' @include Rfunctions.R
#' @usage Initialise(InitType = 0, FreeType = 1, SpType = 0, NrCells, InitIndsFile = "NULL",
#'            InitDens = 1, IndsHaCell, PropStages = 0, InitAge = 2, minX, minY, maxX, maxY,
#'            InitFreezeYear = 0, RestrictRows = 0, RestrictFreq = 0, FinalFreezeYear = 0)
#' @param InitType Type of initialisation:\cr
#' \eqn{0} = Free initialisation according to habitat map (default) (set \code{FreeType}), \cr
#' \eqn{1} = From loaded species distribution map (set \code{SpType}),\cr
#' \eqn{2} = From initial individuals list file (set \code{InitIndsFile}).
#' \cr Must to be \eqn{0} for an \code{\link[RangeShiftR]{ArtificialLandscape}}.
#' @param FreeType,NrCells Option for \emph{free initialisation}, i.e. required only if \code{InitType}\eqn{ = 0}:\cr
#' \code{FreeType} = \eqn{0}: Random; provide number of cells/patches to initialise in \code{NrCells}. \cr
#' \code{FreeType} = \eqn{1}: All suitable cells/patches (default).
#' @param SpType,NrCells Option for \emph{initialisation from species distribution map}, i.e. required only if \code{InitType}\eqn{ = 1}:\cr
#' \code{SpType} = \eqn{0}: All suitable cells within all distribution presence cells (default), \cr
#' \code{SpType} = \eqn{1}: All suitable cells within some randomly chosen presence cells; set number of cells to initialise in \code{NrCells}.
#' @param InitIndsFile Name of \emph{initial individuals list file}, required only if \code{InitType}\eqn{ = 2}.\cr
#' For informaton on the required file format see the Details below.
#' @param InitDens,IndsHaCell Number of individuals to be seeded in each cell/patch:\cr
#' \code{InitDens} = \eqn{0}: At \code{K_or_DensDep},\cr
#' \code{InitDens} = \eqn{1}: At half \code{K_or_DensDep} (default),\cr
#' \code{InitDens} = \eqn{2}: Set the number of individuals per cell/hectare to initialise in \code{IndsHaCell}.
#' @param PropStages For \code{\link[RangeShiftR]{StageStructure}}d models only: Proportion of individuals initialised in each stage.
#' Requires a vector of length equal to the number of stages; its entries must be \eqn{\ge 0} and sum to \eqn{1.0}. However, juveniles
#' can't be initialised and thus stage \eqn{0} (first entry) must have a value of \eqn{0.0}.
#' @param InitAge For \code{StageStructure}d models only: Initial age distribution within each stage:\cr
#' \eqn{0} = Minimum age for the respective stage.\cr
#' \eqn{1} = Age randomly sampled between the minimum and the maximum age for the respective stage.\cr
#' \eqn{2} = According to a quasi-equilibrium distribution (default).
#' @param minX,maxX,minY,maxY Option for \emph{free initialisation} (\code{InitType}\eqn{ = 0}): Restrict initial range in X- and/or Y- coordinates, given in number of cells.
#' All must be \eqn{\ge 0} and (\code{maxX,maxY})\eqn{>}(\code{minX,minY}). (Integer)
#' @param InitFreezeYear Option for \emph{free initialisation} (\code{InitType}\eqn{ = 0}): The year until which species is confined to
#' its initial range limits.
#' @param RestrictRows Option for \emph{free initialisation} (\code{InitType}\eqn{ = 0}): The number of rows at northern front to restrict range.
#' If set to \eqn{0}, the range restriction tured off.
#' @param RestrictFreq Option for \emph{free initialisation} (\code{InitType}\eqn{ = 0}): Frequency in years at which range is restricted to northern front.
#' @param FinalFreezeYear Option for \emph{free initialisation} (\code{InitType}\eqn{ = 0}): The year after which species is confined to its new, current range limits, after a
#' period of range expansion. Will be ignored if set to \eqn{0}, otherwise must be \eqn{>} \code{InitFreezeYear}.
#' @details ## Initialisation Types
#' \emph{Initialisation Types}\cr
#' \itemize{
#'     \item \emph{Free Initialisation.} (\code{InitType}\eqn{ = 0})\cr The population is initialised according to suitable habitat in the landscape.
#' Either all (\code{FreeType}\eqn{=1}), or a specified number (\code{NrCells}) of randomly selected (\code{FreeType}\eqn{=0}) cells/patches
#' will be seeded. When using an artificial landscape, this is the only option available.
#'     \item \emph{From loaded species distribution map.} (\code{InitType}\eqn{ = 1})\cr The population is initialised according to a loaded
#' species distribution map, which needs to be provided through the module \code{\link[RangeShiftR]{ImportedLandscape}}, option \code{SpDistFile}.
#' All habitat cells/patches within either all (\code{SpType}\eqn{=0}), or a specified number (\code{NrCells}) of randomly selected
#' (\code{SpType}\eqn{=1}) presence cells (which can have a lower resolution) specified by this distribution map will seeded.
#'     \item \emph{From initial individuals list file.} (\code{InitType}\eqn{ = 2})\cr The population is initialised according to a list of specific
#' individuals (of given sex, age and stage, if appropriate) in specified cells/patches. This option allows simulation of a reintroduction
#' scenario.\cr The list has to be loaded from a file in the path given by \code{InitIndsFile}. It must be a tab-seperated list with
#' explicit column headers and one row for each individual to be initialized. The expected column headers depend on the model settings and
#' must match the following order exactly: 'Year', 'Species', for cell-/patch-based: 'X', 'Y' / 'PatchID', 'Ninds', for sexual model: 'Sex',
#' for stage-structured population: 'Age', 'Stage'.
#' }
#'
#' \emph{Initialial density, stage, and age}\cr
#' For \code{InitType}\eqn{ = {0,1}}, the number of individuals that should be seeded in each cell or patch has to be set.
#' There are three options:
#' \itemize{
#'     \item \emph{At} \code{K_or_DensDep}. (\code{InitDens}\eqn{=0})\cr The cell/patch will be saturated at its respective \eqn{K} or \eqn{1/b}.
#'     \item \emph{At half} \code{K_or_DensDep}. (\code{InitDens}\eqn{=1})\cr The cell/patch will be saturated at half its \eqn{K} or \eqn{1/b}.
#'     \item \emph{Set value} \code{IndsHaCell}. (\code{InitDens}\eqn{=2})\cr Set the number of individuals to be seeded in each cell or the density
#' in each patch (in units of individuals per hectare).
#' }
#' (These settings have no effect for \code{InitType}\eqn{ = 2}.)
#'
#' In the case of \code{\link[RangeShiftR]{StageStructure}}d models, the initial stage and age distributions must be specified.
#' If \code{InitType}\eqn{ = 2}, this is done via the \code{InitIndsFile}, whereas for \code{InitType}\eqn{ = {0,1}},
#' the proportion of individuals that should be initialised at each stage class is set via the numeric vector \code{PropStages}. It needs
#' to have as many entries as number of stages, starting from the juvenile stage (\eqn{0}). Note that these proportions must sum up to \eqn{1.0},
#' however the proportion of juveniles must be \eqn{0.0}.
#'
#' Options for initial age distributions:
#' \itemize{
#'     \item \emph{Minimal age} of respective stage. (\code{InitAge}=0)
#'     \item \emph{Randomise.} (\code{InitAge}=1)\cr Individuals initialised in each stage will get an age randomly sampled between
#' the minimum and the maximum age for their respective stage.
#'     \item \emph{Quasi-equilibrium.} (\code{InitAge}=2)\cr Initial age distribution is set approximately in accordance with the number of years
#' taken to pass through the stage and the (female) survival rate of the respective stage.
#' }
#' (These settings have no effect for \code{InitType}\eqn{ = 2}.)
#'
#' \emph{Additional options for free initialisation}\cr
#' In the case of free initialisation, either all or a specified number of randomly selected cells/patches will be seeded.
#' It is possible to restrict the landscape area available for initialisation by setting limits to the x- and y- coordinates using
#' \code{minX,maxX,minY} and \code{maxY}).
#'
#' Additionally, for all free initialisation methods, it is possible to restrict the available landscape to the initialised range for
#' a specified period, specifically until year \code{InitFreezeYear}. This option may be particularly useful in an evolutionary model,
#' to allow a population to reach equilibrium in its original range before being allowed to expand into the remainder of the landscape.
#'
#' \emph{Restrict range to northern front} is a further option provided as an additional feature for models of range expansion. Whilst not
#' strictly an initialisation option, it is provided in this context as it will typically be applied together with restriction of the
#' initial range in theoretical eco-evolutional models, where the population is firstly allowed to evolve certain traits in a stable range,
#' then allowed to expand its range during a period of environmental change (when dispersal traits can come under strong selection), and
#' finally further restricted to a stable range (set \code{FinalFreezeYear}). During the period of range expansion, the population can build up to very large numbers
#' of individuals in the core of the range, but the area of interest, where strong selection occurs, lies close to the range front.
#' Therefore, individuals in the core can be discarded, and this is achieved by progressively restricting the suitable range to a number
#' of rows behind the current range front (set \code{RestrictRows}) so that local populations in the core go extinct. This occurs at
#' a frequency \code{RestrictFreq} specified in years. The parameter \code{FinalFreezeYear} sets the year after which the range will be
#' frozen to its new, current limits. It must either be zero, in which case range expansion continues until either the end of the simulation or the northern edge of the
#' landscape is reached, or set to a year after the year specified in \code{InitFreezeYear}.
#'
#' This feature is provided only for expansion in a northerly direction (increasing Y); the initial population should therefore be
#' established in the southern part of the landscape (low Y values). To turn it off, set \code{RestrictRows = 0}.
#
# For Instance a population is firstly allowed to evolve certain traits in a stable range of 1000 years (\code{InitFreezeYear}=1000).
# Then allowed to expand its range during a period of environmental change and is finally further restricted
# to a stable range in the year 5000 (\code{FinalFreezeYear}=5000). During the period of range expansion, the population
# can build up to very large numbers of individuals in the core of the range, but the area of interest, lies close
# to the range front (i.e. 50 rows). Therefore, individuals in the core can be discarded. This is achieved
# by progressively restricting the suitable range to a number of rows behind the front (\code{RestrictRows}=50),
# so that local populations outside of the area of interest go extinct. Note that \code{RestrictRows} is provided
# only for expansion in a northerly direction (increasing \eqn{y}-axis). Thus the initial population should be established in the
# southern part of the landscape (low value on the \eqn{y}-axis).
#
#' @examples init <- Initialise() # use all defaults (a stage-structured model requires setting PropStages)
#' @return a parameter object of class "InitialisationParams"
#' @author Anne-Kathleen Malchow
#' @name Initialise
#' @export Initialise
Initialise <- setClass("InitialisationParams", slots = c(InitType = "integer_OR_numeric",
                                                         FreeType = "integer_OR_numeric",
                                                         SpType = "integer_OR_numeric",
                                                         NrCells = "integer_OR_numeric",
                                                         InitIndsFile = "character",
                                                         InitDens = "integer_OR_numeric",
                                                         IndsHaCell = "integer_OR_numeric",
                                                         PropStages = "numeric",
                                                         InitAge = "integer_OR_numeric",
                                                         minX = "integer_OR_numeric",
                                                         minY = "integer_OR_numeric",
                                                         maxX = "integer_OR_numeric",
                                                         maxY = "integer_OR_numeric",
                                                         InitFreezeYear = "integer_OR_numeric",
                                                         RestrictRows = "integer_OR_numeric",
                                                         RestrictFreq = "integer_OR_numeric",
                                                         FinalFreezeYear = "integer_OR_numeric")
                       , prototype = list(InitType = 0L, #free initialisation
                                          FreeType = 1L, #all
                                          SpType = 0L,   #all
                                          #NrCells,
                                          InitIndsFile = "NULL",
                                          InitDens = 1L, #K/2
                                          #IndsHaCell,
                                          PropStages = 0.0,
                                          InitAge = 2L, #quasi-eq
                                          #minX,
                                          #minY,
                                          #maxX,
                                          #maxY,
                                          InitFreezeYear = 0L,
                                          RestrictRows = 0L,
                                          RestrictFreq = 0L,
                                          FinalFreezeYear = 0L)
)
# finish initialize

setValidity('InitialisationParams', function(object){
    msg <- NULL
    if (is.na(object@InitType) || length(object@InitType) == 0){
        msg <- c(msg, 'Type of initialisation (InitType) must be set!')
    }
    else {
        if (!object@InitType %in% c(0,1,2)){
            msg <- c(msg, 'Type of initialisation (InitType) must be 0, 1 or 2!')
        }
        else { # Valid InitType
            if (object@InitType == 0){  # Free Init
                if (is.na(object@FreeType) || length(object@FreeType) == 0){
                    msg <- c(msg, 'FreeType is required if InitType = 0 (Free initialisation).')
                }
                else{
                    if (!object@FreeType %in% c(0,1)){
                        msg <-  c(msg, 'FreeType must be 0 or 1.')
                    }
                    else { # Valid FreeType
                        if (object@FreeType == 0){ #random
                            if (is.na(object@NrCells) || length(object@NrCells) == 0){
                                msg <- c(msg, 'NrCells is required if FreeType = 0 (Random Free initialisation).')
                            }
                            else { # Valid NrCells
                                if (object@NrCells < 1){
                                    msg <- c(msg, 'NrCells must be positive.')
                                }
                            }
                        }
                    }
                }
                min.set = FALSE
                max.set = FALSE
                if (length(object@minX) != 0 && object@minX != -9){
                    if (object@minX < 0){
                        msg <- c(msg, 'Minimum X coordinate (minX) has to be greater or equal 0!')
                    }
                    else min.set = TRUE
                }
                if (length(object@maxX) != 0 && object@maxX != -9){
                    if (object@maxX < 0){
                        msg <- c(msg, 'Maximum X coordinate (maxX) has to be greater or equal 0!')
                    }
                    else max.set = TRUE
                }
                if(min.set && max.set){
                    Xextend = object@maxX - object@minX
                    if (object@maxX <= object@minX){
                        msg <- c(msg, 'maxX has to be larger than minX!')
                        Xextend = FALSE
                    }
                }
                else {Xextend = FALSE}
                min.set = FALSE
                max.set = FALSE
                if (length(object@minY) != 0 && object@minY != -9){
                    if (object@minY < 0){
                        msg <- c(msg, 'Minimum Y coordinate (minY) has to be greater or equal 0!')
                    }
                    else min.set = TRUE
                }
                if (length(object@maxY) != 0 && object@maxY != -9){
                    if (object@maxY < 0){
                        msg <- c(msg, 'Maximum Y coordinate (maxY) has to be greater or equal 0!')
                    }
                    else max.set = TRUE
                }
                if(min.set && max.set){
                    Yextend = object@maxY - object@minY
                    if (object@maxY <= object@minY){
                        msg <- c(msg, 'maxY has to be larger than minY!')
                        Yextend = FALSE
                    }
                }
                else {Yextend = FALSE}
                if(Xextend && Yextend && length(object@NrCells) > 0){
                    if (object@NrCells > Xextend*Yextend){
                        msg <- c(msg, 'The area bounded by (minX, maxX, minY, maxY) contains less cells than required by NrCells.')
                    }
                }
            }
            if (object@InitType == 1){  # Init Distmap
                if (is.na(object@SpType) || length(object@SpType) == 0){
                    msg <- c(msg, 'SpType is required if InitType = 1 (from loaded species distribution map).')
                }
                else{
                    if (!object@SpType %in% c(0,1)){
                        msg <-  c(msg, 'SpType must be 0 or 1.')
                    }
                    else { # Valid SpType
                        if (object@SpType == 1){ #random
                            if (is.na(object@NrCells) || length(object@NrCells) == 0){
                                msg <- c(msg, 'NrCells is required if SpType = 1 (Random initialisation from loaded species distribution map).')
                            }
                            else { # Valid NrCells
                                if (object@NrCells < 1){
                                    msg <- c(msg, 'NrCells must be positive.')
                                }
                            }
                        }
                    }
                }
            }
            if (object@InitType == 2){  # Init IndsList
                if (object@InitIndsFile == "NULL"){
                    msg <- c(msg, 'InitIndsFile is required if InitType = 2 (from loaded initial individuals list).')
                }
            }
        }
    }
    if (is.na(object@InitDens) || length(object@InitDens) == 0) {
        msg <- c(msg, 'InitDens must be set!')
    }
    else {
        if (!object@InitDens %in% c(0,1,2)) {
            msg <- c(msg, 'InitDens must be 0, 1 or 2!')
        }
        else {
            if (object@InitDens == 2){
                if(is.na(object@IndsHaCell) || length(object@IndsHaCell) == 0) {
                    msg <- c(msg, 'IndsHaCell is required if InitDens = 2')
                }
                else {
                    if (object@IndsHaCell <= 0){
                        msg <- c(msg, 'IndsHaCell must be positive.')
                    }
                }
            }
        }
    }
    if (object@InitType != 2) {
        if (length(object@InitAge) != 0){
            if (!object@InitAge %in% c(0,1,2)){
                msg <- c(msg, 'InitAge must be 0, 1 or 2!')
            }
        }
        if (length(object@PropStages) != 0){
            if(any(object@PropStages < 0.0) || any(object@PropStages > 1.0)) {
                msg <- c(msg, 'All elements of PropStages must be in the closed interval [0,1]!')
            }
            else{
                if (object@PropStages[1] != 0.0) {
                    msg <- c(msg, 'Initial proportion of the juvenile stage (PropStages[1]) must be 0.0!')
                }
                else{
                    if (length(object@PropStages)>1 && sum(object@PropStages) != 1.0) {
                        msg <- c(msg, 'The elements of PropStages must sum to 1!')
                    }
                }
            }
        }
    }
    if (is.na(object@InitFreezeYear) || length(object@InitFreezeYear) == 0) {
        msg <- c(msg, 'InitFreezeYear must be set!')
    }
    else {
        if(object@InitFreezeYear < 0) {
            msg <- c(msg, 'InitFreezeYear must be greater or equal 0!')
        }
    }
    if (is.na(object@RestrictRows) || length(object@RestrictRows) == 0) {
        msg <- c(msg, 'RestrictRows must be set!')
    }
    else {
        if(object@RestrictRows < 0) {
            msg <- c(msg, 'RestrictRows must be greater or equal 0!')
        }
    }
    if (is.na(object@RestrictFreq) || length(object@RestrictFreq) == 0) {
        msg <- c(msg, 'RestrictFreq must be set!')
    }
    else {
        if(object@RestrictFreq <= 0 && object@RestrictRows > 0) {
            msg <- c(msg, 'RestrictFreq must be strictly positive!')
        }
    }
    if (is.na(object@FinalFreezeYear) || length(object@FinalFreezeYear) == 0) {
        msg <- c(msg, 'FinalFreezeYear must be set!')
    }
    else {
        if(object@FinalFreezeYear < 0) {
            msg <- c(msg, 'FinalFreezeYear must be greater or equal 0!')
        }
        else {
            if(object@FinalFreezeYear > 0 && object@FinalFreezeYear <= object@InitFreezeYear) {
                msg <- c(msg, 'FinalFreezeYear must be greater than InitFreezeYear!')
            }
        }
    }
    if (is.null(msg)) TRUE else msg}
)

setMethod('initialize', 'InitialisationParams', function(.Object, ...) {
    this_func = "Initialise(): "
    args <- list(...)
    .Object <- callNextMethod()
    if (.Object@InitType != 0) {
        .Object@FreeType = -9L
        if (!is.null(args$FreeType)) {
            warning(this_func, "FreeType", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
    }
    if (.Object@InitType != 1) {
        .Object@SpType = -9L
        if (!is.null(args$SpType)) {
            warning(this_func, "SpType", warn_msg_ignored, "since InitType != 1.", call. = FALSE)
        }
    }
    if (.Object@InitType == 2) {
        .Object@InitDens <- 0
        .Object@PropStages <- -9L
        .Object@InitAge <- -9L
        if (!is.null(args$InitDens)) {
            warning(this_func, "InitDens", warn_msg_ignored, "since InitType = 2.", call. = FALSE)
        }
        if (!is.null(args$PropStages)) {
            warning(this_func, "PropStages", warn_msg_ignored, "since InitType = 2.", call. = FALSE)
        }
        if (!is.null(args$InitAge)) {
            warning(this_func, "InitAge", warn_msg_ignored, "since InitType = 2.", call. = FALSE)
        }
    }
    else{
        .Object@InitIndsFile = "NULL"
        if (!is.null(args$InitIndsFile)) {
            warning(this_func, "InitIndsFile", warn_msg_ignored, "since InitType != 2.", call. = FALSE)
        }
    }
    if (!((.Object@InitType == 0 && .Object@FreeType == 0) || (.Object@InitType == 1 && .Object@SpType == 1))) {
        .Object@NrCells = -9L
        if (!is.null(args$NrCells)) {
            warning(this_func, "NrCells", warn_msg_ignored, "since its not used.", call. = FALSE)
        }
    }
    if (.Object@InitDens != 2) {
        .Object@IndsHaCell = -9L
        if (!is.null(args$IndsHaCell)) {
            warning(this_func, "IndsHaCell", warn_msg_ignored, "since InitDens != 2.", call. = FALSE)
        }
    }
    if (.Object@InitType == 0) {
        # if (minX, maxX, minY, maxY) are negative, they will be set to their default values on C++ level
        if (is.null(args$minX)) {
            .Object@minX = -9L
        }
        if (is.null(args$minY)) {
            .Object@minY = -9L
        }
        if (is.null(args$maxX)) {
            .Object@maxX = -9L
        }
        if (is.null(args$maxY)) {
            .Object@maxY = -9L
        }
        if (!.Object@RestrictRows) {
            .Object@RestrictFreq = -9L
            if (!is.null(args$RestrictFreq)) {
                warning(this_func, "RestrictFreq", warn_msg_ignored, "since RestrictRows = 0.", call. = FALSE)
            }
        }
    }
    else {
        .Object@minX = -9L
        .Object@minY = -9L
        .Object@maxX = -9L
        .Object@maxY = -9L
        if (!is.null(args$minX)) {
            warning(this_func, "minX", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
        if (!is.null(args$minY)) {
            warning(this_func, "minY", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
        if (!is.null(args$maxX)) {
            warning(this_func, "maxX", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
        if (!is.null(args$maxY)) {
            warning(this_func, "maxY", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
        if (!is.null(args$InitFreezeYear)) {
            .Object@InitFreezeYear = 0L
            warning(this_func, "InitFreezeYear", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
        if (!is.null(args$RestrictRows)) {
            .Object@RestrictRows = 0L
            warning(this_func, "RestrictRows", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
        if (!is.null(args$RestrictFreq)) {
            .Object@RestrictFreq = 0L
            warning(this_func, "RestrictFreq", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
        if (!is.null(args$FinalFreezeYear)) {
            .Object@FinalFreezeYear = 0L
            warning(this_func, "FinalFreezeYear", warn_msg_ignored, "since InitType != 0.", call. = FALSE)
        }
    }
    return(.Object)
    }
)
setMethod("show", "InitialisationParams", function(object){
    cat(" Initialisation: \n")
    cat("   InitType =", object@InitType, ": ")
    if (object@InitType == 0) {
        cat("Free initialisation \n")
        if (object@FreeType == 0) {cat("   of", object@NrCells, "random suitable cells/patches.")}
        else{cat("                  of all suitable cells/patches.")}
    }
    if (object@InitType == 1) {
        cat("Initialisation from loaded species distribution map\n")
        if (object@FreeType == 0) {cat("   ", object@NrCells, "random presence cells/patches.")}
        else{cat("                  all presence cells/patches.")}
    }
    if (object@InitType == 2) {
        cat("Initialisation from initial individuals list\n                  from file:",object@InitIndsFile)
    }
    cat("\n")

    if (object@InitType != 2) {
        cat("   InitDens =", object@InitDens, ": ")
        if (object@InitDens == 0) {
            cat("At K_or_DensDep \n")
        }
        if (object@InitDens == 1) {
            cat("At half K_or_DensDep \n")
        }
        if (object@InitDens == 2) {
            cat(object@IndsHaCell," individuals per cell/hectare \n")
        }
        if (length(object@PropStages) > 2) {
            cat("   PropStages =", object@PropStages, "\n")
            cat("   InitAge =", object@InitAge, " : ")
            if (object@InitAge == 0) { cat("Minimum age \n") }
            if (object@InitAge == 1) { cat("Random between min and max \n") }
            if (object@InitAge == 2) { cat("Quasi-equilibrium \n") }
        }
    }
    if (object@InitType == 0) {
        if (object@minX >= 0 || object@minY >= 0 || object@maxX >= 0 || object@maxY >= 0) {cat("   ")}
        if (object@minX >= 0) { cat("MinX =",object@minX, "  ") }
        if (object@maxX >= 0) { cat("MaxX =",object@maxX, "  ") }
        if (object@minY >= 0) { cat("MinY =",object@minY, "  ") }
        if (object@maxY >= 0) { cat("MaxY =",object@maxY, "  ") }
        if (object@minX >= 0 || object@minY >= 0 || object@maxX >= 0 || object@maxY >= 0) {cat("\n")}
        if (object@InitFreezeYear > 0) { cat("InitFreezeYear =",object@InitFreezeYear,"\n") }
        if (object@RestrictRows > 0) { cat("RestrictRows =",object@RestrictRows,", RestrictFreq =",object@RestrictFreq,"\n") }
        if (object@FinalFreezeYear > 0) { cat("FinalFreezeYear =",object@FinalFreezeYear,"\n") }
    }
    cat("\n")}
)
