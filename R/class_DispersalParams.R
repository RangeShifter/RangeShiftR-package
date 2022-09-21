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



#### DISPERSAL ####


### SUBCLASS EMIGRATIONPARAMS

# from RS 'Emigration' file

#' Set Emigration Parameters
#'
#' Emigration - the first phase of dispersal - is modelled as the probability that an individual leaves its natal patch during the present year (or season).
#' It is constant by default, but can be set to be density-dependent (\code{DensDep}) and/or to vary for each individual (\code{IndVar}). In case of a stage-structured/sexual
#' population model, the emigration probabilities can also vary with stage/sex (\code{StageDep/SexDep}). If inter-individual variability is
#' enabled, the emigration traits can also be allowed to evolve (also set \code{TraitScaleFactor}).
#'
#' @usage Emigration(EmigProb = 0.0,
#'            SexDep = FALSE, StageDep = FALSE,
#'            DensDep = FALSE,
#'            IndVar = FALSE,
#'            TraitScaleFactor,
#'            EmigStage,
#'            UseFullKern = FALSE)
#' @param EmigProb Matrix containing all parameters (#columns) to determine emigration probabilities for each stage/sex (#rows). Its structure depends on the other parameters, see the Details.
#' If the emigration probability is constant (i.e. \code{DensDep, IndVar, StageDep, SexDep = FALSE}), \code{EmigProb} can take a single numeric. Defaults to \eqn{0.0}.
#' @param SexDep Sex-dependent emigration probability? (default: \code{FALSE})
#' @param StageDep Stage-dependent emigration probability? (default: \code{FALSE}) Must be \code{FALSE} if \code{IndVar=TRUE}.
#' @param DensDep Density-dependent emigration probability? (default: \code{FALSE})
#' @param IndVar Individual variability in emigration traits? (default: \code{FALSE}) Must be \code{FALSE} if \code{StageDep=TRUE}.
#' @param TraitScaleFactor Required if \code{IndVar=TRUE}: The scaling factor(s) for emigration traits. A numeric of length \eqn{1} (if \code{DensDep=FALSE}) or \eqn{3} (if \code{DensDep=TRUE}).
#' @param EmigStage Required for stage-structured populations with \code{IndVar=TRUE}: Stage which emigrates. (\code{StageDep} must be \code{FALSE})
#' @param UseFullKern Applicable only if transfer phase is modelled by a \code{\link[RangeShiftR]{DispersalKernel}} and \code{DensDep=FALSE}: Shall the emigration probability be derived from dispersal kernel?
#' @details Emigration is modelled as the probability \eqn{d} that an individual leaves its natal patch during the present year or season (every reproductive event
#' is always followed by dispersal).
#' Populations with non-overlapping generations have only a single opportunity to emigrate, whereas in a stage-structured population the realised overall emigration
#' rate can be larger than \eqn{d}, if a stage with \eqn{d>0} lasts for more than one year so that an individual has multiple opportunities to emigrate (each with probability \eqn{d}).
#'
#' The emigration probability \eqn{d} can be density-dependent (set \code{DensDep=TRUE}), in which case it is given by the following function, introduced by \insertCite{kun2006evolution;textual}{RangeShiftR}:
#'
#' \ifelse{html}{\out{&emsp;&emsp; d(i,t) = D<sub>0</sub> / ( 1 + e<sup>-&alpha;sub>E</sub> (N(i,t) / K(i,t) - &beta;sub>E</sub>) </sup> ) } }{\deqn{ d(i,t) = D_0 / ( 1 + exp[-α_E (N(i,t)/K(i,t) - β_E) ] ) } }
#'
#' In the case of stage-structured models this equation is modified to:
#'
#' \ifelse{html}{\out{&emsp;&emsp; d(i,t) = D<sub>0</sub> / ( 1 + e<sup>-&alpha;sub>E</sub> (b(i,t) * N(i,t) - &beta;sub>E</sub>) </sup> ) } }{\deqn{ d(i,t) = D_0 / ( 1 + exp[-α_E (b(i,t) N(i,t) - β_E) ] ) } }
#'
#' In the first case, \eqn{K(i,t)} is the carrying capacity of the cell/patch \eqn{i} at time \eqn{t} given by \code{K_or_DensDep}.
#' In the latter case, \eqn{b(i,t)} represents the strength of density dependence and is given by the inverse of \code{K_or_DensDep}.\cr
#' Further, \ifelse{html}{\out{D<sub>0</sub>}}{\eqn{D_0}} is the maximum emigration probability,
#' \eqn{N(i,t)} is the number of individuals in the cell/patch \eqn{i} at time \eqn{t},
#' \ifelse{html}{\out{&beta;<sub>E</sub>}}{\eqn{β_S}} is the inflection point of the function and
#' \ifelse{html}{\out{&alpha;<sub>E</sub>}}{\eqn{α_S}} is the slope at the inflection point.\cr
#'
#' Various functions have been proposed for density dependent emigration \insertCite{hovestadt2010information,poethke2011ability}{RangeShiftR}.
#' This one was chosen here because it is a flexible function that
#' allows for modelling a range of different reaction norms, as well as their emergence through evolution. In the case of density-dependent
#' emigration, we assume individuals to have full knowledge of the population density and carrying capacity (or 1/b, respectively) in their natal patch.
#' Information acquisition is not explicitly modelled.
#'
#' The emigration probability can be allowed to vary between individuals (set \code{IndVar=TRUE}) and to evolve. In the this case, individuals exhibit either one trait
#' determining the density-independent \eqn{d} (when \code{DensDep=FALSE}), or the three traits \ifelse{html}{\out{D<sub>0</sub>}}{\eqn{D_0}}, \eqn{α} and
#' \eqn{β} determining the density-dependent emigration probability (when \code{DensDep=TRUE}).\cr
#' For each trait the initial distribution in the population (as mean and standard variation) must be set in \code{EmigProb} (instead of only one constant value),
#' as well as their scaling factors in \code{TraitScaleFactor} (see \code{\link[RangeShiftR]{Genetics}}).
#' Also, if \code{IndVar=TRUE} is set for a stage-structured population, it is required to specify the stage which emigrates via \code{EmigStage}.
#'
#' It is possible to model sex-specific emigration strategies (set \code{SexDep=TRUE}) \insertCite{greenwood1980mating,lawson2007advances}{RangeShiftR}.
#' In this case the number of traits is doubled; one set coding for the trait(s) in females and the other for the trait(s) in males.
#' As well as being sex-biased, emigration parameters can be stage-biased (set \code{StageDep=TRUE}) when modelling stage-structured populations.
#' However, the current version does not accommodate inter-individual variation in emigration strategies when they are stage-dependent.
#'
#' The parameters that determine the emigration probabilities have to be provided via \code{EmigProb}, which generally takes a matrix, or - if only a single constant probability is
#' used (i.e. \code{DensDep, IndVar, StageDep, SexDep = FALSE}) - a single numeric. The format of the matrix is defined as follows: The number of columns depend on the options \code{DensDep} and \code{IndVar}. If \code{DensDep=FALSE}, the
#' density-independent probability \eqn{d} must be specified. If \code{DensDep=TRUE}, the functional parameters \ifelse{html}{\out{D<sub>0</sub>}}{\eqn{D_0}}, \eqn{α} and \eqn{β} (cp. equation above) must be specified.
#' Additionally, if \code{IndVar=FALSE}, these parameters are fixed, but if \code{IndVar=TRUE} each of them is replaced by two parameters: their respective mean and
#' standard deviation. They are used to normally distribute the traits values among the individuals of the initial population.
#'
#' All parameters have to be given for each stage/sex if the respective dependence is enabled. If \code{StageDep=TRUE}, state the corresponding stage in the first column.
#' If \code{SexDep=TRUE}, state the corresponding stage in the next (i.e. first/second) column, with \eqn{0} for \emph{female} and \eqn{1} for \emph{male}. The following table lists the required columns and their correct order for different settings:
#'
#' \tabular{ccccc}{DensDep \tab IndVar \tab StageDep \tab SexDep \tab columns \cr
#'  F \tab F \tab F \tab F \tab \eqn{d} \cr
#'  F \tab F \tab T \tab F \tab stage, \eqn{d} \cr
#'  F \tab F \tab F \tab T \tab sex, \eqn{d} \cr
#'  F \tab F \tab T \tab T \tab stage, sex, \eqn{d} \cr
#'  T \tab F \tab F \tab F \tab \ifelse{html}{\out{D<sub>0</sub>}}{\eqn{D_0}}, \eqn{α}, \eqn{β} \cr
#'  F \tab T \tab F \tab F \tab mean\eqn{(d)}, sd\eqn{(d)} \cr
#'  T \tab T \tab F \tab F \tab \ifelse{html}{\out{mean(D<sub>0</sub>)}}{mean\eqn{(D_0)}}, \ifelse{html}{\out{sd(D<sub>0</sub>)}}{sd\eqn{(D_0)}}, mean\eqn{(α)}, sd\eqn{(α)}, mean\eqn{(β)}, sd\eqn{(β)} \cr
#'  \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \cr
#'  T \tab T \tab F \tab T \tab sex, \ifelse{html}{\out{mean(D<sub>0</sub>)}}{mean\eqn{(D_0)}}, \ifelse{html}{\out{sd(D<sub>0</sub>)}}{sd\eqn{(D_0)}}, mean\eqn{(α)}, sd\eqn{(α)}, mean\eqn{(β)}, sd\eqn{(β)}
#'  }
#'
#' The column headings need not be included, only the numeric matrix is required. The rows require no particular order, but there must be exactly one row for each stage/sex combination. For example, in the case of density-, stage- and sex-dependent emigration with no individual variability:
#' \tabular{ccccc}{ \out{&emsp;} 0 \tab \out{&emsp;} 0 \tab \out{&emsp;} 1.0 \tab \out{&emsp;} 20 \tab \out{&emsp;} 0.2 \cr
#'  \out{&emsp;} 0 \tab \out{&emsp;} 1 \tab \out{&emsp;} 1.0 \tab \out{&emsp;} 20 \tab \out{&emsp;} 0.1 \cr
#'  \out{&emsp;} 1 \tab \out{&emsp;} 0 \tab \out{&emsp;} 0.7 \tab \out{&emsp;} 25 \tab \out{&emsp;} 0.5 \cr
#'  \out{&emsp;} 1 \tab \out{&emsp;} 1 \tab \out{&emsp;} 0.8 \tab \out{&emsp;} 50 \tab \out{&emsp;} 0.5 \cr
#'  \out{&emsp;} 2 \tab \out{&emsp;} 0 \tab \out{&emsp;} 0.4 \tab \out{&emsp;} 10 \tab \out{&emsp;} 1.0 \cr
#'  \out{&emsp;} 2 \tab \out{&emsp;} 1 \tab \out{&emsp;} 0.5 \tab \out{&emsp;} 20 \tab \out{&emsp;} 1.0
#' }
#'
#' In the special case that \code{DensDep=FALSE} and transfer is realised by \code{\link[RangeShiftR]{DispersalKernel}}, then the option \code{UseFullKern} may be switched on. It
#' will prevent re-sampling from the kernel if the distance sampled does not move the individual out of its natal cell/patch. Such individuals
#' are treated as philopatric recruits, and hence the kernel determines the probability of emigration. In this case, the emigration probability
#' for all stages/sexes which potentially disperse should be set to \eqn{1.0}.
#' @examples # stage- and sex-dependent constant emigration probabilities:
#' emigmat_1 <- matrix(c(0,0,1,0,1,1,1,0,.7,1,1,.8,2,0,.4,2,1,.5), byrow = TRUE, ncol = 3)
#' emig_1 <- Emigration(StageDep = TRUE, SexDep = TRUE, EmigProb = emigmat_1)
#' plotProbs(emig_1)
#'
#' # stage- and sex- and density-dependent emigration:
#' emigmat_2 <- matrix(c(0,0,1,20,.2,0,1,1,20,.1,1,0,.7,25,.5,1,1,.8,50,.5,2,0,.4,10,1,2,1,.5,20,1), byrow = TRUE, ncol = 5)
#' emig_2 <- Emigration(DensDep = TRUE, StageDep = TRUE, SexDep = TRUE, EmigProb = emigmat_2)
#' plotProbs(emig_2)
#'
#' # inter-individual variation with density- and sex-dependence :
#' emigmat_3 <- matrix(c(0,.7,.1,20,7,.2,.01,1,.9,.05,40,4,.5,.05), byrow = TRUE, ncol = 7)
#' emig_3 <- Emigration(DensDep = TRUE, IndVar = TRUE, SexDep = TRUE, EmigProb = emigmat_3, TraitScaleFactor = c(.1,7,.05))
#' plotProbs(emig_3)
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "EmigrationParams"
#' @author Anne-Kathleen Malchow
#' @name Emigration
#' @export Emigration
Emigration <- setClass("EmigrationParams", slots = c(DensDep = "logical",
                                                     IndVar = "logical",
                                                     StageDep = "logical",
                                                     SexDep = "logical",
                                                     EmigProb = "matrix_OR_numeric",
                                                     TraitScaleFactor = "numeric",
                                                     EmigStage = "integer_OR_numeric",
                                                     UseFullKern = "logical")
                       , prototype = list(DensDep = FALSE,
                                          IndVar = FALSE,
                                          StageDep = FALSE,
                                          SexDep = FALSE,
                                          EmigProb = matrix(data = 0.0, nrow = 1, ncol = 1),
                                          #TraitScaleFactor = NA_real_,
                                          #EmigStage = 0L,
                                          UseFullKern = FALSE)
)
        # discuss TraitScaleFactor in the Details


setValidity("EmigrationParams", function(object) {
    msg <- NULL
    if (anyNA(object@DensDep) || length(object@DensDep)!=1) {
        msg <- c(msg, "DensDep must be set and of length 1!")
    }
    if (anyNA(object@IndVar) || length(object@IndVar)!=1) {
        msg <- c(msg, "IndVar must be set and of length 1!")
    }
    else {
        if (anyNA(object@StageDep) || length(object@StageDep)!=1) {
            msg <- c(msg, "StageDep must be set and of length 1!")
        }
        else{
            if (object@IndVar && object@StageDep) {
                msg <- c(msg, "Inter-individual variability (IndVar=TRUE) in stage-dependent (StageDep=TRUE) emigration traits is not implemented!")
            }
        }
    }
    if (anyNA(object@SexDep) || length(object@SexDep)!=1) {
        msg <- c(msg, "SexDep must be set and of length 1!")
    }
    if (anyNA(object@EmigProb) || length(object@EmigProb)==0) {
        msg <- c(msg, "EmigProb must be set!")
    }
    else{
        if (!object@IndVar && !object@DensDep && !object@StageDep && !object@SexDep) {
            if (length(object@EmigProb)==1) {
                if (object@EmigProb < 0 || object@EmigProb > 1) {
                    msg <- c(msg, "EmigProb must be in the closed interval [0,1]!")
                }
            }
            else {
                msg <- c(msg, "EmigProb must have only one entry since it is set to be constant!")
            }
        }
        else {
            if (class(object@EmigProb)[1]!="matrix" && length(object@EmigProb)!=1) {
                msg <- c(msg, "EmigProb must be a matrix!")
            }
            else {
                if (object@IndVar && !object@DensDep && !object@StageDep && !object@SexDep && any(dim(object@EmigProb)!=c(1,2)) ) {
                    msg <- c(msg, "EmigProb must be a 1x2 matrix if DensDep,StageDep,SexDep = FALSE and IndVar = TRUE!")
                }
                if (!object@IndVar && object@DensDep && !object@StageDep && !object@SexDep && any(dim(object@EmigProb)!=c(1,3)) ) {
                    msg <- c(msg, "EmigProb must be a 1x3 matrix if IndVar,StageDep,SexDep = FALSE and DensDep = TRUE!")
                }
                if (!object@IndVar && !object@DensDep && object@StageDep && !object@SexDep && dim(object@EmigProb)[2]!=2 ) {
                    msg <- c(msg, "EmigProb must have 2 columns if IndVar,DensDep,SexDep = FALSE and StageDep = TRUE!")
                }
                if (!object@IndVar && !object@DensDep && !object@StageDep && object@SexDep && any(dim(object@EmigProb)!=c(2,2)) ) {
                    msg <- c(msg, "EmigProb must be a 2x2 matrix if IndVar,DensDep,StageDep = FALSE and SexDep = TRUE!")
                }
                if (!object@IndVar && !object@DensDep && object@StageDep && object@SexDep && dim(object@EmigProb)[2]!=3 ) {
                    msg <- c(msg, "EmigProb must have 3 columns if IndVar,DensDep = FALSE and StageDep,SexDep = TRUE!")
                }
                if (!object@IndVar && object@DensDep && object@StageDep && !object@SexDep && dim(object@EmigProb)[2]!=4 ) {
                    msg <- c(msg, "EmigProb must have 4 columns if IndVar,SexDep = FALSE and DensDep,StageDep = TRUE!")
                }
                if (!object@IndVar && object@DensDep && !object@StageDep && object@SexDep && any(dim(object@EmigProb)!=c(2,4)) ) {
                    msg <- c(msg, "EmigProb must be a 2x4 matrix if IndVar,StageDep = FALSE and DensDep,SexDep = TRUE!")
                }
                if (object@IndVar && !object@DensDep && !object@StageDep && object@SexDep && any(dim(object@EmigProb)!=c(2,3)) ) {
                    msg <- c(msg, "EmigProb must be a 2x3 matrix if DensDep,StageDep = FALSE and IndVar,SexDep = TRUE!")
                }
                if (object@IndVar && object@DensDep && !object@StageDep && !object@SexDep && any(dim(object@EmigProb)!=c(1,6)) ) {
                    msg <- c(msg, "EmigProb must be a 1x6 matrix if StageDep,SexDep = FALSE and IndVar,DensDep = TRUE!")
                }
                if (object@IndVar && object@DensDep && !object@StageDep && object@SexDep && any(dim(object@EmigProb)!=c(2,7)) ) {
                    msg <- c(msg, "EmigProb must be a 2x7 matrix if StageDep = FALSE and IndVar,DensDep,SexDep = TRUE!")
                }
                if (!object@IndVar && object@DensDep && object@StageDep && object@SexDep && dim(object@EmigProb)[2]!=5 ) {
                    msg <- c(msg, "EmigProb must have 5 columns if IndVar = FALSE and DensDep,StageDep,SexDep = TRUE!")
                }
            }
        }
    }
    if (object@IndVar) {
        if (anyNA(object@TraitScaleFactor) || (length(object@TraitScaleFactor)==0)) {
            msg <- c(msg, "TraitScaleFactor must be set!")
        }
        else{
            if (object@DensDep) {
                if (length(object@TraitScaleFactor)!=3) {
                    msg <- c(msg, "TraitScaleFactor must have length 3 if DensDep=TRUE!")
                }
                else {
                    if (object@TraitScaleFactor[1] <= 0.0 || object@TraitScaleFactor[1] > 1.0 ) {
                        msg <- c(msg, "TraitScaleFactor μ(D0) must be in the half-open interval (0,1] !")
                    }
                    if (any(object@TraitScaleFactor[2:3] <= 0.0 )) {
                        msg <- c(msg, "TraitScaleFactor μ(α) and μ(β) must be strictly positive !")
                    }
                }
            }
            else {
                if (length(object@TraitScaleFactor)!=1) {
                    msg <- c(msg, "TraitScaleFactor must have length 1 if DensDep=FALSE!")
                }
                else {
                    if (object@TraitScaleFactor <= 0 || object@TraitScaleFactor > 1 ) {
                        msg <- c(msg, "TraitScaleFactor μ(D0) must be in the half-open interval (0,1] !")
                    }
                }
            }
        }
    }
    if (anyNA(object@UseFullKern) || length(object@UseFullKern)!=1) {
        msg <- c(msg, "UseFullKern must be set and of length 1!")
    }
    else {
        if (object@DensDep && object@UseFullKern) {
                msg <- c(msg, "UseFullKern is not applicable if DensDep = TRUE!")
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "EmigrationParams", function(.Object, ...) {
    this_func = "Emigration(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (!is.null(args$EmigProb)) {
        if (class(args$EmigProb)[1] =="numeric" && length(args$EmigProb)==1) {
            .Object@EmigProb <- as.matrix(args$EmigProb)
        }
    }
    # else {
    #     warning(this_func, "Using default EmigProb = 0.0", call. = FALSE)
    # }
    if (!.Object@IndVar) {
        .Object@TraitScaleFactor = -9L
        if (!is.null(args$TraitScaleFactor)) {
            warning(this_func, "TraitScaleFactor", warn_msg_ignored, "since IndVar = FALSE.", call. = FALSE)
        }
        .Object@EmigStage = -9L
        if (!is.null(args$EmigStage)) {
            warning(this_func, "EmigStage", warn_msg_ignored, "since IndVar = FALSE.", call. = FALSE)
        }
    }
    if (.Object@UseFullKern) {
        if (any(.Object@EmigProb!=0.0 & .Object@EmigProb!=1.0)) {
            warning(this_func, "If UseFullKern = TRUE, it is recommended (in most cases) to set all elements of EmigProb to either 1.0 or 0.0, since the emigration probability will be evaluated using the dispersal kernel.", call. = FALSE)
        }
    }
    .Object}
)
setMethod("show", "EmigrationParams", function(object){
    if (object@DensDep) {
        cat("   DensDep =", object@DensDep, "\n")
    }
    if (object@IndVar) {
        cat("   IndVar =", object@IndVar, "\n")
    }
    if (object@StageDep) {
        cat("   StageDep =", object@StageDep, "\n")
    }
    if (object@SexDep) {
        cat("   SexDep =", object@SexDep, "\n")
    }
    cat("   Emigration probabilities:\n")
    print(object@EmigProb)
    if (object@IndVar) {
        cat("   TraitScaleFactor =", object@TraitScaleFactor, "\n")
        if (!anyNA(object@EmigStage) && length(object@EmigStage)!=0) {
            cat("   EmigStage =", object@EmigStage, "\n")
        }
    }
    if (!object@DensDep && object@UseFullKern) {
        cat("   UseFullKern =", object@UseFullKern, "\n")
    }
})
setMethod("plotProbs", "EmigrationParams", function(x, stage = NULL, sex = NULL, xmax = NULL, ymax = NULL){
    emig <- x@EmigProb
    # error messages
    if (!is.null(stage)){
        if (x@StageDep) {
            emig <- subset(emig, emig[,1] %in% stage)
        }
        else{ print("This emigration module has no stage-dependency.\n") }
    }
    if (!is.null(sex)){
        if (x@SexDep) {
            if (x@StageDep) emig <- subset(emig, emig[,2] %in% sex)
            else emig <- subset(emig, emig[,1] %in% sex)
        }
        else{ print("This emigration module has no sex-dependency.\n") }
    }
    # get column indices
    if (x@StageDep) {
        if (x@SexDep) {ind_D0 <- 3} else {ind_D0 <- 2}
    }else{
        if (x@SexDep) {ind_D0 <- 2} else {ind_D0 <- 1}
    }
    if (x@IndVar) {IV <- 2} else {IV <- 1}
    # New plot
    if (x@DensDep) {
        if (is.null(xmax)) {
            ind_max <- which.max(emig[,ind_D0+2*IV])
            xmax = min(2*emig[ind_max,ind_D0+2*IV], emig[ind_max,ind_D0+2*IV] + 4.0/emig[ind_max,ind_D0+IV])
        }
        xvals = seq(0, xmax, length.out = 100)
    }
    else {if(is.null(xmax)) xmax <- 1} # !DensDep

    if (is.null(ymax)) {ymax = 1}
    plot(NULL, type = "n", ylab = "Emigration probability", xlab = "relative population density (N/K or bN)", xlim = c(0,xmax), ylim = c(0,ymax))
    leg.txt <- c()
    # Go through lines of distances matrix and add curves to plot
    for(line in 1:nrow(emig)){
        if (x@IndVar) {
            if (x@DensDep) {
                res <- matrix(ncol = 8, nrow = length(xvals))
                res[,1] <- densdep(xvals, A0 = emig[line,ind_D0]-emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]-emig[line,ind_D0+3], beta = emig[line,ind_D0+4]-emig[line,ind_D0+5])
                res[,2] <- densdep(xvals, A0 = emig[line,ind_D0]-emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]-emig[line,ind_D0+3], beta = emig[line,ind_D0+4]+emig[line,ind_D0+5])
                res[,3] <- densdep(xvals, A0 = emig[line,ind_D0]-emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]+emig[line,ind_D0+3], beta = emig[line,ind_D0+4]-emig[line,ind_D0+5])
                res[,4] <- densdep(xvals, A0 = emig[line,ind_D0]-emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]+emig[line,ind_D0+3], beta = emig[line,ind_D0+4]+emig[line,ind_D0+5])
                res[,5] <- densdep(xvals, A0 = emig[line,ind_D0]+emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]-emig[line,ind_D0+3], beta = emig[line,ind_D0+4]-emig[line,ind_D0+5])
                res[,6] <- densdep(xvals, A0 = emig[line,ind_D0]+emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]-emig[line,ind_D0+3], beta = emig[line,ind_D0+4]+emig[line,ind_D0+5])
                res[,7] <- densdep(xvals, A0 = emig[line,ind_D0]+emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]+emig[line,ind_D0+3], beta = emig[line,ind_D0+4]-emig[line,ind_D0+5])
                res[,8] <- densdep(xvals, A0 = emig[line,ind_D0]+emig[line,ind_D0+1], alpha = emig[line,ind_D0+2]+emig[line,ind_D0+3], beta = emig[line,ind_D0+4]+emig[line,ind_D0+5])
                polygon(c(xvals,rev(xvals)), c(apply(res, 1, min), rev(apply(res, 1, max))), border=NA, col='grey80')
            }
            else {#constant
                polygon(c(0,xmax,xmax,0), c(rep(emig[line,ind_D0]-emig[line,ind_D0+1],2),rep(emig[line,ind_D0]+emig[line,ind_D0+1],2)), border=NA, col='grey80')
            }
        }
        if (x@DensDep) {
            lines(xvals, densdep(xvals, A0 = emig[line,ind_D0], alpha = emig[line,ind_D0+IV], beta = emig[line,ind_D0+2*IV]), type = "l", lty = 1, col = line)
        }
        else {#constant
            lines(x=c(0,xmax), y=rep(emig[line,ind_D0],2), type = "l", lty = 1, col = line)
        }
        if (x@StageDep) {
            if (x@SexDep) {leg.txt <- c(leg.txt, paste0("Stage ",emig[line,1], ifelse(emig[line,2]," male"," female")))} else {leg.txt <- c(leg.txt, paste0("Stage ",emig[line,1]))}
        }
        else {
            if (x@SexDep) {leg.txt <- c(leg.txt, ifelse(emig[line,1],"male","female"))}
        }
    }
    if (length(leg.txt)>0) {
        legend("topleft", leg.txt, col = 1:nrow(emig), lwd = 1.5)
    }
})


### SUBCLASS TRANSFERPARAMS

# from RS 'Transfer' file
# vitual class acting as superclass for: 'DispersalKernel', 'StochMove', 'CorrRW'

#' Set a Transfer method
#'
#' Transfer (or transience) is the second phase of dispersal. It consists of the movement of an individual departing from its natal patch towards
#' a potential new patch, ending with settlement or mortality. This movement can be modelled by one of three alternative processes:\cr
#' - Dispersal kernel: use \code{\link[RangeShiftR]{DispersalKernel}}\cr
#' - Stochastic movement simulator (SMS): use \code{\link[RangeShiftR]{SMS}}\cr
#' - Correlated random walk (CRW): use \code{\link[RangeShiftR]{CorrRW}}
#'
#' @details
#' The choice between the two main methods to model the transfer phase, i.e. phenomenological dispersal kernels or mechanistic movement processes
#' (SMS and CRW) depends on the information available for a given species and on the level of detail that is considered important to represent in
#' the model (which will depend on its aim and the scale).
#' \cr
#' Dispersal is often a costly process for an organism \insertCite{bonte2012costs}{RangeShiftR} and, in some cases, a dispersing individual may
#' suffer mortality. Obtaining a sensible
#' representation of dispersal requires that these mortality costs are described appropriately. The total dispersal mortality experienced will
#' be the sum of two main sources of mortality: First, as the result of individuals failing to reach suitable habitat, and second, as the result
#' of increased energetic, time or attritional costs that longer-distance dispersers will experience.
#' For more details, refer also to the respective Details sections of the different transfer methods.
#' In parameterising the model, it will be important to recognise this such that dispersal mortality is not double-accounted.
#' @references
#'         \insertAllCited{}
#' @author Anne-Kathleen Malchow
#' @name Transfer
TransferParams <- setClass("TransferParams")
setMethod("show", "TransferParams", function(object){
    if (class(object)[1] == "DispersalKernel") cat("   Dispersal Kernel\n")
    if (class(object)[1] == "StochMove") cat("   Stochastic Movement Simulator\n")
    if (class(object)[1] == "CorrRW") cat("   Correlated Random Walk\n")
    }
)



## Transfer-class DISPERSALKERNEL

#' Set up a Dispersal Kernel
#'
#' A method to describe \code{\link[RangeShiftR]{Transfer}}: Dispersal kernels are statistical distributions that are largely used to describe dispersal distances. The main assumption behind them
#' is that the principal determinant of the probability of an individual dispersing to a particular site is the distance from the starting location.\cr
#' As for the other dispersal phases, movement abilities and strategies are under multiple selective pressures and can evolve separately.
#' As a result, the realised dispersal kernels will themselves evolve.
#'
#' @usage DispersalKernel(Distances = 100, DoubleKernel = FALSE,
#'                 SexDep = FALSE, StageDep = FALSE,
#'                 IndVar = FALSE,
#'                 TraitScaleFactor,
#'                 DistMort = FALSE,
#'                 MortProb = 0.0, Slope, InflPoint)
#' @param Distances Matrix containing all dispersal kernel parameters (#columns) for each stage/sex (#rows). Its structure depends on the other parameters, see the Details.
#' If the mean dispersal distance is constant (i.e. \code{DensDep, IndVar, StageDep, SexDep = FALSE}), \code{Distances} can take a single numeric. Defaults to \eqn{100}.
#' @param DoubleKernel Use a mixed (i.e. double negative exponential) kernel? (default: \code{FALSE}) Set probability for using Kernel-1 in matrix \code{Distances}.
#' @param SexDep Sex-dependent dispersal kernel? (default: \code{FALSE})
#' @param StageDep Stage-dependent dispersal kernel? (default: \code{FALSE}) Must be \code{FALSE} if \code{IndVar=TRUE}.
#' @param IndVar Individual variability in dispersal kernel traits? (default: \code{FALSE}) Must be \code{FALSE}, if \code{StageDep=TRUE}.
#' @param TraitScaleFactor Required if \code{IndVar=TRUE}: The scaling factor(s) for dispersal kernel traits. A numeric of length \eqn{1} (if \code{DoubleKernel=FALSE}) or \eqn{3} (if \code{DoubleKernel=TRUE}).
#' @param DistMort Distance-dependent mortality probability? (default: \code{FALSE})
#' @param MortProb Constant mortality probability. Required if \code{DistMort=FALSE}, defaults to \eqn{0.0}.
#' @param InflPoint Required if \code{DistMort=TRUE}: Inflection point for the mortality distance dependence function.
#' @param Slope Required if \code{DistMort=TRUE}: Slope at inflection point for the mortality distance dependence function.
#' @details
#' Two types of dispersal kernels are implemented: negative exponential and a mixed kernel given by two different negative exponentials.
#' Here, kernels are considered as ‘distance kernels’, i.e. the statistical distribution of the probability that an individual will move a
#' certain distance \insertCite{hovestadt2012evolution,nathan2012dispersal}{RangeShiftR}. These kernels are specifically used for the \code{\link[RangeShiftR]{Transfer}} phase, meaning that they do
#' not incorporate information on the \code{\link[RangeShiftR]{Emigration}} or \code{\link[RangeShiftR]{Settlement}} probabilities, which are modelled independently. Therefore, dispersal kernels are
#' applied only to dispersing individuals and not normally to the entire population. However, the program allows a particular setting where
#' emigration and transfer are not explicitly separated but are both modelled through the kernel (see the parameter \code{UseFullKern} in
#' \code{\link[RangeShiftR]{Emigration}} and the Details there).
#'
#' There are many possible statistical distributions that have been fitted to dispersal data, which in many cases perform better
#' than the negative exponential \insertCite{nathan2012dispersal}{RangeShiftR}. However, the negative exponential is still commonly used, has been found useful for
#' describing dispersal patterns of certain organisms and the combination of two different negative exponentials has been demonstrated to be a
#' valuable method for discerning between common short-distance and rare long-distance dispersal \insertCite{hovestadt2011all}{RangeShiftR}.
#'
#' \emph{Negative exponential} \cr
#' If the individual disperses, the distance and the movement direction are determined in continuous space.
#' The distance is drawn from a negative exponential distribution with a given mean \eqn{δ}, and the direction is selected randomly from a uniform
#' distribution between \eqn{0} and \eqn{2π} radians.
#'
#' \ifelse{html}{\out{&emsp;&emsp; p(d;&delta;) = &delta;<sup>-1</sup> e<sup>- d / &delta;</sup>}}{\deqn{ p(d;δ) = 1/δ exp(-d/δ) } }
#'
#' If the arrival point lies beyond the boundary of the landscape, distance and direction are re-drawn.\cr
#' The individual is displaced from a random point (using continuous coordinates) inside the natal cell to the arrival cell where the model
#' switches back to discrete space \insertCite{bocedi2012projecting}{RangeShiftR}. If the arrival point is inside the natal cell, individual starting position, distance and direction are
#' re-sampled until the individual leaves the natal cell. In the case of patch-based models, the individual is assumed to disperse from a
#' random point in the patch and this position, the dispersal distance and direction are drawn until the individual leaves the patch. In order
#' to separate emigration and transfer explicitly, and to avoid potential infinite re-sampling, the program requires the mean of the kernel
#' to be greater or equal the cell resolution. This condition is relaxed only in the special case where emigration probability is set to be
#' density-independent and the kernel is applied to the entire population without re-sampling (the \code{UseFullKern} option in \code{\link[RangeShiftR]{Emigration}}). Individuals
#' which draw a short movement distance do not leave the natal cell/patch and implicitly become sedentary, and therefore the kernel itself
#' defines the proportion of individuals which emigrate. When this option is selected, the emigration probability for those stages/sexes
#' which disperse should be set to \eqn{1.0}; otherwise, only a proportion of such individuals would use the kernel to determine whether or
#' not they emigrate.
#'
#' \emph{Mixed kernel} \cr
#' The distance an individual moves is sampled from a mixed kernel given by the combination of two negative exponentials
#' with different means \ifelse{html}{\out{&delta;<sub>1</sub>}}{\eqn{δ_1}} and \ifelse{html}{\out{&delta;<sub>2</sub>}}{\eqn{δ_2}},
#' occurring with probability \ifelse{html}{\out{p<sub>I</sub>}}{\eqn{p_I}} and \eqn{1-}\ifelse{html}{\out{p<sub>I</sub>}}{\eqn{p_I}} respectively \insertCite{hovestadt2011all}{RangeShiftR}.
#' Otherwise, the conditions for the single kernel apply.
#'
#' \ifelse{html}{\out{&emsp;&emsp; p(d; &delta;<sub>1</sub>,&delta;<sub>2</sub>) = p<sub>I</sub> p(d;&delta;<sub>1</sub>) + (1-p<sub>I</sub>) p(d;&delta;<sub>1</sub>)}}{\deqn{ p(d; δ_1,δ_2) = p_I p(d;δ_1) + (1-p_I) p(d;δ_2)}}
#'
#' For both types of kernel, inter-individual variability of the kernel traits is possible (set \code{IndVar=TRUE}). Individuals will
#' carry either one trait for \eqn{δ} or three traits for \ifelse{html}{\out{&delta;<sub>1</sub>}}{\eqn{δ_1}}, \ifelse{html}{\out{&delta;<sub>2</sub>}}{\eqn{δ_2}} and
#' \ifelse{html}{\out{p<sub>I</sub>}}{\eqn{p_I}}, which they inherit from their parents.\cr
#' Dispersal kernels can also be sex-dependent (set \code{SexDep=TRUE}). In the case of inter-individual variability, the number of traits is doubled to two trait (female \eqn{δ}
#' and male δ) or six traits (female and male \ifelse{html}{\out{&delta;<sub>1</sub>}}{\eqn{δ_1}}, \ifelse{html}{\out{&delta;<sub>2</sub>}}{\eqn{δ_2}} and \ifelse{html}{\out{p<sub>I</sub>}}{\eqn{p_I}}).\cr
#' For each trait the initial distribution in the population (as mean and standard variation) must be set in \code{Distances} (instead of only one constant value),
#' as well as their scaling factors in \code{TraitScaleFactor} (see \code{\link[RangeShiftR]{Genetics}}).\cr
#'
#' Further, dispersal kernels can be stage-specific (set \code{StageDep=TRUE}). For this case, inter-individual variability is not implemented.
#'
#' All dispersal kernel parameters have to be provided via \code{Distances}, which generally takes a matrix, or - if only a single constant mean distance is
#' used (i.e. \code{DensDep, IndVar, StageDep, SexDep = FALSE}) - a single numeric. The format of the matrix is defined as follows: The number of columns depend on the options \code{IndVar} and \code{DoubleKernel}.
#' If \code{DoubleKernel=FALSE}, the mean dispersal distance \eqn{δ} must be specified (in meters). If \code{DoubleKernel=TRUE}, the mean dispersal distances
#' \ifelse{html}{\out{&delta;<sub>1</sub>}}{\eqn{δ_1}} and \ifelse{html}{\out{&delta;<sub>2</sub>}}{\eqn{δ_2}} (in meters), as well as the probability \ifelse{html}{\out{p<sub>I</sub>}}{\eqn{p_I}} of using Kernel-1 must be specified.
#' Additionally, if \code{IndVar=FALSE}, these parameters are fixed, but if \code{IndVar=TRUE} each of them is replaced by two parameters: their respective mean and
#' standard deviation. They are used to normally distribute the traits values among the individuals of the initial population.
#'
#' All parameters have to be given for each stage/sex if the respective dependence is enabled. If \code{StageDep=TRUE}, state the corresponding stage in the first column.
#' If \code{SexDep=TRUE}, state the corresponding stage in the next (i.e. first/second) column, with \eqn{0} for \emph{female} and \eqn{1} for \emph{male}. The following
#' table lists the required columns and their correct order for different settings:
#'
#' \tabular{ccccc}{IndVar \tab DoubleKernel \tab StageDep \tab SexDep \tab columns \cr
#'  F \tab F \tab F \tab F \tab \eqn{δ} \cr
#'  F \tab F \tab T \tab F \tab stage, \eqn{δ} \cr
#'  F \tab F \tab F \tab T \tab sex, \eqn{δ} \cr
#'  F \tab F \tab T \tab T \tab stage, sex, \eqn{δ} \cr
#'  F \tab T \tab F \tab F \tab \ifelse{html}{\out{&delta;<sub>1</sub>, &delta;<sub>2</sub>, p<sub>I</sub>}}{\eqn{δ_1, δ_2, p_I}} \cr
#'  T \tab F \tab F \tab F \tab mean\eqn{(δ)}, sd\eqn{(δ)} \cr
#'  T \tab T \tab F \tab F \tab \ifelse{html}{\out{mean(&delta;<sub>1</sub>)}}{mean\eqn{(δ_1)}}, \ifelse{html}{\out{sd(&delta;<sub>1</sub>)}}{sd\eqn{(δ_1)}}, \ifelse{html}{\out{mean(&delta;<sub>2</sub>)}}{mean\eqn{(δ_2)}}, \ifelse{html}{\out{sd(&delta;<sub>2</sub>)}}{sd\eqn{(δ_2)}}, mean\ifelse{html}{\out{(p<sub>I</sub>)}}{\eqn{(p_I)}}, sd\ifelse{html}{\out{(p<sub>I</sub>)}}{\eqn{(p_I)}} \cr
#'  \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \cr
#'  T \tab T \tab F \tab T \tab sex, \ifelse{html}{\out{mean(&delta;<sub>1</sub>)}}{mean\eqn{(δ_1)}}, \ifelse{html}{\out{sd(&delta;<sub>1</sub>)}}{sd\eqn{(δ_1)}}, \ifelse{html}{\out{mean(&delta;<sub>2</sub>)}}{mean\eqn{(δ_2)}}, \ifelse{html}{\out{sd(&delta;<sub>2</sub>)}}{sd\eqn{(δ_2)}}, mean\ifelse{html}{\out{(p<sub>I</sub>)}}{\eqn{(p_I)}}, sd\ifelse{html}{\out{(p<sub>I</sub>)}}{\eqn{(p_I)}}
#'  }
#'
#' The column headings need not be included, only the numeric matrix is required. The rows require no particular order, but there must be exactly
#' one row for each stage/sex combination. For example, in the case of a mixed kernel with stage- and sex-dependent distances and no individual variability:
#' \tabular{ccccc}{ \out{&emsp;} 0 \tab \out{&emsp;} 0 \tab \out{&emsp;} 1000 \tab \out{&emsp;} 4500 \tab \out{&emsp;} 0.92 \cr
#'  \out{&emsp;} 0 \tab \out{&emsp;} 1 \tab \out{&emsp;} 1400 \tab \out{&emsp;} 6000 \tab \out{&emsp;} 0.95 \cr
#'  \out{&emsp;} 1 \tab \out{&emsp;} 0 \tab \out{&emsp;} 700 \tab \out{&emsp;} 500 \tab \out{&emsp;} 0.50 \cr
#'  \out{&emsp;} 1 \tab \out{&emsp;} 1 \tab \out{&emsp;} 500 \tab \out{&emsp;} 600 \tab \out{&emsp;} 0.55 \cr
#'  \out{&emsp;} 2 \tab \out{&emsp;} 0 \tab \out{&emsp;} 100 \tab \out{&emsp;} 0 \tab \out{&emsp;} 1.0 \cr
#'  \out{&emsp;} 2 \tab \out{&emsp;} 1 \tab \out{&emsp;} 100 \tab \out{&emsp;} 0 \tab \out{&emsp;} 1.0
#' }
#'
#' In the case that the dispersal kernel is applied to the entire
#' population (i.e. density-independent emigration probability of \eqn{1.0}), the mean dispersal distance can evolve down to zero (i.e.
#' evolution for no dispersal). In all other cases where emigration and transfer are modelled separately, the mean dispersal distance has a
#' lower limit which can evolve equal to the landscape resolution.
#'
#' \emph{Mortality}\cr
#' There are two main sources of mortality: First, dispersal mortality can arise as a result of individuals failing to reach suitable habitat. When a simple dispersal
#' kernel is used with no possibility for individuals to search for locally-suitable habitat (see \code{\link[RangeShiftR]{Settlement}}),
#' mortality occurs to all individuals that arrive in unsuitable habitat. In this first case, dispersal mortality clearly depends upon the
#' proportion of suitable habitat in the landscape and will increase as the availability of habitat declines.\cr
#' A second source of dispersal mortality can be specified via the option \code{DistMort}: The probability of mortality is either a constant
#' (\eqn{m=}\code{MortProb}) or a function of distance \eqn{d} (i.e. individuals that travel further are more likely to die):
#'
#' \ifelse{html}{\out{&emsp;&emsp; m(d) = 1 / ( 1 + e<sup>-a (d- b)</sup> ) } }{\deqn{ m(d) = 1 / ( 1 + exp[-α (d-b) ] ) } }
#'
#' with the inflection point \eqn{b=}\code{InflPoint} at which \eqn{m(d=b)=0.5} and the slope \eqn{a=}\code{Slope}.This option may be thought
#' to represent the increased energetic, time or attritional costs that longer-distance dispersers will experience \insertCite{bonte2012costs}{RangeShiftR}.
#'
#' Note that the total dispersal mortality experienced will be the sum of the mortalities due to the two sources identified above and,
#' in parameterising the model, it will be important to recognize this such that dispersal mortality is not double-accounted.
#' @examples # stage- and sex-dependent mixed kernel
#' dists_1 <- matrix(c(0,0,1000,4500,0.92,0,1,1400,6000,0.95,1,0,700,500,0.50,1,1,500,600,0.55,2,0,100,0,1.0,2,1,100,0,1.0), byrow = TRUE, ncol = 5)
#' disp_1 <- DispersalKernel(Distances = dists_1, SexDep = TRUE, StageDep = TRUE, DoubleKernel = TRUE, DistMort = TRUE, Slope = 0.001, InflPoint = 4000)
#' plotProbs(disp_1, sex = 0, ymax = 0.002)
#' plotProbs(disp_1, mortality = TRUE)
#'
#' # mixed kernel with inter-individual variation
#' dists_2 <- matrix(c(0,1000,300,2500,500,0.62,0.13,1,3400,860,8000,2800,0.72,0.12), byrow = TRUE, ncol = 7)
#' disp_2 <- DispersalKernel(Distances = dists_2, SexDep = TRUE, DoubleKernel = TRUE, TraitScaleFactor = c(900,2800,0.14), IndVar = TRUE)
#' plotProbs(disp_2, xmax = 10000, combinekernels = TRUE)
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "DispersalKernel"
#' @author Anne-Kathleen Malchow
#' @name DispersalKernel
#' @export DispersalKernel
DispersalKernel <- setClass("DispersalKernel", slots = c(IndVar = "logical",
                                                         DoubleKernel = "logical",
                                                         StageDep = "logical",
                                                         SexDep = "logical",
                                                         Distances = "matrix_OR_numeric",
                                                         TraitScaleFactor = "numeric",
                                                         DistMort = "logical",
                                                         MortProb = "numeric",
                                                         Slope = "numeric",
                                                         InflPoint = "numeric")
                       , prototype = list(IndVar = FALSE,
                                          DoubleKernel = FALSE,
                                          StageDep = FALSE,
                                          SexDep = FALSE,
                                          Distances = matrix(data = 100L, nrow = 1, ncol = 1),
                                          #TraitScaleFactor,
                                          DistMort = FALSE,
                                          MortProb = 0.0
                                          #,Slope = -9,
                                          #InflPoint = -9
                                          )
                       , contains = "TransferParams"
)

setValidity("DispersalKernel", function(object) {
    msg <- NULL
    if (anyNA(object@DoubleKernel) || length(object@DoubleKernel)!=1) {
        msg <- c(msg, "DoubleKernel must be set and of length 1!")
    }
    if (anyNA(object@IndVar) || length(object@IndVar)!=1) {
        msg <- c(msg, "IndVar must be set and of length 1!")
    }
    else {
        if (anyNA(object@StageDep) || length(object@StageDep)!=1) {
            msg <- c(msg, "StageDep must be set and of length 1!")
        }
        else{
            if (object@IndVar && object@StageDep) {
                msg <- c(msg, "Inter-individual variability (IndVar=TRUE) in stage-dependent (StageDep=TRUE) dispersal kernel traits is not implemented!")
            }
        }
    }
    if (anyNA(object@SexDep) || length(object@SexDep)!=1) {
        msg <- c(msg, "SexDep must be set and of length 1!")
    }
    if (anyNA(object@Distances) || length(object@Distances)==0) {
        msg <- c(msg, "Distances must be set!")
    }
    else{
        if (!object@IndVar && !object@DoubleKernel && !object@StageDep && !object@SexDep) {
            if (length(object@Distances)!=1) {
                msg <- c(msg, "Distances must be a single value if IndVar,DoubleKernel,StageDep,SexDep = FALSE!")
            }
        }
        else {
            if (class(object@Distances)[1] !="matrix") {
                msg <- c(msg, "Distances must be a matrix!")
            }
        }
    }
    if (is.null(msg)) {
        if (object@IndVar && !object@DoubleKernel && !object@StageDep && !object@SexDep && any(dim(object@Distances)!=c(1,2)) ) {
            msg <- c(msg, "Distances must be a 1x2 matrix if DoubleKernel,StageDep,SexDep = FALSE and IndVar = TRUE!")
        }
        if (!object@IndVar && object@DoubleKernel && !object@StageDep && !object@SexDep && any(dim(object@Distances)!=c(1,3)) ) {
            msg <- c(msg, "Distances must be a 1x3 matrix if IndVar,StageDep,SexDep = FALSE and DoubleKernel = TRUE!")
        }
        if (!object@IndVar && !object@DoubleKernel && object@StageDep && !object@SexDep && dim(object@Distances)[2]!=2 ) {
            msg <- c(msg, "Distances must have 2 columns if IndVar,DoubleKernel,SexDep = FALSE and StageDep = TRUE!")
        }
        if (!object@IndVar && !object@DoubleKernel && !object@StageDep && object@SexDep && any(dim(object@Distances)!=c(2,2)) ) {
            msg <- c(msg, "Distances must be a 2x2 matrix if IndVar,DoubleKernel,StageDep = FALSE and SexDep = TRUE!")
        }
        if (!object@IndVar && !object@DoubleKernel && object@StageDep && object@SexDep && dim(object@Distances)[2]!=3 ) {
            msg <- c(msg, "Distances must have 3 columns if IndVar,DoubleKernel = FALSE and StageDep,SexDep = TRUE!")
        }
        if (!object@IndVar && object@DoubleKernel && object@StageDep && !object@SexDep && dim(object@Distances)[2]!=4 ) {
            msg <- c(msg, "Distances must have 4 columns if IndVar,SexDep = FALSE and DoubleKernel,StageDep = TRUE!")
        }
        if (!object@IndVar && object@DoubleKernel && !object@StageDep && object@SexDep && any(dim(object@Distances)!=c(2,4)) ) {
            msg <- c(msg, "Distances must be a 2x4 matrix if IndVar,StageDep = FALSE and DoubleKernel,SexDep = TRUE!")
        }
        if (object@IndVar && !object@DoubleKernel && !object@StageDep && object@SexDep && any(dim(object@Distances)!=c(2,3)) ) {
            msg <- c(msg, "Distances must be a 2x3 matrix if DoubleKernel,StageDep = FALSE and IndVar,SexDep = TRUE!")
        }
        if (object@IndVar && object@DoubleKernel && !object@StageDep && !object@SexDep && any(dim(object@Distances)!=c(1,6)) ) {
            msg <- c(msg, "Distances must be a 1x6 matrix if StageDep,SexDep = FALSE and IndVar,DoubleKernel = TRUE!")
        }
        if (object@IndVar && object@DoubleKernel && !object@StageDep && object@SexDep && any(dim(object@Distances)!=c(2,7)) ) {
            msg <- c(msg, "Distances must be a 2x7 matrix if StageDep = FALSE and IndVar,DoubleKernel,SexDep = TRUE!")
        }
        if (!object@IndVar && object@DoubleKernel && object@StageDep && object@SexDep && dim(object@Distances)[2]!=5 ) {
            msg <- c(msg, "Distances must have 5 columns if IndVar = FALSE and DoubleKernel,StageDep,SexDep = TRUE!")
        }
    }
    if (object@IndVar) {
        if (anyNA(object@TraitScaleFactor) || (length(object@TraitScaleFactor)==0)) {
            msg <- c(msg, "TraitScaleFactor must be set!")
        }
        else{
            if (object@DoubleKernel) {
                if (length(object@TraitScaleFactor)!=3) {
                    msg <- c(msg, "TraitScaleFactor must have length 3 if DoubleKernel=TRUE!")
                }
                else {
                    if (object@TraitScaleFactor[3] <= 0.0 || object@TraitScaleFactor[3] > 1.0 ) {
                        msg <- c(msg, "TraitScaleFactor μ(p) must be in the half-open interval (0,1] !")
                    }
                    if (any(object@TraitScaleFactor[1:2] <= 0.0 )) {
                        msg <- c(msg, "TraitScaleFactor μ(δ1) and μ(δ2) must be strictly positive !")
                    }
                }
            }
            else {
                if (length(object@TraitScaleFactor)!=1) {
                    msg <- c(msg, "TraitScaleFactor must have length 1 if DoubleKernel=FALSE!")
                }
                else {
                    if (object@TraitScaleFactor <= 0.0) {
                        msg <- c(msg, "TraitScaleFactor μ(δ) must be strictly positive !")
                    }
                }
            }
        }
    }
    if (anyNA(object@DistMort) || length(object@DistMort)!=1) {
        msg <- c(msg, "DistMort must be set and of length 1!")
    }
    else {
        if (object@DistMort) {
            # atm no conditions for Slope, InflPoint
            if (anyNA(object@Slope) || length(object@Slope)!=1) {
                msg <- c(msg, "Slope must be set and of length 1!")
            }
            if (anyNA(object@InflPoint) || length(object@InflPoint)!=1) {
                msg <- c(msg, "InflPoint must be set and of length 1!")
            }
        }
        else {
            if (anyNA(object@MortProb) || length(object@MortProb)!=1) {
                msg <- c(msg, "MortProb must be set and of length 1!")
            }
            else {
                if (object@MortProb < 0 | object@MortProb >= 1 ) {
                    msg <- c(msg, "MortProb must be in the half-open interval [0,1) !")
                }
            }
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "DispersalKernel", function(.Object, ...) {
    this_func = "DispersalKernel(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (class(args$Distances)[1]=="numeric" && length(args$Distances)==1) {
        .Object@Distances <- as.matrix(args$Distances)
    }
    if (!.Object@IndVar) {
        .Object@TraitScaleFactor = -9L
        if (!is.null(args$TraitScaleFactor)) {
            warning(this_func, "TraitScaleFactor", warn_msg_ignored, "since IndVar = FALSE.", call. = FALSE)
        }
    }
    if (.Object@DistMort) {
        .Object@MortProb = -9L
        if (!is.null(args$MortProb)) {
            warning(this_func, "MortProb", warn_msg_ignored, "since DistMort = TRUE.", call. = FALSE)
        }
    }
    else {
        .Object@Slope = -9L
        if (!is.null(args$Slope)) {
            warning(this_func, "Slope", warn_msg_ignored, "since DistMort = FALSE.", call. = FALSE)
        }
        .Object@InflPoint = -9L
        if (!is.null(args$InflPoint)) {
            warning(this_func, "InflPoint", warn_msg_ignored, "since DistMort = FALSE.", call. = FALSE)
        }
    }
    .Object}
)
setMethod("show", "DispersalKernel", function(object){
    callNextMethod()
    if (object@IndVar) {
        cat("   IndVar =", object@IndVar, "\n")
    }
    if (object@DoubleKernel) {
        cat("   DoubleKernel =", object@DoubleKernel, "\n")
    }
    if (object@StageDep) {
        cat("   StageDep =", object@StageDep, "\n")
    }
    if (object@SexDep) {
        cat("   SexDep =", object@SexDep, "\n")
    }
    cat("   Dispersal kernel traits:\n")
    print(object@Distances)
    if (object@IndVar) {
        cat("   TraitScaleFactor =", object@TraitScaleFactor, "\n")
    }
    if (object@DistMort) {
        cat("   Distance-dependent mortality prob with:\n   Inflection point =", object@InflPoint, "\n   Slope =", object@Slope, "\n")
    }
    else {
        cat("   Constant mortality prob =", object@MortProb, "\n")
    }
})
setMethod("plotProbs", "DispersalKernel", function(x, mortality = FALSE, combinekernels = FALSE, stage = NULL, sex = NULL, xmax = NULL, ymax = NULL){
    if (mortality) {
        # New plot
        if (x@DistMort) {
            if (is.null(xmax)) {xmax = min(2*x@InflPoint, x@InflPoint + 3.0/x@Slope)}
            xvals = seq(0, xmax, length.out = 100)
            plot(xvals, densdep(xvals, alpha = x@Slope, beta = x@InflPoint), type = "l", lwd = 1.5, ylab = "Mortality Probability", xlab = "Dispersal distance", ylim = c(0,1))
        }
        else {
            plot(x=c(0,10000), y=rep(x@MortProb,2), type = "l", lwd = 1.5, ylab = "Mortality Probability", xlab = "Dispersal distance", ylim = c(0,1))
        }
    }
    else { # !mortality
        dists <- x@Distances
        # error messages
        if (!is.null(stage)){
            if (x@StageDep) {
                dists <- subset(dists, dists[,1] %in% stage)
            }
            else{ print("This dispersal kernel has no stage-dependency.\n") }
        }
        if (!is.null(sex)){
            if (x@SexDep) {
                if (x@StageDep) dists <- subset(dists, dists[,2] %in% sex)
                else dists <- subset(dists, dists[,1] %in% sex)
            }
            else{ print("This dispersal kernel has no sex-dependency.\n") }
        }
        if (combinekernels){
            if (!x@DoubleKernel) {
                print("Kernels can only be combined if a mixed kernel (DoubleKernel = TRUE) is given.\n")
                combinekernels = FALSE
            }
        }
        # get column indices
        if (x@StageDep) {
            if (x@SexDep) {ind_kernel1 <- 3} else {ind_kernel1 <- 2}
        }
        else{
            if (x@SexDep) {ind_kernel1 <- 2} else {ind_kernel1 <- 1}
        }
        if (x@DoubleKernel) {
            if (x@IndVar) {ind_kernel2 <- ind_kernel1 + 2} else {ind_kernel2 <- ind_kernel1 + 1}
        }
        if (combinekernels){
            if (x@IndVar) {ind_pI <- ind_kernel2 + 2} else {ind_pI <- ind_kernel2 + 1}
        }
        # New plot
        xlab = "Dispersal distance"
        if (x@DoubleKernel && !combinekernels) xlab <- paste(xlab, "(Kernel-1 solid, Kernel-2 dashed)")
        if (x@IndVar) xlab <- paste(xlab, ", initial dist.")
        if (is.null(xmax)) {xmax = 3*max(dists)}
        xvals = seq(0, xmax, length.out = 100)
        if (is.null(ymax)) {
            if (x@DoubleKernel){
                ymax = 1/(min(dists[,c(ind_kernel1,ind_kernel2)][dists[,c(ind_kernel1,ind_kernel2)]>0]))
            }
            else{
                ymax = 1/(min(dists[,ind_kernel1][dists[,ind_kernel1]>0]))
            }
        }
        plot(NULL, type = "n", ylab = "Probability Density", xlab = xlab, xlim = c(0,xmax), ylim = c(0,ymax))
        leg.txt <- c()
        # Go through lines of distances matrix and add curves to plot
        for(line in 1:nrow(dists)){
            if(dists[line,ind_kernel1]>0){
                if (!combinekernels){
                    if (x@IndVar) {
                        res <- matrix(ncol = 3, nrow = length(xvals))
                        res[,1] <- dexp(xvals, rate = 1/(dists[line,ind_kernel1]))
                        res[,2] <- dexp(xvals, rate = 1/(dists[line,ind_kernel1]+dists[line,ind_kernel1+1]))
                        res[,3] <- dexp(xvals, rate = 1/(dists[line,ind_kernel1]-dists[line,ind_kernel1+1]))
                        polygon(c(xvals,rev(xvals)), c(apply(res, 1, min), rev(apply(res, 1, max))), border=NA, col='grey80')
                        }
                    lines(xvals,dexp(xvals,rate = 1/dists[line,ind_kernel1]), type = "l", lty = 1, col = line)
                }
            }
            else {lines(xvals, rep(0, length(xvals)), type = "l", lty = 1, col = line)}
            if (x@DoubleKernel){
                if(dists[line,ind_kernel2]>0){
                    if (combinekernels){
                        pI <- dists[line,ind_pI]
                        if (x@IndVar) {
                            pI_sd <- c(pI+dists[line,ind_pI+1],pI-dists[line,ind_pI+1])
                            res <- matrix(ncol = 8, nrow = length(xvals))
                            res[,1] <- pI_sd[1] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]+dists[line,ind_kernel1+1])) + (1-pI_sd[1]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]+dists[line,ind_kernel2+1]))
                            res[,2] <- pI_sd[2] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]+dists[line,ind_kernel1+1])) + (1-pI_sd[2]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]+dists[line,ind_kernel2+1]))
                            res[,3] <- pI_sd[1] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]-dists[line,ind_kernel1+1])) + (1-pI_sd[1]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]+dists[line,ind_kernel2+1]))
                            res[,4] <- pI_sd[2] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]-dists[line,ind_kernel1+1])) + (1-pI_sd[2]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]+dists[line,ind_kernel2+1]))
                            res[,5] <- pI_sd[1] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]+dists[line,ind_kernel1+1])) + (1-pI_sd[1]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]-dists[line,ind_kernel2+1]))
                            res[,6] <- pI_sd[2] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]+dists[line,ind_kernel1+1])) + (1-pI_sd[2]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]-dists[line,ind_kernel2+1]))
                            res[,7] <- pI_sd[1] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]-dists[line,ind_kernel1+1])) + (1-pI_sd[1]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]-dists[line,ind_kernel2+1]))
                            res[,8] <- pI_sd[2] * dexp(xvals, rate = 1/(dists[line,ind_kernel1]-dists[line,ind_kernel1+1])) + (1-pI_sd[2]) * dexp(xvals, rate = 1/(dists[line,ind_kernel2]-dists[line,ind_kernel2+1]))
                            polygon(c(xvals,rev(xvals)), c(apply(res, 1, min), rev(apply(res, 1, max))), border=NA, col='grey80')
                        }
                        yvals <- pI * dexp(xvals, rate = 1/dists[line,ind_kernel1]) + (1-pI) * dexp(xvals, rate = 1/dists[line,ind_kernel2])
                        lines(xvals, yvals , type = "l", lty = 1, col = line)
                    }
                    else {
                        if (x@IndVar) {
                            res[,1] <- dexp(xvals, rate = 1/(dists[line,ind_kernel2]))
                            res[,2] <- dexp(xvals, rate = 1/(dists[line,ind_kernel2]+dists[line,ind_kernel2+1]))
                            res[,3] <- dexp(xvals, rate = 1/(dists[line,ind_kernel2]-dists[line,ind_kernel2+1]))
                            polygon(c(xvals,rev(xvals)), c(apply(res, 1, min), rev(apply(res, 1, max))), border=NA, col='grey80')
                        }
                        lines(xvals,dexp(xvals,rate = 1/dists[line,ind_kernel2]), type = "l", lty = 2, col = line)
                    }
                }
                else{
                    lines(xvals, rep(0, length(xvals)), type = "l", lty = 2, col = line)
                }
            }
            if (x@StageDep) {
                if (x@SexDep) {leg.txt <- c(leg.txt, paste0("Stage ",dists[line,1], ifelse(dists[line,2]," male"," female")))} else {leg.txt <- c(leg.txt, paste0("Stage ",dists[line,1]))}
            }
            else {
                if (x@SexDep) {leg.txt <- c(leg.txt, ifelse(dists[line,1],"male","female"))}
            }
        }
        if (length(leg.txt)>0) {
            legend("topright", leg.txt, col = 1:nrow(dists), lwd = 1.5)
        }
    }
})



## Transfer-class STOCHMOVE

#' Set up a Stochastic Movement Simulator
#'
#' A method to describe \code{\link[RangeShiftR]{Transfer}}:
#' SMS is a stochastic individual-based movement model where organisms move through grid-based, heterogeneous landscapes. The model uses similar
#' cost surfaces as the least cost path (LCP) method, but it relaxes two of its main assumptions: Firstly, individuals are not assumed to be
#' omniscient, but move according to what they can perceive of the landscape within their perceptual range (\code{PR}). Secondly, individuals
#' do not know a priori their final destination, which is a reasonable assumption for dispersing individuals. For a complete description of the
#' method, see the Details below or refer to \insertCite{palmer2011introducing;textual}{RangeShiftR}.
#'
#' @usage SMS(PR = 1, PRMethod = 1, MemSize = 1,
#'     DP = 1.0,
#'     GoalType = 0,
#'     GoalBias = 1.0, AlphaDB, BetaDB,
#'     IndVar = FALSE,
#'     Costs, StepMort = 0.0,
#'     StraightenPath = TRUE)
#' @param PR Perceptual range. Given in number of cells, defaults to \eqn{1}. (integer)
#' @param PRMethod Method to evaluate the effective cost of a particular step from the landscape within the perceptual range:\cr \eqn{1 = }Arithmetic mean (default)\cr \eqn{2 = }Harmonic
#' mean\cr \eqn{3 = }Weighted arithmetic mean
#' @param MemSize Size of memory, given as the number of previous steps over which to calculate current direction to apply directional persistence
#' (\code{DP}). A maximum of \eqn{14} steps is supported, default is \eqn{1}. (integer)
#' @param DP Directional persistence. Corresponds to the tendency to follow a correlated random walk, must be \eqn{\ge 1.0}, defaults to \eqn{1.0}.\cr
#' If \code{IndVar=TRUE}, expects a vector of length three specifying (Mean, SD, TraitScaleFactor) of \code{DP}.
#' @param GoalType Goal bias type: \eqn{0 = } None (default), \eqn{2 = } Dispersal bias.
#' @param GoalBias Only if \code{GoalType=2}: Goal bias strength. Must be must be \eqn{\ge 1.0}, defaults to \eqn{1.0}. \cr If \code{IndVar=TRUE}, expects a vector of length three
#' specifying (Mean, SD, TraitScaleFactor) of \code{GoalBias}.
#' @param AlphaDB Required if \code{GoalType=2}: Dispersal bias decay rate. Must be must be \eqn{> 0.0}.\cr If \code{IndVar=TRUE}, expects a vector of length three
#' specifying (Mean, SD, TraitScaleFactor) of \code{AlphaDB}.
#' @param BetaDB Required if \code{GoalType=2}: Dispersal bias decay inflection point (given in number of steps). Must be must be \eqn{> 0.0}.\cr If \code{IndVar=TRUE},
#' expects a vector of length three specifying (Mean, SD, TraitScaleFactor) of \code{BetaDB}.
#' @param IndVar Individual variability in SMS traits (i.e. \code{DP}, \code{GoalBias}, \code{AlphaDB} and \code{BetaDB})? Defaults to \code{FALSE}.
#' @param Costs Describes the landscapes resistance to movement. Set either: \cr
#'  - \emph{habitat-specific} costs for each habitat type, or\cr
#'  - \code{"file"}, to indictae to use the \emph{cost raster} map(s) specified in the landscape module and import the cost values from them.\cr
#' In the first case of \emph{habitat-specific} costs a numeric vector is expected with, respectively, length \code{Nhabitats} for an \code{\link[RangeShiftR]{ImportedLandscape}}
#' with habitat codes (i.e. \code{HabPercent=FALSE})) or length \eqn{2} for an \code{\link[RangeShiftR]{ArtificialLandscape}} (matrix and habitat costs).\cr
#' In the second case of importing a \emph{cost raster} file, specify the file name(s) in the \code{\link[RangeShiftR]{ImportedLandscape}} module.
#' The specified map has to match the landscape raster in extent, coordinates and resolution, and each cell contains a cost value (\eqn{\ge 1}).
#' This is the only option for an imported landscape with habitat qualities / cover percentage (i.e. \code{HabPercent=TRUE}).
#' @param StepMort Per-step mortality probability. Can be either \emph{constant}, in which case a single numeric is expected (the default, with
#' value \eqn{0.0}) or \emph{habitat-specific}, in which case a numeric vector (with a length as described above for \code{Costs}) is expected.
#' All values must be within the half-open interval \eqn{[0,1)}. \cr
#' For an imported habitat quality landscape (\code{HabPercent=TRUE}), only constant per-step mortality is allowed.
#' @param StraightenPath Straighten path after decision not to settle in a patch? Defaults to \code{TRUE}, see Details below.
#' @details
#' SMS is a stochastic individual-based model where organisms move through grid-based, heterogeneous landscapes. The model uses similar cost
#' surfaces as the least cost path (LCP) \insertCite{adriaensen2003application,chardon2003incorporating,stevens2006gene,driezen2007evaluating}{RangeShiftR},
#' but it relaxes two of the main assumptions/limitations of the latter. Firstly, individuals are not assumed to be omniscient, but move according to what they
#' can perceive of the landscape within their perceptual range. Secondly, individuals do not know a priori their final destination, which is a
#' reasonable assumption for dispersing individuals. Here, the core components of SMS are briefly described;
#' see \insertCite{palmer2011introducing;textual}{RangeShiftR} for a complete description of the method.
#'
#' SMS uses cost maps where a relative cost to movement is assigned to each habitat type. Costs are integer numbers and represent the cost of
#' moving through a particular land cover relative to the cost of moving through breeding habitat (which is conventionally set to a cost of \eqn{1}).
#' Individuals take single cell steps basing their decisions on three parameters: their perceptual range (\code{PR}) (given in number of cells),
#' the method used to evaluate the landscape within their perceptual range (\code{PRMethod}), and their directional persistence (\code{DP}), which corresponds to
#' their tendency to follow a correlated random walk. At each step, the individual evaluates the surrounding habitat in order to determine
#' the effective cost of taking a particular step to each of the eight neighbouring cells. The effective cost is a mean of the cost of the
#' neighbouring cell and the surrounding cells beyond it within the \code{PR}, and is calculated by one of three possible methods:
#' - \emph{Arithmetic mean}\cr
#' - \emph{Harmonic mean} - The reciprocal of the arithmetic mean of the reciprocals of the observations (cell costs). This method increases the
#' detectability of low cost cells but performs less well than the arithmetic mean in detecting high cost cells. Therefore, the choice between
#' the two depends on whether the main driver of the animal movement is selecting for good habitat or avoiding costly habitat.\cr
#' - \emph{Weighted arithmetic mean} - The cost of each cell is weighted by its inverse distance from the individual (which is assumed to be in
#' the centre of the current cell).
#'
#' The effective cost of each neighbouring cell is weighted by the \code{DP}, which is lowest in the direction of travel. \code{DP} can be
#' calculated over more steps than just the previous one (up to a maximum of \eqn{14}), which is controlled by the memory size parameter (\code{MemSize})
#' \insertCite{palmer2014inter,aben2014simple}{RangeShiftR}. Increasing the memory size means that an individual retains for longer its tendency to move in a
#' certain direction, and hence paths tend to become somewhat smoother.
#'
#' There is an option to include goal bias, i.e. a tendency to move towards a particular destination  \insertCite{aben2014simple}{RangeShiftR},
#' which is implemented in a similar way to \code{DP}. However, as dispersers in \emph{RangeShiftR} are naïve and have no goal, it may be applied only in the ‘negative’
#' sense of moving away from the natal location (\code{GoalType=2}), i.e. as a dispersal bias, which is subject to a decay in strength as a
#' function of the total number of steps taken (set the decay rate \code{AlphaDB} and inflection point \code{BetaDB}).
#' This enables a dispersal path to follow a straighter trajectory initially, and later become more responsive to perceived landscape costs.
#'
#' The product of the reciprocals of effective cost, \code{DP} and dispersal bias, scaled to sum to one,
#' give the probabilities that the individual will move to each neighbouring cell. All the dispersing individuals move simultaneously, i.e. at
#' each time-step they all make one move. In the case of patch-based models, the individual is forced to leave the natal patch by increasing its
#' \code{DP} ten-fold until it has taken a number of steps (\eqn{=2*}\code{PR}) outside the natal patch.
#'
#' When an individual arrives in a non-natal patch and decides not to settle there (as a result of a density-dependent or mate-finding settlement
#' rule; see \code{\link[RangeShiftR]{Settlement}}), then there is the option that its path is straightened (\code{StraightenPath=TRUE}). This means
#' that it leaves the patch as soon as possible in order to search for another patch. This is achieved by increasing its \code{DP} ten-fold. However, in
#' certain types of model, e.g. when arbitrary contiguous patches have been defined for what is basically a continuous population, this can lead
#' to the path always being straightened, as an individual enters a new patch as soon as it has left the one it has rejected. In such cases, it is
#' advisable to disable the feature (\code{StraightenPath=FALSE}), although care must be taken that individuals do not become trapped in patches
#' surrounded by very high cost matrix.
#'
#' When inter-individual variability is activated (set \code{IndVar=TRUE}), the four SMS movement traits \code{DP}, as well as - if dispersal goal bias
#' is enabled - \code{GB}, \code{AlphaDB} and \code{BetaDB} can evolve.
#' For each trait the population initial mean and standard deviations must be set (instead of the constant value), as well as its scaling factor
#' (see \code{\link[RangeShiftR]{Genetics}}).\cr
#'
#' As currently implemented, there is no sex- or stage- dependence of SMS traits.
#'
#' The production of the output ‘dispersal heat maps’, which show the total number of times each cell is visited by dispersing individuals within a
#' single replicate simulation, can be enabled with the option \code{SMSHeatMap} in \code{\link[RangeShiftR]{Simulation}}.
#'
#' \emph{Costs layer} \cr
#' Critical for the outcomes of SMS are the relative costs assigned to the different habitats (as it is also the case for
#' the LCP approach). Habitat costs or resistance to movement can be set manually for each habitat code or imported as a raster map, which
#' allows for costs to be a function of multiple variables instead of a simple value associated to the habitat type.
#' In the latter case, the file names are provided in the \code{\link[RangeShiftR]{ImportedLandscape}} module.
#' The specified map has to match the landscape raster in extent, coordinates and resolution, and each cell contains a cost value, with the minimal possible cost being \eqn{1}.
#' Importing a cost layer is the only option when the landscape comprises habitat coverage or quality.
#'
#' \emph{Mortality} \cr
#' There are two main sources of mortality: First, dispersal mortality can arise as a result of individuals failing to reach suitable habitat. Some individuals may fail
#' to find suitable habitat before they use up their maximum number of movement steps (\code{MaxSteps} in \code{\link[RangeShiftR]{Settlement}}).
#' In this first case, dispersal mortality clearly depends upon the proportion of suitable habitat in the landscape and will increase as the
#' availability of habitat declines.\cr
#' A second source of dispersal mortality can be specified by the user in form of a per-step probability of mortality (\code{StepMort}.
#' This can be useful for representing mortality risks that increase with distance or time spent travelling \insertCite{bonte2012costs}{RangeShiftR}.
#' Additionally, it is possible that the per-step mortality varies according to the nature of the local environment by providing the \code{CostMap}.
#'
#' Note that the total dispersal mortality
#' experienced will be the sum of the mortalities due to the two sources identified above and, in parameterising the model, it will be important
#' to recognize this such that dispersal mortality is not double-accounted.
#'
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "StochMove"
#' @author Anne-Kathleen Malchow
#' @name SMS
#' @export SMS
SMS <- setClass("StochMove", slots = c(PR = "integer_OR_numeric",
                                       PRMethod = "integer_OR_numeric", # Perceptual range method: 1 = arithmetic mean; 2 = harmonic mean; 3 = weighted arithmtic mean
                                       MemSize = "integer_OR_numeric",
                                       GoalType = "integer_OR_numeric", # 0 (none) or 2 (dispersal bias)
                                       IndVar = "logical",
                                       DP = "numeric",
                                       GoalBias = "numeric",
                                       AlphaDB = "numeric",
                                       BetaDB = "integer_OR_numeric",
                                       StraightenPath = "logical",
                                       Costs = "numeric_OR_character", # will determine on C++ level the values of CostHab1 ... CostHabN, CostHabitat, CostMatrix, CostMap, CostMapFile
                                       CostMap = "logical",
                                       StepMort = "numeric")           # will determine on C++ level the values of HabMort, MortConst, MortHab1 ... MortHabN, MortHabitat, MortMatrix
                , prototype = list(PR = 1L,
                                   PRMethod = 1L,
                                   MemSize = 1L,
                                   GoalType = 0L,
                                   IndVar = FALSE,
                                   DP = 1.0,
                                   GoalBias = 1.0,
                                   #AlphaDB = 1.0,
                                   #BetaDB = 100000,
                                   StraightenPath = TRUE,
                                   #Costs,
                                   StepMort = 0.0)
                , contains = "TransferParams"
)
setValidity("StochMove", function(object) {
    msg <- NULL
    if (anyNA(object@PR) || length(object@PR)!=1) {
        msg <- c(msg, "PR must be set and of length 1!")
    }
    else{
        if (object@PR < 1.0) {
            msg <- c(msg, "PR must be >= 1!")
        }
    }
    if (anyNA(object@PRMethod) || length(object@PRMethod)!=1) {
        msg <- c(msg, "PRMethod must be set and of length 1!")
    }
    else{
        if (object@PRMethod != 1 && object@PRMethod != 2 && object@PRMethod != 3) {
            msg <- c(msg, "PRMethod must be either 1, 2 or 3!")
        }
    }
    if (anyNA(object@MemSize) || length(object@MemSize)!=1) {
        msg <- c(msg, "MemSize must be set and of length 1!")
    }
    else{
        if (object@MemSize < 1 || object@MemSize > 14) {
            msg <- c(msg, "MemSize must be between 1 and 14 !")
        }
    }
    if (anyNA(object@GoalType) || length(object@GoalType)!=1) {
        msg <- c(msg, "GoalType must be set and of length 1!")
    }
    else{
        if (object@GoalType != 0 && object@GoalType != 2) {
            msg <- c(msg, "GoalType must be either 0 or 2!")
        }
    }
    if (anyNA(object@IndVar) || length(object@IndVar)!=1) {
        msg <- c(msg, "IndVar must be set and of length 1!")
    }
    if(is.null(msg)){
        if (anyNA(object@DP) || length(object@DP)==0) {
            msg <- c(msg, "DP must be set!")
        }
        else{
            if(object@IndVar){
                if(length(object@DP)==3){
                    if (object@DP[1] < 1.0) {
                        msg <- c(msg, "DP (mean) must be >= 1.0!")
                    }
                    if (object@DP[2] <= 0.0) {
                        msg <- c(msg, "DP (SD) must be strictly positive!")
                    }
                    if (object@DP[3] <= 0.0) {
                        msg <- c(msg, "DP (scaling factor) must be strictly positive !")
                    }
                    if(is.null(msg)){
                        if (object@DP[3] < object@DP[2]) {
                            msg <- c(msg, "DP scaling factor must be greater than or equal to its SD !")
                        }
                    }
                }
                else{
                    msg <- c(msg, "DP must have length 3 if IndVar=TRUE!")
                }
            }
            else {
                if(length(object@DP)==1){
                    if (object@DP < 1.0) {
                        msg <- c(msg, "DP must be >= 1.0 !")
                    }
                }
                else{
                    msg <- c(msg, "DP must have length 1 if IndVar=FALSE!")
                }
            }
        }
        if (object@GoalType) { # GoalType = 2
            if (anyNA(object@GoalBias) || length(object@GoalBias)==0) {
                msg <- c(msg, "GoalBias strength must be set!")
            }
            else{
                if(object@IndVar){
                    if(length(object@GoalBias)==3){
                        if (object@GoalBias[1] < 1.0) {
                            msg <- c(msg, "GoalBias strength (mean) must be >= 1.0!")
                        }
                        if (object@GoalBias[2] <= 0.0) {
                            msg <- c(msg, "GoalBias strength (SD) must be strictly positive!")
                        }
                        if (object@GoalBias[3] <= 0.0) {
                            msg <- c(msg, "GoalBias strength (scaling factor) must be strictly positive !")
                        }
                        if(is.null(msg)){
                            if (object@GoalBias[3] < object@GoalBias[2]) {
                                msg <- c(msg, "GoalBias strength scaling factor must be greater than or equal to its SD !")
                            }
                        }
                    }
                    else{
                        msg <- c(msg, "GoalBias strength must have length 3 if IndVar=TRUE!")
                    }
                }
                else {
                    if(length(object@GoalBias)==1){
                        if (object@GoalBias < 1.0 ) {
                            msg <- c(msg, "GoalBias strength must be >= 1.0 !")
                        }
                    }
                    else{
                        msg <- c(msg, "GoalBias strength must have length 1 if IndVar=FALSE!")
                    }
                }
            }
            if (anyNA(object@AlphaDB) || length(object@AlphaDB)==0) {
                msg <- c(msg, "AlphaDB must be set!")
            }
            else{
                if(object@IndVar){
                    if(length(object@AlphaDB)==3){
                        if (object@AlphaDB[1] <= 0.0) {
                            msg <- c(msg, "AlphaDB (mean) must be strictly positive!")
                        }
                        if (object@AlphaDB[2] <= 0.0) {
                            msg <- c(msg, "AlphaDB (SD) must be strictly positive!")
                        }
                        if (object@AlphaDB[3] <= 0.0) {
                            msg <- c(msg, "AlphaDB (scaling factor) must be strictly positive !")
                        }
                        if(is.null(msg)){
                            if (object@AlphaDB[3] < object@AlphaDB[2]) {
                                msg <- c(msg, "AlphaDB scaling factor must be greater than or equal to its SD !")
                            }
                        }
                    }
                    else{
                        msg <- c(msg, "AlphaDB must have length 3 if IndVar=TRUE!")
                    }
                }
                else {
                    if(length(object@AlphaDB)==1){
                        if (object@AlphaDB <= 0.0) {
                            msg <- c(msg, "AlphaDB must be strictly positive!")
                        }
                    }
                    else{
                        msg <- c(msg, "AlphaDB must have length 1 if IndVar=FALSE!")
                    }
                }
            }
            if (anyNA(object@BetaDB) || length(object@BetaDB)==0) {
                msg <- c(msg, "BetaDB must be set!")
            }
            else{
                if(object@IndVar){
                    if(length(object@BetaDB)==3){
                        if (object@BetaDB[1] < 1.0) {
                            msg <- c(msg, "BetaDB (mean) must be >= 1 steps!")
                        }
                        if (object@BetaDB[2] <= 0.0) {
                            msg <- c(msg, "BetaDB (SD) must be strictly positive!")
                        }
                        if (object@BetaDB[3] <= 0.0) {
                            msg <- c(msg, "BetaDB (scaling factor) must be strictly positive !")
                        }
                        if(is.null(msg)){
                            if (object@BetaDB[3] < object@BetaDB[2]) {
                                msg <- c(msg, "BetaDB scaling factor must be greater than or equal to its SD !")
                            }
                        }
                    }
                    else{
                        msg <- c(msg, "BetaDB must have length 3 if IndVar=TRUE!")
                    }
                }
                else {
                    if(length(object@BetaDB)==1){
                        if (object@BetaDB < 1.0) {
                            msg <- c(msg, "BetaDB must be >= 1 steps!")
                        }
                    }
                    else{
                        msg <- c(msg, "BetaDB must have length 1 if IndVar=FALSE!")
                    }
                }
            }
        }
    }
    if (anyNA(object@StraightenPath) || length(object@StraightenPath)!=1) {
        msg <- c(msg, "StraightenPath must be set and of length 1!")
    }
    if (anyNA(object@Costs) || length(object@Costs)==0) {
        msg <- c(msg, "Costs must be set!")
    }
    else{
        if (class(object@Costs)=="numeric") {
            if (any(object@Costs < 1) ) {
                msg <- c(msg, "Costs must equal 1 (minimum cost) or be larger!")
            }
        }
        else {
            if (class(object@Costs)=="character") {
                if (object@Costs != "file") {
                    msg <- c(msg, "Costs has a wrong format! Must be either numeric or the keyword \"file\".")
                }
            }
            else{ # neither numeric nor character
                msg <- c(msg, "Costs has a wrong format!")
            }
        }
    }
    if (anyNA(object@StepMort) || length(object@StepMort)==0) {
        msg <- c(msg, "StepMort must be set!")
    }
    else{
        if (any(object@StepMort < 0.0 | object@StepMort >= 1.0) ) {
            msg <- c(msg, "StepMort probabilities must be within the half-open interval [0,1) !")
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "StochMove", function(.Object,...) {
    this_func = "SMS(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (!.Object@GoalType) { # GoalType = 0
        .Object@GoalBias = 1.0
        if (!is.null(args$GoalBias)) {
            warning(this_func, "GoalBias", warn_msg_ignored, "since GoalType = 0.", call. = FALSE)
        }
        .Object@AlphaDB = -9L
        if (!is.null(args$AlphaDB)) {
            warning(this_func, "AlphaDB", warn_msg_ignored, "since GoalType = 0.", call. = FALSE)
        }
        .Object@BetaDB = -9L
        if (!is.null(args$BetaDB)) {
            warning(this_func, "BetaDB", warn_msg_ignored, "since GoalType = 0.", call. = FALSE)
        }
    }
    if (class(.Object@Costs)=="character") {
        if (.Object@Costs == "file") .Object@CostMap = TRUE
    }
    else {
        if (class(.Object@Costs)=="numeric") {
            .Object@CostMap = FALSE
        }
    }
    .Object}
)
setMethod("show", "StochMove", function(object){
    callNextMethod()
    cat("   PR =", object@PR, ", MemSize =", object@MemSize, "\n")
    if (object@PRMethod == 1) cat("   Method: Arithmetic mean \n")
    if (object@PRMethod == 2) cat("   Method: Harmonic mean \n")
    if (object@PRMethod == 3) cat("   Method: Weighted arithmetic mean \n")

    if (object@IndVar) {
        cat("   DP       =", object@DP[1], "\u00B1" , object@DP[2], ", scale \u03bc =", object@DP[3], "\n")
    }
    else {
        cat("   DP       =", object@DP, "\n")
    }

    if (object@GoalType) {
        cat("   GoalType: Dispersal bias:\n")
        if (object@IndVar) {
            cat("   GoalBias =", object@GoalBias[1], "\u00B1" , object@GoalBias[2], ", scale \u03bc =", object@GoalBias[3], "\n")
            cat("   AlphaDB =", object@AlphaDB[1], "\u00B1" , object@AlphaDB[2], ", scale \u03bc =", object@AlphaDB[3], "\n")
            cat("   BetaDB  =", object@BetaDB[1], "\u00B1" , object@BetaDB[2], ", scale \u03bc =", object@BetaDB[3], "\n")
        }
        else {
            cat("   GoalBias =", object@GoalBias, "\n")
            cat("   AlphaDB =", object@AlphaDB, ", BetaDB =", object@BetaDB, "\n")
        }
    }
    cat("   StraightenPath =", object@StraightenPath, "\n")
    cat("   Costs ")
    if (object@CostMap) cat(" read from file \n")
    else cat("=", object@Costs, "\n")
    cat("   StepMort =", object@StepMort, "\n")
    }
)
setMethod("plotProbs", "StochMove", function(x, xmax = NULL, ymax = NULL){
    # get parameters
    gb  <- x@GoalBias
    if (x@GoalType == 2) {
        alp <- x@AlphaDB
        bet <- x@BetaDB
        main = "Dispersal Bias"
    } else main = "Goal Bias is disabled"
    # New plot
    if (x@GoalType == 2) {
        if (is.null(xmax)) xmax = 2*bet[1]
        xvals = seq(0, xmax)
    }
    else {if(is.null(xmax)) xmax <- 100}
    if (is.null(ymax)) {ymax = gb[1]*1.1}
    plot(NULL, type = "n", main = main, xlab = "Nr of steps", ylab = "Bias strength", xlim = c(0,xmax), ylim = c(1.0,ymax))
    leg.txt <- c()
    # add to plot:
    if (x@IndVar) {
    # plot shaded sd interval
        if (x@GoalType == 2) {
            res <- matrix(ncol = 8, nrow = length(xvals))
            res[,1] <- 1+densdep(xvals, A0 = (gb[1]-gb[2]-1), alpha = -(alp[1]-alp[2]), beta = (bet[1]-bet[2]))
            res[,2] <- 1+densdep(xvals, A0 = (gb[1]-gb[2]-1), alpha = -(alp[1]-alp[2]), beta = (bet[1]+bet[2]))
            res[,3] <- 1+densdep(xvals, A0 = (gb[1]-gb[2]-1), alpha = -(alp[1]+alp[2]), beta = (bet[1]-bet[2]))
            res[,4] <- 1+densdep(xvals, A0 = (gb[1]-gb[2]-1), alpha = -(alp[1]+alp[2]), beta = (bet[1]+bet[2]))
            res[,5] <- 1+densdep(xvals, A0 = (gb[1]+gb[2]-1), alpha = -(alp[1]-alp[2]), beta = (bet[1]-bet[2]))
            res[,6] <- 1+densdep(xvals, A0 = (gb[1]+gb[2]-1), alpha = -(alp[1]-alp[2]), beta = (bet[1]+bet[2]))
            res[,7] <- 1+densdep(xvals, A0 = (gb[1]+gb[2]-1), alpha = -(alp[1]+alp[2]), beta = (bet[1]-bet[2]))
            res[,8] <- 1+densdep(xvals, A0 = (gb[1]+gb[2]-1), alpha = -(alp[1]+alp[2]), beta = (bet[1]+bet[2]))
            polygon(c(xvals,rev(xvals)), c(apply(res, 1, min), rev(apply(res, 1, max))), border=NA, col='grey80')
        }
        else {#constant
            polygon(c(0,xmax,xmax,0), c(rep(gb[1]-gb[2],2),rep(gb[1]+gb[2],2)), border=NA, col='grey80')
        }
    }
    # plot lines
    if (x@GoalType == 2) {
        lines(xvals, 1+densdep(xvals, A0 = (gb[1]-1), alpha = -alp[1], beta = bet[1]), type = "b", lty = 1, col = "blue")
    }else { # constant
        lines(x=c(0,xmax), y=rep(gb[1],2), type = "b", lty = 1, col = "blue")
    }
})


## Transfer-class CORRRW

#' Set up a Correlated Random Walk
#'
#' A method to describe \code{\link[RangeShiftR]{Transfer}}:
#' A simple correlated random walk without any bias; implemented in continuous space on the top of the landscape grid.
#'
#' @usage CorrRW(StepLength = 1, Rho = 0.5,
#'       IndVar = FALSE,
#'       StraightenPath = FALSE,
#'       StepMort = 0.0)
#' @param StepLength Step length given in meters, defaults to \eqn{1}.\cr If \code{IndVar=TRUE}, expects a vector of length three
#' specifying (Mean, SD, TraitScaleFactor) of \code{StepLength}.
#' @param Rho Correlation parameter \eqn{ρ}, defaults to \eqn{0.5}. Must be in the open interval \eqn{(0,1)}.\cr If \code{IndVar=TRUE},
#' expects a vector of length three specifying (Mean, SD, TraitScaleFactor) of \code{Rho}.
#' @param IndVar Individual variability in CorrRW traits (i.e. \code{StepLength} and \code{Rho})? Defaults to \code{FALSE}.
#' @param StraightenPath Straighten path after decision not to settle in a patch? Defaults to \code{TRUE}, see Details below.
#' @param StepMort Per-step mortality probability. Can be either \emph{constant}, in which case a single numeric is expected (the default, with
#' value \eqn{0.0}) or \emph{habitat-specific}, in which case a numeric vector is expected with a length of, respectively, \code{Nhabitats} for an
#' \code{\link[RangeShiftR]{ImportedLandscape}} with habitat codes (i.e. \code{HabPercent=FALSE})) or length \eqn{2} for an
#' \code{\link[RangeShiftR]{ArtificialLandscape}} (mortality probabilities for matrix and habitat cells).\cr
#' All values must be within the half-open interval \eqn{[0,1)}.\cr
#' For an imported habitat quality landscape (\code{HabPercent=TRUE}), only constant per-step mortality is allowed.
#' @details
#' Individuals take steps of a constant \code{StepLength}; the direction is sampled from a wrapped Cauchy distribution having a
#' correlation parameter \eqn{Rho} in the range \eqn{0} to \eqn{1} \insertCite{barton2009evolution,zollner1999search}{RangeShiftR}.
#' As for \code{\link[RangeShiftR]{SMS}}, all individuals take each step
#' simultaneously. In the case of patch-based models,
#' \eqn{Rho} is automatically set to \eqn{0.99} until the individual steps outside its natal patch, after which the value of
#' \eqn{Rho} set by the user is restored. \cr
#' The \code{StepLength} and \eqn{Rho} can be set to vary between individuals and evolve (set \code{IndVar=TRUE}).
#' In this case, each individual exhibits two traits for these two parameters.
#' For each trait the initial mean and standard deviations must be set, as well as the TraitScaleFactor (see \code{\link[RangeShiftR]{Settlement}}),
#' instead of only one constant value each.\cr
#' Note that the step length may not evolve below one fifth of
#' the landscape resolution, and correlation may not evolve above \eqn{0.999}. \cr
#' Per-step mortality is not allowed to vary between individuals or to evolve. \cr
#' There is no implementation of sex- or stage-specific CRW.
#'
#' When an individual arrives in a non-natal patch and decides not to settle there (as a result of a density-dependent or mate-finding settlement
#' rule, see \code{\link[RangeShiftR]{Settlement}}), then there is the option that its path is straightened (\code{StraightenPath=TRUE}). This means
#' that it leaves the patch as soon as possible in order to search for another patch. This is achieved by increasing its path correlation to
#' \code{Rho}\eqn{=0.999}. However, in certain types of model, e.g. when arbitrary contiguous patches have been defined for what is basically a continuous
#' population, this can lead to the path always being straightened, as an individual enters a new patch as soon as it has left the one it has
#' rejected. In such cases, it is advisable to disable the feature (\code{StraightenPath=FALSE}), although care must be taken that individuals
#' do not become trapped in patches surrounded by very high cost matrix.
#'
#' \emph{Mortality}\cr
#' There are two main sources of mortality: First, dispersal mortality can arise as a result of individuals failing to reach suitable habitat. Some individuals may fail
#' to find suitable habitat before they use up their maximum number of movement steps (\code{MaxSteps} in \code{\link[RangeShiftR]{Settlement}}).
#' In this first case, dispersal mortality clearly depends upon the proportion of suitable habitat in the landscape and will increase as the
#' availability of habitat declines.\cr
#' A second source of dispersal mortality can be specified by the user in form of a per-step probability of mortality (\code{StepMort}.
#' This can be useful for representing mortality risks that increase with distance or time spent travelling \insertCite{bonte2012costs}{RangeShiftR}.
#'
#' Note that the total dispersal mortality experienced will be the sum of the mortalities due to the two sources identified above and, in parameterising the model,
#' it will be important to recognize this such that dispersal mortality is not double-accounted.
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "CorrRW"
#' @author Anne-Kathleen Malchow
#' @name CorrRW
#' @export CorrRW
CorrRW <- setClass("CorrRW", slots = c(IndVar = "logical",
                                       StepLength = "numeric",
                                       Rho = "numeric",
                                       StraightenPath = "logical",
                                       StepMort = "numeric")  # will determine on C++ level the values of HabMort, MortConst, MortHab1 ... MortHabN, MortHabitat, MortMatrix
                , prototype = list(IndVar = FALSE,
                                   StepLength = 1L,
                                   Rho = 0.5,
                                   StraightenPath = FALSE,
                                   StepMort = c(0.0))
                , contains = "TransferParams"
)
setValidity("CorrRW", function(object) {
    msg <- NULL
    if (anyNA(object@IndVar) || length(object@IndVar)!=1) {
        msg <- c(msg, "IndVar must be set and of length 1!")
    }
    if (anyNA(object@StepLength) || length(object@StepLength)==0) {
        msg <- c(msg, "StepLength must be set!")
    }
    else{
        if(object@IndVar){
            if(length(object@StepLength)==3){
                if (object@StepLength[1] <= 0.0) {
                    msg <- c(msg, "Step Length (mean) must be strictly positive!")
                }
                if (object@StepLength[2] <= 0.0) {
                    msg <- c(msg, "Step Length (SD) must be strictly positive!")
                }
                if (object@StepLength[3] <= 0.0) {
                    msg <- c(msg, "Step Length (scaling factor) must be strictly positive !")
                }
                if(is.null(msg)){
                    if (object@StepLength[3] < object@StepLength[2]) {
                        msg <- c(msg, "Step Length scaling factor must be greater than or equal to its SD !")
                    }
                }
            }
            else{
                msg <- c(msg, "StepLength must have length 3 if IndVar=TRUE!")
            }
        }
        else {
            if(length(object@StepLength)==1){
                if (object@StepLength <= 0.0) {
                    msg <- c(msg, "StepLength must be strictly positive!")
                }
            }
            else{
                msg <- c(msg, "StepLength must have length 1 if IndVar=FALSE!")
            }
        }
    }
    if (anyNA(object@Rho) || length(object@Rho)==0) {
        msg <- c(msg, "Rho must be set!")
    }
    else{
        if(object@IndVar){
            if(length(object@Rho)==3){
                if (object@Rho[1] <= 0.0 || object@Rho[1] >= 1.0) {
                    msg <- c(msg, "Rho (mean) must be within the open interval (1,0)!")
                }
                if (object@Rho[2] <= 0.0 || object@Rho[2] >= 1.0) {
                    msg <- c(msg, "Rho (SD) must be within the open interval (1,0)!")
                }
                if (object@Rho[3] <= 0.0 || object@Rho[3] >= 1.0) {
                    msg <- c(msg, "Rho (scaling factor) must be within the open interval (1,0)!")
                }
                if(is.null(msg)){
                    if (object@Rho[3] < object@Rho[2]) {
                        msg <- c(msg, "Rho scaling factor must be greater than or equal to its SD !")
                    }
                }
            }
            else{
                msg <- c(msg, "Rho must have length 3 if IndVar=TRUE!")
            }
        }
        else {
            if(length(object@Rho)==1){
                if (object@Rho <= 0.0 || object@Rho >= 1.0) {
                    msg <- c(msg, "Correlation coefficient (Rho) must be within the open interval (1,0)!")
                }
            }
            else{
                msg <- c(msg, "Rho must have length 1 if IndVar=FALSE!")
            }
        }
    }
    if (anyNA(object@StraightenPath) || length(object@StraightenPath)!=1) {
        msg <- c(msg, "StraightenPath must be set and of length 1!")
    }
    if (anyNA(object@StepMort) || length(object@StepMort)==0) {
        msg <- c(msg, "StepMort must be set!")
    }
    else{
        if (any(object@StepMort < 0.0 | object@StepMort >= 1.0) ) {
            msg <- c(msg, "StepMort probabilities must be within the half-open interval [0,1) !")
        }
    }
    if (is.null(msg)) TRUE else msg}
)
# setMethod("initialize", "CorrRW", function(.Object,...) {
#     this_func = "CorrRW(): "
#     args <- list(...)
#     .Object <- callNextMethod()
#     if ( length(args) == 0 ) {
#         validObject(.Object)
#     }
#     .Object}
# )
setMethod("show", "CorrRW", function(object){
    callNextMethod()
    if (object@IndVar) {
        cat("   StepLength =", object@StepLength[1], "\u00B1" , object@StepLength[2], ", scale \u03bc =", object@StepLength[3], "\n")
        cat("   Rho =", object@Rho[1], "\u00B1" , object@Rho[2], ", scale \u03bc =", object@Rho[3], "\n")
    }
    else {
        cat("   StepLength =", object@StepLength, "\n")
        cat("   Rho =", object@Rho, "\n")
    }
    cat("   StraightenPath =", object@StraightenPath, "\n")
    cat("   StepMort =", object@StepMort, "\n")
    }
)


### SUBCLASS SETTLEMENTPARAMS

# from RS 'Settlement' file

#' Set Settlement Parameters
#'
#' Settlement, or immigration, is the last phase of dispersal, when the organism stops in a new cell or patch of breeding habitat. The
#' available settlement conditions vary depending on the used \code{\link[RangeShiftR]{Transfer}} type. In any case, dispersing individuals
#' are not allowed to settle in their natal cell or patch and can only settle in suitable habitat.\cr
#' \emph{RangeShiftR} incorporates some basic settlement rules that can be stage- or sex-specific or both (set \code{StageDep}, \code{SexDep}).
#' If a movement process is used, density-dependence (\code{DensDep}) and/or inter-individual variability (\code{IndVar}) are available.
#' @usage Settlement(StageDep = FALSE, SexDep = FALSE,
#'           Settle = 0, FindMate = FALSE,
#'           DensDep = FALSE,
#'           IndVar = FALSE, TraitScaleFactor,
#'           MinSteps = 0, MaxSteps = 0, MaxStepsYear = 0)
#' @param StageDep Stage-dependent settlement requirements? (default: \code{FALSE})
#' @param SexDep Sex-dependent settlement requirements? (default: \code{FALSE})
#' @param Settle Settlement codes (for \code{DispersalKernel}) or settlement probability parameters (for \emph{Movement process}) for all
#' stages/sexes, defaults to \eqn{0} (i.e. 'die when unsuitable' for \emph{DispersalKernel} and
#' 'always settle when suitable' for \emph{Movement process}). See the Details below.
#' @param FindMate Mating requirements to settle? Set for all stages/sexes. Must be \code{FALSE} (default) in a female-only model.
#' @param DensDep Movement process only: Density-dependent settlement probability? (default: \code{FALSE})
#' @param IndVar Movement process only: Individual variability in settlement probability traits? Must be \code{FALSE} (default),
#' if \code{DensDep=FALSE} or \code{StageDep=TRUE}.
#' @param TraitScaleFactor Movement process only, required if \code{IndVar=TRUE}: Scaling factors for the three traits of density-dependent settlement,
#' a numeric of length \eqn{3}.
#' @param MinSteps Movement process only: Minimum number of steps. Defaults to \eqn{0}.
#' @param MaxSteps Movement process only: Maximum number of steps. Must be \eqn{0} or more; set to \eqn{0} (default) for “per-step mortality only”.
#' @param MaxStepsYear Movement process and stage-structured population only: Maximum number of steps per year, if there are more than \eqn{1} reproductive seasons (option \code{RepSeasons} in \code{\link[RangeShiftR]{StageStructure}}).
#' Must be \eqn{0} or more. If \eqn{0}, every individual completes the dispersal phase in one year, i.e. between two successive reproduction phases.
#' @details
#' The settlement phase is determined by a suite of strategies, behaviours and reaction norms that lead individuals to the decision to stop
#' in a particular place. Habitat selection, mate finding and density dependence are probably three of the main processes involved,
#' but not the only ones. Like emigration, settlement is a complex process affected by multiple criteria including inter-individual
#' variability and context dependencies.
#'
#' The type of implemented settlement rules depends on the movement model utilized for the \code{\link[RangeShiftR]{Transfer}}.
#' In any case, dispersing individuals are not allowed to settle in their natal cell or patch.\cr
#' \emph{RangeShiftR} incorporates some basic settlement rules that can be stage- or sex-specific or both (set \code{StageDep}, \code{SexDep}).
#' Inter-individual variability (\code{IndVar}) is implemented only for movement processes and then for the three traits
#' determining density-dependent settlement (\ifelse{html}{\out{S<sub>0</sub>}}{\eqn{S_0}}, \ifelse{html}{\out{&alpha;<sub>S</sub>}}{\eqn{α_S}},
#' \ifelse{html}{\out{&beta;<sub>S</sub>}}{\eqn{β_S}}; see below). In this case, settlement may not be stage-dependent.\cr
#'
#' \emph{Settlement with dispersal kernels}\cr
#' When using a \code{\link[RangeShiftR]{DispersalKernel}}, individuals are displaced directly from the starting location to the arrival location. The suitability
#' of the arrival cell or patch determines whether the disperser is successful or not.\cr For species with \emph{non-overlapping generations},
#' where individuals have only one chance to disperse and reproduce, the model has two options if the arrival cell is unsuitable: the
#' individual either dies (\code{Settle} code \eqn{0}) or it can move to one of the eight neighbouring cells in the case that at least one
#' of them is suitable (\code{Settle} code \eqn{2}). In
#' the latter case, if more than one of the neighbouring cells is suitable, the individual is placed in one of them chosen randomly.
#' For patch-based models, if the arrival patch is unsuitable, the individual either dies or can move to a randomly chosen neighbouring
#' suitable patch, provided that the new patch is only one cell apart from the arrival patch.\cr For species with \emph{overlapping generations}
#' (i.e. a \code{\link[RangeShiftR]{StageStructure}}d population),
#' where individuals can disperse over multiple seasons, there are two additional options: First, if the arrival cell/patch is unsuitable,
#' the individual can stay there waiting until the next dispersal event when it will disperse again according to the set kernel (\code{Settle} code \eqn{1}). Second,
#' if both the arrival cell/patch and all eight neighbouring cells, or all eventual neighbouring patches, are unsuitable the individual
#' can wait in the arrival cell/patch before moving again at the next dispersal event (\code{Settle} code \eqn{3}). The arrival cell/patch is considered suitable if
#' it contains breeding habitat.
#'
#' If the settlement condition is the same for the entire population (\code{StageDep=FALSE} and \code{SexDep=FALSE}), then the parameter
#' \code{Settle} takes a single integer, specifying the corresponding settlement condition code (\eqn{0,1,2} or \eqn{3}).
#' In case of stage- and/or sex-specific settlement conditions, \code{Settle} must be an integer matrix with the first one/two columns
#' stating 'stage' and/or 'sex' (\eqn{0} for \emph{female} and \eqn{1} for \emph{male}), and the last column giving the respective settlement condition.
#' Exactly one row for each stage/sex-combination is required, in no particular order.
#'
#' Settlement condition codes: If the individuals current step ends in unsuitable habitat, it:\cr
#' \eqn{0} = die (default),\cr
#' \eqn{1} = wait (stage-structured models only),\cr
#' \eqn{2} = randomly choose a suitable neighbouring cell or die,\cr
#' \eqn{3} = randomly choose a suitable neighbouring cell or wait (stage-structured models only).\cr
#'
#' Simple example for sex-dependence only: Females choose a neighbouring cell or wait, males wait: \tabular{cc}{\out{&emsp;} 0 \tab \out{&emsp;} 3 \cr \out{&emsp;} 1 \tab \out{&emsp;} 0 }
#'
#' \emph{Settlement with movement processes}\cr
#' If individuals are dispersing by one of the two movement processes implemented (\code{\link[RangeShiftR]{SMS}} or
#' \code{\link[RangeShiftR]{CorrRW}}), at each step (made simultaneously) they each evaluate their current cell or patch for the
#' possibility of settling. This allows for the implementation of more complex settlement rules. The simplest one is that the individual
#' decides to stop if there is suitable habitat; this is in any case a necessary condition (set \code{DensDep=FALSE}).\cr
#' If a Settlement module with a constant \eqn{S_0=0} (the default) is used in a model with \emph{movement process},
#' it gets converted to \eqn{S_0=1.0}, i.e. 'always settle when habitat is suitable'. \cr
#' \cr
#' Furthermore, the settlement decision can be density-dependent (set \code{DensDep=TRUE}). In this case, the individual has a probability \ifelse{html}{\out{p<sub>S</sub>}}{\eqn{p_S}}
#' of settling in the cell or patch \eqn{i}, given by:
#'
#' \ifelse{html}{\out{&emsp;&emsp; p<sub>S</sub>(i,t) = S<sub>0</sub> / ( 1 + e<sup>-&alpha;<sub>S</sub> (N(i,t) / K(i,t) - &beta;<sub>S</sub>) </sup> ) } }{\deqn{ p_S(i,t) = S_0 / ( 1 + exp[-α_S (N(i,t)/K(i,t) - β_S) ] ) } }
#'
#' In the case of stage-structured models the above equation is modified to:
#'
#' \ifelse{html}{\out{&emsp;&emsp; p<sub>S</sub>(i,t) = S<sub>0</sub> / ( 1 + e<sup>-&alpha;<sub>S</sub> (b(i,t) * N(i,t) - &beta;<sub>S</sub>) </sup> ) } }{\deqn{ p_S(i,t) = S_0 / ( 1 + exp[-α_S (b(i,t) N(i,t) - β_S) ] ) } }
#'
#' In the first case, \eqn{K(i,t)} is the carrying capacity of the cell/patch \eqn{i} at time \eqn{t} given by \code{K_or_DensDep}.
#' In the latter case, \eqn{b(i,t)} represents the strength of density dependence that is given by the inverse of \code{K_or_DensDep}.\cr
#' Further, \ifelse{html}{\out{S<sub>0</sub>}}{\eqn{S_0}} is the maximum settlement probability,
#' \eqn{N(i,t)} is the number of individuals in the cell/patch \eqn{i} at time \eqn{t},
#' \ifelse{html}{\out{&beta;<sub>S</sub>}}{\eqn{β_S}} is the inflection point of the function and
#' \ifelse{html}{\out{&alpha;<sub>S</sub>}}{\eqn{α_S}} is the slope at the inflection point.\cr
#'
#' Inter-individual variability \code{IndVar=TRUE} and thus evolution is implemented only for the three traits determining density-dependent settlement
#' (\code{DensDep=TRUE}), and if so, it may not be stage-dependent (\code{StageDep=FALSE}).
#' For each trait the initial distribution in the population (as mean and standard variation) must be set in \code{Settle} (instead of only one constant value),
#' as well as their scaling factors in \code{TraitScaleFactor} (see \code{\link[RangeShiftR]{Genetics}}).
#'
#' The parameters that determine the settlement probabilities have to be provided via the parameter \code{Settle}, which generally takes a numeric matrix, or - if only a single constant probability is
#' used (i.e. \code{DensDep, IndVar, StageDep, SexDep = FALSE}) - a single numeric.
#' The format of the matrix is defined as follows: The number of columns depend on the options \code{DensDep} and \code{IndVar}. If \code{DensDep=FALSE}, the
#' density-independent probability \ifelse{html}{\out{p<sub>S</sub>}}{\eqn{p_S}} must be specified. If \code{DensDep=TRUE}, the functional parameters \ifelse{html}{\out{S<sub>0</sub>}}{\eqn{S_0}},
#' \ifelse{html}{\out{&alpha;<sub>S</sub>}}{\eqn{α_S}} and \ifelse{html}{\out{&beta;<sub>S</sub>}}{\eqn{β_S}} (cf. equation above) must be specified.
#' Additionally, if \code{IndVar=FALSE}, these traits are fixed, but if \code{IndVar=TRUE} each of them is replaced by two parameters: their respective initial mean and
#' standard deviation. They are used to normally distribute the traits values among the individuals of the initial population. Additionally, the \code{TraitScaleFactor} of
#' these traits have to be set.
#'
#' All parameters have to be given for each stage/sex if the respective dependency is enabled. If \code{StageDep=TRUE}, state the corresponding stage in the first column.
#' If \code{SexDep=TRUE}, state the corresponding stage in the next (i.e. first/second) column, with \eqn{0} for \emph{female} and \eqn{1} for \emph{male}.
#' The following table lists the required columns and their correct order for different settings:
#'
#' \tabular{ccccc}{DensDep \tab IndVar \tab StageDep \tab SexDep \tab columns \cr
#  F \tab F \tab F \tab F \tab \ifelse{html}{\out{p<sub>S</sub>}}{\eqn{p_S}} \cr
#'  F \tab F \tab F \tab F \tab - not applicable - \cr
#'  F \tab F \tab T \tab F \tab stage \cr
#'  F \tab F \tab F \tab T \tab sex \cr
#'  F \tab F \tab T \tab T \tab stage, sex \cr
#'  T \tab F \tab F \tab F \tab \ifelse{html}{\out{S<sub>0</sub>}}{\eqn{S_0}}, \ifelse{html}{\out{&alpha;<sub>S</sub>}}{\eqn{α_S}}, \ifelse{html}{\out{&beta;<sub>S</sub>}}{\eqn{β_S}} \cr
#'  \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \tab \out{&#8942;} \cr
#'  T \tab F \tab T \tab T \tab stage, sex, \ifelse{html}{\out{S<sub>0</sub>}}{\eqn{S_0}}, \ifelse{html}{\out{&alpha;<sub>S</sub>}}{\eqn{α_S}}, \ifelse{html}{\out{&beta;<sub>S</sub>}}{\eqn{β_S}} \cr
#'  T \tab T \tab F \tab F \tab mean\ifelse{html}{\out{(S<sub>0</sub>)}}{\eqn{(S_0)}}, sd\ifelse{html}{\out{(S<sub>0</sub>)}}{\eqn{(S_0)}}, mean\ifelse{html}{\out{(&alpha;<sub>S</sub>)}}{(\eqn{α_S})}, sd\ifelse{html}{\out{(&alpha;<sub>S</sub>)}}{(\eqn{α_S})}, mean\ifelse{html}{\out{(&beta;<sub>S</sub>)}}{(\eqn{β_S})}, sd\ifelse{html}{\out{(&beta;<sub>S</sub>)}}{(\eqn{β_S})} \cr
#'  T \tab T \tab F \tab T \tab sex, mean\ifelse{html}{\out{(S<sub>0</sub>)}}{\eqn{(S_0)}}, sd\ifelse{html}{\out{(S<sub>0</sub>)}}{\eqn{(S_0)}}, mean\ifelse{html}{\out{(&alpha;<sub>S</sub>)}}{(\eqn{α_S})}, sd\ifelse{html}{\out{(&alpha;<sub>S</sub>)}}{(\eqn{α_S})}, mean\ifelse{html}{\out{(&beta;<sub>S</sub>)}}{(\eqn{β_S})}, sd\ifelse{html}{\out{(&beta;<sub>S</sub>)}}{(\eqn{β_S})}
#'  }
#'
#' The column headings need not be included, only the numeric matrix is required. The rows require no particular order, but there must be exactly one row for each stage/sex combination.
#' For example, in the case of density-, stage- and sex-dependent settlement with no individual variability:
#' \tabular{ccccc}{ \out{&emsp;} 0 \tab \out{&emsp;} 0 \tab \out{&emsp;} 1.0 \tab \out{&emsp;} 0.2 \tab \out{&emsp;} 4.0 \cr
#'  \out{&emsp;} 0 \tab \out{&emsp;} 1 \tab \out{&emsp;} 1.0 \tab \out{&emsp;} 0.1 \tab \out{&emsp;} 6.0 \cr
#'  \out{&emsp;} 1 \tab \out{&emsp;} 0 \tab \out{&emsp;} 0.7 \tab \out{&emsp;} 0.5 \tab \out{&emsp;} 2.0 \cr
#'  \out{&emsp;} 1 \tab \out{&emsp;} 1 \tab \out{&emsp;} 0.5 \tab \out{&emsp;} 0.5 \tab \out{&emsp;} 2.0 \cr
#'  \out{&emsp;} 2 \tab \out{&emsp;} 0 \tab \out{&emsp;} 0.05 \tab \out{&emsp;} 1.0 \tab \out{&emsp;} 1.0 \cr
#'  \out{&emsp;} 2 \tab \out{&emsp;} 1 \tab \out{&emsp;} 0.05 \tab \out{&emsp;} 1.0 \tab \out{&emsp;} 1.0
#' }
#'
#' To avoid having individuals moving perpetually because they cannot find suitable conditions to settle, the model requires a maximum number
#' of steps (\code{MaxSteps}) or a per-step mortality (within the \code{Transfer} object), or both, to be set.
#' The maximum number of steps defines the maximum time length of the transfer period.
#' When an individual reaches the maximum number of steps, it stops where it is regardless of the suitability of the location. In
#' the case of \emph{non-overlapping generations} this results in automatic death if the individual stops in unsuitable habitat. For species that can disperse
#' over multiple seasons, i.e. a \emph{stage-structured population}, the model additionally requires a maximum number of steps per dispersal event (\code{MaxStepsYear}); on reaching that limit, the individual will
#' stop where it is and the next season, if still alive, it will move again.
#'
#' An additional rule that can be set constitutes a minimum number of steps (\code{MinSteps}) that each individual must take before settlement can take
#' place. This is useful for simulating situations where animals, in a ‘dispersal mode’, will keep moving and not consider settling even
#' if suitable conditions are available \insertCite{@e.g. @barton2012risky}{RangeShiftR}.
#'
#' \emph{Mating requirements}\cr
#' Sexual species may be required to find a mate, i.e. there has to be at least one individual of the opposite sex present for the cell/patch to be considered suitable for settlement.
#' Density-dependence and mating requirements can also be combined together to determine the settlement decision.
#'
#' The application of this mating condition can be switched on or off using the parameter \code{FindMate}, which takes either a single
#' logical if \code{StageDep=FALSE} and \code{SexDep=FALSE} or, otherwise, a logical vector of same length as the number rows in
#' \code{Settle}, using same order.
#'
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "SettlementParams"
#' @author Anne-Kathleen Malchow
#' @name Settlement
#' @export Settlement
Settlement <- setClass("SettlementParams", slots = c(StageDep = "logical",
                                                     SexDep = "logical",
                                                     Settle = "matrix_OR_numeric",      # Settlement conditions for all sexes/stages. Settlement rule if the arrival cell/patch is unsuitable: 0 = die, 1 = wait, 2 = randomly choose a suitable cell/patch or die, 3 = randomly choose a suitable cell/patch or wait
                                                     FindMate = "logical",
                                                     DensDep = "logical",               # For MovementProcess only!
                                                     IndVar = "logical",                # For MovementProcess only!
                                                     TraitScaleFactor = "numeric",        # For MovementProcess only! For IndVar only. Set for stage=0 only and 1 bzw 2 sexes
                                                     MinSteps = "integer_OR_numeric",   # For MovementProcess only!
                                                     MaxSteps = "integer_OR_numeric",   # For MovementProcess only!
                                                     MaxStepsYear = "integer_OR_numeric")   # For MovementProcess only!
                       , prototype = list(StageDep = FALSE,
                                          SexDep = FALSE,
                                          Settle = matrix(0L),  # for DispKernel this is "die"
                                          FindMate = FALSE,
                                          DensDep = FALSE,
                                          IndVar = FALSE,
                                          #TraitScaleFactor,
                                          MinSteps = 0L,
                                          MaxSteps = 0L,
                                          MaxStepsYear = 0L)
)

setValidity("SettlementParams", function(object) {
    msg <- NULL
    if (anyNA(object@StageDep) || length(object@StageDep)!=1) {
        msg <- c(msg, "StageDep must be set and of length 1!")
    }
    if (anyNA(object@SexDep) || length(object@SexDep)!=1) {
        msg <- c(msg, "SexDep must be set and of length 1!")
    }
    if (anyNA(object@Settle) || length(object@Settle)==0) {
        msg <- c(msg, "Settle must be set!")
    }
    else{
        if (class(object@Settle)[1]=="numeric" && length(object@Settle)!=1) {
            msg <- c(msg, "Settle must be a matrix!")
        }
    }
    if (anyNA(object@FindMate) || length(object@FindMate)==0) {
        msg <- c(msg, "FindMate must be set!")
    }
    else {
        if (is.null(msg)) {
            if( length(object@FindMate)!=nrow(object@Settle) && length(object@FindMate)!=1 ) {
                msg <- c(msg, "FindMate must have either 1 entry or as many as rows in the Settle matrix!")
            }
        }
    }
    if (anyNA(object@DensDep) || length(object@DensDep)!=1) {
        msg <- c(msg, "DensDep must be set and of length 1!")
    }
    if (anyNA(object@IndVar) || length(object@IndVar)!=1) {
        msg <- c(msg, "IndVar must be set and of length 1!")
    }
    if (is.null(msg)) {
        if ( object@IndVar && !object@DensDep ) {
            msg <- c(msg, "Inter-individual variability (IndVar=TRUE) in settlement traits requires density-dependence (DensDep=TRUE) !")
        }
        if ( object@IndVar && object@StageDep ) {
            msg <- c(msg, "Inter-individual variability (IndVar=TRUE) in stage-dependent (StageDep=TRUE) settlement traits is not implemented!")
        }
    }
    if (is.null(msg)) {
        if (object@IndVar) {
            if (anyNA(object@TraitScaleFactor) || length(object@TraitScaleFactor)==0) {
                msg <- c(msg, "TraitScaleFactor must be set!")
            }
            else{
                if (object@DensDep) {
                    if (length(object@TraitScaleFactor)!=3) {
                        msg <- c(msg, "TraitScaleFactor must have length 3 if DensDep=TRUE!")
                    }
                    else {
                        if (object@TraitScaleFactor[1] <= 0.0 || object@TraitScaleFactor[1] > 1.0 ) {
                            msg <- c(msg, "TraitScaleFactor μ(S_0) must be in the half-open interval (0,1] !")
                        }
                        if (any(object@TraitScaleFactor[2:3] <= 0.0 )) {
                            msg <- c(msg, "TraitScaleFactor μ(α_s) and μ(β_s) must be strictly positive !")
                        }
                    }
                }
                else {
                    if (length(object@TraitScaleFactor)!=1) {
                        msg <- c(msg, "TraitScaleFactor must have length 1 if DensDep=FALSE!")
                    }
                    else {
                        if (object@TraitScaleFactor <= 0 || object@TraitScaleFactor > 1 ) {
                            msg <- c(msg, "TraitScaleFactor μ(S_0) must be in the half-open interval (0,1] !")
                        }
                    }
                }
            }
        }
    }
    if (anyNA(object@MinSteps) || length(object@MinSteps)==0) {
        msg <- c(msg, "MinSteps must be set!")
    }
    else {
        if (any(object@MinSteps < 0)){
            msg <- c(msg, "MinSteps can't be negative!")
        }
    }
    if (anyNA(object@MaxSteps) || length(object@MaxSteps)==0) {
        msg <- c(msg, "MaxSteps must be set!")
    }
    else {
        if (any(object@MaxSteps < 0)){
            msg <- c(msg, "MaxSteps can't be negative!")
        }
    }
    if (is.null(msg)) {
        if( length(object@MinSteps)!=nrow(object@Settle) && length(object@MinSteps)!=1 ) {
            msg <- c(msg, "MinSteps must have either 1 entry or as many as rows in the Settle matrix!")
        }
        if( length(object@MaxSteps)!=nrow(object@Settle) && length(object@MaxSteps)!=1 ) {
            msg <- c(msg, "MaxSteps must have either 1 entry or as many as rows in the Settle matrix!")
        }
    }
    if (object@StageDep) {
        if (anyNA(object@MaxStepsYear) || length(object@MaxStepsYear)==0) {
            msg <- c(msg, "MaxStepsYear must be set!")
        }
        else {
            if (any(object@MaxStepsYear < 0)){
                msg <- c(msg, "MaxStepsYear can't be negative!")
            }
            else{
                if( length(object@MaxStepsYear)!=nrow(object@Settle) && length(object@MaxStepsYear)!=1 ) {
                    msg <- c(msg, "MaxStepsYear must have either 1 entry or as many as rows in the Settle matrix!")
                }
                else {
                    if (any(object@MaxStepsYear>object@MaxSteps)) {
                        msg <- c(msg, "MaxStepsYear can't be larger than MaxSteps!")
                    }
                }
            }
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "SettlementParams", function(.Object,...) {
    this_func = "Settlement(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (!is.null(args$Settle)) {
        if (class(args$Settle)[1]=="numeric" && length(args$Settle)==1) {
            .Object@Settle <- as.matrix(args$Settle)
        }
    }
    if (!.Object@IndVar) {
        .Object@TraitScaleFactor = -9L
        if (!is.null(args$TraitScaleFactor)) {
            warning(this_func, "TraitScaleFactor", warn_msg_ignored, "since IndVar = FALSE.", call. = FALSE)
        }
    }
    .Object}
)
setMethod("show", "SettlementParams", function(object){
    if (object@DensDep) {
        cat("   DensDep =", object@DensDep, "\n")
    }
    if (object@IndVar) {
        cat("   IndVar =", object@IndVar, "\n")
    }
    if (object@StageDep) {
        cat("   StageDep =", object@StageDep, "\n")
    }
    if (object@SexDep) {
        cat("   SexDep =", object@SexDep, "\n")
    }
    cat("   Settlement conditions:\n")
    print(object@Settle)
    if (object@IndVar) {
        cat("   TraitScaleFactor =", object@TraitScaleFactor, "\n")
    }
    cat("   FindMate =", object@FindMate, "\n")
    }
)
setMethod("plotProbs", "SettlementParams", function(x, stage = NULL, sex = NULL, xmax = NULL, ymax = NULL){
    sett <- x@Settle
    if (x@DensDep) {
        # error messages
        if (!is.null(stage)){
            if (x@StageDep) {
                sett <- subset(sett, sett[,1] %in% stage)
            }
            else{ print("This settlement module has no stage-dependency.\n") }
        }
        if (!is.null(sex)){
            if (x@SexDep) {
                if (x@StageDep) sett <- subset(sett, sett[,2] %in% sex)
                else sett <- subset(sett, sett[,1] %in% sex)
            }
            else{ print("This settlement module has no sex-dependency.\n") }
        }
        # get column indices
        if (x@StageDep) {
            if (x@SexDep) {ind_D0 <- 3} else {ind_D0 <- 2}
        }else{
            if (x@SexDep) {ind_D0 <- 2} else {ind_D0 <- 1}
        }
        if (x@IndVar) {IV <- 2} else {IV <- 1}
        # New plot
        if (is.null(xmax)) {
            ind_max <- which.max(sett[,ind_D0+2*IV])
            xmax = min(2*sett[ind_max,ind_D0+2*IV], sett[ind_max,ind_D0+2*IV] + 6.0/abs(sett[ind_max,ind_D0+IV]))
        }
        xvals = seq(0, xmax, length.out = 100)

        if (is.null(ymax)) {ymax = 1}

        plot(NULL, type = "n", ylab = "Settlement probability", xlab = "relative population density (N/K or bN)", xlim = c(0,xmax), ylim = c(0,ymax))
        leg.txt <- c()
        # Go through lines of distances matrix and add curves to plot
        for(line in 1:nrow(sett)){
            if (x@IndVar) {
                res <- matrix(ncol = 8, nrow = length(xvals))
                res[,1] <- densdep(xvals, A0 = sett[line,ind_D0]-sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]-sett[line,ind_D0+3], beta = sett[line,ind_D0+4]-sett[line,ind_D0+5])
                res[,2] <- densdep(xvals, A0 = sett[line,ind_D0]-sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]-sett[line,ind_D0+3], beta = sett[line,ind_D0+4]+sett[line,ind_D0+5])
                res[,3] <- densdep(xvals, A0 = sett[line,ind_D0]-sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]+sett[line,ind_D0+3], beta = sett[line,ind_D0+4]-sett[line,ind_D0+5])
                res[,4] <- densdep(xvals, A0 = sett[line,ind_D0]-sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]+sett[line,ind_D0+3], beta = sett[line,ind_D0+4]+sett[line,ind_D0+5])
                res[,5] <- densdep(xvals, A0 = sett[line,ind_D0]+sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]-sett[line,ind_D0+3], beta = sett[line,ind_D0+4]-sett[line,ind_D0+5])
                res[,6] <- densdep(xvals, A0 = sett[line,ind_D0]+sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]-sett[line,ind_D0+3], beta = sett[line,ind_D0+4]+sett[line,ind_D0+5])
                res[,7] <- densdep(xvals, A0 = sett[line,ind_D0]+sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]+sett[line,ind_D0+3], beta = sett[line,ind_D0+4]-sett[line,ind_D0+5])
                res[,8] <- densdep(xvals, A0 = sett[line,ind_D0]+sett[line,ind_D0+1], alpha = sett[line,ind_D0+2]+sett[line,ind_D0+3], beta = sett[line,ind_D0+4]+sett[line,ind_D0+5])
                polygon(c(xvals,rev(xvals)), c(apply(res, 1, min), rev(apply(res, 1, max))), border=NA, col='grey80')
            }
            lines(xvals, densdep(xvals, A0 = sett[line,ind_D0], alpha = sett[line,ind_D0+IV], beta = sett[line,ind_D0+2*IV]), type = "l", lty = 1, col = line)

            if (x@StageDep) {
                if (x@SexDep) {leg.txt <- c(leg.txt, paste0("Stage ",sett[line,1], ifelse(sett[line,2]," male"," female")))} else {leg.txt <- c(leg.txt, paste0("Stage ",sett[line,1]))}
            }
            else {
                if (x@SexDep) {leg.txt <- c(leg.txt, ifelse(sett[line,1],"male","female"))}
            }
        }
        if (length(leg.txt)>0) {
            legend("topright", leg.txt, col = 1:nrow(sett), lwd = 1.5)
        }
    }
    else{ print("Plotting is only implemented for density-dependent settlement (in a movement process).\n") }
})


#### CLASS DISPERSALPARAMS

# Superclass holding the sub-classes 'EmigrationParams', 'TransferParams' and 'SettlementParams'

#' Set Dispersal Parameters
#'
#' Dispersal is defined as movement leading to spatial gene flow. It typically involves three phases, which are all modelled explicitly:
#' \code{\link[RangeShiftR]{Emigration}}, \code{\link[RangeShiftR]{Transfer}} and \code{\link[RangeShiftR]{Settlement}}.\cr
#' The specific parameters of each phase are set through their respective functions. For more details, see their documentation.\cr
#' \cr
#' The potential evolution of several key dispersal traits is implemented through the possibility of inter-individual variability and heritability.
#' This option (called \code{IndVar}) can be enabled for each dispersal phase in their respective modules. See the Details below for information on how
#' to set the associated parameters. Additionally, the \code{\link[RangeShiftR]{Genetics}} module must be defined.
#'
#' @usage Dispersal(Emigration = Emigration(),
#'           Transfer   = DispersalKernel(),
#'           Settlement = Settlement())
#' @param Emigration The first phase of dispersal; determines if an individual leaves its natal patch.
#' @param Transfer (or transience) The second phase of dispersal; consists of the movement of an individual departing from its natal patch towards
#' a potential new patch, ending with settlement or mortality. This movement can be modelled by one of three alternative processes:\cr
#' - Dispersal kernel: use \code{\link[RangeShiftR]{DispersalKernel}}\cr
#' - Stochastic movement simulator (SMS): use \code{\link[RangeShiftR]{SMS}}\cr
#' - Correlated random walk (CRW): use \code{\link[RangeShiftR]{CorrRW}}
#' @param Settlement (or immigration) The last phase of dispersal; determines when the individual stops in a new cell or patch of
#' breeding habitat.
#' @details
#' Dispersal is defined as movement leading to spatial gene flow, and it typically involves three phases:
#' emigration, transfer and settlement
#'  \insertCite{stenseth1992,clobert2001,clobert2009,clobert2012,bowler2005,ronce2007}{RangeShiftR}.\cr
#'
#' The key role of dispersal in species persistence and responses to environmental change is increasingly recognized \insertCite{travis2013dispersal}{RangeShiftR}.
#' Moreover, the importance of modelling dispersal as a complex process, explicitly considering its three phases, each of which has
#' ´its own mechanisms and costs, has been recently highlighted \insertCite{bonte2012costs,travis2012modelling,travis2013dispersal}{RangeShiftR}.
#' The implementation of the dispersal process in \emph{RangeShiftR} is based on these recent frameworks and the substantial dispersal
#' theory that has been developed so far \insertCite{clobert2012}{RangeShiftR}.
#'
#' It is possible to model all three phases - emigration, transfer and settlement - with stage- and/or sex-specific parameters via the
#' switches \code{StageDep} and \code{SexDep} in the respective functions. In the case of sex-dependence, the number of traits is doubled, with
#' one set coding for the trait in females and the other for the trait in males. As well as being sex-biased, all dispersal phases can be
#' stage-biased, meaning that parameters can vary for different stage or age classes.
#'
#' The options of inter-individual variability in the various dispersal traits can be enabled by setting \code{IndVar=TRUE} for the respective dispersal phase
#' module and defining the genetic module to simulate heritability and evolution of traits (\code{\link[RangeShiftR]{Genetics}}).\cr
#' For each such heritable dispersal trait, this requires to set the mean and standard deviation of the distribution of trait values
#' in the initial population (instead of a constant value each, as in the case of \code{IndVar=FALSE}),
#' as well as the \code{TraitScaleFactor} for each heritable trait (see \code{\link[RangeShiftR]{Genetics}} documentaion).\cr
#'
#' Note, however, that not all combinations of sex-/stage-dependencies with inter-individual variability are implemented.\cr
#'
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "DispersalParams"
#' @author Anne-Kathleen Malchow
#' @name Dispersal
#' @export Dispersal
Dispersal <- setClass("DispersalParams", slots = c(Emigration = "EmigrationParams",
                                                   Transfer   = "TransferParams",
                                                   Settlement = "SettlementParams")
                      , prototype = list(Emigration = Emigration(),
                                         Transfer   = DispersalKernel(),
                                         Settlement = Settlement())
)
    ## add references to all sub-classes to this explanation of IndVar and TraitScaleFactor

setValidity("DispersalParams", function(object) {
    msg <- NULL
    validObject(object@Emigration)
    if (object@Emigration@UseFullKern) {
        if (!class(object@Transfer)[1] == "DispersalKernel") {
            msg <- c(msg, "Dispersal(): The emigration option \"UseFullKern\" can only be used if a dispersal kernel is used as transfer method!")
        }
    }
    validObject(object@Transfer)
    validObject(object@Settlement)
    if (class(object@Transfer)[1] == "DispersalKernel") {
        if (object@Settlement@DensDep) {
            msg <- c(msg, "Dispersal(): Settlement can only be density-dependent (DensDep = TRUE) if a movement process is used as transfer method!")
        }
    }else{
        if (object@Settlement@MaxSteps<=0 && all(object@Transfer@StepMort<=0) ) {
            msg <- c(msg, "Dispersal(): At least one of the two options MaxSteps and StepMort must be set (>0) if a movement process is used as transfer method!")
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "DispersalParams", function(.Object,...) {
    this_func = "Dispersal(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    .Object}
)
setMethod("show", "DispersalParams", function(object){
    cat(" Dispersal: \n  Emigration:\n")
    print(object@Emigration)
    cat("\n  Transfer:\n")
    print(object@Transfer)
    cat("\n  Settlement:\n")
    print(object@Settlement)
    if (!class(object@Transfer)[1] == "DispersalKernel") {
        cat("   MinSteps =", object@Settlement@MinSteps, "\n")
        cat("   MaxSteps =", object@Settlement@MaxSteps, "\n")
        cat("   MaxStepsYear =", object@Settlement@MaxStepsYear, "\n")
    }
})

# RS manual 2.5 (page 35) - 2.4.3 (page 34)
