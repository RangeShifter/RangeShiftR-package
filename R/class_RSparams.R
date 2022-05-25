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



#### PARAMETER MASTER CLASS ####

#' @include class_ControlParams.R
#' @include class_SimulationParams.R
#' @include class_LandParams.R
#' @include class_DemogParams.R
#' @include class_DispersalParams.R
#' @include class_GeneticsParams.R
#' @include class_InitialisationParams.R
RSparams <- setClass("RSparams", slots = c(control = "ControlParams",
                                           simul = "SimulationParams",
                                           land = "LandParams",
                                           demog = "DemogParams",
                                           dispersal = "DispersalParams",
                                           gene = "GeneticsParams",
                                           init = "InitialisationParams")
)
setValidity("RSparams", function(object) {
    msg <- NULL
    #CONTROL
    validObject(object@control)
    #SIMULATION
    validObject(object@simul)
    if (object@control@patchmodel) {
        if (object@simul@Gradient) {
            msg <- c(msg, "Environmental gradients are not implemented for patch-based models!")
        }
        if (object@simul@LocalExt) {
            msg <- c(msg, "Local extinction is not implemented for patch-based models!")
        }
        if (object@simul@EnvStoch==2) {
            msg <- c(msg, "Local environmental stochasticity (EnvStoch=2) is not implemented for patch-based models!")
        }
        if (object@simul@OutIntTraitRow) { # or rather both:  (object@simul@OutIntTraitCell || object@simul@OutIntTraitRow) ???
            msg <- c(msg, "Traits output by row is only applicable for a cell-based model!")
        }
    }
    else {
        if (object@simul@OutIntConn) {
            msg <- c(msg, " Connectivity output is only applicable for a patch-based model!")
        }
    }
    if (object@simul@EnvStochType==1) {
        if (object@control@landtype==0 || object@control@landtype==2) {
            msg <- c(msg, "Environmental stochasticity in carrying capacity (EnvStochType=1) is implemented for artificial landscapes only!")
        }
    }
    if (object@simul@OutIntTraitCell || object@simul@OutIntTraitRow) {
        if (!object@dispersal@Emigration@IndVar && !object@dispersal@Transfer@IndVar && !object@dispersal@Settlement@IndVar) {
            msg <- c(msg, "Traits output is only applicable for a model with inter-individual variability in (at least one) dispersal trait(s)!")
        }
    }
    if (class(object@dispersal@Transfer)[1] != "StochMove"){
        if (object@simul@OutIntPaths) {
            msg <- c(msg, "SMS paths output is only applicable for a model with the SMS transfer method!")
        }
        if (object@simul@SMSHeatMap){
            msg <- c(msg, "SMS heat map output is only applicable for a model with the SMS transfer method!")
        }
    }
    #LAND
    validObject(object@land)
    if (any(object@control@landtype==c(0,2))){
        if (any(object@land@DynamicLandYears>object@simul@Years)) {
            warning("ImportedLandscape(): Dynamic landscape contains years that exceed the simulated years, so that some land changes will not apply.", call. = FALSE)
        }
        if (object@land@CostsFile[1] !="NULL") {
            if (class(object@dispersal@Transfer)[1] == "StochMove") {
                if (object@dispersal@Transfer@Costs[1] != "file") {
                    warning("ImportedLandscape(): Landscape module contains SMS cost layers, but SMS module does not use them.", call. = FALSE)
                }
            }
            else{
                warning("ImportedLandscape(): Landscape module contains SMS cost layers, but Transfer module does not use SMS().", call. = FALSE)
            }
        }
    }
    #DEMOGRAPHY
    validObject(object@demog)
    #DISPERSAL
    validObject(object@dispersal)
    ## Emigration: check dimensions and values of EmigProb and EmigStage:
    # check StageDep and SexDep
    dim_ok = TRUE
    if (object@dispersal@Emigration@StageDep) {
        if (object@control@stagestruct) {
            rows = object@control@stages
            offset = 1
            if (object@dispersal@Emigration@IndVar) {
                if (is.na(object@dispersal@Emigration@EmigStage) || length(object@dispersal@Emigration@EmigStage)!=1) {
                    msg <- c(msg, "Emigration(): EmigStage (exactly 1) must be set!")
                }
                else {
                    if (object@dispersal@Emigration@EmigStage<0 || object@dispersal@Emigration@EmigStage>rows) {
                        msg <- c(msg, "Emigration(): EmigStage must be within [0, #stages]!")
                    }
                }
            }
        }
        else {
            dim_ok = FALSE
            msg <- c(msg, "Emigration can only be stage-dependent in a stage-structured model (StageStruct = StageStructure()) !")
        }
    }
    else {
        rows = 1
        offset = 0
    }
    if (object@dispersal@Emigration@SexDep) {
        if (object@control@reproductn) {
            rows = 2*rows
            offset = 1+offset
        }
        else {
            dim_ok = FALSE
            msg <- c(msg, "Emigration can only be sex-dependent in a sexually explicit model (ReproductionType = {1,2}) !")
        }
    }
    # check dimensions of EmigProb
    if (dim_ok) {
        if (object@dispersal@Emigration@DensDep) {
            cols = 3
        }
        else {
            cols = 1
        }
        if (object@dispersal@Emigration@IndVar) {
            cols = 2*cols
        }
        if (dim(object@dispersal@Emigration@EmigProb)[1]!=rows) {
            dim_ok = FALSE
            msg <- c(msg, paste0("Matrix of emigration probability traits (EmigProb) must have ", rows ," rows (with the current settings)!"))
        }
        if (dim(object@dispersal@Emigration@EmigProb)[2]!=(offset+cols)) {
            dim_ok = FALSE
            msg <- c(msg, paste0("Matrix of emigration probability traits (EmigProb) must have ", (offset+cols) ," columns (with the current settings)!"))
        }
    }
    if (dim_ok) {
        # check stage column of EmigProb
        if (object@dispersal@Emigration@StageDep) {
            if(any(object@dispersal@Emigration@EmigProb[,1]%%1!=0)){
                msg <- c(msg, "First column of emigration probability traits matrix (EmigProb) must contain the stage numbers (but non-integer number(s) found)!")
            }
            else {
                if(any(object@dispersal@Emigration@EmigProb[,1]<0 | object@dispersal@Emigration@EmigProb[,1]>object@control@stages)){
                    msg <- c(msg, "First column of emigration probability traits matrix (EmigProb) must contain the stage numbers (found <0 or >#stages)!")
                }
                else{
                    if(length(unique(object@dispersal@Emigration@EmigProb[,1])) != object@control@stages){
                        msg <- c(msg, "First column of emigration probability traits matrix (EmigProb) must contain the stage numbers (but found incorrect stage numbers)!")
                    }
                }
            }
        }
        # check sex column of EmigProb
        if (object@dispersal@Emigration@SexDep) {
            if(any(object@dispersal@Emigration@EmigProb[,offset]!=0 && object@dispersal@Emigration@EmigProb[,offset]!=1)){
                msg <- c(msg, paste0(offset,". column of emigration probability traits matrix (EmigProb) must contain the sex numbers (0 for female, 1 for male)!"))
            }
            else {
                if (object@dispersal@Emigration@StageDep) {
                    comb_ok = TRUE
                    for(i in 0:(object@control@stages-1)) {
                        if (length(unique(object@dispersal@Emigration@EmigProb[object@dispersal@Emigration@EmigProb[,1]==i, offset])) != 2) {
                            comb_ok = FALSE
                        }
                    }
                    if (!comb_ok) {
                        dim_ok = FALSE
                        msg <- c(msg, "The emigration probability traits matrix (EmigProb) must contain exactly one row for each stage (1. column) and sex (2. column) combination!")
                    }
                }
                else {
                    if (length(unique(object@dispersal@Emigration@EmigProb[, offset])) != 2) {
                        dim_ok = FALSE
                        msg <- c(msg, "The emigration probability traits matrix (EmigProb) must contain exactly one row for each sex (1. column)!")
                    }
                }
            }
        }
    }
    if (dim_ok) {
        # check value columns of EmigProb
        if (object@dispersal@Emigration@IndVar) {
            if(any(object@dispersal@Emigration@EmigProb[,(offset+1):(offset+2)]<=0 | object@dispersal@Emigration@EmigProb[,(offset+1):(offset+2)]>1)){
                msg <- c(msg, paste0("Columns ", (offset+1), " and ", (offset+2), " of emigration traits matrix (EmigProb) must contain mean(D0) and sd(D0), with values in the half-open inerval (0,1] !"))
            }
            else {
                if(any(object@dispersal@Emigration@EmigProb[,(offset+2)] > object@dispersal@Emigration@TraitScaleFactor[1])) {
                    msg <- c(msg, paste0("Column ", (offset+2), " of emigration traits matrix (EmigProb) must contain sd(D0), with values less than or equal to TraitScaleFactor μ(D0)!"))
                }
            }
            if (object@dispersal@Emigration@DensDep) {
                if(any(object@dispersal@Emigration@EmigProb[,c((offset+4),(offset+6))]<=0 )){
                    msg <- c(msg, paste0("Columns ", (offset+4), " and ", (offset+6), " of emigration traits matrix (EmigProb) must contain sd(α) and sd(β), with strictly positive values!"))
                }
                else {
                    if(any(object@dispersal@Emigration@EmigProb[,(offset+4)] > object@dispersal@Emigration@TraitScaleFactor[2])) {
                        msg <- c(msg, paste0("Column ", (offset+4), " of emigration traits matrix (EmigProb) must contain sd(α), with values less than or equal to TraitScaleFactor μ(α)!"))
                    }
                    if(any(object@dispersal@Emigration@EmigProb[,(offset+6)] > object@dispersal@Emigration@TraitScaleFactor[3])) {
                        msg <- c(msg, paste0("Column ", (offset+6), " of emigration traits matrix (EmigProb) must contain sd(β), with values less than or equal to TraitScaleFactor μ(β)!"))
                    }
                }
            }
        }
        else { # !IndVar
            if (object@dispersal@Emigration@DensDep) {
                if(any(object@dispersal@Emigration@EmigProb[,(offset+1)]<0 | object@dispersal@Emigration@EmigProb[,(offset+1)]>1)){
                    msg <- c(msg, paste0("Column ", (offset+1), " of emigration traits matrix (EmigProb) must contain the maximum emigration probability D0, with values in the closed inerval [0,1] !"))
                }
                #// NB alpha and beta may take any value
            }
            else {
                if(any(object@dispersal@Emigration@EmigProb[,(offset+1)]<0 | object@dispersal@Emigration@EmigProb[,(offset+1)]>1)){
                    msg <- c(msg, paste0("Column ", (offset+1), " of emigration traits matrix (EmigProb) must contain the constant emigration probability d, with values in the closed inerval [0,1] !"))
                }
            }
        }
    }
    ## Transfer / DispersalKernel:
    # check StageDep and SexDep
    dim_ok = TRUE
    if (class(object@dispersal@Transfer)[1] == "DispersalKernel") {
        if (object@control@transfer != 0) {
            msg <- c(msg, "DispersalKernel(): Something went wrong adding a DispersalKernel to the Parameter Master!")
        }
        else {
            if (object@dispersal@Transfer@StageDep) {
                if (object@control@stagestruct) {
                    rows = object@control@stages
                    offset = 1
                }
                else {
                    dim_ok = FALSE
                    msg <- c(msg, "Transfer (DispersalKernel) can only be stage-dependent in a stage-structured model (use Demography(StageStruct = StageStructure()) ) !")
                }
            }
            else {
                rows = 1
                offset = 0
            }
            if (object@dispersal@Transfer@SexDep) {
                if (object@control@reproductn) {
                    rows = 2*rows
                    offset = 1+offset
                }
                else {
                    dim_ok = FALSE
                    msg <- c(msg, "Transfer (DispersalKernel) can only be sex-dependent in a sexually explicit model (ReproductionType = {1,2}) !")
                }
            }
            if (dim_ok) {
                if (object@dispersal@Transfer@DoubleKernel) {
                    cols = 3
                }
                else {
                    cols = 1
                }
                if (object@dispersal@Transfer@IndVar) {
                    cols = 2*cols
                }
                if (dim(object@dispersal@Transfer@Distances)[1]!=rows) {
                    dim_ok = FALSE
                    msg <- c(msg, paste0("Matrix of dispersal kernel traits (Distances) must have ", rows ," rows (with the current settings)!"))
                }
                if (dim(object@dispersal@Transfer@Distances)[2]!=(offset+cols)) {
                    dim_ok = FALSE
                    msg <- c(msg, paste0("Matrix of dispersal kernel traits (Distances) must have ", (offset+cols) ," columns (with the current settings)!"))
                }
            }
            if (dim_ok) {
                # check stage column of Distances matrix
                if (object@dispersal@Transfer@StageDep) {
                    if(any(object@dispersal@Transfer@Distances[,1]%%1!=0)){
                        msg <- c(msg, "First column of dispersal kernel traits (Distances) matrix must contain the stage numbers (but non-integer number(s) found)!")
                    }
                    else {
                        if(any(object@dispersal@Transfer@Distances[,1]<0 | object@dispersal@Transfer@Distances[,1]>object@control@stages)){
                            msg <- c(msg, "First column of dispersal kernel traits (Distances) matrix must contain the stage numbers (found <0 or >#stages)!")
                        }
                        else{
                            if(length(unique(object@dispersal@Transfer@Distances[,1])) != object@control@stages){
                                msg <- c(msg, "First column of dispersal kernel traits (Distances) matrix must contain the stage numbers (but found incorrect stage numbers)!")
                            }
                        }
                    }
                }
                if (object@dispersal@Transfer@SexDep) {
                    # check sex column of Distances matrix
                    if(any(object@dispersal@Transfer@Distances[,offset]!=0 && object@dispersal@Transfer@Distances[,offset]!=1)){
                        msg <- c(msg, paste0(offset,". column of dispersal kernel traits (Distances) matrix must contain the sex numbers (0 for female, 1 for male)!"))
                    }
                    else {
                        if (object@dispersal@Transfer@StageDep) {
                            comb_ok = TRUE
                            for(i in 0:(object@control@stages-1)) {
                                if (length(unique(object@dispersal@Transfer@Distances[object@dispersal@Transfer@Distances[,1]==i, offset])) != 2) {
                                    comb_ok = FALSE
                                }
                            }
                            if (!comb_ok) {
                                dim_ok = FALSE
                                msg <- c(msg, "The dispersal kernel traits (Distances) matrix must contain exactly one row for each stage (1. column) and sex (2. column) combination!")
                            }
                        }
                        else {
                            if (length(unique(object@dispersal@Transfer@Distances[, offset])) != 2) {
                                dim_ok = FALSE
                                msg <- c(msg, "The emigration dispersal kernel traits (Distances) matrix must contain exactly one row for each sex (1. column) !")
                            }
                        }
                    }
                }
            }
            if (dim_ok) {
                # check value columns of Distances matrix
                if (object@dispersal@Emigration@UseFullKern) {
                    resol = 0.0000000000001
                }
                else {
                    resol = object@control@resolution
                }
                if (object@dispersal@Transfer@IndVar) {
                    if (object@dispersal@Transfer@DoubleKernel) {
                        if(any(object@dispersal@Transfer@Distances[,(offset+1)]<resol)){
                            msg <- c(msg, paste0("Column ", (offset+1), " of dispersal kernel traits (Distances) matrix must contain the mean of the mean distance of Kernel 1, mean(δ1), with values >= ", (resol), " (=landscape resolution (unless UseFullKernel=TRUE)) !"))
                        }
                        if(any(object@dispersal@Transfer@Distances[,(offset+2)]<=0.0)){
                            msg <- c(msg, paste0("Column ", (offset+2), " of dispersal kernel traits (Distances) matrix must contain the std. dev. of the mean distance of Kernel 1, sd(δ1), with strictly positive values!"))
                        }
                        else {
                            if(any(object@dispersal@Transfer@Distances[,(offset+2)] > object@dispersal@Transfer@TraitScaleFactor[1])) {
                                msg <- c(msg, paste0("Column ", (offset+2), " of dispersal kernel traits (Distances) matrix must contain sd(δ1), with values less than or equal to TraitScaleFactor μ(δ1)!"))
                            }
                        }
                        if(any(object@dispersal@Transfer@Distances[,(offset+3)]<resol)){
                            msg <- c(msg, paste0("Column ", (offset+3), " of dispersal kernel traits (Distances) matrix must contain the mean of the mean distance of Kernel 2, mean(δ2), with values >= ", (resol), " (=landscape resolution (unless UseFullKernel=TRUE)) !"))
                        }
                        if(any(object@dispersal@Transfer@Distances[,(offset+4)]<=0.0)){
                            msg <- c(msg, paste0("Column ", (offset+4), " of dispersal kernel traits (Distances) matrix must contain the std. dev. of the mean distance of Kernel 2, sd(δ2), with strictly positive values!"))
                        }
                        else {
                            if(any(object@dispersal@Transfer@Distances[,(offset+4)] > object@dispersal@Transfer@TraitScaleFactor[2])) {
                                msg <- c(msg, paste0("Column ", (offset+4), " of dispersal kernel traits (Distances) matrix must contain sd(δ2), with values less than or equal to TraitScaleFactor μ(δ2)!"))
                            }
                        }
                        if(any(object@dispersal@Transfer@Distances[,(offset+5)]<=0.0 | object@dispersal@Transfer@Distances[,(offset+5)]>=1.0)){
                            msg <- c(msg, paste0("Column ", (offset+5), " of dispersal kernel traits (Distances) matrix must contain the mean of the probability of using Kernel 1, mean(p), with values in the open interval (0,1) !"))
                        }
                        if(any(object@dispersal@Transfer@Distances[,(offset+6)]<=0.0)){
                            msg <- c(msg, paste0("Column ", (offset+6), " of dispersal kernel traits (Distances) matrix must contain the std. dev. of the probability of using Kernel 1, sd(p), with strictly positive values!"))
                        }
                        else {
                            if(any(object@dispersal@Transfer@Distances[,(offset+6)] > object@dispersal@Transfer@TraitScaleFactor[3])) {
                                msg <- c(msg, paste0("Column ", (offset+6), " of dispersal kernel traits (Distances) matrix must contain sd(p), with values less than or equal to TraitScaleFactor μ(p)!"))
                            }
                        }
                    }
                    else { # !DoubleKernel
                        if(any(object@dispersal@Transfer@Distances[,(offset+1)]<resol)){
                            msg <- c(msg, paste0("Column ", (offset+1), " of dispersal kernel traits (Distances) matrix must contain the mean of the mean distance, mean(δ), with values >= ", (resol), " (=landscape resolution (unless UseFullKernel=TRUE)) !"))
                        }
                        if(any(object@dispersal@Transfer@Distances[,(offset+2)]<=0.0)){
                            msg <- c(msg, paste0("Column ", (offset+2), " of dispersal kernel traits (Distances) matrix must contain the std. dev. of the mean distance, sd(δ), with strictly positive values!"))
                        }
                        else {
                            if(any(object@dispersal@Transfer@Distances[,(offset+2)] > object@dispersal@Transfer@TraitScaleFactor[1])) {
                                msg <- c(msg, paste0("Column ", (offset+2), " of dispersal kernel traits (Distances) matrix must contain sd(δ), with values less than or equal to TraitScaleFactor μ(δ)!"))
                            }
                        }
                    }
                }
                else { # !IndVar
                    if (object@dispersal@Transfer@DoubleKernel) {
                        if(any(object@dispersal@Transfer@Distances[,(offset+1):(offset+2)]<resol)){
                            msg <- c(msg, paste0("Columns ", (offset+1), " and ", (offset+2), " of dispersal kernel traits (Distances) matrix must contain the two mean distance δ1 and δ2, with values >= ", (resol), " (=landscape resolution (unless UseFullKernel=TRUE)) !"))
                        }
                        if(any(object@dispersal@Transfer@Distances[,(offset+3)]<=0 | object@dispersal@Transfer@Distances[,(offset+3)]>=1)) {
                            msg <- c(msg, paste0("Column ", (offset+3), " of dispersal kernel traits (Distances) matrix must contain the probability p of using kernel 1, with values in the open interval (0,1) !"))
                        }
                    }
                    else{
                        if(any(object@dispersal@Transfer@Distances[,(offset+1)]<resol)){
                            msg <- c(msg, paste0("Column ", (offset+1), " of dispersal kernel traits (Distances) matrix must contain the mean distance δ, with values >= ", (resol), " (=landscape resolution (unless UseFullKernel=TRUE)) !"))
                        }
                    }
                }
            }
        }
    }
    ## Transfer / StochMove:
    if (class(object@dispersal@Transfer)[1] == "StochMove") {
        if (object@control@transfer != 1) {
            msg <- c(msg, "SMS(): Something went wrong adding a SMS Transfer method to the Parameter Master!")
        }
        else {
            # check transfer cost parameters
            if (object@control@landtype==0) { # imported land with habitat codes
                if (!any(length(object@dispersal@Transfer@StepMort) == c(1,object@control@maxNhab) )) {
                    msg <- c(msg, "SMS(): Per-step mortality probability must have either 1 or Nhabitat entries for an imported landscape with habitat codes!")
                }
                if (class(object@dispersal@Transfer@Costs)=="numeric") {
                    if (length(object@dispersal@Transfer@Costs) != object@control@maxNhab) {
                        msg <- c(msg, "SMS(): Costs must have Nhabitat entries for an imported landscape with habitat codes!")
                    }
                }
                if (class(object@dispersal@Transfer@Costs)=="character") {
                    if (object@dispersal@Transfer@Costs == "file") {
                        if (object@land@CostsFile[1] == "NULL") {
                            msg <- c(msg, "SMS(): No cost map filenames found in the landscape module!")
                        }
                    }
                    else{
                        msg <- c(msg, "SMS(): Costs has a wrong format! Must be either numeric or the keyword \"file\".")
                    }
                }
            }
            else {
                if (object@control@landtype==2) { # imported land with habitat quality
                    if (length(object@dispersal@Transfer@StepMort)!=1) {
                        msg <- c(msg, "SMS(): Per-step mortality probability must be a constant for an imported habitat percentage landscape!")
                    }
                    if (class(object@dispersal@Transfer@Costs)=="character") {
                        if (object@dispersal@Transfer@Costs == "file") {
                            if (object@land@CostsFile[1] == "NULL") {
                                msg <- c(msg, "SMS(): No cost map filenames found in the landscape module!")
                            }
                        }
                        else{
                            msg <- c(msg, "SMS(): Costs has a wrong format! Must be either numeric or the keyword \"file\".")
                        }
                    }
                    else{
                        msg <- c(msg, "SMS(): Costs must be imported from a raster map for an imported habitat percentage landscape!")
                    }
                }
                else {
                    if (object@control@landtype==9) { # artificial landscape
                        if (!any(length(object@dispersal@Transfer@StepMort) == c(1,2)) ) {
                            msg <- c(msg, "SMS(): Per-step mortality probability must have one or two (i.e. for matrix and habitat) entries for an artificial landscape!")
                        }
                        if (class(object@dispersal@Transfer@Costs)=="numeric") {
                            if (length(object@dispersal@Transfer@Costs) != 2) {
                                msg <- c(msg, "SMS(): Costs must have two entries (i.e. for matrix and habitat) for an artificial landscape!")
                            }
                            else {
                                if (class(object@dispersal@Transfer@Costs)=="character") {
                                    msg <- c(msg, "SMS(): Costs can not be imported from a raster map for an artificial landscape!")
                                }
                            }
                        }
                    }
                    else { # no valid land code
                        msg <- c(msg, "SMS(): Something went wrong adding a SMS Transfer method to the Parameter Master: No valid landscape code.")
                    }
                }
            }
            if(length(object@dispersal@Transfer@CostMap)==0) {
                msg <- c(msg, "SMS(): Something went wrong adding a SMS Transfer method to the Parameter Master: No set CostMap switch.")
            }
            else{
                if(object@dispersal@Transfer@CostMap){
                    if(object@dispersal@Transfer@Costs[1] != "file") {
                        msg <- c(msg, "SMS(): Something went wrong adding a SMS Transfer method to the Parameter Master: CostMap switch incorrect.")
                    }
                }
                else{
                    if(object@dispersal@Transfer@Costs[1] == "file") {
                        msg <- c(msg, "SMS(): Something went wrong adding a SMS Transfer method to the Parameter Master: CostMap switch incorrect.")
                    }

                }
            }
        }
    }
    ## Transfer / CorrRW:
    if (class(object@dispersal@Transfer)[1] == "CorrRW") {
        if (object@control@transfer != 2) {
            msg <- c(msg, "CorrRW(): Something went wrong adding a Correlated RW Transfer method to the Parameter Master!")
        }
    }

    ## Settlement:
    # check StageDep and SexDep
    dim_ok = TRUE
    if (object@dispersal@Settlement@StageDep) {
        if (object@control@stagestruct) {
            rows = object@control@stages
            offset = 1
        }
        else {
            dim_ok = FALSE
            msg <- c(msg, "Settlement can only be stage-dependent in a stage-structured model (StageStruct = StageStructure()) !")
        }
    }
    else {
        rows = 1
        offset = 0
    }
    if (object@dispersal@Settlement@SexDep) {
        if (object@control@reproductn) {
            rows = 2*rows
            offset = 1+offset
        }
        else {
            dim_ok = FALSE
            msg <- c(msg, "Settlement can only be sex-dependent in a sexually explicit model (ReproductionType = {1,2}) !")
            if (any(FindMate)){
                msg <- c(msg, "FindMate can only be TRUE in a sexually explicit model (ReproductionType = {1,2}) !")
            }
        }
    }
    #check dimensions and values of matrix Settle:
    if (dim_ok) {
        if (object@dispersal@Settlement@DensDep) {
            cols = 3
        }
        else {
            if (object@control@transfer) cols = 0
            else  cols = 1
        }
        if (object@dispersal@Settlement@IndVar) {
            cols = 2*cols
        }
        if (dim(object@dispersal@Settlement@Settle)[1]!=rows) {
            dim_ok = FALSE
            msg <- c(msg, paste0("Matrix of Settlement traits (Settle) must have ", rows ," rows (with the current settings)!"))
        }
        if((offset+cols) == 0){
            if (any(dim(object@dispersal@Settlement@Settle)!=c(1,1) )) {
                msg <- c(msg, paste0("Matrix of Settlement traits is not needed with the current settings!"))
            }
        }
        else{
            if (dim(object@dispersal@Settlement@Settle)[2]!=(offset+cols)) {
                dim_ok = FALSE
                msg <- c(msg, paste0("Matrix of Settlement traits (Settle) must have ", (offset+cols) ," columns (with the current settings)!"))
            }
        }
    }
    if (dim_ok) {
        # check stage column of Settle
        if (object@dispersal@Settlement@StageDep) {
            if(any(object@dispersal@Settlement@Settle[,1]%%1!=0)){
                msg <- c(msg, "First column of Settlement traits matrix (Settle) must contain the stage numbers (but non-integer number(s) found)!")
            }
            else {
                if(any(object@dispersal@Settlement@Settle[,1]<0 | object@dispersal@Settlement@Settle[,1]>object@control@stages)){
                    msg <- c(msg, "First column of Settlement traits matrix (Settle) must contain the stage numbers (found <0 or >#stages)!")
                }
                else{
                    if(length(unique(object@dispersal@Settlement@Settle[,1])) != object@control@stages){
                        msg <- c(msg, "First column of Settlement traits matrix (Settle) must contain the stage numbers (but found incorrect stage numbers)!")
                    }
                }
            }
        }
        # check sex column of Settle
        if (object@dispersal@Settlement@SexDep) {
            if(!all(object@dispersal@Settlement@Settle[,offset] %in% c(0,1))){
                msg <- c(msg, paste0(offset,". column of Settlement traits matrix (Settle) must contain the sex numbers (0 for female, 1 for male)!"))
            }
            else {
                if (object@dispersal@Settlement@StageDep) {
                    comb_ok = TRUE
                    for(i in 0:(object@control@stages-1)) {
                        if (length(unique(object@dispersal@Settlement@Settle[object@dispersal@Settlement@Settle[,1]==i, offset])) != 2) {
                            comb_ok = FALSE
                        }
                    }
                    if (!comb_ok) {
                        dim_ok = FALSE
                        msg <- c(msg, "The Settlement traits matrix (Settle) must contain exactly one row for each stage (1. column) and sex (2. column) combination!")
                    }
                }
                else {
                    if (length(unique(object@dispersal@Settlement@Settle[, offset])) != 2) {
                        dim_ok = FALSE
                        msg <- c(msg, "The Settlement traits matrix (Settle) must contain exactly one row for each sex (1. column)!")
                    }
                }
            }
        }
    }
    if (dim_ok) {
        # check value columns of Settle
        if (object@control@transfer) {          #if Movemodel:
            if (object@dispersal@Settlement@IndVar) {
                if (object@dispersal@Settlement@DensDep) {
                    if(any(object@dispersal@Settlement@Settle[,(offset+1):(offset+2)]<=0 | object@dispersal@Settlement@Settle[,(offset+1):(offset+2)]>1)){
                        msg <- c(msg, paste0("Columns ", (offset+1), " and ", (offset+2), " of settlement traits matrix (Settle) must contain mean(S0) and sd(S0), with values in the half-open inerval (0,1] !"))
                    }
                    else {
                        if(any(object@dispersal@Settlement@Settle[,(offset+2)] > object@dispersal@Settlement@TraitScaleFactor[1])) {
                            msg <- c(msg, paste0("Column ", (offset+2), " of settlement traits matrix (Settle) must contain sd(S0), with values less than or equal to TraitScaleFactor μ(S0)!"))
                        }
                    }
                    if(any(object@dispersal@Settlement@Settle[,c((offset+4),(offset+6))]<=0 )){
                        msg <- c(msg, paste0("Columns ", (offset+4), " and ", (offset+6), " of settlement traits matrix (Settle) must contain sd(α_s) and sd(β_s), with strictly positive values!"))
                    }
                    else {
                        if(any(object@dispersal@Settlement@Settle[,(offset+4)] > object@dispersal@Settlement@TraitScaleFactor[2])) {
                            msg <- c(msg, paste0("Column ", (offset+4), " of settlement traits matrix (Settle) must contain sd(α_s), with values less than or equal to TraitScaleFactor μ(α_s)!"))
                        }
                        if(any(object@dispersal@Settlement@Settle[,(offset+6)] > object@dispersal@Settlement@TraitScaleFactor[3])) {
                            msg <- c(msg, paste0("Column ", (offset+6), " of settlement traits matrix (Settle) must contain sd(β_s), with values less than or equal to TraitScaleFactor μ(β_s)!"))
                        }
                    }
                }
                else {
                    msg <- c(msg, paste0("Settlement(): Inter-individual variability (IndVar=TRUE) in settlement traits requires density-dependence (DensDep=TRUE) !"))
                }
            }
            else { # !IndVar
                if (object@dispersal@Settlement@DensDep) {
                    if(any(object@dispersal@Settlement@Settle[,(offset+1)]<=0 | object@dispersal@Settlement@Settle[,(offset+1)]>1)){
                        msg <- c(msg, paste0("Column ", (offset+1), " of settlement traits matrix (Settle) must contain the maximum settlement probability S0, with values in the half-open inerval (0,1] !"))
                    }
                    #// NB alpha and beta may take any value
                }
                #else {}  // no more columns required other than stage and sex if applicable
            }
        }
        else {      # DispersalKernel
            if (object@control@stagestruct) {
                if(!all(object@dispersal@Settlement@Settle[,(offset+1)] %in% c(0,1,2,3))){
                    msg <- c(msg, paste0("Column ", (offset+1), " of settlement traits matrix (Settle) must contain the settlement condition codes, with valid values: 0, 1, 2 or 3."))
                }
            }
            else{
                if(!all(object@dispersal@Settlement@Settle[,(offset+1)] %in% c(0,2))){
                    msg <- c(msg, paste0("Column ", (offset+1), " of settlement traits matrix (Settle) must contain the settlement condition codes; for a non-StageStructured population the only valid values are 0 and 2."))
                }
            }
        }
        #if Movemodel, check its additional parameters:
        if (object@control@transfer) {
            if (!length(object@dispersal@Settlement@MinSteps) %in% c(1,rows) ){
                msg <- c(msg, paste0("Settlement(): MinSteps must have either 1 or ", rows ," rows (with the current settings)!"))
            }
            if (!length(object@dispersal@Settlement@MaxSteps) %in% c(1,rows) ){
                msg <- c(msg, paste0("Settlement(): MaxSteps must have either 1 or ", rows ," rows (with the current settings)!"))
            }
            if (object@dispersal@Settlement@StageDep) {
                if (!length(object@dispersal@Settlement@MaxStepsYear) %in% c(1,rows) ){
                    msg <- c(msg, paste0("Settlement(): MaxStepsYear must have either 1 or ", rows ," rows (with the current settings)!"))
                }
            }
        }
    }

    #GENETICS
    validObject(object@gene)
    if(any(object@dispersal@Emigration@IndVar,object@dispersal@Transfer@IndVar,object@dispersal@Settlement@IndVar)) anyIndVar <- TRUE
    else anyIndVar <- FALSE
    # if( ...(object@gene) & !anyIndVar) # TODO: check if genetics module is set but no variable traits -> give warning in this case

    #INITIALISATION
    validObject(object@init)
    if (object@control@landtype == 9) { # artificial land
        if (object@init@InitType) {
            msg <- c(msg, 'Initialise(): InitType must be set to 0 (Free initialisation) for an artificial landscape!')
        }
        else {
            if (object@init@minX > object@land@dimX) {
                msg <- c(msg, 'Initialise(): minX exceeds the landscape dimensions!')
            }
            else {
                if (object@init@maxX > object@land@dimX) {
                    msg <- c(msg, 'Initialise(): maxX exceeds the landscape dimensions!')
                }
            }
            if (object@init@minY > object@land@dimY) {
                msg <- c(msg, 'Initialise(): minY exceeds the landscape dimensions!')
            }
            else {
                if (object@init@maxY > object@land@dimY) {
                    msg <- c(msg, 'Initialise(): maxY exceeds the landscape dimensions!')
                }
            }
        }
    }
    else {  # imported land
        if (object@init@InitType == 1 && !object@control@speciesdist) {
            msg <- c(msg, 'Initialise(): A species distribution map has to be loaded via the \'land\' module if InitType = 1 (initialisation from loaded species distribution map) !')
        }
    }
    if (object@control@stagestruct) {
        if (is.na(object@init@InitAge) || length(object@init@InitAge) == 0){
            msg <- c(msg, 'Initialise(): InitAge must be set!')
        }
        else {
            if (!object@init@InitAge %in% c(0,1,2) && object@init@InitType != 2){
                msg <- c(msg, 'Initialise(): InitAge must be 0, 1 or 2!')
            }
        }
        if (object@init@InitType != 2) {
            if ( is.na(object@init@PropStages) || length(object@init@PropStages) == 0 ){
                msg <- c(msg, 'Initialise(): PropStages must be set!')
            }
            else {
                if(length(object@init@PropStages) != object@control@stages) {
                    msg <- c(msg, 'Initialise(): PropStages must have a length equal to number of stages!')
                }
                else{
                    if(any(object@init@PropStages < 0.0) | any(object@init@PropStages > 1.0)) {
                        msg <- c(msg, 'Initialise(): All elements of PropStages must be in the closed interval [0,1]!')
                    }
                    else{
                        if (object@init@PropStages[1] != 0.0) {
                            msg <- c(msg, 'Initialise(): Initial proportion of the juvenile stage (PropStages[1]) must be 0.0!')
                        }
                        else{
                            if (length(object@init@PropStages)>1 && sum(object@init@PropStages) != 1.0) {
                                msg <- c(msg, 'Initialise(): The elements of PropStages must sum to 1!')
                            }
                        }
                    }
                }
            }
        }
    }
    else{
        if (length(object@init@PropStages)>1 || !(object@init@PropStages[1] %in% c(0,-9))) {
            msg <- c(msg, 'PropStages is not used for a population without stage structure.')
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("show", "RSparams", function(object){
    print(object@control)
    cat("\n")
    print(object@simul)
    cat("\n")
    print(object@land)
    cat("\n")
    print(object@demog)
    cat("\n")
    print(object@dispersal)
    cat("\n")
    if(any(object@dispersal@Emigration@IndVar,object@dispersal@Transfer@IndVar,object@dispersal@Settlement@IndVar)
       || object@gene@Architecture == 1 ){
        print(object@gene)
        cat("\n")
    }
    print(object@init)}
)
