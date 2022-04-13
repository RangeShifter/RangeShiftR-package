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



# -----
#
# Define Addition operations between modules
#
# -----


setGeneric("+")

setMethod("+", signature(e1 = "RSparams", e2 = "SimulationParams"), function(e1, e2) {
    validObject(e2)
    e1@simul <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "LandParams"), function(e1, e2) {
    validObject(e2)
    if (class(e2)[1] == "ImportedLandscape") {
        if (any(e2@PatchFile=="NULL")) {
            e1@control@patchmodel = FALSE
        }
        else {
            e1@control@patchmodel = TRUE
        }
        e1@control@resolution = e2@Resolution
        if (e2@HabPercent) {
            e1@control@landtype = 2L
            e1@control@maxNhab = 1L
        }
        else { # habitat codes
            e1@control@landtype = 0L
            e1@control@maxNhab = e2@Nhabitats
        }
        if (e2@SpDistFile=="NULL") {
            e1@control@speciesdist = FALSE
            e1@control@distresolution = -9L
        }
        else {
            e1@control@speciesdist = TRUE
            e1@control@distresolution = e2@SpDistResolution
        }
    }
    if (class(e2)[1] == "ArtificialLandscape") {
        e1@control@patchmodel = FALSE
        e1@control@resolution = e2@Resolution
        e1@control@landtype = 9L
        e1@control@maxNhab = 1L
        e1@control@speciesdist = FALSE
        e1@control@distresolution = 0L
    }
    e1@land = e2
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "DemogParams"), function(e1, e2) {
    validObject(e2)
    e1@control@reproductn = e2@ReproductionType
    if (class(e2@StageStruct)[1] == "StagesParams") {
        e1@control@repseasons = e2@StageStruct@RepSeasons
        e1@control@stagestruct = TRUE
        e1@control@stages = e2@StageStruct@Stages
        if (length(e2@StageStruct@MinAge)==1 && e2@StageStruct@MinAge==0) {
            if (e2@ReproductionType==2) {
                e2@StageStruct@MinAge <- rep(0L, 2*e2@StageStruct@Stages-1)
            }
            else {
                e2@StageStruct@MinAge <- rep(0L, e2@StageStruct@Stages)
            }
        }
    }
    else {
        e1@control@repseasons = 1L
        e1@control@stagestruct = FALSE
        e1@control@stages = 1L
    }
    e1@demog <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "DemogParams", e2 = "StagesParams"), function(e1, e2) {
    validObject(e2)
    e1@StageStruct <- e2
    e1@Rmax <- -9L
    e1@bc <- -9L
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "DispersalParams"), function(e1, e2) {
    validObject(e2)
    if (class(e2@Transfer)[1] == "DispersalKernel") {
        e1@control@transfer = 0
    }
    else {
        if (class(e2@Transfer)[1] == "StochMove") {
            e1@control@transfer = 1
        }
        else {
            if (class(e2@Transfer)[1] == "CorrRW") {
                e1@control@transfer = 2
            }
            else {
                e1@control@transfer = NA_integer_
            }
        }
    }
    e1@dispersal <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "DispersalParams", e2 = "EmigrationParams"), function(e1, e2) {
    validObject(e2)
    e1@Emigration <- e2
    validObject(e1)
    return(e1)}
)

setMethod("+", signature(e1 = "DispersalParams", e2 = "TransferParams"), function(e1, e2) {
    validObject(e2)
    e1@Transfer <- e2
    validObject(e1)
    return(e1)}
)

setMethod("+", signature(e1 = "DispersalParams", e2 = "SettlementParams"), function(e1, e2) {
    validObject(e2)
    e1@Settlement <- e2
    validObject(e1)
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "GeneticsParams"), function(e1, e2) {
    validObject(e2)
    e1@gene <- e2
    return(e1)}
)

setMethod("+", signature(e1 = "RSparams", e2 = "InitialisationParams"), function(e1, e2) {
    validObject(e2)
    e1@init <- e2
    return(e1)}
)
