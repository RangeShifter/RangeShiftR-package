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
 
 

#### CLASS CONTROLPARAMS

# from RS 'CONTROL.txt' file

ControlParams <- setClass("ControlParams", slots = c(
                                   #nSimuls = "integer_OR_numeric",         # Not yet used by R version, in C++ version its read from parameterfile
                                   #nLandscapes = "integer_OR_numeric",     # Not yet used by R version, in C++ version its read from landfile
                                   batchnum = "integer_OR_numeric",         # only variable to set from RSsim(), optional
                                   patchmodel = "logical",                  # set via +Land
                                   resolution = "integer_OR_numeric",       # set via +Land
                                   landtype = "integer_OR_numeric",         # set via +Land
                                   maxNhab = "integer_OR_numeric",          # set via +Land
                                   speciesdist = "logical",                 # set via +Land
                                   distresolution = "integer_OR_numeric",   # set via +Land
                                   reproductn = "integer_OR_numeric",       # set via +Demography
                                   repseasons = "integer_OR_numeric",       # set via +Demography if (StageStruct is Type "StagesParams") {repseasons = 1} else {demo@StageStruct@repseasons}
                                   stagestruct = "logical",                 # set via +Demography
                                   stages = "integer_OR_numeric",           # set via +Demography@StageStruct
                                   transfer = "integer_OR_numeric",         # set via +Dispersal     Transfer method: 0 = dispersal kernels, 1 = SMS, 2 = CRW)
                                   seed = "integer_OR_numeric")
                         ,prototype = list(#nSimuls = 1L,
                                  #nLandscapes = 1L,
                                  batchnum = 0L,
                                  patchmodel = FALSE,
                                  resolution = 100L,
                                  landtype = 9L,
                                  maxNhab = 1L,
                                  speciesdist = FALSE,
                                  distresolution = NA_integer_,
                                  reproductn = 0L,
                                  repseasons = 1L,
                                  stagestruct = FALSE,
                                  stages = NA_integer_,
                                  transfer = 0L,
                                  seed = 0L)
)
setValidity("ControlParams", function(object){
    msg <- NULL
    if(object@stagestruct) {
        stg <- as.integer(object@stages)
        if(stg < 2 || stg > 10) msg <- c(msg, "Number of Stages must be in the interval [2; 10]")
                                                          #paste("Unequal x,y lengths: ", length(object@x), ", ", length(object@y), sep="")
    }
    if(object@seed > 4000000000) {
        warning("RSsim(): Seeds greater than 4e9 should be avoided, since they might not be handled correctly.", call. = FALSE)
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("show", "ControlParams", function(object){
    cat(" Batch #", object@batchnum, "\n")
    if(object@seed){
        cat(" Seed =", object@seed)
        if(object@seed<0) cat("  (generate random seed)")
        if(object@seed>0) cat("  (fixed seed)")
        cat("\n")
    }
})
