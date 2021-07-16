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



#-------------------
# Output handling and plotting functions
#-------------------


#---------------------------------------------------------

### READING OUTPUT FILES


#' Read 'range' file
#'
#' Read the RangeShiftR output file 'range' into a data.frame, if it was generated.
#' @param s RSmaster parameter object
#' @param dirpath RS directory path
#' @return a data.frame
#' @export
setGeneric("readRange", function(s,dirpath) standardGeneric("readRange") )

setMethod("readRange", c(s="RSparams", dirpath="character"), function(s,dirpath) {
    path <- paste0(dirpath, "Outputs/Batch", s@control@batchnum, "_Sim", s@simul@Simulation, "_Land", s@land@LandNum, "_Range.txt")
    if(file.exists(path)){
        return(read.table(path, h = T, sep = "\t"))
    }else {
        warning("The 'range' output file for this simulation does not exist.", call. = FALSE)
    }
})


#' Read 'pop' file
#'
#' Read the RangeShiftR output file 'pop' into a data.frame, if it was generated.
#' @param s RSmaster parameter object
#' @param dirpath RS directory path
#' @param center In a cell-based model, the x- and y- coordinates will be converted from *nrow* and *ncol* (as in the 'pop' output file)
#' to true coordinates of the center of the cell measured in metres.
#' @return a data.frame
#' @export
setGeneric("readPop", function(s,dirpath,...) standardGeneric("readPop") )

setMethod("readPop", c(s="RSparams", dirpath="character"), function(s,dirpath,center=TRUE) {
    path <- paste0(dirpath, "Outputs/Batch", s@control@batchnum, "_Sim", s@simul@Simulation, "_Land", s@land@LandNum, "_Pop.txt")
    if(file.exists(path)){
        pop_table <- read.table(path, h = T, sep = "\t")
        if (sum(grepl('x',names(pop_table))) & center){
            pop_table$x <- (pop_table$x+0.5)*s@land@Resolution
            pop_table$y <- (pop_table$y+0.5)*s@land@Resolution
        }
        return(pop_table)
    }else {
        warning("The 'pop' output file for this simulation does not exist.", call. = FALSE)
    }
})


#---------------------------------------------------------

### PROCESS OUTPUT


#' ColonisationStats
#'
#' This function produces patch statistics and maps on occupancy probability and mean time to colonisation.
#'
#' It uses the RangeShiftR 'population' output data.
#' @usage ColonisationStats(x, y = getwd(), years = numeric(0), maps = FALSE)
#' ColonisationStats(x, y = NULL,    years = numeric(0))
#' @param x,y Either the parameter master (\code{x}) and RS directory path (\code{y}) of the simulation or, alternatively,\cr
#' the population output as a dataframe (\code{x}). In this case \code{y} is an optional parameter taking the patch map(s) as a raster layer (stack),
#' which will then be used to produce maps of the patch statistics. For more info on this, see the Details.
#' @param years Years at which the probabilty of occupancy over all replicates will be calculated.
#' @param maps Only in the parameter master (=\code{x}) notation: For each given \code{year}, uses the current patch map (in case of a dynamic landscape, otherwise the static patch map)
#' to produce a raster of occupancy probabilties. For time to colonisation, uses the current patch map of last year of population record.
#' @details In the population dataframe (=\code{x}) notation, there are 3 options on how many maps to produce:
#' (1) y = NULL: no raster maps produced.
#' (2) y is a RasterLayer: All statistics projected onto the same patch map.
#' (3) y is a RasterStack with 2 layers: The first will be used for probabilty of occupancy for all years, the second will be used for time to colonisation.
#' (4) y is a RasterStack with length(years)+1 layers: The first ones will be used for probabilty of occupancy on each year, the last for time to colonisation.
#' @return a list with the elements
#' \code{occ_prob} (a dataframe with mean (over all replicates) probabilty of occupancy for each patch (in rows) at the given year(s) in columns),\cr
#' \code{col_time} (a dataframe with first recorded year of colonisation for each patch; with replicates in columns and rows named by patch-IDs),\cr
#' Optional:\cr\cr
#' \code{patch_occ_prob} (a raster (stack) with the data from \code{occ_prob} stored in the respective patch cells per given year),\cr
#' \code{patch_col_time} (a raster with the data from \code{col_time_mean} stored in their respective patch cells)
#' @export
setGeneric("ColonisationStats", function(x, ...) standardGeneric("ColonisationStats") )

setMethod("ColonisationStats", "data.frame", function(x, y = NULL, years = numeric(0)) {
#test <- function(x, y = NULL, years = numeric(0)) {
        if(class(years) %in% c("integer","numeric") ){
            if(length(years)==0) {
                years <- max(x$Year)
            }else{
                # check that requested years are in dataset
                if(!all(years %in% unique(x$Year))){
                    years <- sort(years[years %in% unique(x$Year)])
                    warning("ColonisationStats(): Some of the given years are not included in the population dataframe. They will be ignored.", call. = FALSE)
                }
            }
        }else{
            warning("ColonisationStats(): Years must be of class numeric or integer.", call. = FALSE)
            return(NULL)
        }

        if("PatchID" %in% names(x)) {
            patchbased <- TRUE
            x <- x[x$PatchID!=0,] # exclude matrix patch
        } else { # cell-based -> create PatchIDs from (x,y)-coordinates:
            patchbased <- FALSE
            x[c("x","y")] <- x[c("x","y")]+1  # have counting start at 1
            maxX <- max(x$x)
            maxY <- max(x$y)
            digitsY <- ifelse(maxY < 1, 0, floor(log10(maxY)) + 1)
            x$PatchID <- (x$x*(10^digitsY)) + (x$y) #+ (maxY+1-x$y)
        }

        x <- subset(x, NInd>0, select = c("Rep","PatchID","Year","NInd") )
        reps <- sort(unique(x$Rep))
        n=length(reps)

        # calculate per patch
        patchstats <- sapply(split(x, f=x$PatchID),
                             FUN=function(xs){
                                 return(list(occprob = sapply(years,
                                                              FUN=function(year){
                                                                  sum(xs$Year==year)/n
                                                              }),
                                             coltime = sapply(reps,
                                                              FUN=function(rep){
                                                                  xxs <- xs[xs$Rep==rep,"Year"]
                                                                  ifelse(length(xxs),min(xxs),NA)
                                                              }
                                             )
                                 ))
                             }
        )
        patches <- as.integer(colnames(patchstats))

        # Occupancy Probability
        occ_prob <- data.frame(patchstats["occprob",])
        rownames(occ_prob) <- years
        if(patchbased) occ_prob <- cbind(data.frame(patch=patches), t(occ_prob))
        else occ_prob <- cbind(data.frame(x=floor(patches/(10^digitsY)),y=patches%%(10^digitsY)), t(occ_prob))
        rownames(occ_prob) <- NULL

        # Time to Colonisation
        col_time <- data.frame(patchstats["coltime",])
        if(patchbased){
            col_time <- cbind(data.frame(patch=patches), t(col_time))
            names(col_time) <- c("patch",paste0("rep.",as.character(reps)))
        }else{
            col_time <- cbind(data.frame(x=occ_prob$x,y=occ_prob$y), t(col_time))
            names(col_time) <- c("x","y",paste0("rep.",as.character(reps)))
        }
        rownames(col_time) <- NULL


        # create maps
        if(!is.null(y)){
            if(patchbased){
                col_time_mean <- rowMeans(col_time[,-1], na.rm = TRUE)
                names(col_time_mean) <- patches
            }else{
                col_time_mean <- rowMeans(col_time[,-(1:2)], na.rm = TRUE)
            }
            #require('raster')

            if(class(y) == "RasterLayer" || class(y) == "RasterStack" ){
                onelayer = FALSE
                if(class(y) == "RasterLayer") onelayer = TRUE
                if(class(y) == "RasterStack") if(length(y@layers)==1 ) onelayer = TRUE

                if(patchbased){
                    # non-dynamic landscape
                    if(onelayer){
                        # initialise output rasters
                        if(class(y) == "RasterStack") y <- y[[1]]
                        patch_occ_prob <- patch_col_time <- y
                        # denote matrix with NA
                        raster::values(patch_occ_prob)[raster::values(y)==0] <- raster::values(patch_col_time)[raster::values(y)==0] <- NA
                        # init all habitat patches to also address those that never had a population (these don't occur in RS output)
                        raster::values(patch_occ_prob)[raster::values(y) >0] <- 0
                        raster::values(patch_col_time)[raster::values(y) >0] <- -9

                        # fill output rasters
                        #if(length(years)>1){
                            patch_outstack <- raster::stack()
                            for (j in 1:length(years)){
                                patch_outstack <- raster::addLayer(patch_outstack, patch_occ_prob)
                                for (i in patches){
                                    raster::values(patch_outstack[[j]])[raster::values(y)==i] <- occ_prob[occ_prob$patch==i,paste(years[j])]
                                }
                            }
                        #}
                        #else{
                        #    for (i in patches){
                        #        raster::values(patch_occ_prob)[raster::values(y)==i] <- occ_prob[occ_prob$patch==i,"occ_prob"]
                        #    }
                        #    patch_outstack <- patch_occ_prob
                        #}

                        for (i in patches){
                            raster::values(patch_col_time)[raster::values(y)==i] <- ifelse(is.na(col_time_mean[paste(i)]),-9,col_time_mean[paste(i)])
                        }

                        return(list(occ_prob=occ_prob, col_time=col_time, map_occ_prob=patch_outstack, map_col_time=patch_col_time))

                    }#else warning("ColonisationStats(): The given map raster stack must have only one layer for a non-dynamic landscape.", call. = FALSE)

                    # dynamic landscape
                    if(!onelayer){

                        N_layers <- length(years)+1
                        if( length(y@layers) != N_layers ){
                            warning("ColonisationStats(): Number of raster layers must be either 1 or number of years plus 1.", call. = FALSE)
                        }
                        else{
                            # initialise output rasters
                            patch_outstack <- y
                            # denote matrix with NA
                            raster::values(patch_outstack)[raster::values(y)==0] <- NA
                            # all habitat patches to address those that never had a population
                            raster::values(patch_outstack)[raster::values(y)>0] <- 0.0
                            raster::values(patch_outstack[[N_layers]])[raster::values(y[[N_layers]])>0] <- -9

                            # fill output rasters
                            for (i in patches){
                                #if(N_layers==2){
                                #    raster::values(patch_outstack[[1]])[raster::values(y[[1]])==i] <- occ_prob[occ_prob$patch==i,"occ_prob"]
                                #}else {
                                    for (j in 1:length(years)){
                                        raster::values(patch_outstack[[j]])[raster::values(y[[j]])==i] <- occ_prob[occ_prob$patch==i,paste(years[j])]
                                    }
                                #}
                                raster::values(patch_outstack[[N_layers]])[raster::values(y[[N_layers]])==i] <- ifelse(is.na(col_time_mean[paste(i)]),-9,col_time_mean[paste(i)])
                            }
                            return(list(occ_prob=occ_prob, col_time=col_time, map_occ_prob=patch_outstack[[-N_layers]], map_col_time=patch_outstack[[N_layers]]))

                        }
                    }#else warning("ColonisationStats(): Given map must be a raster stack for a dynamic landscape.", call. = FALSE)

                }else{ # cell-based
                    if(onelayer){

                        # initialise output rasters
                        if(class(y) == "RasterStack") y <- y[[1]]
                        patch_occ_prob <- patch_col_time <- y
                        # init all habitat patches to also address those that never had a population (these don't occur in RS output)
                        raster::values(patch_occ_prob)[!is.na(raster::values(y))] <- 0
                        raster::values(patch_col_time)[!is.na(raster::values(y))] <- -9
                        # make value index from patchIDs
                        value_ix <- floor(patches/(10^digitsY))+(nrow(y)-patches%%(10^digitsY))*ncol(y)

                        # fill output rasters
                        #if(length(years)>1){
                            patch_outstack <- raster::stack()
                            for (j in 1:length(years)){
                                patch_outstack <- raster::addLayer(patch_outstack, patch_occ_prob)
                                raster::values(patch_outstack[[j]])[value_ix] <- occ_prob[,j+2]
                            }
                        #} else{
                        #    raster::values(patch_occ_prob)[value_ix] <- occ_prob[,1]
                        #    patch_outstack <- patch_occ_prob
                        #}
                        raster::values(patch_col_time)[value_ix] <- ifelse(is.na(col_time_mean[]),-9,col_time_mean[])

                        return(list(occ_prob=occ_prob, col_time=col_time, map_occ_prob=patch_outstack, map_col_time=patch_col_time))

                    }else warning("ColonisationStats(): The given map raster stack must have only one layer for a cell-based landscape.", call. = FALSE)
                }

            }else warning("ColonisationStats(): Given map is not a raster layer or raster stack.", call. = FALSE)
        }
        return(list(occ_prob=occ_prob, col_time=col_time))
})

setMethod("ColonisationStats", "RSparams", function(x, y = getwd(), years = numeric(0), maps = FALSE) {
    res <- NULL
    if(class(x@land)=="ImportedLandscape" || class(maps)=="logical") {
        if(x@simul@OutIntPop>0){
            if(!is.null(y) & class(y)=="character" ){
                if(class(years) %in% c("integer","numeric") ){

                    # read population output
                    pop_df <- try(readPop(x, y, center=FALSE))
                    if ( class(pop_df) == "try-error" ) {            # try() returned a warning
                        warning("ColonisationStats(): Couldn't read population output for this simulation.", call. = FALSE)
                        return(NULL)
                    }
                    if ( class(pop_df) == "character" ) return(NULL) # readPop() returned a warning

                    # set year to last recorded year if none is given
                    if(length(years)==0) years <- max(pop_df$Year)

                    # read patch rasters if needed
                    if(maps){
                        #require('raster')
                        if(x@control@patchmodel){ # for patch-based model, read all relevant patch-maps
                            # non-dynamic landscape
                            if(length(x@land@LandscapeFile)==1){

                                patch_r <- try(raster::raster(paste0(dirpath, "Inputs/", x@land@PatchFile)))
                                if( class(patch_r) == "try-error" ) warning("ColonisationStats(): Couldn't read patch raster file ", x@land@PatchFile , call. = FALSE)

                                if( class(pop_df) == "data.frame" & class(patch_r) == "RasterLayer" ) res <- ColonisationStats(pop_df,patch_r,years)

                            }
                            # dynamic landscape
                            else{
                                patch_r <- raster::stack()
                                # rasters for occ_prob output
                                for(year in years){
                                    current <- which(x@land@DynamicLandYears == max(x@land@DynamicLandYears[x@land@DynamicLandYears<=year]) )
                                    patch_curr <- try(raster::raster(paste0(dirpath, "Inputs/", x@land@PatchFile[current])))
                                    if ( class(patch_curr) == "try-error" ) warning("ColonisationStats(): Couldn't read patch raster file nr ", current , " for this simulation.", call. = FALSE)
                                    else patch_r <- raster::addLayer(patch_r , patch_curr)
                                }
                                # rasters for col_time output
                                year <- max(pop_df$Year)
                                current <- which(x@land@DynamicLandYears == max(x@land@DynamicLandYears[x@land@DynamicLandYears<=year]) )
                                patch_curr <- try(raster::raster(paste0(dirpath, "Inputs/", x@land@PatchFile[current])))
                                if ( class(patch_curr) == "try-error" ) warning("ColonisationStats(): Couldn't read patch raster file nr ", current , " for this simulation.", call. = FALSE)
                                else patch_r <- raster::addLayer(patch_r , patch_curr)

                                if(class(pop_df) == "data.frame" & length(patch_r@layers)==(length(years)+1) ) res <- ColonisationStats(pop_df,patch_r,years)
                            }
                        }else{
                            # for cell-based model, read only main habitat maps to use as raster template
                            patch_r <- try(raster::raster(paste0(dirpath, "Inputs/", x@land@LandscapeFile[1])))
                            if ( class(patch_r) == "try-error" ) warning("ColonisationStats(): Couldn't read patch raster file nr ", current , " for this simulation.", call. = FALSE)
                            if(class(pop_df) == "data.frame" & class(patch_r) == "RasterLayer" ) res <- ColonisationStats(pop_df,patch_r,years)
                        }
                    }else { # no maps requested
                        if(class(pop_df) == "data.frame") res <- ColonisationStats(pop_df,NULL,years)
                    }

                    if(is.null(res)) warning("ColonisationStats(): Couldn't get population output or patch raster file.", call. = FALSE)
                    else return(res)

                }else {
                    warning("ColonisationStats(): Years must be of class numeric or integer.", call. = FALSE)}
            }else {
                warning("ColonisationStats(): dirpath must be set", call. = FALSE)}
        }else{
            warning("ColonisationStats(): This simulation has population output turned off (OutIntPop=0), but that is needed to calculate the colonisation statistics.", call. = FALSE)}
    }else {
        warning("ColonisationStats(): For Artificial Landscape models, maps must be FALSE", call. = FALSE)}
})




#---------------------------------------------------------

### PLOTTING


#' Plot Abundance
#'
#' Uses the RangeShiftR output data 'range' to generate abundance time series.
#' Plots the mean abundance over all replicates, and optionally the standard deviation and/or the single replicates.
#' @param s RSparams object or a data.frame in the 'range' file format
#' @param dirpath RS directory path; required if \code{s} is a \code{RSparams}
#' @param sd plot standard deviation? (default is \code{FALSE})
#' @param replicates plot the replicates? (default is \code{TRUE})
#' @param ylim upper limit to the y-axis
#' @export
setGeneric("plotAbundance", function(s,...) standardGeneric("plotAbundance") )

setMethod("plotAbundance", "data.frame", function(s, sd = FALSE, replicates = TRUE, ylim=NULL, ...) {
    # Calculate means
    rep_means <- aggregate(NInds~Year, data = s, FUN = "mean")
    # Calculate standard deviation
    if (sd) {
        rep_sd <- aggregate(NInds~Year, data = s, FUN = "sd")
        rep_sd[is.na(rep_sd)] <- 0
    }
    # Set y-limits
    if (is.null(ylim)) {
        if (replicates) {
            ylim <- c(min(s$NInds),max(s$NInds))
        }else {
            if (sd) {
                ylim <- c(max(0,min(rep_means$NInds-rep_sd$NInds)), max(rep_means$NInds+rep_sd$NInds))
            }else{
                ylim <- c(min(rep_means$NInds),max(rep_means$NInds))
            }
        }
    }
    # New plot
    plot(NULL, type = "n", ylab = "Abundance", xlab = "Year", xlim=c(0, max(s$Year)), ylim=ylim, ...)
    # Plot standard deviation
    if (sd) {
        polygon(c(rep_sd$Year,rev(rep_sd$Year)), c(rep_means$NInds+rep_sd$NInds, rev(pmax(0,rep_means$NInds-rep_sd$NInds))), border=NA, col='grey80')
    }
    # Plot replicates
    if (replicates) {
        for (i in 0:max(s$Rep)) {
            lines(s$Year[s$Rep==i], s$NInds[s$Rep==i], type = "l", lwd = 0.5)
        }
    }
    # Plot abundance
    lines(rep_means$Year, rep_means$NInds, type = "l", lwd = 3, col = "red")
})
setMethod("plotAbundance", "RSparams", function(s, dirpath, ...) {
    if (class(dirpath)=="character"){
        range_table <- try(readRange(s, dirpath))
        if ( class(range_table) == "data.frame") {
            plotAbundance(range_table, ...)
        }
    }else{
        warning("plotAbundance(): dirpath must be of type character.", call. = TRUE)
    }
})


#' Plot Occupancy
#'
#' Uses the RangeShiftR output data 'range' to generate occupancy time series.
#' Plots the mean occupancy over all replicates, and optionally the standard deviation and/or the single replicates.
#' @param s RSparams object or a data.frame in the 'range' file format
#' @param dirpath RS directory path; required if \code{s} is a \code{RSparams}
#' @param sd plot standard deviation? (default is \code{FALSE})
#' @param replicates plot the replicates? (default is \code{TRUE})
#' @param ylim upper limit to the y-axis
#' @export
setGeneric("plotOccupancy", function(s,...) standardGeneric("plotOccupancy") )

setMethod("plotOccupancy", "data.frame", function(s, sd = FALSE, replicates = TRUE, ylim=NULL, ...) {
    names(s)[grep('NOccup',names(s))] <- 'NOccup'
    # Calculate means
    rep_means <- aggregate(NOccup~Year, data = s, FUN = "mean")
    # Calculate standard deviation
    if (sd) {
        rep_sd <- aggregate(NOccup~Year, data = s, FUN = "sd")
        rep_sd[is.na(rep_sd)] <- 0
    }
    # Set y-limits
    if (is.null(ylim)) {
        if (replicates) {
            ylim <- c(min(s$NOccup),max(s$NOccup))
        } else {
            if (sd) {
                ylim <- c(max(0,min(rep_means$NOccup-rep_sd$NOccup)), max(rep_means$NOccup+rep_sd$NOccup))
            } else {
                ylim <- c(min(rep_means$NOccup),max(rep_means$NOccup))
            }
        }
    }
    # New plot
    plot(NULL, type = "n", ylab = "Occupancy", xlab = "Year", xlim=c(0, max(s$Year)), ylim=ylim, ...)
    # Plot standard deviation
    if (sd) {
        polygon(c(rep_sd$Year,rev(rep_sd$Year)), c(rep_means$NOccup+rep_sd$NOccup, rev(pmax(0,rep_means$NOccup-rep_sd$NOccup))), border=NA, col='grey80')
    }
    # plot replicates
    if (replicates) {
        for (i in 0:max(s$Rep)) {
            lines(s$Year[s$Rep==i], s$NOccup[s$Rep==i], type = "l", lwd = 0.5)
        }
    }
    # Plot occupancy
    lines(rep_means$Year, rep_means$NOccup, type = "l", lwd = 3, col = "blue")
})
setMethod("plotOccupancy", "RSparams", function(s, dirpath, ...) {
    if (class(dirpath)=="character"){
        range_table <- try(readRange(s, dirpath))
        if ( class(range_table) == "data.frame") {
            plotOccupancy(range_table, ...)
        }
    }else{
        warning("plotOccupancy(): dirpath must be of type character.", call. = TRUE)
    }
})


#---------------------------------------------------------

#--- SMS PathLenghts

#' Get the distribution of SMS Path Lengths
#'
#' Reads the RangeShiftR output files 'MovePaths' (if they were generated) to get the distribution of lengths of SMS
#' paths taken in all recorded years.
#' @param s RSmaster parameter object
#' @param dirpath RS directory path
#' @return a data.frame that contains the mean (over replicates) number of SNS paths of a given length (rows) per simulated map (columns)
#' @export
setGeneric("SMSpathLengths", function(s,dirpath,...) standardGeneric("SMSpathLengths") )

setMethod("SMSpathLengths", c(s="RSparams", dirpath="character"), function(s,dirpath) {
    if(class(s@dispersal@Transfer)!="StochMove"){
        warning("SMSpathLengths(): Transfer must be of type SMS.", call. = TRUE)
        return(NULL)
    }
    if(s@simul@OutIntPaths<1){
        warning("SMSpathLengths(): Output of SMS paths ('OutIntPaths') is disabled in this RS simulation.", call. = TRUE)
        return(NULL)
    }

    maxLength = s@dispersal@Settlement@MaxSteps
    if(maxLength<1) maxLength = as.integer(log(.015)/log(1-s@dispersal@Transfer@StepMort)) + 1
    year_blocks = s@land@DynamicLandYears
    nr_maps = length(year_blocks)

    steps_sum <- matrix(0, nrow = maxLength, ncol = nr_maps)
    rownames(steps_sum) <- sapply(1:maxLength, FUN=paste0)
    null_tb <- rep(0,maxLength)
    names(null_tb) <- sapply(1:maxLength, FUN=paste0)
    for(i in 0:(s@simul@Replicates-1)){
        steps <- try(read.table(paste0(dirpath,
                                   "Outputs/Batch",s@control@batchnum,
                                   "_Sim",s@simul@Simulation,
                                   "_Land",s@land@LandNum,
                                   "_Rep",i,
                                   "_MovePaths.txt"),
                            header = T))
        if ( class(steps) == "try-error" ) {
            warning(cat("SMSpathLengths(): Couldn't read MovePaths output for replicate",i,", skipping it."), call. = FALSE)
            #return(NULL)
        }else{
            steps_y <- sapply(seq(nr_maps), FUN = function(y){
                tb <- null_tb
                lo <- year_blocks[y]
                if (y==nr_maps) tb_n <- table(aggregate(Step~IndID, data=subset(steps, Year >= lo), FUN = max)$Step)
                else {
                    up <- year_blocks[y+1]
                    tb_n <- table(aggregate(Step~IndID, data=subset(steps, Year >= lo & Year < up), FUN = max)$Step)
                }
                ix <- names(tb) %in% names(tb_n)
                tb[ix] <- tb_n[1:sum(ix)]
                tb
            })
            steps_sum <- steps_sum + steps_y
        }
    }
    colnames(steps_sum) <- sapply(1:nr_maps, function(m){paste0('Map_',m)})
    # mean over replicates
    steps_sum <- data.frame(steps_sum / s@simul@Replicates)

    return(steps_sum)
})




#---------------------------------------------------------

### B FUNCTION


## ---- Backend matrix model simulation function -----

get_eq_pop <- function(b, demog, N_0 = NULL, t_max = 1000, t_rec = 1, delta = .1, diagnostics = FALSE, rm.stage0 = TRUE){
    if(class(demog)!="DemogParams"){
        warning("This function expects an object of class DemogParams.", call. = FALSE)
        return(NULL)
    }
    else{
        if(class(demog@StageStruct)!="StagesParams"){
            warning("demog needs to have a StageStructure", call. = FALSE)
            return(NULL)
        }
    }

    TraMa <- demog@StageStruct@TransMatrix

    if(demog@ReproductionType>1){
        # for sex-specific (ReproductionType = 2) need to include Sex-ratio to weigh fecundities and make transition matrix quadratic
        SexDep <- TRUE
        TraMa <- rbind(TraMa[1,]*demog@PropMales,TraMa)
        TraMa[2,] <- TraMa[2,]*(1-demog@PropMales)
    }
    else SexDep <- FALSE
    TraMa_t <- TraMa
    lines <- nrow(TraMa)
    N_rec <- matrix(0, ncol = t_rec, nrow = lines-rm.stage0*(1+SexDep))

    if(demog@StageStruct@FecDensDep){
        FecDensDep <- TRUE
        if(demog@StageStruct@FecStageWts) {
            FecStageDep <- TRUE
            FecStageWts <- demog@StageStruct@FecStageWtsMatrix
        }
        else FecStageDep <- FALSE
    }
    else FecDensDep <- FALSE

    if(demog@StageStruct@DevDensDep){
        DevDensDep <- TRUE
        C_dev <- demog@StageStruct@DevDensCoeff
        if(demog@StageStruct@DevStageWts) {
            DevStageDep <- TRUE
            DevStageWts <- demog@StageStruct@DevStageWtsMatrix
        }
        else DevStageDep <- FALSE
        # disentangle survival and development rates:
        surv <- devs <- rep(0,lines)
        for(s in 1:lines ){
            ss <- TraMa[s,s]
            if (SexDep) if((s+2)>lines) dd <- 0 else dd <- TraMa[s+2,s]
            else if(s==lines) dd <- 0 else dd <- TraMa[s+1,s]
            surv[s] <- ss+dd
            devs[s] <- dd/(ss+dd)
        }
        surv_t <- surv
        devs_t <- devs
    }
    else DevDensDep <- FALSE

    if(demog@StageStruct@SurvDensDep){
        SurvDensDep <- TRUE
        C_surv <- demog@StageStruct@SurvDensCoeff
        if(demog@StageStruct@SurvStageWts) {
            SurvStageDep <- TRUE
            SurvStageWts <- demog@StageStruct@SurvStageWtsMatrix
        }
        else SurvStageDep <- FALSE
    }
    else SurvDensDep <- FALSE

    if(is.null(N_0)) N_0 <- matrix( rep(1/b/lines, lines), ncol = 1)
    N_t2 <- N_0
    t <- 0
    repeat{
        # reset abundance vector and time-dependent transition matrix
        N_t1 <- N_t2
        TraMa_t <- TraMa
        # calculate density-dependent rates
        abund <- sum(N_t1)
        if(SexDep){
            if(FecDensDep) {
                if(FecStageDep) {
                    TraMa_t[1:2,] <- TraMa[1:2,]*rep(exp(-b*FecStageWts%*%N_t1),each=2)
                }
                else TraMa_t[1:2,] <- TraMa[1:2,]*exp(-b*abund)
            }
            if(DevDensDep) {
                if(DevStageDep) {
                    # dev diagonals
                    devs_t <- devs*exp(-C_dev*b*DevStageWts%*%N_t1)
                    diag(TraMa_t) <- surv*(1-devs_t)
                    # dev off-diagonals
                    TraMa_t[row(TraMa)==col(TraMa)+2] <- TraMa[row(TraMa)==col(TraMa)+2]*exp(-C_dev*b*(DevStageWts%*%N_t1)[-c(lines-1,lines)])
                }
                else {
                    # dev diagonals
                    devs_t <- devs*exp(-C_dev*b*abund)
                    diag(TraMa_t) <- surv*(1-devs_t)
                    # dev off-diagonals
                    TraMa_t[row(TraMa)==col(TraMa)+2] <- TraMa[row(TraMa)==col(TraMa)+2]*exp(-C_dev*b*abund)
                }
                if(SurvDensDep) {
                    if(SurvStageDep){
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa_t)*exp(-C_surv*b*SurvStageWts%*%N_t1)
                        # surv off-diagonals
                        TraMa_t[row(TraMa_t)==col(TraMa_t)+2] <- TraMa_t[row(TraMa_t)==col(TraMa_t)+2]*exp(-C_surv*b*(SurvStageWts%*%N_t1)[-c(lines-1,lines)])
                    }
                    else {
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa_t)*exp(-C_surv*b*abund)
                        # surv off-diagonals
                        TraMa_t[row(TraMa_t)==col(TraMa_t)+2] <- TraMa_t[row(TraMa_t)==col(TraMa_t)+2]*exp(-C_surv*b*abund)
                    }
                }
            }
            else{
                if(SurvDensDep) {
                    if(SurvStageDep){
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa)*exp(-C_surv*b*SurvStageWts%*%N_t1)
                        # surv off-diagonals
                        TraMa_t[row(TraMa)==col(TraMa)+2] <- TraMa[row(TraMa)==col(TraMa)+2]*exp(-C_surv*b*(SurvStageWts%*%N_t1)[-c(lines-1,lines)])
                    }
                    else {
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa)*exp(-C_surv*b*abund)
                        # surv off-diagonals
                        TraMa_t[row(TraMa)==col(TraMa)+2] <- TraMa[row(TraMa)==col(TraMa)+2]*exp(-C_surv*b*abund)
                    }
                }
            }
        }
        else{ # !SexDep
            if(FecDensDep) {
                if(FecStageDep) TraMa_t[1,] <- TraMa[1,]*exp(-b*FecStageWts%*%N_t1)
                else TraMa_t[1,] <- TraMa[1,]*exp(-b*abund)
            }
            if(DevDensDep) {
                if(DevStageDep) {
                    # dev diagonals
                    devs_t <- devs*exp(-C_dev*b*DevStageWts%*%N_t1)
                    diag(TraMa_t) <- surv*(1-devs_t)
                    # dev off-diagonals
                    TraMa_t[row(TraMa)==col(TraMa)+1] <- TraMa[row(TraMa)==col(TraMa)+1]*exp(-C_dev*b*(DevStageWts%*%N_t1)[-lines])
                }
                else {
                    # dev diagonals
                    devs_t <- devs*exp(-C_dev*b*abund)
                    diag(TraMa_t) <- surv*(1-devs_t)
                    # dev off-diagonals
                    TraMa_t[row(TraMa)==col(TraMa)+1] <- TraMa[row(TraMa)==col(TraMa)+1]*exp(-C_dev*b*abund)
                }
                if(SurvDensDep) {
                    if(SurvStageDep){
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa_t)*exp(-C_surv*b*SurvStageWts%*%N_t1)
                        # surv off-diagonals
                        TraMa_t[row(TraMa_t)==col(TraMa_t)+1] <- TraMa_t[row(TraMa_t)==col(TraMa_t)+1]*exp(-C_surv*b*(SurvStageWts%*%N_t1)[-lines])
                    }
                    else {
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa_t)*exp(-C_surv*b*abund)
                        # surv off-diagonals
                        TraMa_t[row(TraMa_t)==col(TraMa_t)+1] <- TraMa_t[row(TraMa_t)==col(TraMa_t)+1]*exp(-C_surv*b*abund)
                    }
                }
            }
            else{
                if(SurvDensDep) {
                    if(SurvStageDep){
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa)*exp(-C_surv*b*SurvStageWts%*%N_t1)
                        # surv off-diagonals
                        TraMa_t[row(TraMa)==col(TraMa)+1] <- TraMa[row(TraMa)==col(TraMa)+1]*exp(-C_surv*b*(SurvStageWts%*%N_t1)[-lines])
                    }
                    else {
                        # surv diagonals
                        diag(TraMa_t) <- diag(TraMa)*exp(-C_surv*b*abund)
                        # surv off-diagonals
                        TraMa_t[row(TraMa)==col(TraMa)+1] <- TraMa[row(TraMa)==col(TraMa)+1]*exp(-C_surv*b*abund)
                    }
                }
            }
        }
        # turn RS transition matrix into actual Leslie matrix (i.e. remove zeroth stage)
        if(rm.stage0){
            if(SexDep){
                TraMa_t[1:2,] <- TraMa_t[1:2,] * diag(TraMa_t[c(3,4),c(1,2)]) # put survival rate into fecundity
                TraMa_t[3:4,] <- TraMa_t[3:4,] + TraMa_t[1:2,]
                TraMa_t <- TraMa_t[3:lines,3:lines]                          # get rid of stage 0
                # compress state vector as well
                N_t1[3] <- sum(N_t1[c(1,3)])
                N_t1[4] <- sum(N_t1[c(2,4)])
                N_t1 <- matrix(N_t1[3:lines], ncol = 1)
            }
            else{ # !SexDep
                TraMa_t[1,] <- TraMa_t[1,] * TraMa_t[2,1]                   # -"-
                TraMa_t[2,] <- TraMa_t[2,] + TraMa_t[1,]
                TraMa_t <- TraMa_t[2:lines,2:lines]
                # compress state vector as well
                N_t1[2] <- sum(N_t1[1:2])
                N_t1 <- matrix(N_t1[2:lines], ncol = 1)
            }
        }
        # take a time step
        N_t2 <- TraMa_t %*% N_t1
        # evaluate convergence
        del <- sqrt(sum((N_t2-N_t1)^2))
        N_rec[,1+(t%%t_rec)] <- N_t2
        if(del < delta || t >= t_max ){
            if(t >= t_max) warning("No convergence")
            if(t < t_rec )  N_rec <- N_rec[,1:t]
            if(diagnostics) return(list(N=N_rec, TraMa_t=TraMa_t, t=t, del=del))
            else return(N_rec)
            break
        }
        # end time step
        t <- t+1
        # expand state vector again
        if(rm.stage0){
            if(SexDep){
                N_t2[1:2] <- N_t2[1:2]/2
                N_t2 <- matrix(c(N_t2[1:2],N_t2), ncol = 1)
            }
            else{ # !SexDep
                N_t2[1] <- N_t2[1]/2
                N_t2 <- matrix(c(N_t2[1],N_t2), ncol = 1)
            }
        }
    }
}


## ---- Frontend Plot function -----

#' Calculates the equilibrium density and stage distribution for a localised (i.e. non-spatial) closed population
#'
#' Uses the \emph{RangeShiftR} Demography module to create the corresponding matrix model and runs it until equilibrium is reached.
#' This corresponds to a localised population in a single cell without dispersal (i.e. no immigration or emigration).
#' Since the matrix model representation is used, some options (e.g. maximum age) of the \code{\link[RangeShiftR]{Demography}} module can not be taken into account.
#' @param demog DemogParams object with a \code{StageStructure}
#' @param DensDep_values values of \eqn{1/b} to run the matrix model for
#' @param plot plot the equilibrium densities? (default is \code{TRUE})
#' @param stages_out which stages to plot? (defaults to all)
#' @param juv.stage use explicit juvenile (zeroth) stage? (default is \code{TRUE})
#' @param t_rec number of time steps to record at the end (defaults to \eqn{1}); if \code{t_rec}\eqn{>1}, the mean over all recorded time steps is returned
#' @param N_0 initial condition, i.e. population density at time zero; must include stage zero regardless of the value of \code{juv.stage}
#' @param t_max allowed number of time steps to reach equilibrium (default is \eqn{1000})
#' @param delta tolerance to check for equilibrium (default is \eqn{.1}); the utilised measure is euclidian distance of current to previous time steps
#' @param diagnostics in addition to recorded population vectors, returns the number of steps taken as well as the transition matrix and the value if delta at the last step (default is \code{FALSE})
#' @details \emph{RangeShiftR} requires an additional juvenile stage to be added to the common transition matrix as stage \eqn{0} (in order
#' to allow for juvenile dispersal). For the simulation with \code{RunMatrixModel()}, this stage can be kept (\code{juv.stage=TRUE})
#' or removed to yield the corresponding Lefkovitch matrix (\code{juv.stage=FALSE}).\cr
#' The default initial state \code{N_0} is a population at its respective density \eqn{1/b} with unpopulated juvenile stage and
#' all higher stages equally populated.\cr
#' @return a matrix of population densities with a named row for each stage and a column for each given value of \eqn{1/b}
#' @export
setGeneric("getLocalisedEquilPop", function(demog,...) standardGeneric("getLocalisedEquilPop") )

setMethod("getLocalisedEquilPop", "DemogParams", function(demog, DensDep_values, plot=TRUE, stages_out=NULL, juv.stage=TRUE, t_rec=1,
                                                    t_max = 1000, N_0 = NULL, delta=.1, diagnostics=FALSE){
    # make b from 1/b
    b_vector <- 1/DensDep_values
    # calculate equilibrium population
    if(t_rec==1) res <- sapply(b_vector, get_eq_pop, demog=demog, t_rec = t_rec, rm.stage0 = !juv.stage)
    if(t_rec>1) {res <- lapply(b_vector, get_eq_pop, demog=demog, t_rec = t_rec, rm.stage0 = !juv.stage)
                 res <- sapply(res, rowMeans)}
    # name rows
    if(demog@ReproductionType <2){row.names(res) <- seq((!juv.stage),demog@StageStruct@Stages-1)}
    if(demog@ReproductionType==2){row.names(res) <- rep(seq((!juv.stage),demog@StageStruct@Stages-1), each = 2)}
    # plot stages ?
    if(plot) {
        # x-axis labels
        xlabs <- print(DensDep_values, digits = 2)
        # which stages to plot?
        if(is.null(stages_out)) stages_out = seq((!juv.stage),demog@StageStruct@Stages-1)
        colors <- hcl.colors(length(stages_out), palette = "Harmonic")
        if(demog@ReproductionType <2){
            barplot(res[as.character(stages_out),], names.arg = xlabs, beside = F, col = colors,
                          main = "Localised equilibrium densities", xlab = "1/b", ylab = "Population density")
        }
        if(demog@ReproductionType==2){
            if(length(stages_out)<2) warning("getLocalisedEquilPop(): Please specify more than one stage when plotting a sex-explicit model.", call. = TRUE)
            else {
                res_2 <- res[which(rownames(res) %in% stages_out),]
                mal <- seq.int(1,length(stages_out)*2,2)
                fem <- seq.int(2,length(stages_out)*2,2)
                res_2 <- cbind(res_2[mal,1:length(DensDep_values)],res_2[fem,1:length(DensDep_values)])
                res_2 <- cbind(res_2[,c(sapply(1:length(DensDep_values), function(i){c(0,1)*length(DensDep_values)+i}))])
                barplot(res_2, space=c(0.3,0.1), names.arg = c(rbind(xlabs,NA)), beside = F, col = rep(colors, 2),
                        main = "Localised equilibrium densities", xlab = "1/b", ylab = "Population density")
                text(seq(0.5,length(DensDep_values)*2,2)*1.2, colSums(res_2[,seq(1,length(DensDep_values)*2,2)])*1.1, "m", cex=1, col="black")
                text(seq(1.5,length(DensDep_values)*2,2)*1.2, colSums(res_2[,seq(2,length(DensDep_values)*2,2)])*1.1, "f", cex=1, col="black")
            }
        }
        legend("topleft", legend = rev(sapply(stages_out, function(s){paste("Stage",s)})), col = rev(colors), pch = 16)
    }
    # return equilibrium population densities
    return(res)
})

