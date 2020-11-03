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
 
 

### CLASS SIMULATIONPARAMS

# from RS 'Parameter' file

#' Set Simulation parameters
#'
#' @description Set basic simulation parameters and control output types.\cr
#' Furthermore, optionally define a (\code{Shifting}) Environmental \code{Gradient}, Environmental Stochasticity (\code{EnvStoch}) and/or Local extinction (\code{LocalExt}).
#' (These options are to be moved to separate classes in future versions.)
#'
#' @author Anne-Kathleen Malchow
#' @usage Simulation(Simulation = 1, Replicates = 2, Years = 50, Absorbing = FALSE,
#'            Gradient = 0, GradSteep, Optimum, f, ExtinctOptim,
#'            Shifting = FALSE, ShiftRate, ShiftStart, ShiftEnd,
#'            LocalExt = FALSE, LocalExtProb,
#'            EnvStoch = 0, EnvStochType, std, ac, minR, maxR, minK, maxK,
#'            OutIntRange = 1, OutIntOcc = 0,
#'            OutIntPop = 1, OutIntInd = 0,
#'            OutIntTraitCell = 0, OutIntTraitRow = 0,
#'            OutIntConn = 0, OutIntPaths = 0, OutIntGenetic = 0,
#'            OutGenType = 0, OutGenCrossTab = FALSE,
#'            OutStartPop = 0, OutStartInd = 0,
#'            OutStartTraitCell = 0, OutStartTraitRow = 0,
#'            OutStartConn = 0, OutStartPaths = 0, OutStartGenetic = 0,
# #'            SaveMaps = FALSE, MapsInterval, DrawLoadedSp = FALSE,,
# #'            ReturnPopRaster = FALSE, CreatePopFile = TRUE
#'            SMSHeatMap = FALSE)
#' @param Simulation ID number of current simulation, defaults to \eqn{1}. (integer)
#' @param Replicates Number of simulation iterations, defaults to \eqn{2}. (integer)
#' @param Years The number of simulated years, defaults to \eqn{50}. (integer)
#' @param Absorbing If \code{FALSE} (default), every move in the \code{\link[RangeShiftR]{Transfer}} process will be
#' repeated until a valid cell is met.\cr
#' If \code{TRUE}, an individual which hits a non-valid cell or
#' transgresses the landscape boundary during the dispersal act is eliminated from the simulation.
#' @param Gradient Whether to apply north-south gradient:
#' \eqn{0} = None (default), \cr
#' \eqn{1} = decreasing \code{K_or_DensDep}\cr
#' \eqn{2} = decreasing growth rate \eqn{r} or, for a \code{\link[RangeShiftR]{StageStructure}}d population,
#' fecundity \ifelse{html}{\out{&phi;}}{\eqn{\phi}}, \cr
#' \eqn{3} = increasing local extinction probability \eqn{e}. \cr
#' If activated, a gradient will be imposed along the north-south (\eqn{y}-) axis in which the selected parameter varies linearly with distance from the
#' optimum \eqn{y}-value \code{Optimum}. Note that a \code{Gradient} can not be applied in for patch-based models (must be \code{Gradient}\eqn{=0}).
#' @param GradSteep Required if \code{Gradient} \ifelse{html}{\out{&ne; 0}}{\eqn{> 0}}: gradient steepness in units of (fraction of local value) per cell. Must be \eqn{\ge 0}.
#' @param Optimum Required if \code{Gradient} \ifelse{html}{\out{&ne; 0}}{\eqn{> 0}}: \eqn{y}-value at which the extremum is obtained. Must be \eqn{\ge 0}.
#' @param f Required if \code{Gradient} \ifelse{html}{\out{&ne; 0}}{\eqn{> 0}}: local scaling factor that determines the magnitude of stochastic local heterogeneity relative to the optimal value. Must be \eqn{\ge 0}.
#' @param ExtinctOptim Required if \code{Gradient} \eqn{= 3}: optimum (i.e. minimal) local extinction probability at \code{Optimum}. Must be between \eqn{0} and \eqn{1}.
#' @param Shifting Only applicable if \code{Gradient} \ifelse{html}{\out{&ne; 0}}{\eqn{> 0}}:\cr
#' If \code{FALSE} (default), the gradient is stationary.\cr
#' If \code{TRUE}, the \code{Gradient} shifts along the \eqn{y}-axis towards increasing \eqn{y} (northwards). Requires to set \code{ShiftRate}, \code{ShiftStart} and \code{ShiftEnd},
#' @param ShiftRate Required if \code{Shifting=TRUE}: shift rate of the gradient in units of rows per year. (integer)
#' @param ShiftStart Required if \code{Shifting=TRUE}: year in which the gradient shifting starts (integer)
#' @param ShiftEnd Required if \code{Shifting=TRUE}: year in which the gradient shifting stops. (integer)
#' @param LocalExt If \code{FALSE} (default), no additional extinction probability is applied.\cr
#' If \code{TRUE}, an independent constant extinction probability \code{LocalExtProb} is applied, defined as
#' the probability that each population goes extinct at each year.
#' Note that \code{LocalExt} must be \code{FALSE} for a patch-based model or if \code{Gradient}\eqn{=3}.
#' @param LocalExtProb Required if \code{LocalExt=TRUE}: independent yearly extinction probability of populations.
#' @param EnvStoch Scale of environmental stochasticity:\cr
#' \eqn{0} = none (default),\cr
#' \eqn{1} = global (a single time series for entire landscape),\cr
#' \eqn{2} = local (each cell fluctuates independently, only permitted for cell-based model).\cr
#' Environmental stochasticity is always applied on a yearly basis.
#' @param EnvStochType Required if \code{EnvStoch} \ifelse{html}{\out{&ne; 0}}{\eqn{> 0}}: Parameter to which environmental stochasticity is applied:\cr
#' \eqn{0} = growth rate \eqn{r} or, for a \code{\link[RangeShiftR]{StageStructure}}d population, fecundity (\eqn{\phi}).\cr
#' \eqn{1} = demographic density dependence \code{K_or_DensDep} (carrying capacity or 1/b) (allowed for artificial landscapes only!).
#' @param std Required if \code{EnvStoch} \ifelse{html}{\out{&ne; 0}}{\eqn{> 0}}: magnitude of stochastic fluctuations. Must be \eqn{> 0.0} and \eqn{\le 1.0}.
#' @param ac Required if \code{EnvStoch} \ifelse{html}{\out{&ne; 0}}{\eqn{> 0}}: temporal autocorrelation coefficient. Must be \eqn{\ge 0.0} and \eqn{<1.0}.
#' @param minR,maxR Required if \code{EnvStochType}\eqn{=0}: minimum and maximum growth rates.
#' @param minK,maxK Required if \code{EnvStochType}\eqn{=1}: minimum and maximum value of \eqn{K} or \eqn{1/b}, respectively.
#' @param OutIntRange,OutIntOcc,OutIntPop,OutIntInd,OutIntGenetic,OutIntTraitCell,OutIntTraitRow,OutIntConn,OutIntPaths Control the various types
#' of Output files, i.e. range, occupancy, populations, individuals, traits (by cell or by row), connectivity, SMS paths and genetics:\cr
#'  \eqn{=0 }: Output disabled.\cr
#'  \eqn{>0 }: Output enabled; sets interval (in years) in which output is generated.\cr
#' If the output is enabled, start values are required. By default, only the output of Range and Population are enabled.\cr
#' Occupancy output is only applicable if \code{Replicates>1}.
#' Traits output is only applicable for a cell-based model with inter-individual variability.
#' Connectivity output is only applicable for a patch-based model.
#' SMS paths is only applicable for a model with SMS transfer method.
#' @param OutStartPop,OutStartInd,OutStartGenetic,OutStartTraitCell,OutStartTraitRow,OutStartConn,OutStartPaths
#' Starting years for output generation. Note that the first year is year \eqn{0}. Defaults to \eqn{0} for all output types. (integer)
#' @param OutGenType Required if \code{OutIntGenetic}\eqn{>0}: Genetics output will be generated for:\cr
#' \eqn{0} = juveniles only (default), \eqn{1} =  all individuals, \eqn{2} = adults only.
#' @param OutGenCrossTab Required if \code{OutIntGenetic}\eqn{>0}:\cr
#' If \code{FALSE} (default), Genetics output will be written to several files.
#' \cr If \code{TRUE} Genetics output is generated as a cross table.
# #' @param SaveMaps If \code{FALSE} (default), no maps will be created.\cr If \code{TRUE}, maps will be generated.
# #' beginning in the first year in accordance to \code{MapsInterval}.
# #' @param MapsInterval Required if \code{SaveMaps=TRUE}: save maps every \eqn{n} reproductive seasons. (integer)
# #' @param DrawLoadedSp If \code{FALSE} (default), only the simulated distribution is drawn into the output map.\cr
# #' If \code{TRUE}, the initial species distribution is drawn additionally.
#' @param SMSHeatMap Produce SMS heat map raster as output? Defaults to \code{FALSE}.
# #' @param ReturnPopRaster Return population data to R (as data frame)? Defaults to \code{TRUE}.
# #' @param CreatePopFile Create population output file? Defaults to \code{TRUE}.
#' @details \emph{Environmental Gradient}\cr
#' In \emph{RangeShiftR}, it is possible to superimpose an artificial gradient on top of the landscape map (real or artificial).
#' Gradients are implemented for cell-based models only.\cr
#' An environmental gradient can be superimposed on the habitat map to describe gradual change in abiotic factors through space. Use the option \code{Gradient}
#' and choose one of three implemented parameter gradients. These are, respectively for non-structured / stage-structured population models:\cr
#' \itemize{
#'    \item{Decreasing values of \code{K_or_DensDep}, that mediates demographic density dependence (stronger with lower values) and is interpreted as
#'          the \emph{carrying capacity} \eqn{K} / the \emph{strength of density dependence} \ifelse{html}{\out{b<sup>-1</sup>}}{\eqn{1/b}} (set \code{Gradient=1});}
#'    \item{Decreasing growth rate \eqn{r} / fecundity (\eqn{\phi}) (set \code{Gradient=2});}
#'    \item {Increasing local extinction probability \eqn{e} (set \code{Gradient=3}).}
#' }
#' The gradient is restrictively implemented along the north-south (\eqn{y})-axis and the selected parameter declines linearly with (\eqn{y})-distance from
#' an optimum location (\code{Optimum}).
#'
#' Gradients are implemented following the method of \insertCite{travis2004;textual}{RangeShiftR}, which combines linear variability with local heterogeneity.
#' If \eqn{Z} is one of the gradient variables listed above, \eqn{Z={K, 1/b, r, \phi, e}}, the value of \eqn{Z(x,y)} for a cell with \eqn{x} and \eqn{y}-coordinates
#' is given by the following equation:
#'
#' \ifelse{html}{\out{&nbsp; &nbsp; Z(x,y) = Z<sub>0</sub> * z(x,y) &nbsp; &nbsp; &nbsp; for K and r,}}{\deqn{Z(x,y)=Z_0 * z(x,y)} for K and r,}
#' or\cr
#' \ifelse{html}{\out{&nbsp; &nbsp; Z(x,y) = 1 - z(x,y) + e<sub>opt</sub> &nbsp; &nbsp; for e}}{\deqn{Z(x,y)= 1 - z(x,y) + e_{opt}} for e}
#' \cr
#'
#' with \cr
#' \ifelse{html}{\out{&nbsp; &nbsp; z(x,y) = 1 - |y - y<sub>opt</sub>| G + U(-1,1) f }}{\deqn{z(x,y) = 1 - |y - y_{opt}| G + U(-1,1) f }}
#' , constrained to \eqn{\ge 0}
#'
#' where \ifelse{html}{\out{Z<sub>0</sub>}}{Z_0} is the original parameter value at \eqn{(x,y)},
#' \ifelse{html}{\out{e<sub>opt</sub>}}{e_{opt}} (\code{ExtinctOptim}) the minimum extinction probability.
#' The linear variability is specified by
#' \ifelse{html}{\out{y<sub>opt</sub>}}{y_{opt}} (\code{Optimum}), the \eqn{y}-value at which the extremum (i.e. \ifelse{html}{\out{Z<sub>0</sub>}}{Z_0} or \ifelse{html}{\out{e<sub>opt</sub>}}{e_{opt}}) is obtained,
#' and \eqn{G} (\code{GradSteep}), the gradient steepness in units of fraction of the local value per cell.
#' The local heterogeneity is determined by a random number drawn from a uniform distribution between \eqn{-1} and \eqn{1} for each cell
#' and \code{f}, the local scaling factor that determines the magnitude of this stochastic local variation relative to the extremal value.
#'
#' The gradient in fecundity φ applies to the fecundity of each stage. Negative local values in \eqn{z(x,y)} are set to \eqn{0}.
#'
#' It is also possible to simulate the shifting of the gradient by setting the option \code{Shifting}. Here the position \eqn{y} of the species’
#' optimum is shifted northwards (increasing \eqn{y}) at a given rate \code{ShiftRate} (in units of rows per year),
#' starting from year \code{ShiftStart} to year \code{ShiftEnd}.
#'
#' Environmental gradients are available for cell-based models only, due to the cell-based character of operations and therefore \code{Gradient}
#' has to be \eqn{0} for patch-based models.
#'
#' \emph{Local Extinctions}\cr
#' An additional, independent extinction probability can be added using the option \code{LocalExt}. If set,
#' in each year, every population has an identical probability \code{LocalExtProb} of going extinct.
#' This does not affect any demographic parameters but simply kills off the local population.
#'
#' \emph{Environmental Stochasticity}\cr
#' It is possible to model environmental stochasticity via the option \code{EnvStoch} acting at a global or local scale
#' and can be applied to \code{K_or_DensDep}, the demographic density dependence (\code{EnvStoch=1}), or to growth rate / fecundity (\code{EnvStoch=0}).
#' It is implemented using a first order autoregressive process to generate time series of the noise value \eqn{ε}
#' \insertCite{ruokolainen2009}{RangeShiftR}:
#'
#' \deqn{ε(t+1) = κ ε(t) + \omega(t) \sqrt(1-\kappa^2)}
#'
#' where κ is the autocorrelation coefficient (\code{ac}) and ω is a random normal variable drawn from \eqn{N(0,σ)}.
#' Changing σ (\code{std}) changes the magnitude of the fluctuations. The spatial scale of the variation can either be global (a single time series
#' for the entire landscape) or local (each cell fluctuates independently), and is always applied on a yearly basis.
#' Different degrees of spatial autocorrelation are not implemented in the current version.
#'
#' The noise affects the species' selected parameter \eqn{Z} as follows:
#' \deqn{Z(x,y,t) = Z(x,y,0) + Z ε(t)}
#' where \eqn{x} and \eqn{y} are the cell coordinates and \eqn{Z} is the original parameter value in absence of stochasticity and gradients.
#' In the presence of an environmental gradient, \eqn{Z(x,y,0)} is the gradient value at the cell location, otherwise its equal to \eqn{Z}.
#' The resulting values \eqn{Z(x,y,t)} are limited to the maximum and minimum values \code{minR,maxR} or \code{minK,maxK}, respectively.
#'
#' \emph{Output files}\cr
#' Seven different types of outputs can be produced, plus one more for patch-based models.
#' All the output files will be named with a standard name reporting the simulation ID number and
#' the type of output. The file name will start with the batch number, and also indicate the number
#' of the landscape to which the output refers. Additionally, for each simulation all the set parameters
#' will be automatically written to a text file named \"Sim0_Parameters.txt\" in the case of simulation\eqn{#=0}.
#'
#' - \emph{Species range} (\code{Sim0_Range.txt}) \cr
#' contains the following general information regarding the species’ range:\cr
#' Replicate number (Rep), Year (Year), Reproductive season within the year (RepSeason), Total number of individuals (NInds),
#' Total number of individuals in each stage (NInd_stageX; only in case of stage-structured models),
#' Total number of juveniles born (NJuvs; only in case of stage-structured models),
#' Total number of occupied cells (NOccupCells) or total number of occupied patches (NOccupPatches),
#' Ratio between occupied and suitable cells or patches (OccupSuit),
#' Species’ range, in terms of maximum and minimum coordinates (min_X, max_X, min_Y, max_Y).\cr
#' Data are written before reproduction at each reproductive season at the specified yearly interval. An extra line is written
#' at the end of the simulation.
#'
#' - \emph{Occupancy} \cr
#' reports the cell/patch probability of occupancy. This is only possible if the number of replicates is greater than \eqn{1}.
#' Two files will be produced:\cr
#'
#'      1) \code{Sim0_Occupancy.txt}:  contains a list of all the cells in the landscape (\eqn{x-} and \eqn{y-}coordinates) or
#' of all the patches (PatchID). The remaining columns give the occupancy probability of the cell/patch at defined time steps.
#' The occupancy probability is obtained by dividing the number of times (replicates) that the cell/patch has been occupied in
#' a given year, by the total number of replicates.
#'
#'      2) \code{Sim0_Occupancy_Stats.txt}: Summary occupancy statistics, i.e. the mean ratio between occupied and suitable cells
#' (Mean_OccupSuit) and its standard error (Std_error) at the set time interval.
#'
#'     Data will be recorded at the beginning of the year before any other process (and only once a year no matter the number
#' of reproductive seasons per year).
#'
#' - \emph{Populations} (\code{Sim0_Pop.txt}) \cr
#' contains statistics regarding each population present in the landscape at a given time interval:\cr
#' Replicate number (Rep), Year (Year), Reproductive season within the year (RepSeason), Cell location (\eqn{x-} and \eqn{y-}coordinates) or
#' patch ID (PatchID), Species number (Species; not yet used, always \eqn{0}), Number of individuals in the population (NInd),
#' Number of individuals in each stage (NInd_stageX; only in case of stage-structured models). If the reproduction is sexual,
#' these columns will be replaced by the number of females (Nfemales_stageX) and of males (Nmales_stageX) in each stage. In the case
#' of sexual model without stage structure, two columns will indicate the number of females (Nfemales) and of males (Nmales) in
#' the population. In the case of a stage-structured population, the number of juveniles born (NJuvs). If the reproduction is sexual,
#' these columns will be replaced by the number of females juveniles (NJuvFemales) and males (NJuvMales).\cr
#' As for the species’ range output, data are collected before reproduction at each reproductive
#' season at the specified yearly interval and at the end of the simulation.
#'
#' - \emph{Individuals} (\code{Sim0_Rep0_Inds.txt}) \cr
#' contains information regarding each individual at a given time step. To avoid the production of huge files, a separate file is
#' saved for each replicate. Data are recorded after settlement and before aging (in the case of overlapping generations). For each
#' individual the following data are saved: \cr
#' Replicate number (Rep), Year, Reproductive season within the year (RepSeason), Species ID (always \eqn{0}),
#' Individual ID (IndID), the individual’s Status (Status), Natal cell (Natal_X and Natal_Y) and current cell (\eqn{x} and
#' \eqn{y}) coordinates or natal and current patch IDs (Natal_patch and PatchID), Sex (0 = female, 1 = male), Age in years
#' (in case of overlapping generations), Stage (in case of stage structure), Emigration traits, Transfer traits (depending on
#' transfer method).
#'
#' - \emph{Genetics} (\code{Sim0_Rep0_Genetics.txt}) \cr
#' lists the full genome of each individual selected for output (i.e. all individuals if the population is not structured) during the reporting year
#' (or present in the initial population at year \eqn{0}) for the current replicate. This file can therefore be \emph{extremely large}, and should be
#' produced only for temporally short simulations, small populations or at infrequent reporting time intervals. It comprises:\cr
#'     - Replicate number (Rep), Year, Species ID (always \eqn{0}), Individual ID (IndID), \cr
#' and then \emph{either} one or more lines listing:
#'     - Chromosome number (starting from 0), Locus on this chromosome (starting from 0), value of the only allele at the locus for a
#'     haploid species (Allele0) or the values of both alleles at the locus for a diploid species (Allele0,Allele1) \cr
#' \emph{or} a single line of:
#'     - a set of columns having compound headings of the form \code{Chr0Loc0Allele0} derived from each chromosome, locus and allele (as above).
#'
#' - \emph{Traits} \cr
#' In the case of inter-individual variability and evolution of the dispersal traits, it is possible to output the mean traits of
#' the population. There are two types of traits output:\cr
#'
#'     1) \code{Sim0_TraitsXcell.txt / Sim0_TraitsXpatch.txt} reports mean and standard deviation of the varying/evolving traits for each
#' cell/patch, for each replicate and reproductive season at the set year interval.\cr
#'
#'     2) \code{Sim0_TraitsXrow.txt} mean and standard deviation of the varying/evolving traits computed at the row (\eqn{y}) level,
#' pulling together all the populations occupying cells in \eqn{y}. Values are reported for each replicate and reproductive season
#' at the specified yearly interval. This is particularly useful for analyzing the structuring of traits along latitudinal gradients.
#' It is possible to compute this output only for cell-based models. \cr
#'
#'     Data for these outputs are collected at the same time as for the
#' range and population outputs, i.e. before reproduction at each reproductive season at the set year interval and at the end of the
#' simulation. For sexual models, the standard deviation relates to the variation between all alleles in the local population (which
#' is greater than the variation in phenotypic expression; if the phenotypic s.d. is required, it must be calculated from
#' individual-level output data).
#'
#' - \emph{Connectivity matrix} (\code{Sim0_Connect.txt}) \cr
#' is available for a patch-based model only. It presents counts of the number of individuals successfully dispersing from each patch
#' to each other patch for each year specified by \code{OutIntConn}, starting from \code{OutStartConn}. If there is more than one
#' reproductive season during the year, cumulative year-end totals are reported. Although the file contains the data required for
#' true \eqn{NxN} matrices, the data are presented in list format:\cr
#' Replicate number (Rep), Year (Year), ID number of natal patch (StartPatch), ID number of settlement patch (EndPatch), Number of
#' individuals dispersing from StartPatch to EndPatch (NInds).
#'
#' - \emph{SMS paths} (\code{Sim0_Rep0_SMSpaths.txt}) \cr
#' is available for a model with transfer method SMS only. It lists the cell-based trajectories of all (successfully or unsuccessfully)
#' dispersed individuals from the natal to the final (settlement or fatal) patch for each year specified by \code{OutIntPaths}, starting
#' from \code{OutStartPaths}. The data are presented in list format with the columns:\cr
#' Year (Year), Individual ID (IndID), consecutive step number (Step), coordinates of cell at this step (\eqn{x} and \eqn{y}),
#' status of individual (Status).
#' The status is an integer number that codes for the following possible states:\cr
#'    0 = natal patch,\cr
#'    1 = disperser,\cr
#'    2 = disperser awaiting settlement in possible suitable patch,\cr
#'    3 = waiting between dispersal events,\cr
#'    4 = completed settlement,\cr
#'    5 = completed settlement in a suitable neighbouring cell,\cr
#'    6 = died during transfer by failing to find a suitable patch (includes exceeding maximum number of steps or crossing absorbing boundary),\cr
#'    7 = died during transfer by constant, step-dependent, habitat-dependent or distance-dependent mortality,\cr
#'    8 = failed to survive annual (demographic) mortality,\cr
#'    9 = exceeded maximum age.\cr\cr
#'
#' - \emph{SMS Heat map} (\code{OutputMaps/Sim0_Land0_Rep0_Visits.txt}) \cr
#' When the transfer model is \emph{SMS}, an additional optional output is a series of maps in ASCII raster format, showing how many times each
#' cell has been visited by a dispersing individual across the whole time period of the simulation. These heat maps may be useful, for example,
#' for identifying corridors which are heavily used during the dispersal phase. One raster map is created in the \emph{Output_Maps} folder for
#' each replicate simulation, and is in the same format as the input habitat file.
#'
#' - \emph{Log file} (\code{Batch1_RS_log.csv}) \cr
#' An additional log file will be created automatically. In it is listed the time taken (in seconds) to run the simulation.
#' It may also possibly include error codes, which can occur in rare occasions when the batch input files are in themselves valid,
#' but there is an inconsistency between files or an invalid habitat code or patch number occurs in an input map file.
#' Error codes are listed in the \emph{Batch_error_codes.xlsx} file.
#'
#' @references
#'         \insertAllCited{}
#' @name Simulation
#' @export Simulation
Simulation <- setClass("SimulationParams", slots = c(Simulation = "integer_OR_numeric",
                                                 Replicates = "integer_OR_numeric",
                                                 Years = "integer_OR_numeric",
                                                 Absorbing = "logical",
                                                 Gradient = "integer_OR_numeric",  #Environmental gradient: 0 = none, 1 = K or 1/b, 2 = growth rate (or fecundity), 3 = local extinction probability.
                                                 GradSteep = "numeric",
                                                 Optimum = "numeric",
                                                 f = "numeric",
                                                 ExtinctOptim = "numeric",
                                                 Shifting = "logical",
                                                 ShiftRate = "numeric",
                                                 ShiftStart = "integer_OR_numeric",
                                                 ShiftEnd = "integer_OR_numeric",
                                                 LocalExt = "logical",
                                                 LocalExtProb = "numeric",
                                                 EnvStoch = "integer_OR_numeric", #Environmental stochasticity: 0 = none, 1 = global, 2 = local
                                                 EnvStochType = "integer_OR_numeric",   #Environmental stochasticity type: FALSE = in growth rate, TRUE = in K or 1/b
                                                 std = "numeric",
                                                 ac = "numeric",
                                                 minR = "numeric",
                                                 maxR = "numeric",
                                                 minK = "numeric",
                                                 maxK = "numeric",
                                                 OutIntRange = "integer_OR_numeric",
                                                 OutIntOcc = "integer_OR_numeric",
                                                 OutIntPop = "integer_OR_numeric",
                                                 OutIntInd = "integer_OR_numeric",
                                                 OutIntGenetic = "integer_OR_numeric",
                                                 OutGenType = "integer_OR_numeric",  #Output genetics for: 0 = juveniles only, 1 =  all individuals, 2 = adults only
                                                 OutGenCrossTab = "logical",
                                                 OutIntTraitCell = "integer_OR_numeric",
                                                 OutIntTraitRow = "integer_OR_numeric",
                                                 OutIntConn = "integer_OR_numeric",
                                                 OutIntPaths = "integer_OR_numeric",
                                                 OutStartPop = "integer_OR_numeric",
                                                 OutStartInd = "integer_OR_numeric",
                                                 OutStartGenetic = "integer_OR_numeric",
                                                 OutStartTraitCell = "integer_OR_numeric",
                                                 OutStartTraitRow = "integer_OR_numeric",
                                                 OutStartConn = "integer_OR_numeric",
                                                 OutStartPaths = "integer_OR_numeric",
                                                 SaveMaps = "logical",
                                                 MapsInterval = "integer_OR_numeric",
                                                 DrawLoadedSp = "logical",
                                                 SMSHeatMap = "logical",
                                                 ReturnPopRaster = "logical",
                                                 CreatePopFile = "logical"
                                                 #moved! PropMales = "integer_OR_numeric", #move to Demography
                                                 #moved! Harem = "integer_OR_numeric",     #move to Demography
                                                 #moved! bc = "integer_OR_numeric",        #move to Demography - determines density dependence
                                                 #moved! Rmax = "integer_OR_numeric",      #move to Demography -> dem.lambda
                                                 #moved! K = "integer_OR_numeric",         #move to Land
                                                 )
                       , prototype = list(Simulation = 1L,
                                          Replicates = 2L,
                                          Years = 50L,
                                          Absorbing = FALSE,
                                          Gradient = 0L,
                                          #GradSteep,
                                          #Optimum,
                                          #f,
                                          #ExtinctOptim,
                                          Shifting = FALSE,
                                          #ShiftRate,
                                          #ShiftStart,
                                          #ShiftEnd,
                                          LocalExt = FALSE,
                                          #LocalExtProb,
                                          EnvStoch = 0L,
                                          #EnvStochType,
                                          #std,
                                          #ac,
                                          #minR,
                                          #maxR,
                                          #minK,
                                          #maxK,
                                          OutIntRange = 1L,
                                          OutIntOcc = 0L,
                                          OutIntPop = 1L,
                                          OutIntInd = 0L,
                                          OutIntGenetic = 0L,
                                          OutGenType = 0L,
                                          OutGenCrossTab = FALSE,
                                          OutIntTraitCell = 0L,
                                          OutIntTraitRow = 0L,
                                          OutIntConn = 0L,
                                          OutIntPaths = 0L,
                                          OutStartPop = 0L,
                                          OutStartInd = 0L,
                                          OutStartGenetic = 0L,
                                          OutStartTraitCell = 0L,
                                          OutStartTraitRow = 0L,
                                          OutStartConn = 0L,
                                          OutStartPaths = 0L,
                                          SaveMaps = FALSE,
                                          #MapsInterval,
                                          DrawLoadedSp = FALSE,
                                          SMSHeatMap = FALSE,
                                          ReturnPopRaster = FALSE,
                                          CreatePopFile = FALSE
                                          #moved! PropMales,
                                          #moved! Harem,
                                          #moved! bc,
                                          #moved! Rmax,
                                          #moved! K,
                                          )
)
setValidity('SimulationParams', function(object){
    msg <- NULL
    if (is.na(object@Simulation || length(object@Simulation)==0 )){
        msg <- c(msg, 'ID Number of current simulation must be set!')
    }
    else {
        if (object@Simulation<0){
            msg <- c(msg, 'Simulation ID Number must be positive or 0.')
        }
    }
    if (is.na(object@Replicates || length(object@Replicates)==0 )){
        msg <- c(msg, 'Number of replicates must be set!')
    }
    else {
        if (object@Replicates<=0) {
            msg <- c(msg, 'Number of replicates must be positive.')
        }
    }
    if (is.na(object@Years) || length(object@Years)==0 ){
        msg <- c(msg,'Number of years must be set!')
    }
    else {
        if (object@Years<=0){
            msg <- c(msg, 'Number of year must be positive.')
        }
    }
    if (is.na(object@Absorbing) || length(object@Absorbing)==0 ){
        msg <- c(msg, 'Absorbing must be set!')
    }
    if (is.na(object@Gradient) || length(object@Gradient)==0 ){
        msg <- c(msg, 'Gradient option must be set!')
    }
    else{
        if (object@Gradient!=0 && object@Gradient!=1 && object@Gradient!=2 && object@Gradient!=3){
            msg <- c(msg, 'Gradient must be set to 0, 1, 2 or 3!')
        }
        else{
            if (object@Gradient){   # Gradient = {1,2,3}
                if (is.na(object@GradSteep) || length(object@GradSteep)==0 ){
                    msg <- c(msg, 'GradSteep is required if Gradient is > 0.')
                }
                else {
                    if (object@GradSteep<0){
                        msg <- c(msg, 'GradSteep has to be >= 0.')
                    }
                }
                if (is.na(object@Optimum) || length(object@Optimum)==0 ){
                    msg <-  c(msg, 'Optimum is required if Gradient is > 0.')
                }
                else {
                    if (object@Optimum<0){
                        msg <- c(msg, 'Optimum has to be >= 0.')
                    }
                }
                if (is.na(object@f) || length(object@f)==0){
                    msg <- c(msg, 'Local scaling factor f is required if Gradient is > 0.')
                }
                else {
                    if (object@f<0){
                        msg <- c(msg, 'Local scaling factor f has to be >= 0.')
                    }
                }
                if (object@Gradient == 3){
                    if (is.na(object@ExtinctOptim) || length(object@ExtinctOptim)==0){
                        msg <- c(msg, 'ExtinctOptim has to be set.')
                    }
                    else{
                        if (object@ExtinctOptim<0 || object@ExtinctOptim>=1){
                            msg <- c(msg, 'Value of ExtinctOptim mustlie in the half-open interval [0,1).')
                        }
                    }
                }
                if (is.na(object@Shifting) || length(object@Shifting)==0 ){
                    msg <- c(msg, 'Shifting must be set.')
                }
                else{
                    if (object@Shifting){
                        if (is.na(object@ShiftRate) || length(object@ShiftRate)==0 ) {
                            msg <- c(msg, 'ShiftRate must be set.')
                        }
                        else{
                            if (object@ShiftRate <= 0){
                                msg <- c(msg, 'ShiftRate has to be > 0, if Shifting is = TRUE.')
                            }
                        }
                        if (is.na(object@ShiftStart)|| length(object@ShiftStart)==0 ){
                            msg <- c(msg, 'ShiftStart must be set.')
                        }
                        else{
                            if (object@ShiftStart <= 0){
                                msg <- c(msg, 'ShiftStart has to be > 0, if Shifting is = TRUE.')
                            }
                        }
                        if (is.na(object@ShiftEnd)|| length(object@ShiftEnd)==0 ){
                            msg <- c(msg, 'ShiftEnd has to be set.')
                        }
                        else{
                            if(object@ShiftEnd <= 0)
                                msg <- c(msg, 'ShiftEnd has to be > 0, if Shifting is = TRUE.')
                        }
                        if (is.null(msg)) {
                            if (object@ShiftEnd <= object@ShiftStart ){
                                msg <- c(msg, 'ShiftEnd must be greater than ShiftStart.')
                            }
                        }
                    }
                }
            }
            else { # no gradient
                if (is.na(object@Shifting) || length(object@Shifting)==0 ){
                    msg <- c(msg, 'Shifting must be set.')
                }
                else{
                    if (object@Shifting){
                        msg <- c(msg, 'Shifting is not applicable if there is no gradient! (Gradient = 0).')

                    }
                }
            }
        }
    }
    if (is.na(object@LocalExt || length(object@LocalExt)==0)){
        msg <- c(msg, 'LocalExt must be set!')
    }
    else{
        if (object@LocalExt){ # LocalExt = TRUE
            if (object@Gradient == 3){
                msg <- c(msg, 'LocalExt has to be FALSE if Gradient = 3 (i.e. environmental gradient in extinction probability)')
            }
            else {
                if (is.na(object@LocalExtProb) || length(object@LocalExtProb)==0 ){
                    msg <- c(msg, 'LocalExtProb has to be set.')
                }
                else{
                    if (object@LocalExtProb<=0 || object@LocalExtProb>=1){
                        msg <- c(msg, 'Value of LocalExtProb must be between 0 and 1.')
                    }
                }
            }
        }
    }
    if (is.na(object@EnvStoch) || length(object@EnvStoch)==0){
        msg <- c(msg, 'Environmental stochasticity has to be set.')
    }
    else{
        if (object@EnvStoch!=0 & object@EnvStoch!=1 & object@EnvStoch!=2) {
            msg <- c(msg, 'Environmental stochasticity option (EnvStoch) has to be 0, 1 or 2.')
        }
        else{
            if (object@EnvStoch){ #EnvStoch = 1 or 2
                if (is.na(object@ac) || length(object@ac)==0 ){
                    msg <- c(msg, 'Autocorrelation coefficient (ac) has to be set if environmental stochasticity is set.')
                }
                else {
                    if (object@ac < 0.0 || object@ac >= 1.0){
                        msg <- c(msg, 'Autocorrelation coefficient (ac) must be in the half-open interval [0,1).')
                    }
                }
                if (is.na(object@std) || length(object@std)==0 ){
                    msg <- c(msg, 'std has to be set if environmental stochasticity is set.')
                }
                else {
                    if (object@std <= 0.0 || object@std > 1.0){
                        msg <- c(msg, 'std must be in the half-open interval (0,1].')
                    }
                }
                if (is.na(object@EnvStochType) || length(object@EnvStochType)==0 ){
                    msg <- c(msg, 'Type of environmental stochasticity (EnvStochType) must be set.')
                }
                else {
                    if (object@EnvStochType != 1 & object@EnvStochType != 0 ){
                        msg <- c(msg, 'Type of environmental stochasticity (EnvStochType) must be 0 or 1.')
                    }
                    else{
                        if (object@EnvStochType == 0){
                            if (is.na(object@minR) || length(object@minR)==0){
                                msg <- c(msg, 'Minimum growth rate (minR) has to be set.')
                            }
                            else {
                                if (object@minR <= 0){
                                    msg <- c(msg, 'Minimum growth rate (minR) must be positive.')
                                }
                                else {
                                    if (is.na(object@maxR) || length(object@maxR)==0){
                                        msg <- c(msg, 'Maximum growth rate (maxR) has to be set.')
                                    }
                                    else {
                                        if (object@maxR <= object@minR){
                                            msg <- c(msg, 'Maximum growth rate (maxR) must be greater than minR.')
                                        }
                                    }
                                }
                            }
                        }
                        if (object@EnvStochType == 1){
                            if (is.na(object@minK) || length(object@minK)==0){
                                msg <- c(msg, 'Minimum growth rate (minK) has to be set.')
                            }
                            else {
                                if (object@minK <= 0){
                                    msg <- c(msg, 'Minimum growth rate (minK) must be positive.')
                                }
                                else {
                                    if (is.na(object@maxK) || length(object@maxK)==0){
                                        msg <- c(msg, 'Maximum growth rate (maxK) has to be set.')
                                    }
                                    else {
                                        if (object@maxK < 0){
                                            msg <- c(msg, 'Maximum growth rate (maxK) must be greater than minK.')
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    # Range
    if (is.na(object@OutIntRange) || length(object@OutIntRange)==0 ){
        msg <- c(msg, 'Output interval of range (OutIntRange) has to be set.')
    }
    else{
        if (object@OutIntRange < 0){
            msg <- c(msg, 'Output interval of range (OutIntRange) must be positive or zero.')
        }
    }
    # Occupancy
    if (is.na(object@OutIntOcc) || length(object@OutIntOcc)==0 ){
        msg <- c(msg, 'Output interval of occupancy (OutIntOcc) has to be set.')
    }
    else{
        if (object@OutIntOcc > 0){
            if (object@Replicates < 2){
                msg <- c(msg, 'For Occupancy output, Replicates must be at least 2.')
            }
        }
        else{
            if (object@OutIntOcc){
                msg <- c(msg, 'Output interval of occupancy (OutIntOcc) must be positive or zero.')
            }
        }
    }
    # Population
    if (is.na(object@OutIntPop) || length(object@OutIntPop)==0 ){
        msg <- c(msg, 'Output interval of population (OutIntPop) has to be set.')
    }
    else {
        if (object@OutIntPop > 0){
            if (is.na(object@OutStartPop) || length(object@OutStartPop)==0 ){
                msg <- c(msg, 'Start year of population output (OutStartPop) has to be set.')
            }
            else{
                if (object@OutStartPop < 0 || object@OutStartPop > object@Years){
                    msg <- c(msg, 'Invalid value of output parameter OutStartPop: Value has to be positive and less than simulated Years')
                }
            }
        }
        else {
            if (object@OutIntPop){
                msg <- c(msg, 'Output interval of population (OutIntPop) must be positive or zero.')
            }
        }
    }
    # Individuals
    if (is.na(object@OutIntInd) || length(object@OutIntInd)==0 ){
        msg <- c(msg, 'Output interval of individuals (OutIntInd) has to be set.')
    }
    else {
        if (object@OutIntInd > 0){
            if (is.na(object@OutStartInd) || length(object@OutStartInd)==0 ){
                msg <- c(msg, 'Start year of individuals output (OutStartInd) has to be set.')
            }
            else{
                if (object@OutStartInd < 0 || object@OutStartInd > object@Years){
                    msg <- c(msg, 'Invalid value of output parameter OutStartInd: Value has to be positive and less than simulated Years')
                }
            }
        }
        else {
            if (object@OutIntInd){
                msg <- c(msg, 'Output interval of individuals (OutIntInd) must be positive or zero.')
            }
        }
    }
    # Genetics
    if (is.na(object@OutIntGenetic) || length(object@OutIntGenetic)==0 ){
        msg <- c(msg, 'Output interval of genetics (OutIntGenetic) has to be set.')
    }
    else {
        if (object@OutIntGenetic > 0){
            if (is.na(object@OutStartGenetic) || length(object@OutStartGenetic)==0 ){
                msg <- c(msg, 'Start year of genetics output (OutStartGenetic) has to be set.')
            }
            else{
                if (object@OutStartGenetic < 0 || object@OutStartGenetic > object@Years){
                    msg <- c(msg, 'Invalid value of output parameter OutStartGenetic: Value has to be positive and less than simulated Years')
                }
            }
            if (is.na(object@OutGenType) || length(object@OutGenType)==0 ){
                msg <- c(msg, 'Type of genetics output (OutGenType) has to be specified.')
            }
            else{
                if (object@OutGenType != 0 && object@OutGenType != 1 && object@OutGenType != 2) {
                    msg <- c(msg, 'OutGenType has to be 0 (juveniles only), 1 (all individuals) or 2 (adults only).')
                }
            }
            if (is.na(object@OutGenCrossTab) || length(object@OutGenCrossTab)==0 ){
                msg <- c(msg, 'OutGenCrossTab has to be set')
            }
        }
        else {
            if (object@OutIntGenetic){
                msg <- c(msg, 'Output interval of genetics (OutIntGenetic) must be positive or zero.')
            }
        }
    }
    # TraitCell
    if (is.na(object@OutIntTraitCell) || length(object@OutIntTraitCell)==0 ){
        msg <- c(msg, 'Output interval of traits per cell (OutIntTraitCell) has to be set.')
    }
    else {
        if (object@OutIntTraitCell > 0){
            if (is.na(object@OutStartTraitCell) || length(object@OutStartTraitCell)==0){
                msg <- c(msg, 'Start year of traits (per cell) output (OutStartTraitCell) has to be set.')
            }
            else {
                if (object@OutStartTraitCell < 0 || object@OutStartTraitCell > object@Years){
                    msg <- c(msg, 'Invalid value of output parameter OutStartTraitCell: Value has to be positive and less than simulated Years')
                }
            }
        }
        else{
            if(object@OutIntTraitCell){
                msg <- c(msg, 'Output interval of traits per cell (OutIntTraitCell) must be positive or zero.')
            }
        }
    }
    # TraitRow
    if (is.na(object@OutIntTraitRow) || length(object@OutIntTraitRow)==0 ){
        msg <- c(msg, 'Output interval of traits per row (OutIntTraitRow) has to be set.')
    }
    else {
        if (object@OutIntTraitRow > 0){
            if (is.na(object@OutStartTraitRow) || length(object@OutStartTraitRow)==0){
                msg <- c(msg, 'Start year of traits (per row) output (OutStartTraitRow) has to be set.')
            }
            else {
                if (object@OutStartTraitRow < 0 || object@OutStartTraitRow > object@Years){
                    msg <- c(msg, 'Invalid value of output parameter OutStartTraitRow: Value has to be positive and less than simulated Years')
                }
            }
        }
        else {
            if(object@OutIntTraitRow){
                msg <- c(msg, 'Interval of traits per row output (OutIntTraitRow) must be positive or zero.')
            }
        }
    }
    # Connectivity matrix
    if (is.na(object@OutIntConn || length(object@OutIntConn)==0 )){
        msg <- c(msg, 'Output interval of connectivity matrix (OutIntConn) has to be set.')
    }
    else {
        if (object@OutIntConn > 0){
            if (is.na(object@OutStartConn) || length(object@OutStartConn)==0 ){
                msg <- c(msg, 'Start year of connectivity matrix output (OutStartConn) has to be set.')
            }
            else {
                if (object@OutStartConn < 0 || object@OutStartConn > object@Years){
                    msg <- c(msg, 'Invalid value of output parameter OutStartConn: Value has to be positive and less than simulated Years')
                }
            }
        }
        else {
            if(object@OutIntConn){
                msg <- c(msg, 'Interval of connectivity output (OutIntConn) must be positive or zero.')
            }
        }
    }
    # Paths record
    if (is.na(object@OutIntPaths || length(object@OutIntPaths)==0 )){
        msg <- c(msg, 'Output interval of SMS paths (OutIntPaths) has to be set.')
    }
    else {
        if (object@OutIntPaths > 0){
            if (is.na(object@OutStartPaths) || length(object@OutStartPaths)==0 ){
                msg <- c(msg, 'Start year of SMS paths output (OutStartPaths) has to be set.')
            }
            else {
                if (object@OutStartPaths < 0 || object@OutStartPaths > object@Years){
                    msg <- c(msg, 'Invalid value of output parameter OutStartPaths: Value has to be positive and less than simulated Years')
                }
            }
        }
        else {
            if(object@OutIntPaths){
                msg <- c(msg, 'Interval of SMS paths output (OutIntConn) must be positive or zero.')
            }
        }
    }
    # Maps
    if (is.na(object@SaveMaps) || length(object@SaveMaps)==0 ){
        msg <- c(msg, 'SaveMaps has to be set')
    }
    else {
        if (object@SaveMaps){ # TRUE
            if (is.na(object@MapsInterval) || length(object@MapsInterval)==0 ){
                msg <- c(msg, 'MapsInterval has to be set.')
            }
            else{
                if(object@MapsInterval < 1){
                    msg <- c(msg, 'MapsInterval must be positive.')
                }
            }
            if (is.na(object@DrawLoadedSp) || length(object@DrawLoadedSp)==0 ){
                msg <- c(msg, 'DrawLoadedSp has to be set.')
            }
        }
    }
    if (is.na(object@SMSHeatMap) || length(object@SMSHeatMap)==0 ){
        msg <- c(msg, 'SMSHeatMap has to be set')
    }
    # R Output
    if (is.na(object@ReturnPopRaster) || length(object@ReturnPopRaster)==0 ){
        msg <- c(msg, 'ReturnPopRaster has to be set')
    }
    if (is.na(object@CreatePopFile) || length(object@CreatePopFile)==0 ){
        msg <- c(msg, 'CreatePopFile has to be set')
    }
    if (is.null(msg)) TRUE else msg
})

setMethod("initialize", "SimulationParams", function(.Object, ...) {
    this_func = "Simulation(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if(.Object@Gradient == 0){
        .Object@GradSteep = -9L
        if (!is.null(args$GradSteep)) {
            warning(this_func, 'GradSteep', warn_msg_ignored, 'since Gradient is 0.', call. = FALSE)
        }
        .Object@Optimum = -9L
        if (!is.null(args$Optimum)) {
            warning(this_func, 'Optimum', warn_msg_ignored, 'since Gradient is 0.', call. = FALSE)
        }
        .Object@f = -9L
        if (!is.null(args$f)) {
            warning(this_func, 'f', warn_msg_ignored, 'since Gradient is 0.', call. = FALSE)
        }
    }
    if (.Object@Gradient != 3){
        .Object@ExtinctOptim = -9L
        if (!is.null(args$ExtinctOptim)) {
            warning(this_func, 'ExtinctOptim is only used if Gradient is 3 and thus ', warn_msg_ignored , call. = FALSE)
        }
    }
    if (!.Object@Shifting){
        .Object@ShiftRate = -9L
        if (!is.null(args$ShiftRate)) {
            warning(this_func, 'ShiftRate', warn_msg_ignored, 'since Shifting is = FALSE', call. = FALSE)
        }
        .Object@ShiftStart = -9L
        if (!is.null(args$ShiftStart)) {
            warning(this_func, 'ShiftStart', warn_msg_ignored, 'since Shifting is = FALSE', call. = FALSE)
        }
        .Object@ShiftEnd = -9L
        if (!is.null(args$ShiftEnd)) {
            warning(this_func, 'ShiftEnd', warn_msg_ignored, 'since Shifting is = FALSE', call. = FALSE)
        }
    }
    if (!.Object@LocalExt){
        .Object@LocalExtProb = -9L
        if (!is.null(args$LocalExtProb)) {
            warning(this_func, 'LocalExtProb', warn_msg_ignored, 'since LocalExt is = FALSE', call. = FALSE)
        }
    }
    if (.Object@EnvStoch){ # EnvStoch = {1,2}
        if (.Object@EnvStochType == 1) { # EnvStoch in K
            .Object@minR = -9L
            if (!is.null(args$minR)) {
                warning(this_func, 'minR', warn_msg_ignored, 'since EnvStochType is 1.', call. = FALSE)
            }
            .Object@maxR = -9L
            if (!is.null(args$maxR)) {
                warning(this_func, 'maxR', warn_msg_ignored, 'since EnvStochType is 1.', call. = FALSE)
            }
        }
        if (.Object@EnvStochType == 0){ # EnvStoch in r
            .Object@minK = -9L
            if (!is.null(args$minK)) {
                warning(this_func, 'minK', warn_msg_ignored, 'since EnvStochType is 0.', call. = FALSE)
            }
            .Object@maxK = -9L
            if (!is.null(args$maxK)) {
                warning(this_func, 'maxK', warn_msg_ignored, 'since EnvStochType is 0.', call. = FALSE)
            }
        }
    }
    else{                   # EnvStoch = 0
        .Object@EnvStochType = -9L
        if (!is.null(args$EnvStochType)) {
            warning(this_func, 'EnvStochType', warn_msg_ignored, 'since EnvStoch is 0.', call. = FALSE)
        }
        .Object@ac = -9L
        if (!is.null(args$ac)) {
            warning(this_func, 'ac', warn_msg_ignored, 'since EnvStoch is 0.', call. = FALSE)
        }
        .Object@std = -9L
        if (!is.null(args$std)) {
            warning(this_func, 'std', warn_msg_ignored, 'since EnvStoch is 0.', call. = FALSE)
        }
        .Object@minR = -9L
        if (!is.null(args$minR)) {
            warning(this_func, 'minR', warn_msg_ignored, 'since EnvStoch is 0.', call. = FALSE)
        }
        .Object@maxR = -9L
        if (!is.null(args$maxR)) {
            warning(this_func, 'maxR', warn_msg_ignored, 'since EnvStoch is 0.', call. = FALSE)
        }
        .Object@minK = -9L
        if (!is.null(args$minK)) {
            warning(this_func, 'minK', warn_msg_ignored, 'since EnvStoch is 0.', call. = FALSE)
        }
        .Object@maxK = -9L
        if (!is.null(args$maxK)) {
            warning(this_func, 'maxK', warn_msg_ignored, 'since EnvStoch is 0.', call. = FALSE)
        }
    }
    if (!.Object@OutIntPop){
        .Object@OutStartPop = 0L
        if(!is.null(args$OutStartPop)){
            warning(this_func, 'OutStartPop', warn_msg_ignored, 'since OutIntPop is zero (output disabled).', call. = FALSE)
        }
    }
    if (!.Object@OutIntInd){
        .Object@OutStartInd = 0L
        if(!is.null(args$OutStartInd)){
            warning(this_func, 'OutStartInd', warn_msg_ignored, 'since OutIntInd is zero (output disabled).', call. = FALSE)
        }
    }
    if (!.Object@OutIntTraitCell){
        .Object@OutStartTraitCell = 0L
        if(!is.null(args$OutStartTraitCell)){
            warning(this_func, 'OutStartTraitCell', warn_msg_ignored, 'since OutIntTraitCell is zero (output disabled).', call. = FALSE)
        }
    }
    if (!.Object@OutIntTraitRow) {
        .Object@OutStartTraitRow = 0L
        if(!is.null(args$OutStartTraitRow)){
            warning(this_func, 'OutStartTraitRow', warn_msg_ignored, 'since OutIntTraitRow is zero (output disabled).', call. = FALSE)
        }
    }
    if (!.Object@OutIntConn){
        .Object@OutStartConn = 0L
        if(!is.null(args$OutStartConn)){
            warning(this_func, 'OutStartConn', warn_msg_ignored, 'since OutIntConn is zero (output disabled).', call. = FALSE)
        }
    }
    if (!.Object@OutIntPaths){
        .Object@OutStartPaths = 0L
        if(!is.null(args$OutStartPaths)){
            warning(this_func, 'OutStartPaths', warn_msg_ignored, 'since OutIntPaths is zero (output disabled).', call. = FALSE)
        }
    }
    if (!.Object@OutIntGenetic){
        .Object@OutStartGenetic = 0L
        if(!is.null(args$OutStartGenetic)){
            warning(this_func, 'OutStartGenetic', warn_msg_ignored, 'since OutIntGenetic is zero (output disabled)', call. = FALSE)
        }
        .Object@OutGenType = 0L
        if(!is.null(args$OutGenType)){
            warning(this_func, 'OutGenType', warn_msg_ignored, 'since OutIntGenetic is zero (output disabled)', call. = FALSE)
        }
        .Object@OutGenCrossTab = FALSE
        if(!is.null(args$OutGenCrossTab)){
            warning(this_func, 'OutGenCrossTab', warn_msg_ignored, 'since OutIntGenetic is zero (output disabled)', call. = FALSE)
        }
    }
    if (!.Object@SaveMaps){
        .Object@MapsInterval = -9L
        if(!is.null(args$MapsInterval)){
            warning(this_func, 'MapsInterval', warn_msg_ignored, 'since SaveMaps = FALSE.', call. = FALSE)
        }
        if(!is.null(args$DrawLoadedSp)){
            warning(this_func, 'DrawLoadedSp', warn_msg_ignored, 'since SaveMaps = FALSE.', call. = FALSE)
        }
    }
    return(.Object)
})
setMethod("show", "SimulationParams", function(object) {
    cat(" Simulation #", object@Simulation, "\n")
    cat(" -----------------\n")
    cat("   Replicates = ", object@Replicates, "\n")
    cat("   Years      = ", object@Years, "\n")
    cat("   Absorbing  = ", object@Absorbing, "\n")
    if (object@Gradient) {
        if (object@Shifting) cat(" Shifting Environmental Gradient in")
        else cat(" Environmental Gradient in")
        if (object@Gradient==1) cat(" K_or_DensDep:\n")
        if (object@Gradient==2) cat(" r/phi:\n")
        if (object@Gradient==3) cat(" e:\n")
        cat("   G = ", object@GradSteep, ", y_opt = ", object@Optimum , sep = "")
        if (object@Gradient==3) cat(", e_opt =", object@ExtinctOptim, "\n") else cat("\n")
        cat("   f =", object@f, "\n")
    }
    if (object@Shifting) {
        cat("   ShiftRate =", object@ShiftRate, "rows per year; from year", object@ShiftStart, "to", object@ShiftEnd, "\n")
    }
    if (object@LocalExt) {
        cat("  Local Extinction probalitity =", object@LocalExtProb, "\n")
    }
    if (object@EnvStoch) {
        if (object@EnvStoch==1) cat("  Global")
        if (object@EnvStoch==2) cat("  Local")
        cat(" Environmental Stochasticity in:")
        if (object@EnvStochType) cat(" K_or_DensDep\n") else cat(" r/phi \n")
        cat(" Std = ", object@std, "\n")
        cat(" ac = ", object@ac, "\n Min/Max limits = ")
        if (object@EnvStochType) cat(object@minK, object@maxK, "\n") else cat(object@minR, object@maxR, "\n")
    }
    cat(" File Outputs:\n")
    if (object@OutIntRange) {
        cat("   Range,       every", object@OutIntRange, "years\n")
    }
    if (object@OutIntOcc) {
        cat("   Occupancy,   every", object@OutIntOcc, "years\n")
    }
    if (object@OutIntPop) {
        cat("   Populations, every", object@OutIntPop, "years, starting year", object@OutStartPop)
        if (object@ReturnPopRaster) {
            cat("\n   (R output: RasterStack)")
        }
        cat("\n")
    }
    if (object@OutIntInd) {
        cat("   Individuals, every", object@OutIntInd, "years, starting year", object@OutStartInd, "\n")
    }
    if (object@OutIntTraitCell) {
        cat("   Traits/cell, every", object@OutIntTraitCell, "years, starting year", object@OutStartTraitCell, "\n")
    }
    if (object@OutIntTraitRow) {
        cat("   Traits/row,  every", object@OutIntTraitRow, "years, starting year", object@OutStartTraitRow, "\n")
    }
    if (object@OutIntConn) {
        cat("   Connectivity,  every", object@OutIntConn, "years, starting year", object@OutStartConn, "\n")
    }
    if (object@OutIntPaths) {
        cat("   SMS paths,  every", object@OutIntPaths, "years, starting year", object@OutStartPaths, "\n")
    }
    if (object@OutIntGenetic) {
        cat("   Genetics")
        if (object@OutGenCrossTab) cat(" cross table")
        if (object@OutGenType==0) cat(" of juveniles")
        if (object@OutGenType==1) cat(" of all individuals")
        if (object@OutGenType==2) cat(" of adults")
        cat(",    every", object@OutIntGenetic, "years, starting year", object@OutStartGenetic, "\n")
    }
    if (object@SMSHeatMap) {
        cat("   SMS heat map\n")
    }
})
