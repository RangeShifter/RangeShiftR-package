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
# RangeShiftR parameter master
#
# -----

#' Define a RangeShiftR parameter master object
#'
#' Set up a parameter master that can be used in \code{\link[RangeShiftR]{RunRS}}() to run a simulation.\cr
#' All parameter modules can be added to an existing parameter master via the "+"-functions. However, note that the entire respective module will be overwritten.\cr
#'
#' @usage RSsim(batchnum = 1L,
#'       simul = Simulation(),
#'       land = ArtificialLandscape(),
#'       demog = Demography(Rmax = 1.5),
#'       dispersal = Dispersal(),
#'       gene = Genetics(),
#'       init = Initialise(),
#'       seed = 0)
#' @include class_RSparams.R
#' @param batchnum Batch ID is part of output files names and can be used to prevent overwriting.
#' @param simul Set \code{\link[RangeShiftR]{Simulation}} parameters
#' @param land Set landscape parameters. Can be either \code{\link[RangeShiftR]{ArtificialLandscape}} or \code{\link[RangeShiftR]{ImportedLandscape}}.
#' @param demog Set \code{\link[RangeShiftR]{Demography}} parameters
#' @param dispersal Set \code{\link[RangeShiftR]{Dispersal}} parameters
#' @param gene Set \code{\link[RangeShiftR]{Genetics}} parameters
#' @param init Set \code{\link[RangeShiftR]{Initialise}} parameters
#' @param seed Set seed for random number generator. If non-positive, a random seed will be generated.
#' @return returns a \emph{RangeShiftR} parameter master object (class 'RSparams')
#' @details
#' \emph{Demographic stochasticity} \cr Demographic stochasticity is fundamentally important for the dynamics of populations that are naturally small or have declined to low abundances owing to
#' anthropogenic pressures. Additionally, inter-individual variability within populations can have a major influence on dynamics. Modelling stochastic events
#' that happen to individuals is crucial for avoiding systematic overestimation of population viability or rate of spread
#' \insertCite{clark2001invasion,kendall2003unstructured,robert2003variation,grimm2005individual,jongejans2008dispersal,travis2011improving}{RangeShiftR}.
#' Thus, population dynamics in \emph{RangeShiftR} were constructed to be
#' fully individual-based and stochastic. Each reproductive individual produces a discrete number of offspring sampled from a Poisson distribution with a mean
#' that is influenced by the species’ demographic parameters and the local population density. As \emph{RangeShiftR} has been designed for modelling a variety of
#' species with different life-history traits, a range of different population models can be chosen, depending on the species being modelled and on the
#' available information (see \code{\link[RangeShiftR]{Demography}}). In all cases demographic stochasticity is implemented.
#'
#' \emph{Cell-based vs. patch-based model} \cr
#' \emph{RangeShiftR} can be run as a cell-based or patch-based model \insertCite{bian2003representation}{RangeShiftR}. It should be noted
#' that the selection between cell-based or patch-based model is of fundamental importance for population dynamics calculations because
#' it influences the spatial extent at which density dependence operates. In both cases, the landscape is represented as a grid with cells
#' belonging to a particular habitat type, holding a percentage of habitat cover or being assigned a habitat quality index. However,
#' when \emph{RangeShiftR} is run using the cell-based setting, the cell is the scale at which processes such as population dynamics and dispersal
#' act. The individuals present in a cell define a distinct population, and density-dependencies for reproduction, emigration and settlement
#' all operate at this scale. Even in the case where two habitat cells are adjacent, they still hold separate populations. In contrast, in
#' the patch-based model, population dynamics happen at the patch level, a patch being an assemblage of landscape cells of potentially
#' different habitat types. Patches are not defined automatically by \emph{RangeShiftR} (see \code{\link[RangeShiftR]{ImportedLandscape}}).
#' Rather, the user is required to define which cells belong
#' to which patch, taking into account the ecological understanding of the study species. Density-dependencies regarding reproduction,
#' development, survival, emigration and settlement will depend on the density of individuals in a patch. However, discrete step-wise
#' movements during the transfer phase will always use the cell as the resolution at which steps occur, thus retaining important information
#' about the landscape heterogeneity.\cr
#' The choice between cell- and patch-based modelling can be of crucial importance. While a cell-based model provides an excellent abstraction
#' of space for many theoretical studies, for some applied studies it may be insufficient. This is because the misrepresentation of population
#' dynamics and dispersal (in terms of the scale at which they operate) can lead to substantial biases in projections regarding, for example,
#' rate of range expansion and population persistence \insertCite{bocedi2012projecting}{RangeShiftR}. Ideally, the scales at which population dynamics and dispersal
#' processes are modelled (by choosing the cell resolution or by defining the patches) should be those that are relevant for the species.
#' Importantly, the patch-based implementation allows separating the scales used for population dynamics and movements. In this case, the
#' landscape can be modelled at very fine resolution in order to capture the features that are likely to influence movements (e.g. narrow linear
#' features) without constraining the local population dynamics to operate at too small a scale.
#'
#' \emph{Temporal and spatial scales} \cr
#' It is important to note an essential difference in spatial scale between
#' the cell-based and the patch-based version. In the cell-based model, the cell resolution represents the spatial scale at which the two fundamental
#' processes of population dynamics and dispersal happen. This means that all the density-dependencies in the model (reproduction, survival,
#' emigration, settlement, etc...) act at the cell scale and the same scale is used as a single step unit for discrete movement models. In
#' the patch-based version, two spatial scales are simultaneously present: the cell scale, which in this case is used just for the transfer phase
#' of dispersal (movements) and the patch scale, at which the density-dependences are acting. The choice of type of model and cell resolution (as
#' well as the definition/scale of patches) is of fundamental importance because, depending on the system and on the question being tackled, it can
#' systematically bias the outcomes of the model.\cr
#' The user also defines the temporal scales. There are three distinct temporal scales. The highest-level one has years as units and represents the
#' scale at which variations in the abiotic environment are modelled (\emph{RangeShiftR} does not explicitly model within-year variability in conditions).
#' The intermediate scale is the species’ reproductive season. The model can be used to simulate the case where there is only one reproductive
#' season per year but it is also possible to simulate situations where there more than one per year or only one every \eqn{N} years. A single
#' reproductive event is always followed by dispersal. Finally, the smallest time scale is represented by the number of steps that emigrants take
#' during the movement phase of dispersal. This can be determined by a maximum number of steps, per-step mortality or both.
#'
#' @references
#'         \insertAllCited{}
#' @return if given parameters pass validity check, returns a RangeShiftR parameter master object of class "RSparams", otherwise NULL
#' @author Anne-Kathleen Malchow
#' @export
RSsim <- function(batchnum = 1L,
                  simul = NULL,
                  land = NULL,
                  demog = NULL,
                  dispersal = NULL,
                  gene = NULL,
                  init = NULL,
                  seed = 0L){
    args <- as.list(match.call())
    # filter for names in ... that are also in slots(ControlParams) and pass them on
    s <- RSparams(control = ControlParams(batchnum = batchnum, seed = seed),
                  simul = Simulation(),
                  land = ArtificialLandscape(),
                  demog = Demography(Rmax = 1.5),
                  dispersal = Dispersal(),
                  gene = Genetics(),
                  init = Initialise())
    if (!is.null(args$land))  s <- s + land
    if (!is.null(args$simul)) s <- s + simul
    if (!is.null(args$demog)) s <- s + demog
    if (!is.null(args$dispersal))  s <- s + dispersal
    if (!is.null(args$gene))  s <- s + gene
    if (!is.null(args$init))  s <- s + init
    # check validity
    if(validObject(s)) return(s)
    else return(NULL)
}
