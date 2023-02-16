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



#### DEMOGRAPHY ####


### SUBCLASS STAGESPARAMS

# from RS 'StageStruct' file

#' Define a Stage-structure
#'
#' Set all demographic parameters for a stage-structured population. All elements of the transition matrix, i.e. Fecundity, Development, and Survival can be subjected to a density-dependent scaling.
#' In this case, the base probabilities given via the transition matrix will be reduced exponentially with the strength of density dependence \eqn{b(i,t)} times abundance \eqn{N(i,t)}, in patch \eqn{i} at time \eqn{t}.
#' Additionally, the effect of the abundances of each stage/sex on these parameters can be weighted. For more information, see the details.
#'
#' @usage StageStructure(Stages, TransMatrix, MaxAge = 100, MinAge = 0,
#'                RepSeasons = 1, PRep = 1.0, RepInterval = 0, SurvSched = 0,
#'                FecDensDep = FALSE, DevDensDep = FALSE, SurvDensDep = FALSE,
#'                FecStageWtsMatrix, DevStageWtsMatrix, SurvStageWtsMatrix,
#'                FecLayer, DevLayer, SurvLayer,
#'                PostDestructn = FALSE)
#' @param Stages Number of life-stages (integer). Requires a minimum of \eqn{2}, for "juvenile" (stage 0) and "adult". Maximum is 10.
#' @param TransMatrix Transition matrix. Defines the development probabilities from each stage into the next, as well as the respective survival probabilities and fecundities. See Details for matrix structure and notation.
#' @param MaxAge Maximum age in years, defaults to \eqn{100}. (integer)
#' @param MinAge Integer vector containing the ages which an individual in stage \eqn{i-1} (with sex \eqn{m/f}, if applicable) must already have reached before it can develop into the next stage \eqn{(i)}. The defaults are \eqn{0} for all stages and sexes. The required format is described in the Details.
#' @param RepSeasons Number of potential reproduction events per year. Defaults to \eqn{1}. (integer)
#' @param RepInterval Number of reproductive seasons which must be missed following a reproduction attempt, before another reproduction attempt may occur. Defaults to \eqn{0}. (integer)
#' @param PRep Probability of reproducing in subsequent reproductive seasons. Defaults to \eqn{1.0}.
#' @param SurvSched Scheduling of Survival. When should survival and development occur?\cr 0 = At reproduction\cr 1 = Between reproductive events (default)\cr 2 = Annually (only for \code{RepSeasons>1})
#' @param FecDensDep,DevDensDep,SurvDensDep Density-dependent fecundity / development / survival probabilities? Defaults to \code{FALSE}. See Details for functional form of density-dependence.
#' @param DevDensCoeff,SurvDensCoeff Relative density dependence coefficient for development / survival. Defaults to \eqn{1.0}.
# @param FecStageWts,DevStageWts,SurvStageWts Stage-dependent density dependence in fecundity / development / survival? Defaults to \code{FALSE}.
# @param FecStageWtsMatrix,DevStageWtsMatrix,SurvStageWtsMatrix Required if stage-dependent density dependence, i.e. if respective \code{FecStageWts}/\code{DevStageWts}/\code{SurfStageWts=TRUE}. Takes a quadratic matrix with (#stages * #sexes) rows/columns indicating the weight of the abundance of each sex/stage on the fecundity/development/survival of each other sex/stage.
#' @param FecStageWtsMatrix,DevStageWtsMatrix,SurvStageWtsMatrix Stage-dependent weights in density dependence of fecundity / development / survival. Takes a quadratic matrix with (#sexes * #stages) rows/columns indicating the weight of the abundance of each stage/sex on the fecundity/development/survival of each other stage/sex. If not set, all stages are weighted equally.
#' @param FecLayer,DevLayer,SurvLayer A matrix of layer indices for the three demographic rates (fecundity/development/survival) if they are spatially varying. The indices match the
#' scaling layers given in the \code{\link[RangeShiftR]{ImportedLandscape}} module to each demographic rate. The value NA denotes a spatially constant rate. The row number corresponds
#' to the stage; the first/second column contains the layer index for females/males. In case of a female-only or simple sexual model (\code{ReproductionType={0,1}}) only the first
#' column is needed and a vector is accepted to represent it. In any case, values will be coerced to integer.
#' @param PostDestructn In a dynamic landscape, determine if all individuals of a population should die (\code{FALSE}, default) or disperse (\code{TRUE}) if its patch gets destroyed.
#' @details In stage-structured populations, generations can overlap and individuals can be classified in different stages (e.g. immature vs. breeding individuals) differing in their
#' demographic parameters. Individuals are characterized by their age and stage. Each stage has a certain fecundity, survival probability and probability of developing to the next stage.
#'
#' These parameters are provided through classical transition matrices \insertCite{caswell2001}{RangeShiftR}.
#' \ifelse{html}{\out{&phi;<sub>i</sub>}}{\eqn{φ_i}} is the fecundity of stage \eqn{i} ;
#' \ifelse{html}{\out{&sigma;<sub>i</sub>}}{\eqn{σ_i}} is the survival probability of an individual in stage \eqn{i} ;
#' and \ifelse{html}{\out{&gamma;<sub>i-j</sub>}}{\eqn{γ_(i-j)}} is the probability of developing from stage \eqn{i} to stage \eqn{j}.
#' In this way, the transition matrix describes the effect of each individuals current stage (column) on all stages at the next timestep (rows).
#' Since all offspring is born into the juvenile stage (stage 0), all fecundities are always located in the first row of the matrix.
#'
#' However, in \emph{RangeShiftR}, these parameters are not used deterministically as \emph{rates} (like it is typical for matrix models) but, instead, as \emph{probabilities} which are
#' applied stochastically at the individual level. Hence, each female at stage \eqn{i}, if it reproduces, produces a number of offspring given by \eqn{Poisson}(\ifelse{html}{\out{&phi;<sub>i</sub>}}{\eqn{φ_i}}),
#' while Bernoulli trials \eqn{Bern}(\ifelse{html}{\out{&sigma;<sub>i</sub>}}{\eqn{σ_i}}) and \eqn{Bern}(\ifelse{html}{\out{&gamma;<sub>i,i+1</sub>}}{\eqn{γ_(i,i+1)}}) determine if an individual/female survives or not
#' and - if it survives - if it develops to the next stage or not.
#'
#' An example transition matrix for a 3-staged only-female or simple sexual (\code{ReproductionType={0,1}}) population model:
#' \tabular{ccccc}{0 \tab | \tab \ifelse{html}{\out{&phi;<sub>1</sub>}}{\eqn{φ_1}} \tab | \tab \ifelse{html}{\out{&phi;<sub>2</sub>}}{\eqn{φ_2}} \cr
#' \eqn{1.0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1</sub> (1-&gamma;<sub>1-2</sub>)}}{\eqn{σ_1 (1 − γ_(1-2))}} \tab | \tab \eqn{0} \cr
#' \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1</sub> &gamma;<sub>1-2</sub>}}{\eqn{σ_1 γ_(1-2)}} \tab | \tab \ifelse{html}{\out{&sigma;<sub>2</sub>}}{\eqn{σ_2}} \cr}
#'
#' In a female-only model (\code{ReproductionType=0}), φ represents the number of female offspring per female and reproductive event. \cr
#' Note that for an implicit (simple) sexual model (\code{ReproductionType=1}), the demographic parameters are not sex-specific. However, individuals are defined by their sex, which is acknowledged for example in the dispersal
#' process and transmission of alleles. The fecundities φ refer to the number of all offspring (males and females) per female and reproductive event.
#'
#' In case of an explicit (complex) sexual model (\code{ReproductionType=2}), in contrast, each stage must be represented by two columns and rows to distinguish male and female demographic parameters.
#' Note, however, that in any case the juvenile stage has only one row; it contains the fecundities. Male fecundities should have one of two possible values: set \ifelse{html}{\out{&phi;<sub>m</sub>}}{\eqn{φ_m}} \eqn{=1.0} for reproductive males or \ifelse{html}{\out{&phi;<sub>m</sub>}}{\eqn{φ_m}} \eqn{=0.0} for non-reproductive males.\cr
#' An example transition matrix for a 3-staged explicit sexual population model \insertCite{weeks1986,lindstrom1998}{RangeShiftR}:
#' \tabular{ccccccccccc}{\eqn{0} \tab | \tab \eqn{0} \tab | \tab \ifelse{html}{\out{&phi;<sub>1m</sub>}}{\eqn{φ_1m}} \tab | \tab \ifelse{html}{\out{&phi;<sub>1f</sub>}}{\eqn{φ_1f}} \tab | \tab \ifelse{html}{\out{&phi;<sub>2m</sub>}}{\eqn{φ_2m}} \tab | \tab \ifelse{html}{\out{&phi;<sub>2f</sub>}}{\eqn{φ_2f}} \cr
#' \eqn{1.0} \tab | \tab \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1m</sub> (1-&gamma;<sub>1-2,m</sub>)}}{\eqn{σ_1m (1 − γ_(1-2,m))}} \tab | \tab \eqn{0} \tab | \tab \eqn{0} \tab | \tab \eqn{0} \cr
#' \eqn{0} \tab | \tab \eqn{1.0} \tab | \tab \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1f</sub> &gamma;<sub>1-2,f</sub>}}{\eqn{σ_1f γ_(1-2,f)}} \tab | \tab \eqn{0} \tab | \tab \eqn{0} \cr
#' \eqn{0} \tab | \tab \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1m</sub> &gamma;<sub>1-2,m</sub>}}{\eqn{σ_1m γ_(1-2,m)}} \tab | \tab \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>2m</sub>}}{\eqn{σ_2m}} \tab | \tab \eqn{0} \cr
#' \eqn{0} \tab | \tab \eqn{0} \tab | \tab \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>1f</sub> &gamma;<sub>1-2,f</sub>}}{\eqn{σ_1f γ_(1-2,f)}} \tab | \tab \eqn{0} \tab | \tab \ifelse{html}{\out{&sigma;<sub>2f</sub>}}{\eqn{σ_2f}} \cr}
#' The mating system is explicitly modelled and a female’s probability of reproducing is determined as described in \code{\link[RangeShiftR]{Demography}} \insertCite{weeks1986,caswell2001}{RangeShiftR}.
#'
#' A common mistake in building a transition matrix is made when offspring produced at year \eqn{t} develop to the next stage in the same year \insertCite{@ @caswell2001 pg. 60-62}{RangeShiftR}. To avoid this problem without losing the offspring stage, and hence the chance for simulating post-natal dispersal,
#' we require an additional explicit juvenile stage (stage 0). Juveniles have to develop to stage 1 in the same year they are born. Hence the minimum number of stages of a stage-structured model is always \eqn{2}. It is important to note that juvenile mortality can be accounted for in
#' two ways. Either, as in the examples above, it is included in adult fecundity φ (by appropriately reducing its value), and \ifelse{html}{\out{&sigma;<sub>0</sub> &gamma;<sub>(0-1)</sub>}}{\eqn{σ_0 γ_(0-1)}} is equal to \eqn{1.0}. This is how it is typically accounted for in matrix models. Or, alternatively,
#' φ is equal to the true maximum fecundity and \ifelse{html}{\out{&sigma;<sub>0</sub> &gamma;<sub>(0-1)</sub>}}{\eqn{σ_0 γ_(0-1)}} is less than \eqn{1.0}.
#' Only the first approach allows straightforward direct comparison with standard analytical matrix models.
#'
#' \emph{Minimum ages:} For every row in the transition matrix, a minimum age must be provided through a vector in argument \code{MinAge}. It specifies the age which an individual in stage \eqn{i-1} (with sex \eqn{m/f}, if applicable) must already have reached before it can develop into the next stage \eqn{(i)}. The vector must be in
#' the order of increasing stage, stating first male then female values; i.e. \eqn{0,1m,1f,2m,2f,...}. Note that (analogous to the transition matrix) the juvenile stage (\eqn{i=0}) has only one entry for male and female. The defaults are \eqn{0} for all stages and sexes.
#' Note that the minimum age for juveniles (stage \eqn{0}) is by definition zero, and that the minimum age for stage \eqn{1} is required to also be zero because individuals may not persist as juveniles beyond the breeding season in which they are born.
#'
#' It is possible to have one or more reproductive seasons per year (\code{RepSeasons}), or a reproductive event once every few years (controlled by \code{RepInterval}). At each reproductive season, two parameters
#' control the likelihood that each individual / female reproduces:
#' \enumerate{
#' \item First, it is determined whether a reproductively mature female is a potential reproducer. The user specifies a minimum interval (\code{RepInterval}) before an individual, that has
#' already reproduced, is able to reproduce again. Only those mature individuals that are either yet to reproduce, or last reproduced more than this number of reproductive seasons previously, are potential breeders.
#' \item Potential breeders all reproduce with a set probability (\code{PRep}). Note that this probability is different from the probability of reproducing \ifelse{html}{\out{p<sub>r</sub>}}{\eqn{p_r}} described in \code{\link[RangeShiftR]{Demography}}.
#' The latter will be additionally applied only in the case of more complex modelling of the mating system (\code{ReproductionType=2}) and it is determined by the number of reproductive males and females present in the cell/patch.
#' }
#' Note that in the current implementation, reproductive attempts that result in zero offspring still count in terms of an individual having to wait for the chance to reproduce again.
#'
#' A major difference between transition matrices and this individual-based model is that in the first, the three processes of reproduction, survival and development happen simultaneously while, in the second, they are explicitly modelled in sequence.
#' The sequence of these events and the time of the dispersal phase in relation to them can change the actual dynamics and density-dependencies in both population growth and dispersal. At the beginning of each year, reproduction is always the first
#' process to be modelled. However, there can be multiple reproductive seasons per year (default is one); in this case the year starts with a succession of all reproductive seasons. There are three choices for the scheduling of reproduction, survival, development and dispersal:
#' \itemize{
#' \item \code{SurvSched=0}: For each reproductive season: reproduction; survival and development of all stages (apart from stage 0); dispersal; survival and development of stage 0. Then: aging; end of the year.
#' \item \code{SurvSched=1}: For each reproductive season: reproduction; dispersal; survival and successive development of all stages. Then: aging; end of the year.
#' \item \code{SurvSched=2}: For each reproductive season: reproduction; dispersal. Then: survival and development of all stages; aging; end of the year. This option applies only for species having multiple reproductive seasons in a year, otherwise it is equivalent to \code{SurvSched=1}.
#' }
#' Option \code{SurvSched=0} gives results that are comparable with the deterministic solution of the matrix. The choice will depend on the biology of the species. If the main mortality happens overwinter, option \code{SurvSched=1} might be more appropriate.
#'
#' Note that \code{SurvSched=1} in combination with multiple reproductive seasons (\code{RepSeasons>1}) implies several evaluations of the fecundity and the survival and development probabilities, so that the transition matrix should be set accordingly.
#' If the transition matrix contains the annual survival and development rates, \code{SurvSched=2} is the appropriate option (fecundity, however, is still given per reproductive event).
#'
#' \emph{Density dependence} can act on each of the three demographic phases (i.e. reproduction, survival and development) and is controlled by \code{FecDensDep,DevDensDep,SurvDensDep}.
#' It is implemented as an exponential decay \insertCite{neubert2000}{RangeShiftR}:
#'
#' \ifelse{html}{\out{&emsp;&emsp; &phi;<sub>i</sub>(r,t) = &phi;<sub>0,i</sub> &ast; e<sup> - b(r) N(t) </sup>}}{\deqn{φ_i(r,t)=φ_(0,i) * exp(- b(r) N(t) ) }}
#'
#' \ifelse{html}{\out{&emsp;&emsp; &sigma;<sub>i</sub>(r,t) = &sigma;<sub>0,i</sub> &ast; e<sup> - C<sub>&sigma;</sub> b(r) N(t) </sup>}}{\deqn{σ_i(r,t)=σ_(0,i) * exp(- C_\sigma b(r) N(t) ) }}
#'
#' \ifelse{html}{\out{&emsp;&emsp; &gamma;<sub>i</sub>(r,t) = &gamma;<sub>0,i</sub> &ast; e<sup> - C<sub>&gamma;</sub> b(r) N(t) </sup>}}{\deqn{γ_i(r,t)=γ_(0,i) * exp(- C_γ b(r) N(t) ) }}
#'
#' where \eqn{b(r)} is the strength of density dependence in fecundity at site \eqn{r}, which is given by the argument \code{K_or_DensDep} in the landscape module.
#' Furthermore, \ifelse{html}{\out{C<sub>&sigma;</sub>}}{\eqn{C_\sigma}} and \ifelse{html}{\out{C<sub>&gamma;</sub>}}{\eqn{C_γ}} (\code{DevDensCoeff,SurvDensCoeff})
#' scale the strength of density dependence in survival and development relative to that in fecundity.
#'
#' Moreover, the strength of density-dependence can be uniform for all stages or stage-dependent. Even greater complexity can be incorporated with
#' different stages contributing differently to density-dependence \insertCite{caswell2004}{RangeShiftR}:
#'
#' \ifelse{html}{\out{&emsp; &phi;<sub>i</sub>(r,t) = &phi;<sub>0,i</sub> &ast; e<sup> - b(r) &Sigma;<sub>j</sub><sup>S</sup> &omega;<sub>&phi;,ij</sub> N(j,t)</sup>}}{\deqn{φ_i(r,t)=φ_(0,i) * exp(- b(r) \Sigma_j^S ω_{φ,ij} N_j(t) ) }}
#'
#' \ifelse{html}{\out{&emsp; &sigma;<sub>i</sub>(r,t) = &sigma;<sub>0,i</sub> &ast; e<sup> - C<sub>&sigma;</sub> b(r) &Sigma;<sub>j</sub><sup>S</sup> &omega;<sub>&sigma;,ij</sub> N(j,t) </sup>}}{\deqn{σ_i(r,t)=σ_(0,i) * exp(- C_\sigma b(r) \Sigma_j^S ω_{σ,ij} N_j(t) )}}
#'
#' \ifelse{html}{\out{&emsp; &gamma;<sub>i</sub>(r,t) = &gamma;<sub>0,i</sub> &ast; e<sup> - C<sub>&gamma;</sub> b(r) &Sigma;<sub>j</sub><sup>S</sup> &omega;<sub>&gamma;,ij</sub> N(j,t)</sup>}}{\deqn{γ_i(r,t)=γ_(0,i) * exp(- C_γ b(r) \Sigma_j^S ω_{γ,ij} N_j(t) )}}
#'
#' where \ifelse{html}{\out{&omega;<sub>&phi;</sub>}}{\eqn{ω_φ}}, \ifelse{html}{\out{&omega;<sub>&sigma;</sub>}}{\eqn{ω_σ}}, \ifelse{html}{\out{&omega;<sub>&gamma;</sub>}}{\eqn{ω_γ}} are weight matrices given by \code{FecStageWtsMatrix, DevStageWtsMatrix, SurvStageWtsMatrix}. Their elements \ifelse{html}{\out{&omega;<sub>ij</sub>}}{\eqn{ω_ij}}
#' represent the contributions of the abundance of stage \eqn{j} to the density dependence in the fecundity / survival / development of stage \eqn{i}, thus they are quadratic matrices of size \code{Stages}\eqn{^2}. Note that the row sums are not required to be normalized, therefore they can be used
#' to scale the density-dependence for the different stages. In fact, any real value will be accepted for the single weights, so care should be taken when setting them.
#'
#' The demographic rates (fecundity, development probability and survival probability) can be allowed to vary spatially, if the landscape is an imported habitat quality map. In this case, additional maps with the same resolution and extent as the habitat quality map(s) need to be given
#' in \code{\link[RangeShiftR]{ImportedLandscape}()}. These additional layers contain percentage values between 0 and 100 and the matrices \code{FecLayer}, \code{DevLayer}, and \code{SurvLayer} indicate the mapping of each demographic rate to these layers. The local value
#' of a given demographic rate for a certain stage and sex in a cell or patch is then determined as the maximum value (the value given in the transition matrix \code{TransMatrix}) scaled with the percentage in the respective mapped layer.
#' For a patch-based landscape, the scaling percentage of a patch is given by the average percentage of its constituent cells.
#' Potential density-dependence mediated by the strength 1/b still takes effect also for spatially-varying demographic rates. The respective base values φ_0, σ_0 or γ_0 are then replaced by their locally scaled values.
#' @examples  # Stage-structure for simple sexual model
#' transmat <- matrix(c(0,1,0,0,0,0.3,0.4,0,1.5,0,0.6,0.3,2.5,0,0,0.8), ncol = 4)
#' stg <- StageStructure(Stages = 4, TransMatrix = transmat, FecDensDep = TRUE, SurvDensDep = TRUE, SurvDensCoeff = 1.5)
#' plotProbs(stg, stage=1:3)
#'
#'  # Stage-structure for explicit sexual model
#' transmat_2 <- matrix(c(0,0.5,0,0,0,0,0,0.5,0,0,1.0,0.4,0,0.3,0,3.0,0,0.3,0,0.4,1.0,0,0,0.85,0,5,0,0,0,0.8), ncol = 6)
#' stg_2 <- StageStructure(Stages = 3, TransMatrix = transmat_2, FecDensDep = TRUE, DevDensDep = TRUE, DevDensCoeff = 1.2)
#' plotProbs(stg_2, stage=c(1,2), sex = 1)
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "StagesParams"
#' @author Anne-Kathleen Malchow
#' @name StageStructure
#' @export StageStructure
StageStructure <- setClass("StagesParams", slots = c(Stages = "integer_OR_numeric",
                                                     TransMatrix = "matrix",
                                                     MaxAge = "integer_OR_numeric",
                                                     MinAge = "integer_OR_numeric",
                                                     RepSeasons = "integer_OR_numeric",   # No. of reproductive seasons per year
                                                     RepInterval = "integer_OR_numeric",
                                                     PRep = "numeric",
                                                     SurvSched = "integer_OR_numeric",
                                                     FecDensDep = "logical",
                                                     FecStageWts = "logical",
                                                     FecStageWtsMatrix = "matrix",
                                                     DevDensDep = "logical",
                                                     DevDensCoeff = "integer_OR_numeric",
                                                     DevStageWts = "logical",
                                                     DevStageWtsMatrix = "matrix",
                                                     SurvDensDep = "logical",
                                                     SurvDensCoeff = "integer_OR_numeric",
                                                     SurvStageWts = "logical",
                                                     SurvStageWtsMatrix = "matrix",
                                                     FecLayer = "matrix_OR_numeric",
                                                     DevLayer = "matrix_OR_numeric",
                                                     SurvLayer = "matrix_OR_numeric",
                                                     PostDestructn = "logical")
                           , prototype = list(#Stages = 2L,
                                              #TransMatrix = matrix(data = NA, nrow = 1, ncol = 1),
                                              MaxAge = 100L,
                                              MinAge = 0L,
                                              RepSeasons = 1L,
                                              RepInterval = 0L,
                                              PRep = 1.0,
                                              SurvSched = 1L,
                                              FecDensDep = FALSE,
                                              FecStageWts = FALSE,
                                              #FecStageWtsMatrix = matrix(data = NA, nrow = 1, ncol = 1),
                                              DevDensDep = FALSE,
                                              DevDensCoeff = 1.0,
                                              DevStageWts = FALSE,
                                              #DevStageWtsMatrix = matrix(data = NA, nrow = 1, ncol = 1),
                                              SurvDensDep = FALSE,
                                              SurvDensCoeff = 1.0,
                                              SurvStageWts = FALSE,
                                              #SurvStageWtsMatrix = matrix(data = NA, nrow = 1, ncol = 1),
                                              #FecLayer,DevLayer,SurvLayer,
                                              PostDestructn = FALSE)
)
    # check forbidden transitions in TransMatrix (e.g. between non-successive stages or between sexes)
setValidity("StagesParams", function(object) {
    msg <- NULL
    if (anyNA(object@Stages) || length(object@Stages)!=1) {
        msg <- c(msg, "Number of Stages must be set and of length 1!")
    }
    else {
        if (object@Stages<2) {
            msg <- c(msg, "Number of Stages must be at least 2!")
        }
    }
    if (anyNA(object@TransMatrix)) {
        msg <- c(msg, "Transition matrix must be given!")
    }
    else {
        if (any(object@TransMatrix[-1, ] < 0 ) || any(object@TransMatrix[-1, ] > 1 ) ) {
            msg <- c(msg, "All elements exept the first row of the transition matrix must be between 0 and 1!")
        }
        if (any(object@TransMatrix[ 1, ] < 0 ) ) {
            msg <- c(msg, "All elements in the first row of the transition matrix (i.e. the fecundities) must be positive!")
        }
    }
    if (anyNA(object@MaxAge) || length(object@MaxAge)!=1) {
        msg <- c(msg, "MaxAge must be set and of length 1!")
    }
    else {
        if (object@MaxAge<2) {
            msg <- c(msg, "MaxAge must be greater than 1!")
        }
    }
    if (anyNA(object@MinAge) || length(object@MinAge)==0) {
        msg <- c(msg, "MinAge must be set!")
    }
    else {
        if (any(object@MinAge<0)) {
            msg <- c(msg, "MinAge must be positive!")
        }
        if (any(object@MinAge>=object@MaxAge)) {
            msg <- c(msg, "MinAge must be smaller than MaxAge for all stages and sexes!")
        }
    }
    if (anyNA(object@RepSeasons) || length(object@RepSeasons)!=1) {
        msg <- c(msg, "RepSeasons must be set and of length 1!")
    }
    else {
        if (object@RepSeasons<1) {
            msg <- c(msg, "RepSeasons must be at least 1!")
        }
    }
    if (anyNA(object@RepInterval) || length(object@RepInterval)!=1) {
        msg <- c(msg, "RepInterval must be set and of length 1!")
    }
    else {
        if (object@RepInterval<0) {
            msg <- c(msg, "RepInterval must be positive!")
        }
    }
    if (anyNA(object@PRep) || length(object@PRep)!=1) {
        msg <- c(msg, "PRep must be set and of length 1!")
    }
    else {
        if (object@PRep <= 0 || object@PRep > 1.0) {
            msg <- c(msg, "PRep must in the half-open interval (0;1]!")
        }
    }
    if (anyNA(object@SurvSched) || length(object@SurvSched)!=1) {
        msg <- c(msg, "SurvSched must be set and of length 1!")
    }
    else {
        if (! object@SurvSched %in% (0:2) ) {
            msg <- c(msg, "SurvSched must be 0, 1 or 2!")
        }
    }
    if (anyNA(object@FecDensDep) || length(object@FecDensDep)!=1) {
        msg <- c(msg, "FecDensDep must be set and of length 1!")
    }
    else{
        if(object@FecDensDep) {
            if (anyNA(object@FecStageWts) || length(object@FecStageWts)!=1) {
                msg <- c(msg, "FecStageWts must be set and of length 1!")
            }
            else{
                if (object@FecStageWts) {
                    if (anyNA(object@FecStageWtsMatrix)) {
                        msg <- c(msg, "FecStageWtsMatrix must be set (since FecStageWts = TRUE) !")
                    }
                    # else {
                    #     if(any(object@FecStageWtsMatrix < 0)){
                    #         msg <- c(msg, "All elements of FecStageWtsMatrix must be positive or zero!")
                    #     }
                    # }
                }
            }
        }
    }
    if (anyNA(object@DevDensDep) || length(object@DevDensDep)!=1) {
        msg <- c(msg, "DevDensDep must be set and of length 1!")
    }
    else{
        if(object@DevDensDep) {
            if (anyNA(object@DevDensCoeff) || length(object@DevDensCoeff)!=1) {
                msg <- c(msg, "DevDensCoeff must be set and of length 1!")
            }
            else {
                if (object@DevDensCoeff<=0) {
                    msg <- c(msg, "DevDensCoeff must be greater than zero!")
                }
            }
            if (anyNA(object@DevStageWts) || length(object@DevStageWts)!=1) {
                msg <- c(msg, "DevStageWts must be set and of length 1!")
            }
            else{
                if (object@DevStageWts) {
                    if (anyNA(object@DevStageWtsMatrix)) {
                        msg <- c(msg, "DevStageWtsMatrix must be set (since DevStageWts = TRUE) !")
                    }
                    # else {
                    #     if(any(object@DevStageWtsMatrix < 0)){
                    #         msg <- c(msg, "All elements of DevStageWtsMatrix must be positive or zero!")
                    #     }
                    # }
                }
            }
        }
    }
    if (anyNA(object@SurvDensDep) || length(object@SurvDensDep)!=1) {
        msg <- c(msg, "SurvDensDep must be set and of length 1!")
    }
    else{
        if(object@SurvDensDep){
            if (anyNA(object@SurvDensCoeff) || length(object@SurvDensCoeff)!=1) {
                msg <- c(msg, "SurvDensCoeff must be set and of length 1!")
            }
            else {
                if (object@SurvDensCoeff<=0) {
                    msg <- c(msg, "SurvDensCoeff must be greater than zero!")
                }
            }
            if (anyNA(object@SurvStageWts) || length(object@SurvStageWts)!=1) {
                msg <- c(msg, "SurvStageWts must be set and of length 1!")
            }
            else{
                if (object@SurvStageWts) {
                    if (anyNA(object@SurvStageWtsMatrix)) {
                        msg <- c(msg, "SurvStageWtsMatrix must be set (since SurvStageWts = TRUE) !")
                    }
                    # else{
                    #     if(any(object@SurvStageWtsMatrix < 0)){
                    #         msg <- c(msg, "All elements of SurvStageWtsMatrix must be positive or zero!")
                    #     }
                    # }
                }
            }
        }
    }
    if (length(object@FecLayer)>0) {
        if( any( !is.na( object@FecLayer[object@FecLayer<=0] ))) {
            msg <- c(msg, "Elements of FecLayer must be strictly postive integers or NA for no spatial variation!")
        }
        if( min(nrow(object@FecLayer),length(object@FecLayer)) != object@Stages) {
            msg <- c(msg, "FecLayer must have as many rows as Stages!")
        }
    }
    if (length(object@DevLayer)>0) {
        if( any( !is.na( object@DevLayer[object@DevLayer<=0] ))) {
            msg <- c(msg, "Elements of DevLayer must be strictly postive integers or NA for no spatial variation!")
        }
        if( min(nrow(object@DevLayer),length(object@DevLayer)) != object@Stages) {
            msg <- c(msg, "DevLayer must have as many rows as Stages!")
        }
    }
    if (length(object@SurvLayer)>0) {
        if( any( !is.na( object@SurvLayer[object@SurvLayer<=0] ))) {
            msg <- c(msg, "Elements of SurvLayer must be strictly postive integers or NA for no spatial variation!")
        }
        if( min(nrow(object@SurvLayer),length(object@SurvLayer)) != object@Stages) {
            msg <- c(msg, "SurvLayer must have as many rows as Stages!")
        }
    }
    if (is.na(object@PostDestructn) || length(object@PostDestructn)!=1) {
        msg <- c(msg, "PostDestructn must be set and of length 1!")
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "StagesParams", function(.Object,...) {
    this_func = "StageStructure(): "
    args <- list(...)
    #set weights switches if matrix is given
    if (!is.null(args$FecStageWtsMatrix)) {
        .Object@FecStageWts <- TRUE
    }
    if (!is.null(args$DevStageWtsMatrix)) {
        .Object@DevStageWts <- TRUE
    }
    if (!is.null(args$SurvStageWtsMatrix)) {
        .Object@SurvStageWts <- TRUE
    }
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (!.Object@FecDensDep) {
        if (.Object@FecStageWts) {
            warning(this_func, "FecStageWts stage weigths", warn_msg_ignored, "since FecDensDep = FALSE.", call. = FALSE)
        }
    }
    if (!.Object@DevDensDep) {
        if (!is.null(args$DevDensCoeff)) {
            warning(this_func, "DevDensCoeff", warn_msg_ignored, "since DevDensDep = FALSE.", call. = FALSE)
        }
        if (.Object@DevStageWts) {
            warning(this_func, "DevStageWts stage weigths", warn_msg_ignored, "since DevDensDep = FALSE.", call. = FALSE)
        }
    }
    if (!.Object@SurvDensDep) {
        if (!is.null(args$SurvDensCoeff)) {
            warning(this_func, "SurvDensCoeff", warn_msg_ignored, "since SurvDensDep = FALSE.", call. = FALSE)
        }
        if (.Object@SurvStageWts) {
            warning(this_func, "SurvStageWts stage weigths", warn_msg_ignored, "since SurvDensDep = FALSE.", call. = FALSE)
        }
    }
    # replace NAs & convert to matrix
    if (!is.null(args$FecLayer)) {
        #.Object@FecLayer[is.na(.Object@FecLayer)] <- -9
        if(class(.Object@FecLayer)[1]!="matrix") .Object@FecLayer <- matrix(.Object@FecLayer)
    }
    if (!is.null(args$DevLayer)) {
        #.Object@DevLayer[is.na(.Object@DevLayer)] <- -9
        if(class(.Object@DevLayer)[1]!="matrix") .Object@DevLayer <- matrix(.Object@DevLayer)
    }
    if (!is.null(args$SurvLayer)) {
        #.Object@SurvLayer[is.na(.Object@SurvLayer)] <- -9
        if(class(.Object@SurvLayer)[1]!="matrix") .Object@SurvLayer <- matrix(.Object@SurvLayer)
    }
    # convert to integer
    if(length(.Object@FecLayer)>0) mode(.Object@FecLayer) <- "integer"
    if(length(.Object@DevLayer)>0) mode(.Object@DevLayer) <- "integer"
    if(length(.Object@SurvLayer)>0) mode(.Object@SurvLayer) <- "integer"
    .Object}
)
setMethod("show", "StagesParams", function(object){
    cat("   Stages           :", paste(object@Stages), "\n")
    cat("   Transition matrix:\n")
    print(object@TransMatrix)
    cat("   Maximum age      :", paste(object@MaxAge) , "\n")
    cat("   Minimum age      :", paste(object@MinAge) , "\n")
    cat("   Reprod. seasons  :", paste(object@RepSeasons) , "\n")
    cat("   Prob. Reprod.    :", paste(object@PRep) , "\n")
    cat("   Reprod. Intervals:", paste(object@RepInterval) , "\n")
    cat("   SurvScheduling   :", paste(object@SurvSched) , "\n")
    if (object@FecDensDep || object@DevDensDep || object@SurvDensDep) {
        cat("   Density-dependence in:\n")
    }
    if (object@FecDensDep) {
        cat("    - Fecundity\n")
        if (object@FecStageWts) {
            cat("    with stage-dependent density depedence with weights:\n")
            print(object@FecStageWtsMatrix)
        }
    }
    if (object@DevDensDep) {
        cat("    - Development, with coefficient = ", paste(object@DevDensCoeff) , "\n")
        if (object@DevStageWts) {
            cat("                 and stage-dependence with weights: \n")
            print(object@DevStageWtsMatrix)
        }
    }
    if (object@SurvDensDep) {
        cat("    - Survival,    with coefficient =", paste(object@SurvDensCoeff) , "\n")
        if (object@SurvStageWts) {
            cat("                 and stage-dependence with weights: \n")
            print(object@SurvStageWtsMatrix)
        }
    }
    if(length(object@FecLayer)>0) {
        cat("   Spatially varying fecundity layer indices:\n")
        print(object@FecLayer)
    }
    if(length(object@DevLayer)>0) {
        cat("   Spatially varying development layer indices:\n")
        print(object@DevLayer)
    }
    if(length(object@SurvLayer)>0) {
        cat("   Spatially varying survival layer indices:\n")
        print(object@SurvLayer)
    }
    cat("   PostDestructn    :", paste(object@PostDestructn) )
    if (object@PostDestructn) {
        cat(" (all disperse)\n")
    }
    else {
        cat(" (all die)\n")
    }
    cat ("\n")}
)
setMethod("plotProbs", "StagesParams", function(x, stage = NULL, sex = NULL, xmax = NULL, ymax = NULL){
    lines <- dim(x@TransMatrix)[2]
    if (lines==x@Stages) {SexDep = FALSE}
    else if (lines==(2*x@Stages)) {SexDep = TRUE}
    else {
        print("Wrong format of transition matrix.\n")
        return(NULL)
    }
    fecs <- x@TransMatrix[1,]
    if (SexDep) {
        surv <- devs <- rep(0,lines)
        for(s in 1:lines ){
            if(s==1) ss <- x@TransMatrix[s,s] else ss <- x@TransMatrix[s-1,s]
            if((s+2)>lines) dd <- 0 else dd <- x@TransMatrix[s+1,s]
            surv[s] <- ss+dd
            devs[s] <- dd/(ss+dd)
        }
    }
    else { # !SexDep
        if (!is.null(sex)) print("This stage structure has no sex-dependency.\n")
        surv <- devs <- rep(0,lines)
        for(s in 1:lines ){
            ss <- x@TransMatrix[s,s]
            if(s!=lines) dd <- x@TransMatrix[s+1,s] else dd <- 0
            surv[s] <- ss+dd
            devs[s] <- dd/(ss+dd)
        }
    }
    if (is.null(stage)){stage <- 0:10}
    if (is.null(sex)){sex <- 0:1}
    # New plot
    if (is.null(xmax)) {xmax = 4.0}
    xvals = seq(0, xmax, length.out = 100)
    if (is.null(ymax)) { ymax = max(1.0,fecs) }
    plot(NULL, type = "n", ylab = "Fecundity, Survival, Development", xlab = "relative population density bN", xlim = c(0,xmax), ylim = c(0,ymax))
    leg.txt <- leg.col <- c()
    for(line in 1:lines ){
        if (SexDep) {
            this.stage = as.integer(line/2-0.1)
            this.sex = (line+1)%%2
        }
        else {
            this.stage = line-1
            this.sex = 0
        }
        if (this.stage %in% stage && this.sex %in% sex) {
            if(x@FecDensDep){ lines(xvals, fecs[line]*exp(-xvals), type = "l", lty = 1, col = line) }
            else{ lines(xvals, rep(fecs[line], length(xvals)), type = "l", lty = 1, col = line) }
            if(x@SurvDensDep){ lines(xvals, surv[line]*exp(-x@SurvDensCoeff*xvals), type = "l", lty = 2, col = line) }
            else{ lines(xvals, rep(surv[line], length(xvals)), type = "l", lty = 2, col = line) }
            if(x@DevDensDep){ lines(xvals, devs[line]*exp(-x@DevDensCoeff*xvals), type = "l", lty = 3, col = line) }
            else{ lines(xvals, rep(devs[line], length(xvals)), type = "l", lty = 3, col = line) }
            if(SexDep) {leg.txt <- c(leg.txt, paste0("Stage ",this.stage, ifelse(this.sex," female"," male")))}
            else {leg.txt <- c(leg.txt, paste0("Stage ",this.stage))}
            leg.col <- c(leg.col, line)
        }
    }
    if (length(leg.txt)>0) {
        leg.txt <- c("Fecundity","Survival prob.","Developmt. prob.",leg.txt)
        legend("topright", leg.txt, col = c(rep(1,3),leg.col), lwd = 1.5, lty = c(1:3,rep(1,length(leg.col))) )
    }
})


#### CLASS DEMOGRAPHY

# contains basic demographic parameters (originally in RS 'Rarameters'-file) and optionally the 'StageStruct' object

# define this ClassUnion so that the 'stages' slot in the parameter master class 'RSparams' can be FALSE for option 'population with non-overlapping generations'
setClassUnion("StagesSlot", c("logical", "StagesParams"))

#' Set Demographic Parameters
#'
#' For a simple non-structured population, set its basic demographic parameters here, i.e. the maximum growth rate (\code{Rmax}) and the competition coefficient (\code{bc}).
#' For a stage-structured population, define its corresponding parameters via \code{\link[RangeShiftR]{StageStructure}} and add it to Demography.\cr
#' \cr
#' Choose the Reproduction model that determines if sexes are considered implicitly or explicitly and if a mating system is used. If applicable, set the corresponding parameters, i.e. the proportion of males (\code{PropMales}) and the maximum harem size (\code{Harem}).
#'
#' @usage Demography(Rmax, bc = 1.0, StageStruct = FALSE,
#'            ReproductionType = 0, PropMales = 0.5, Harem = 1)
#' @param Rmax Maximum growth rate. Describes the mean number of offspring per female and reproductive event at very low density. Only required if \code{StageStruct=FALSE}.
#' @param bc Competition coefficient. Describes the type of density regulation, providing the possibility for under-compensatory (\eqn{b_c < 1}), compensatory (\eqn{b_c = 1}) (default) or over-compensatory (\eqn{b_c > 1}) dynamics. Only required if \code{StageStruct=FALSE}.
#' @param StageStruct \code{FALSE} (default) yields a population model with non-overlapping generations.\cr For a stage-structured population, this takes the corresponding parameter object generated by \code{\link[RangeShiftR]{StageStructure}}, which holds all demographic parameters.
#' @param ReproductionType 0 = asexual / only female model (default)\cr1 = simple sexual model\cr2 = sexual model with explicit mating system
#' @param PropMales Required if \code{ReproductionType={1,2}}: Proportion of males in the population, between \eqn{0} and \eqn{1}. Defaults to \eqn{0.5}.
#' @param Harem Required if \code{ReproductionType=2}: Maximum harem size. The maximum number of pair bonds that a male can establish. \eqn{Harem = 1} (default) corresponds to monogamy, \eqn{0<Harem<1} to polyandry and \eqn{Harem>1} to polygyny.
#' @details The following information regards the population dynamics of a \strong{non-structured} (\code{StageStruct=FALSE}) population.\cr
#' For more information on the population dynamics of a \strong{structured} population, see \code{\link[RangeShiftR]{StageStructure}}.
#'
#' Populations with non-overlapping generations, i.e. with \strong{no stage-structure} are the appropriate way to model species that have discrete generations.
#' At each generation the life cycle comprises - in that order - reproduction, death of the adults and offspring dispersal.
#' Two parameters determine the nature of the demographic density dependence:
#' the carrying capacity \eqn{K} (given by the argument \code{K_or_DensDep} in the landscape module) and
#' the competition coefficient \eqn{b_c} (given by the argument \code{bc}).
#' These discrete generation models can be applied to asexual species, species for which it is assumed that females play the dominant role in spatial dynamics
#' and for species for which it is considered crucial to model both sexes explicitly.
#'
#' \emph{Asexual / only-female models:}  (\code{ReproductionType=0})\cr
#' Recruitment is determined by a stochastic, individual-based formulation of the \insertCite{smith1973;textual}{RangeShiftR} population model, where the number of offspring produced by a single individual in the cell/patch \eqn{i} at time \eqn{t}, is drawn from the following distribution:\cr
#'
#' \ifelse{html}{\out{&emsp;&emsp;N<sub>juv</sub>(i,t) = Poisson( R(i,t) / (1+|R(i,t) - 1| &ast; (N(i,t) / K(i,t))<sup>b<sub>c</sub></sup> ) ) } }{\deqn{N_juv(i,t) = Poisson( R(i,t) / (1 + |R(i,t) - 1|*( N(i,t) / K(i,t) )^bc ) ) } }
#'
#' Here, \eqn{R(i,t)} is the maximum growth rate \code{Rmax} (obtained at very low density only) and \eqn{K(i,t)} is the carrying capacity
#' at patch \eqn{i} and time \eqn{t}.
#' Both \eqn{R(i,t)} and \eqn{K(i,t)} can vary in space and time, depending on the model setting. \ifelse{html}{\out{b<sub>c</sub>}}{\eqn{b_c}}
#' is the competition coefficient which describes the type of density regulation, providing the possibility for under-compensatory
#' (\ifelse{html}{\out{b<sub>c</sub> < 1}}{\eqn{b_c < 1}}), compensatory (\ifelse{html}{\out{b<sub>c</sub> = 1}}{\eqn{b_c = 1}}) or
#' over-compensatory (\ifelse{html}{\out{b<sub>c</sub> > 1}}{\eqn{b_c > 1}}) dynamics.\cr
#' \cr
#' \emph{Sexual models}\cr
#' In this type of models, individuals are explicitly characterized by their sex. The proportion of each sex in the population is controlled by setting the proportion of males (\code{PropMales}). There are two types of possible sexual sub-models:\cr
#' \cr
#' \emph{Simple mating system:}  (\code{ReproductionType=1})\cr
#' This is the simplest form of mate limitation. Each female individual is assumed to mate, as long as there is at least one male in the population. As for the asexual case, the Maynard Smith and Slatkin model is used to determine the expected number of
#' offspring produced by each female. To maintain equivalence between the asexual and sexual versions, the expected value of the Poisson distribution is multiplied by \eqn{2} (Lindström & Kokko 1998):\cr
#'
#' \ifelse{html}{\out{&emsp;&emsp;N<sub>juv</sub>(i,t) = Poisson( 2 R(i,t) / (1+|R(i,t) - 1| &ast; (N(i,t) / K(i,t) )<sup>b<sub>c</sub></sup> ) ) } }{\deqn{N_juv(i,t) = Poisson( 2 R(i,t) / (1 + |R(i,t) - 1|*( N(i,t) / K(i,t) )^bc ) ) } }
#'
#' \emph{Complex mating system:}  (\code{ReproductionType=2})\cr
#' More complex and flexible mating system. Mating is explicitly modelled through a mating function
#' \insertCite{lindstrom1998,legendre2004,bessa2010}{RangeShiftR}, where the number of mated females \eqn{c} is given by:
#' \deqn{c = f * min(1, 2hm/(f+hm) )}
#' where \eqn{f} and \eqn{m} are the numbers of potentially reproductive females and males, respectively, and \eqn{h} is the maximum harem size \code{Harem}, i.e. the maximum number of pair bonds that a male can establish. \eqn{h=1} corresponds to monogamy, \eqn{0<h<1} to polyandry
#' and \eqn{h>1} to polygyny. Each potentially reproductive female has a probability of reproducing \ifelse{html}{\out{p<sub>r</sub>}}{\eqn{p_r}} given by:
#' \deqn{p_r=c/f}
#' A Bernoulli trial \eqn{Bern}\ifelse{html}{\out{(p<sub>r</sub>)}}{\eqn{(p_r)}} determines if the female reproduces or not. For females that reproduce, the number of offspring is determined through the Poisson distribution as stated above. Hence, the specification of the mating system determines the probability for each female to reproduce. However, no explicit pair bonds are formed, and in the cases where traits inheritance is
#' involved, the father (of all the offspring produced by a single female in a single reproductive event) is selected randomly from the males in the population.\cr
#' \cr
#' Populations with overlapping generations, i.e. with \strong{with stage-structure} are the appropriate choice for species in which generations can overlap and individuals can be classified in different stages
#' (e.g. immature vs. breeding individuals) differing in their demographic parameters. Individuals are characterized by their age and stage. Each stage has a certain fecundity (which will be used as expected value of the Poisson distribution for reproduction), survival
#' and probability of developing to the next stage. The parameters are provided as classical transition matrices \insertCite{caswell2001}{RangeShiftR} through a \code{\link[RangeShiftR]{StageStructure}}
#' parameter object. For more information on the population dynamics in this case, see the details there.
#' @references
#'         \insertAllCited{}
#' @return a parameter object of class "DemogParams"
#' @author Anne-Kathleen Malchow
#' @name Demography
#' @export Demography
Demography <- setClass("DemogParams", slots = c(Rmax = "integer_OR_numeric",
                                                bc = "numeric",
                                                StageStruct = "StagesSlot",
                                                ReproductionType = "integer_OR_numeric",
                                                PropMales = "numeric",
                                                Harem = "integer_OR_numeric")
                       , prototype = list(#Rmax = 1.5,
                                          bc = 1.0,
                                          StageStruct = FALSE,
                                          ReproductionType = 0L,
                                          PropMales = 0.5,
                                          Harem = 1L)
)
setValidity("DemogParams", function(object) {
    msg <- NULL
    if (anyNA(object@ReproductionType) || length(object@ReproductionType)!=1) {
        msg <- c(msg, "ReproductionType must be set and of length 1!")
    }
    else {
        if (! object@ReproductionType %in% (0:2) ) {
            msg <- c(msg, "ReproductionType must be 0, 1 or 2!")
        }
        else{
            if (object@ReproductionType) {
                if (anyNA(object@PropMales) || length(object@PropMales)!=1) {
                    msg <- c(msg, "Proportion of males must be set and of length 1!")
                }
                else {
                    if (object@PropMales <= 0.0 || object@PropMales >= 1.0) {
                        msg <- c(msg, "Proportion of males must be in the open interval (0,1)!")
                    }
                }
            }
            if (object@ReproductionType==2) {
                if (anyNA(object@Harem) || length(object@Harem)!=1) {
                    msg <- c(msg, "Maximum harem size must be set and of length 1!")
                }
                else {
                    if (object@Harem<=0) {
                        msg <- c(msg, "Maximum harem size must be positive!")
                    }
                }
            }
        }
    }
    validObject(object@StageStruct)
    if (class(object@StageStruct)[1]=="logical") {
        if (object@StageStruct) {                # StageStruct=TRUE
            msg <- c(msg, "StageStruct must either be FALSE or an object of class \"StagesParams\" !")
        }
        else {                                   # StageStruct=FALSE
            if (anyNA(object@Rmax) || length(object@Rmax)!=1) {
                msg <- c(msg, "Maximum growth rate Rmax must be set for a non-structured population!")
            }
            else {
                if (object@Rmax<=0) {
                    msg <- c(msg, "Maximum growth rate Rmax must be positive!")
                }
            }
            if (anyNA(object@bc) || length(object@bc)!=1) {
                msg <- c(msg, "Competition coefficient bc must be set for a non-structured population!")
            }
            else {
                if (object@bc<=0) {
                    msg <- c(msg, "Competition coefficient bc must be positive!")
                }
            }
        }
    }
    else{
        if (class(object@StageStruct)[1]=="StagesParams") { #stage-structured model
            stgs <- object@StageStruct@Stages
            if (object@ReproductionType == 2) { # explicit sexual model
                if ( any(dim(object@StageStruct@TransMatrix)!=c(2*stgs-1,2*stgs) )) {
                    msg <- c(msg, "Transition Matrix must have dimensions 2*\'Stages\'-1 x 2*\'Stages\' !")
                }
                else { # check transition matrix
                    if (any(object@StageStruct@TransMatrix[1,1:2]!=0) ) {
                        msg <- c(msg, "Transition Matrix: Fecundity of juvenile stage must be zero!")
                    }
                    if (any(object@StageStruct@TransMatrix[1,]<0) ) {
                        msg <- c(msg, "Transition Matrix: Fecunditiy must be positive or zero!")
                    }
                    if (any(object@StageStruct@TransMatrix[-1,]<0) || any(object@StageStruct@TransMatrix[-1,]>1) ) {
                        msg <- c(msg, "Transition Matrix: All transition probabilities must be in the closed interval [0,1]")
                    }
                    else{ # passed first line of tests
                        for (ss in 1:2){
                            if (any(object@StageStruct@TransMatrix[-(ss+1),ss]!=0)) {
                                msg <- c(msg, "Transition Matrix: Found non-zero entries (in juvenile) where zeroes are expected!")
                            }
                        }
                        if (stgs>2){
                            for (ss in 3:((stgs-1)*2)){
                                if (sum(object@StageStruct@TransMatrix[c((ss-1),(ss+1)),ss])>1) {
                                    msg <- c(msg, "Transition Matrix: Sum of probabilities to stay in current stage and to develop can't exceed 1!")
                                }
                                if (any(object@StageStruct@TransMatrix[-c(1,(ss-1),(ss+1)),ss]!=0)) {
                                    msg <- c(msg, "Transition Matrix: Found non-zero entries where zeroes are expected!")
                                }
                            }
                        }
                        for (ss in (stgs*2-1):(stgs*2)){
                            if (any(object@StageStruct@TransMatrix[-c(1,(ss-1)),ss]!=0)) {
                                msg <- c(msg, "Transition Matrix: Found non-zero entries where zeroes are expected!")
                            }
                        }
                    }
                }
                if (length(object@StageStruct@MinAge)==1 && object@StageStruct@MinAge==0) {
                    # do nothing and figure out a way to set as default a null vector of the correct length
                }
                else {
                    if (length(object@StageStruct@MinAge)!=(2*stgs-1)) {
                        msg <- c(msg, "Minimum Age vector must have dimensions 2*\'Stages\'-1 !")
                    }
                    else{
                        if (any(object@StageStruct@MinAge[1:2]!=0) ) {
                            msg <- c(msg, "MinAge of stages 0 and 1 must be zero!")
                        }
                        if (any(object@StageStruct@MinAge<0)) {
                            msg <- c(msg, "MinAge must be positive!")
                        }
                    }
                }
                if (object@StageStruct@FecStageWts) {
                    if (!anyNA(object@StageStruct@FecStageWtsMatrix)) {
                        if ( any(dim(object@StageStruct@FecStageWtsMatrix)!=c(2*stgs,2*stgs) )) {
                            msg <- c(msg, "FecStageWtsMatrix must have dimensions 2*\'Stages\' x 2*\'Stages\' !")
                        }
                    }
                }
                if (object@StageStruct@DevStageWts) {
                    if (!anyNA(object@StageStruct@DevStageWtsMatrix)) {
                        if ( any(dim(object@StageStruct@DevStageWtsMatrix)!=c(2*stgs,2*stgs) )) {
                            msg <- c(msg, "DevStageWtsMatrix must have dimensions 2*\'Stages\' x 2*\'Stages\' !")
                        }
                    }
                }
                if (object@StageStruct@SurvStageWts) {
                    if (!anyNA(object@StageStruct@SurvStageWtsMatrix)) {
                        if ( any(dim(object@StageStruct@SurvStageWtsMatrix)!=c(2*stgs,2*stgs) )) {
                            msg <- c(msg, "SurvStageWtsMatrix must have dimensions 2*\'Stages\' x 2*\'Stages\' !")
                        }
                    }
                }
                if (length(object@StageStruct@FecLayer)>0) {
                    if ( any(dim(object@StageStruct@FecLayer)!=c(stgs,2) )) {
                        msg <- c(msg, "FecLayer must have dimensions \'Stages\' x 2 !")
                    }
                }
                if (length(object@StageStruct@DevLayer)>0) {
                    if ( any(dim(object@StageStruct@DevLayer)!=c(stgs,2) )) {
                        msg <- c(msg, "DevLayer must have dimensions \'Stages\' x 2 !")
                    }
                }
                if (length(object@StageStruct@SurvLayer)>0) {
                    if ( any(dim(object@StageStruct@SurvLayer)!=c(stgs,2) )) {
                        msg <- c(msg, "SurvLayer must have dimensions \'Stages\' x 2 !")
                    }
                }
            }
            else{ # asexual model or simple sexual model
                if ( any(dim(object@StageStruct@TransMatrix)!=c(stgs,stgs) )) {
                    msg <- c(msg, "Transition Matrix must have dimensions \'Stages\' x \'Stages\' !")
                }
                else { # check transition matrix
                    if (object@StageStruct@TransMatrix[1,1]!=0) {
                        msg <- c(msg, "Transition Matrix: Fecundity of juvenile stage must be zero!")
                    }
                    if (any(object@StageStruct@TransMatrix[1,]<0) ) {
                        msg <- c(msg, "Transition Matrix: Fecunditiy must be positive or zero!")
                    }
                    if (any(object@StageStruct@TransMatrix[-1,]<0) || any(object@StageStruct@TransMatrix[-1,]>1) ) {
                        msg <- c(msg, "Transition Matrix: All transition probabilities must be in the closed interval [0,1]")
                    }
                    else{ # passed first line of tests
                        for (ss in 1:(stgs-1)){
                            if (sum(object@StageStruct@TransMatrix[ss:(ss+1),ss])>1) {
                                msg <- c(msg, "Transition Matrix: Sum of probabilities to stay in current stage and to develop can't exceed 1!")
                            }
                            if (any(object@StageStruct@TransMatrix[-c(1,ss,(ss+1)),ss]!=0)) {
                                msg <- c(msg, "Transition Matrix: Found non-zero entries where zeroes are expected!")
                            }
                        }
                        if (any(object@StageStruct@TransMatrix[-c(1,stgs),stgs]!=0)) {
                            msg <- c(msg, "Transition Matrix: Found non-zero entries where zeroes are expected!")
                        }
                    }
                }
                if (length(object@StageStruct@MinAge)==1 && object@StageStruct@MinAge==0) {
                    # do nothing and figure out a way to set as default a null vector of the correct length
                }
                else {
                    if (length(object@StageStruct@MinAge)!=stgs) {
                        msg <- c(msg, "Minimum Age vector must have dimensions \'Stages\' !")
                    }
                    else{
                        if (any(object@StageStruct@MinAge[1:2]!=0) ) {
                            msg <- c(msg, "MinAge of stages 0 and 1 must be zero!")
                        }
                        if (any(object@StageStruct@MinAge<0)) {
                            msg <- c(msg, "MinAge must be positive!")
                        }
                    }
                }
                if (object@StageStruct@FecStageWts) {
                    if (!anyNA(object@StageStruct@FecStageWtsMatrix)) {
                        if ( any(dim(object@StageStruct@FecStageWtsMatrix)!=c(stgs,stgs) )) {
                            msg <- c(msg, "FecStageWtsMatrix must have dimensions \'Stages\' x \'Stages\' !")
                        }
                    }
                }
                if (object@StageStruct@DevStageWts) {
                    if (!anyNA(object@StageStruct@DevStageWtsMatrix)) {
                        if ( any(dim(object@StageStruct@DevStageWtsMatrix)!=c(stgs,stgs) )) {
                            msg <- c(msg, "DevStageWtsMatrix must have dimensions \'Stages\' x \'Stages\' !")
                        }
                    }
                }
                if (object@StageStruct@SurvStageWts) {
                    if (!anyNA(object@StageStruct@SurvStageWtsMatrix)) {
                        if ( any(dim(object@StageStruct@SurvStageWtsMatrix)!=c(stgs,stgs) )) {
                            msg <- c(msg, "SurvStageWtsMatrix must have dimensions \'Stages\' x \'Stages\' !")
                        }
                    }
                }
                if (length(object@StageStruct@FecLayer)>0) {
                    if ( any(dim(object@StageStruct@FecLayer)!=c(stgs,1) )) {
                        msg <- c(msg, "FecLayer must have dimensions \'Stages\' x 1 !")
                    }
                }
                if (length(object@StageStruct@DevLayer)>0) {
                    if ( any(dim(object@StageStruct@DevLayer)!=c(stgs,1) )) {
                        msg <- c(msg, "DevLayer must have dimensions \'Stages\' x 1 !")
                    }
                }
                if (length(object@StageStruct@SurvLayer)>0) {
                    if ( any(dim(object@StageStruct@SurvLayer)!=c(stgs,1) )) {
                        msg <- c(msg, "SurvLayer must have dimensions \'Stages\' x 1 !")
                    }
                }
            }
        }
        else{
            msg <- c(msg, "StageStruct must either be FALSE or an object of class \"StagesParams\" !")
        }
    }
    if (is.null(msg)) TRUE else msg}
)
setMethod("initialize", "DemogParams", function(.Object,...) {
    this_func = "Demography(): "
    args <- list(...)
    .Object <- callNextMethod()
    if ( length(args) == 0 ) {
        validObject(.Object)
    }
    if (class(.Object@StageStruct)[1]=="StagesParams") {
        .Object@Rmax = -9L
        if (!is.null(args$Rmax)) {
            warning(this_func, "Maximum growth rate Rmax", warn_msg_ignored, "in a stage-structured model.", call. = FALSE)
        }
        .Object@bc = -9L
        if (!is.null(args$bc)) {
            warning(this_func, "Competition coefficient bc", warn_msg_ignored, "in a stage-structured model.", call. = FALSE)
        }
    }
    if (.Object@ReproductionType==0) {
        .Object@PropMales = -9L
        if (!is.null(args$PropMales)) {
            warning(this_func, "Proportion of males", warn_msg_ignored, "in an asexual / only-female model.", call. = FALSE)
        }
    }
    if (.Object@ReproductionType!=2) {
        .Object@Harem = -9L
        if (!is.null(args$Harem)) {
            warning(this_func, "Maximum harem size", warn_msg_ignored, "in a model without mating system. (Only used when ReproductionType = 2)", call. = FALSE)
        }
    }
    .Object}
)
setMethod("show", "DemogParams", function(object){
    cat(" Demography:\n")
    if (class(object@StageStruct)[1]=="StagesParams") {
        cat("  Stage-structured population:\n")
        print(object@StageStruct)
    }
    else {
        cat("  Unstructured population:\n")
        cat("   Rmax :", paste(object@Rmax), "\n")
        cat("   bc   :", paste(object@bc) , "\n")
    }
    cat("  Reproduction Type :", paste(object@ReproductionType))
    if(object@ReproductionType) {
        cat(" (sexual model)\n")
        cat("   PropMales :", paste(object@PropMales), "\n")
    }
    else{
        cat(" (female only)\n")
    }
    if(object@ReproductionType==2) {
        cat("   Harem :", paste(object@Harem), "\n")
    }
})

# RS manual 2.4.2 (page 23) - 2.4.3 (page 34)
