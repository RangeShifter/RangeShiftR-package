% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_DemogParams.R
\docType{class}
\name{Demography}
\alias{Demography}
\title{Set Demographic Parameters}
\usage{
Demography(Rmax, bc = 1.0, StageStruct = FALSE,
           ReproductionType = 0, PropMales = 0.5, Harem = 1)
}
\arguments{
\item{Rmax}{Maximum growth rate. Describes the mean number of offspring per female and reproductive event at very low density. Only required if \code{StageStruct=FALSE}.}

\item{bc}{Competition coefficient. Describes the type of density regulation, providing the possibility for under-compensatory (\eqn{b_c < 1}), compensatory (\eqn{b_c = 1}) (default) or over-compensatory (\eqn{b_c > 1}) dynamics. Only required if \code{StageStruct=FALSE}.}

\item{StageStruct}{\code{FALSE} (default) yields a population model with non-overlapping generations.\cr For a stage-structured population, this takes the corresponding parameter object generated by \code{\link[RangeShiftR]{StageStructure}}, which holds all demographic parameters.}

\item{ReproductionType}{0 = asexual / only female model (default)\cr1 = simple sexual model\cr2 = sexual model with explicit mating system}

\item{PropMales}{Required if \code{ReproductionType={1,2}}: Proportion of males in the population, between \eqn{0} and \eqn{1}. Defaults to \eqn{0.5}.}

\item{Harem}{Required if \code{ReproductionType=2}: Maximum harem size. The maximum number of pair bonds that a male can establish. \eqn{Harem = 1} (default) corresponds to monogamy, \eqn{0<Harem<1} to polyandry and \eqn{Harem>1} to polygyny.}
}
\value{
a parameter object of class "DemogParams"
}
\description{
For a simple non-structured population, set its basic demographic parameters here, i.e. the maximum growth rate (\code{Rmax}) and the competition coefficient (\code{bc}).
For a stage-structured population, define its corresponding parameters via \code{\link[RangeShiftR]{StageStructure}} and add it to Demography.\cr
\cr
Choose the Reproduction model that determines if sexes are considered implicitly or explicitly and if a mating system is used. If applicable, set the corresponding parameters, i.e. the proportion of males (\code{PropMales}) and the maximum harem size (\code{Harem}).
}
\details{
The following information regards the population dynamics of a \strong{non-structured} (\code{StageStruct=FALSE}) population.\cr
For more information on the population dynamics of a \strong{structured} population, see \code{\link[RangeShiftR]{StageStructure}}.

Populations with non-overlapping generations, i.e. with \strong{no stage-structure} are the appropriate way to model species that have discrete generations.
At each generation the life cycle comprises - in that order - reproduction, death of the adults and offspring dispersal.
Two parameters determine the nature of the demographic density dependence:
the carrying capacity \eqn{K} (given by the argument \code{K_or_DensDep} in the landscape module) and
the competition coefficient \eqn{b_c} (given by the argument \code{bc}).
These discrete generation models can be applied to asexual species, species for which it is assumed that females play the dominant role in spatial dynamics
and for species for which it is considered crucial to model both sexes explicitly.

\emph{Asexual / only-female models:}  (\code{ReproductionType=0})\cr
Recruitment is determined by a stochastic, individual-based formulation of the \insertCite{smith1973;textual}{RangeShiftR} population model, where the number of offspring produced by a single individual in the cell/patch \eqn{i} at time \eqn{t}, is drawn from the following distribution:\cr

\ifelse{html}{\out{&emsp;&emsp;N<sub>juv</sub>(i,t) = Poisson( R(i,t) / (1+|R(i,t) - 1| &ast; (N(i,t) / K(i,t))<sup>b<sub>c</sub></sup> ) ) } }{\deqn{N_juv(i,t) = Poisson( R(i,t) / (1 + |R(i,t) - 1|*( N(i,t) / K(i,t) )^bc ) ) } }

Here, \eqn{R(i,t)} is the maximum growth rate \code{Rmax} (obtained at very low density only) and \eqn{K(i,t)} is the carrying capacity
at patch \eqn{i} and time \eqn{t}.
Both \eqn{R(i,t)} and \eqn{K(i,t)} can vary in space and time, depending on the model setting. \ifelse{html}{\out{b<sub>c</sub>}}{\eqn{b_c}}
is the competition coefficient which describes the type of density regulation, providing the possibility for under-compensatory
(\ifelse{html}{\out{b<sub>c</sub> < 1}}{\eqn{b_c < 1}}), compensatory (\ifelse{html}{\out{b<sub>c</sub> = 1}}{\eqn{b_c = 1}}) or
over-compensatory (\ifelse{html}{\out{b<sub>c</sub> > 1}}{\eqn{b_c > 1}}) dynamics.\cr
\cr
\emph{Sexual models}\cr
In this type of models, individuals are explicitly characterized by their sex. The proportion of each sex in the population is controlled by setting the proportion of males (\code{PropMales}). There are two types of possible sexual sub-models:\cr
\cr
\emph{Simple mating system:}  (\code{ReproductionType=1})\cr
This is the simplest form of mate limitation. Each female individual is assumed to mate, as long as there is at least one male in the population. As for the asexual case, the Maynard Smith and Slatkin model is used to determine the expected number of
offspring produced by each female. To maintain equivalence between the asexual and sexual versions, the expected value of the Poisson distribution is multiplied by \eqn{2} (Lindström & Kokko 1998):\cr

\ifelse{html}{\out{&emsp;&emsp;N<sub>juv</sub>(i,t) = Poisson( 2 R(i,t) / (1+|R(i,t) - 1| &ast; (N(i,t) / K(i,t) )<sup>b<sub>c</sub></sup> ) ) } }{\deqn{N_juv(i,t) = Poisson( 2 R(i,t) / (1 + |R(i,t) - 1|*( N(i,t) / K(i,t) )^bc ) ) } }

\emph{Complex mating system:}  (\code{ReproductionType=2})\cr
More complex and flexible mating system. Mating is explicitly modelled through a mating function
\insertCite{lindstrom1998,legendre2004,bessa2010}{RangeShiftR}, where the number of mated females \eqn{c} is given by:
\deqn{c = f * min(1, 2hm/(f+hm) )}
where \eqn{f} and \eqn{m} are the numbers of potentially reproductive females and males, respectively, and \eqn{h} is the maximum harem size \code{Harem}, i.e. the maximum number of pair bonds that a male can establish. \eqn{h=1} corresponds to monogamy, \eqn{0<h<1} to polyandry
and \eqn{h>1} to polygyny. Each potentially reproductive female has a probability of reproducing \ifelse{html}{\out{p<sub>r</sub>}}{\eqn{p_r}} given by:
\deqn{p_r=c/f}
A Bernoulli trial \eqn{Bern}\ifelse{html}{\out{(p<sub>r</sub>)}}{\eqn{(p_r)}} determines if the female reproduces or not. For females that reproduce, the number of offspring is determined through the Poisson distribution as stated above. Hence, the specification of the mating system determines the probability for each female to reproduce. However, no explicit pair bonds are formed, and in the cases where traits inheritance is
involved, the father (of all the offspring produced by a single female in a single reproductive event) is selected randomly from the males in the population.\cr
\cr
Populations with overlapping generations, i.e. with \strong{with stage-structure} are the appropriate choice for species in which generations can overlap and individuals can be classified in different stages
(e.g. immature vs. breeding individuals) differing in their demographic parameters. Individuals are characterized by their age and stage. Each stage has a certain fecundity (which will be used as expected value of the Poisson distribution for reproduction), survival
and probability of developing to the next stage. The parameters are provided as classical transition matrices \insertCite{caswell2001}{RangeShiftR} through a \code{\link[RangeShiftR]{StageStructure}}
parameter object. For more information on the population dynamics in this case, see the details there.
}
\references{
\insertAllCited{}
}
\author{
Anne-Kathleen Malchow
}
