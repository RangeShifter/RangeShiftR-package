% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_DispersalParams.R
\docType{class}
\name{CorrRW}
\alias{CorrRW}
\title{Set up a Correlated Random Walk}
\usage{
CorrRW(StepLength = 1, Rho = 0.5,
      IndVar = FALSE,
      StraightenPath = FALSE,
      StepMort = 0.0)
}
\arguments{
\item{StepLength}{Step length given in meters, defaults to \eqn{1}.\cr If \code{IndVar=TRUE}, expects a vector of length three
specifying (Mean, SD, TraitScaleFactor) of \code{StepLength}.}

\item{Rho}{Correlation parameter \eqn{ρ}, defaults to \eqn{0.5}. Must be in the open interval \eqn{(0,1)}.\cr If \code{IndVar=TRUE},
expects a vector of length three specifying (Mean, SD, TraitScaleFactor) of \code{Rho}.}

\item{IndVar}{Individual variability in CorrRW traits (i.e. \code{StepLength} and \code{Rho})? Defaults to \code{FALSE}.}

\item{StraightenPath}{Straighten path after decision not to settle in a patch? Defaults to \code{TRUE}, see Details below.}

\item{StepMort}{Per-step mortality probability. Can be either \emph{constant}, in which case a single numeric is expected (the default, with
value \eqn{0.0}) or \emph{habitat-specific}, in which case a numeric vector is expected with a length of, respectively, \code{Nhabitats} for an
\code{\link[RangeShiftR]{ImportedLandscape}} with habitat codes (i.e. \code{HabPercent=FALSE})) or length \eqn{2} for an
\code{\link[RangeShiftR]{ArtificialLandscape}} (mortality probabilities for matrix and habitat cells).\cr
All values must be within the half-open interval \eqn{[0,1)}.\cr
For an imported habitat quality landscape (\code{HabPercent=TRUE}), only constant per-step mortality is allowed.}
}
\value{
a parameter object of class "CorrRW"
}
\description{
A method to describe \code{\link[RangeShiftR]{Transfer}}:
A simple correlated random walk without any bias; implemented in continuous space on the top of the landscape grid.
}
\details{
Individuals take steps of a constant \code{StepLength}; the direction is sampled from a wrapped Cauchy distribution having a
correlation parameter \eqn{Rho} in the range \eqn{0} to \eqn{1} \insertCite{barton2009evolution,zollner1999search}{RangeShiftR}.
As for \code{\link[RangeShiftR]{SMS}}, all individuals take each step
simultaneously. In the case of patch-based models,
\eqn{Rho} is automatically set to \eqn{0.99} until the individual steps outside its natal patch, after which the value of
\eqn{Rho} set by the user is restored. \cr
The \code{StepLength} and \eqn{Rho} can be set to vary between individuals and evolve (set \code{IndVar=TRUE}).
In this case, each individual exhibits two traits for these two parameters.
For each trait the initial mean and standard deviations must be set, as well as the TraitScaleFactor (see \code{\link[RangeShiftR]{Settlement}}),
instead of only one constant value each.\cr
Note that the step length may not evolve below one fifth of
the landscape resolution, and correlation may not evolve above \eqn{0.999}. \cr
Per-step mortality is not allowed to vary between individuals or to evolve. \cr
There is no implementation of sex- or stage-specific CRW.

When an individual arrives in a non-natal patch and decides not to settle there (as a result of a density-dependent or mate-finding settlement
rule, see \code{\link[RangeShiftR]{Settlement}}), then there is the option that its path is straightened (\code{StraightenPath=TRUE}). This means
that it leaves the patch as soon as possible in order to search for another patch. This is achieved by increasing its path correlation to
\code{Rho}\eqn{=0.999}. However, in certain types of model, e.g. when arbitrary contiguous patches have been defined for what is basically a continuous
population, this can lead to the path always being straightened, as an individual enters a new patch as soon as it has left the one it has
rejected. In such cases, it is advisable to disable the feature (\code{StraightenPath=FALSE}), although care must be taken that individuals
do not become trapped in patches surrounded by very high cost matrix.

\emph{Mortality}\cr
There are two main sources of mortality: First, dispersal mortality can arise as a result of individuals failing to reach suitable habitat. Some individuals may fail
to find suitable habitat before they use up their maximum number of movement steps (\code{MaxSteps} in \code{\link[RangeShiftR]{Settlement}}).
In this first case, dispersal mortality clearly depends upon the proportion of suitable habitat in the landscape and will increase as the
availability of habitat declines.\cr
A second source of dispersal mortality can be specified by the user in form of a per-step probability of mortality (\code{StepMort}.
This can be useful for representing mortality risks that increase with distance or time spent travelling \insertCite{bonte2012costs}{RangeShiftR}.

Note that the total dispersal mortality experienced will be the sum of the mortalities due to the two sources identified above and, in parameterising the model,
it will be important to recognize this such that dispersal mortality is not double-accounted.
}
\references{
\insertAllCited{}
}
\author{
Anne-Kathleen Malchow
}
