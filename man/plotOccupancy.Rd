% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_handling.R
\name{plotOccupancy}
\alias{plotOccupancy}
\title{Plot Occupancy}
\usage{
plotOccupancy(s, ...)
}
\arguments{
\item{s}{RSparams object or a data.frame in the 'range' file format}

\item{dirpath}{RS directory path; required if \code{s} is a \code{RSparams}}

\item{sd}{plot standard deviation? (default is \code{FALSE})}

\item{replicates}{plot the replicates? (default is \code{TRUE})}

\item{ylim}{upper limit to the y-axis}
}
\description{
Uses the RangeShiftR output data 'range' to generate occupancy time series.
Plots the mean occupancy over all replicates, and optionally the standard deviation and/or the single replicates.
}
