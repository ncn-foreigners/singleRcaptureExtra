% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate}
\alias{simulate}
\alias{simulate.singleRforeign}
\title{Generating data from foreign objects}
\usage{
\method{simulate}{singleRforeign}(object, nsim = 1, seed = NULL, truncated = FALSE, ...)
}
\arguments{
\item{object}{an object representing a fitted model.}

\item{nsim}{a numeric scalar specifying:
\itemize{
\item number of response vectors to simulate in \code{simulate.singleRStaticCountData}, defaults to \code{1L}.
\item number of units to draw in \code{simulate.singleRfamily}, defaults to \code{NROW(eta)}.
}}

\item{seed}{an object specifying if and how the random number generator should be initialized (‘seeded’).}

\item{truncated}{logical value indicating whether to sample from truncated or
full distribution.}

\item{...}{additional optional arguments.}
}
\value{
a \code{data.frame} with \code{n} rows and \code{nsim} columns.
}
\description{
An S3 method for \code{stats::simulate} to handle \code{singleRforeign} and
\code{singleRfamily} classes.
}
\seealso{
\code{\link[stats:simulate]{stats::simulate()}} \code{\link[singleRcapture:estimatePopsize]{singleRcapture::estimatePopsize()}}
}
