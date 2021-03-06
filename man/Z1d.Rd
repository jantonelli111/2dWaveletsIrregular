% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OneDimensional.R
\name{Z1d}
\alias{Z1d}
\title{Creating design matrix with 1d wavelet basis functions as predictors}
\usage{
Z1d(x, numLevels = 6, filterNumber = 5, resolution = 16384)
}
\arguments{
\item{x}{Vector of locations at which we observe data}

\item{numLevels}{Number of wavelet levels. Should be an integer between 2 and 10}

\item{filterNumber}{Which Daubechies wavelet filter is desired. Default value is 5}

\item{resolution}{The number of points in the grid used to evaluate wavelet basis functions. Should be
a large integer that is a power of 2}
}
\value{
An N x K matrix containing wavelet basis functions evaluated at observed data
}
\description{
Takes in the locations at which data is observed and creates a design matrix where each column corresponds
to a given wavelet basis function evaluated at the observed data locations. This function is taken from the
code provided in Wand and Ormerod (2011).
}
\examples{

n <- 1000
x <- runif(n)

numLevels <- 3
k = 2^numLevels - 1
Zx <- cbind(rep(1, n), ZDaub(x, numLevels=numLevels))
}

