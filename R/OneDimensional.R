#' Creating design matrix with 1d wavelet basis functions as predictors
#'
#' Takes in the locations at which data is observed and creates a design matrix where each column corresponds
#' to a given wavelet basis function evaluated at the observed data locations. This function is taken from the
#' code provided in Wand and Ormerod (2011).
#'
#' @param x              Vector of locations at which we observe data
#' @param numLevels      Number of wavelet levels. Should be an integer between 2 and 10
#' @param filterNumber   Which Daubechies wavelet filter is desired. Default value is 5
#' @param resolution     The number of points in the grid used to evaluate wavelet basis functions. Should be
#'                       a large integer that is a power of 2    
#'
#' @return An N x K matrix containing wavelet basis functions evaluated at observed data
#'
#' @export
#' @examples
#' 
#' n <- 1000
#' x <- runif(n)
#' 
#' numLevels <- 3
#' k = 2^numLevels - 1
#' Zx <- cbind(rep(1, n), ZDaub(x, numLevels=numLevels))

Z1d <- function(x,numLevels=6,filterNumber=5,
                  resolution=16384)
{
  # Load required package:
  
  range.x=range(x)
      
  # Ensure that the number of levels is `allowable'.
  
  if (!any(numLevels==(1:10)))
    stop("Number of levels should be between 2 and 10.")
  
  # Ensure the resolution value is a power of 2 and within a reasonable range.
  
  if (!any(resolution==(2^(10:20))))
    stop("Resolution value should be a power of 2, with the
         power between 10 and 20.")
  
  # Transform x to the unit interval and obtain variables
  # required for linear interpolation:
  
  xUnit <- (x - range.x[1])/(range.x[2] - range.x[1])
  xUres <- xUnit*resolution
  fXuRes <- floor(xUres)
  
  # Set filter and wavelet family  
  
  family <- "DaubExPhase"
  K <- 2^numLevels - 1
  
  # Create a dummy wavelet transform object
  
  wdObj <- wavethresh::wd(rep(0,resolution),filter.number=filterNumber,
              family="DaubExPhase")
  
  Z <- matrix(0,length(x),K)
  for (k in 1:K)
  {
    # Create wobj so that it contains the Kth basis
    # function of the Z matrix with `resolution' regularly 
    # spaced points:
    
    putCobj <- wavethresh::putC(wdObj,level=0,v=0)
    putCobj$D <- putCobj$D*0
    putCobj$D[resolution-k] <- 1
    
    # Obtain kth column of Z via linear interpolation
    # of the wr(putCobj) grid values:
    
    wtVec <- xUres - fXuRes
    wvVec <- wr(putCobj)
    wvVec <- c(wvVec,rep(wvVec[length(wvVec)],2))
    Z[,k] <- sqrt(resolution)*((1 - wtVec)*wvVec[fXuRes+1]
                               + wtVec*wvVec[fXuRes+2])
  }
  
  # Create column indices to impose "left-to-right" ordering
  # within the same level:
  
  newColInds <- 1
  for (ell in 1:(numLevels-1))
    newColInds <- c(newColInds,(2^(ell+1)-1):(2^(ell)))
  
  Z <- Z[,newColInds]
  
  return(Z)
}