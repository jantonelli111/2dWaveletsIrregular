#' Creating design matrix with 2d wavelet basis functions as predictor
#'
#' Takes in two matrices: The N x K matrices generated from 1d wavelet basis functions in each direction
#'
#' @param Zx      An N x K matrix where the columns are the 1d wavelet basis functions in the X direction
#'                evaluated at the N locations of the data
#' @param Zy      An N x K matrix where the columns are the 1d wavelet basis functions in the Y direction
#'                evaluated at the N locations of the data
#'
#' @return An N x K^2 matrix containing products of wavelet basis functions
#'
#' @export
#' @examples
#' 
#' n <- 1000
#' x <- runif(n)
#' y <- runif(n)
#' 
#' numLevels <- 3
#' k = 2^numLevels - 1
#' Zx <- cbind(rep(1, n), ZDaub(x, numLevels=numLevels))
#' Zy <- cbind(rep(1, n), ZDaub(y, numLevels=numLevels))
#' Zxy <- Z2d(Zx,Zy)


Z2d <- function(Zx,Zy) {
  n = dim(Zx)[1]
  k = dim(Zx)[2]
  Zxy <- matrix(NA, nrow=n, ncol=k^2)
  for (i in 1 : n) {
    for (j1 in 1 : k) {
      for (j2 in 1 : k) {
        Zxy[i, ((j1 - 1)*k) + j2] = Zx[i, j1]*Zy[i, j2]
      }
    }
  }
  return(Zxy)
}




#' Setting wavelet coefficients to zero that correspond to certain levels
#'
#' Takes in a vector of coefficients and returns the same vector, but with entries forced to zero that correspond
#' to the desired wavelet levels
#'
#' @param beta           Vector of wavelet coefficients
#' @param numLevels      Number of wavelet levels. Should be an integer between 2 and 10
#' @param remove.x       Vector of wavelet levels in the X direction that you want removed. Should contain
#'                       numbers between 1 and numLevels
#' @param remove.y       Vector of wavelet levels in the Y direction that you want removed. Should contain
#'                       numbers between 1 and numLevels
#'
#' @return A vector of wavelet coefficients with those corresponding to basis functions involving levels in
#'         remove.x or remove.y set to zero
#'
#' @export
#' @examples
#' 
#' n <- 1000
#' x <- runif(n)
#' y <- runif(n)
#' 
#' numLevels <- 3
#' k = 2^numLevels - 1
#' Zx <- cbind(rep(1, n), ZDaub(x, numLevels=numLevels))
#' Zy <- cbind(rep(1, n), ZDaub(y, numLevels=numLevels))
#' Zxy <- Z2d(Zx,Zy)

threshold2d <- function(beta, numLevels, remove.x, remove.y) {
  threshold.beta <- beta
  k <- 2^numLevels
  j1.remove <- c()
  j2.remove <- c()
  for (i in 1 : length(remove.x)) {
    temp.remove <- remove.x[i]
    j1.remove <- c(j1.remove, (2^(temp.remove-1) + 1) : (2^temp.remove))
  }
  
  for (i in 1 : length(remove.y)) {
    temp.remove <- remove.y[i]
    j2.remove <- c(j2.remove, (2^(temp.remove-1) + 1) : (2^temp.remove))
  }
  
  for (j1 in 1 : k) {
    for (j2 in 1 : k) {
      if (j1 %in% j1.remove | j2 %in% j2.remove) threshold.beta[((j1 - 1)*k) + j2] <- 0
    }
  }
  
  return(threshold.beta)
}



#' Setting wavelet coefficients to zero that correspond to certain levels
#'
#' Takes in a vector of coefficients and returns the same vector, but with entries forced to zero that correspond
#' to the desired wavelet levels
#'
#' @param x              Vector of locations of the data in the first direction
#' @param y              Vector of locations of the data in the second direction
#' @param f              Vector of data values observed at the given locations
#' @param numLevels      Number of wavelet levels. Should be an integer between 2 and 10
#'
#' @return A matrix of wavelet basis functions evaluated at the data as well as the estimated wavelet coefficients
#'
#' @export
#' @examples
#' 
#' n <- 1000
#' x <- runif(n)
#' y <- runif(n)
#' 
#' numLevels <- 3
#' k = 2^numLevels - 1
#' Zx <- cbind(rep(1, n), ZDaub(x, numLevels=numLevels))
#' Zy <- cbind(rep(1, n), ZDaub(y, numLevels=numLevels))
#' Zxy <- Z2d(Zx,Zy)

Irregular2dWavelet <- function(x, y, f, numLevels) {
  n = length(x)
  k = 2^numLevels - 1
  Zx <- cbind(rep(1, n), Z1d(x, numLevels=numLevels))
  Zy <- cbind(rep(1, n), Z1d(y, numLevels=numLevels))
  Zxy <- Z2d(Zx,Zy)
  
  fCenter = f - mean(f)
  cvGLMNET <- cv.glmnet(Zxy[,-1], fCenter, intercept=FALSE)
  betauHat <- c(mean(f), as.numeric(coef(cvGLMNET, s="lambda.1se"))[-1])
  
  l = list(Zxy = Zxy, beta=betauHat)
  return(l)
}