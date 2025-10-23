#' Random Generation from a Bivariate Flat-Top Normal Distribution
#'
#' This function generates samples from a bivariate flat-top normal distribution.
#'
#' @param n A numerical value defining the number of observations to be sampled.
#' @param mean A vector of 2 values defining the mean (centre) of the distribution. 
#' @param sd A numerical value defining the standard deviation of the distribution. 
#' @param range A numerical value defining the radius of a circle in which the flat-top of the kernel is defined.
#' @param truncDist A numerical value defining the distance beyond which values are considered too far from the center of the distribution to be considered reasonable to consider.
#' 
#' @details
#' 
#' Rejection sampling is used to carry out the sampling using as a starting point a random sample from a uniform distribution bounded by the truncation distance around the mean. Because the distribution function is continuous and unbounded, conceptually, the truncation could be +/- infinite in both dimensions, but since rejection-sampling is used this would not be computationally efficient. In this respect, we recommend to choose the truncation distance (`truncDist`) using a reasonable value for the sampling to be both computationally efficient and representative (so, not too large but not too small). 
#' 
#' Even if the function samples from a bivariate flat-top normal distribution, it is designed to focus only on a special case where the standard deviation (`sd`) and the radius of a circle in which the flat-top of the kernel is defined (`range`) is unique. In doing so, the function assumes weights are equal in both direction. 
#' 
#' If the `range` is set at 0
#' 
#' @return A two-column matrix with the sample from the distribution 
#' 
#' @author F. Guillaume Blanchet, Département de biologie, mathématiques et des sciences de la santé communautaire, Université de Sherbrooke, Canada. 
#' 
#' @examples
#' rFlatTopNorm(10)
#' rFlatTopNorm(5, mean = c(10, 12), sd = 3, range = 1, truncDist = 10)
#' 
#' @export
rFlatTopNorm<-function(n, mean = c(0, -1), sd = 1, range = 2, truncDist = 5){

  # Result object
  result <- matrix(NA, nrow = 0, ncol = 2)
  colnames(result) <- c("x", "y")

  # Circle defining the flat-top portion of the distribution (Circle built with 100 coordinates)
  xyCircle <- seq(0, 2 * pi, l = 100)
  cx <- mean[1] + range * cos(xyCircle)
  cy <- mean[2] + range * sin(xyCircle)
  
  # Kernel definition
  while(nrow(result) < n){
    # Define sampling location
    x <- runif(1, min = mean[1] - truncDist,  max = mean[1] + truncDist)
    y <- runif(1, min = mean[2] - truncDist,  max = mean[2] + truncDist)
    
    # Find the distance between the mean and the sampling location
    distSmpl <- sqrt((mean[1] - x)^2 + (mean[2] - y)^2)
    
    # If the sampling location is within the flat-top portion of the kernel
    if(distSmpl < range){
      acceptProb <- 1
      # If it is not
    }else{
      # Bivariate normal distribution outside the flat top
      acceptProb <- exp(-(min(sqrt((((x - cx)^2)+((y - cy)^2)) / (sdDistr^2))))^2)
    }
    
    # Rejection sampling
    if(runif(1) < acceptProb){
      result <- rbind(result, c(x, y))
    }
  }
  
  return(result)
}