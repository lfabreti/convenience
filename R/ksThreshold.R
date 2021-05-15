#' KS threshold
#'
#' Calculates the threshold for the Kolmogorov-Smirnov test
#'  
#' @param alpha The chosen level of the test
#' @param ess The ESS threshold
#' 
#' @return The calculated threshold
#' 
#' @export


ksThreshold <- function(alpha, ess){
  c_alpha <- sqrt( -(log(alpha/2)) * 0.5) 
  return( (c_alpha*sqrt((2*ess)/(ess*ess))) )
}