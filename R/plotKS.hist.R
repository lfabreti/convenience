#' Histogram of the KS values
#' 
#' Plots the histogram of the KS values for the combination of all runs. The MCMC must have at least 2 runs
#' 
#' @param x A list of convenience.diag type
#' @param precision The precision of the mean estimates. Default is 0.01
#' 
#' @return Histogram
#' 
#' @export


plotKS.hist <- function(x, precision = 0.01){
  
  minimumESS <- minESS(precision)
  minimumKS <- ksThreshold(1.95,minimumESS, minimumESS)
  
  KS_values <- vector()
  
  for (i in 1:ncol(x$continuous_parameters$compare_runs)) {
    KS_values <- c(KS_values, x$continuous_parameters$compare_runs[[i]])
  }
  plot <- hist(KS_values, 
               xlab = "KS score", 
               main = "KS histogram",
               xlim = c( min(minimumKS, KS_values), max(minimumKS, KS_values) ),
               col = "grey",
               yaxs="i")
  plot <- box("plot", "solid")
  plot <- abline( v = minimumKS, col = "red", lwd= 2)
               
  }

