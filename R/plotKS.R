#' Histogram of the KS values
#' 
#' Plots the histogram of the KS values for the combination of all runs. The MCMC must have at least 2 runs
#' 
#' @param x A list of convenience.diag type
#' @param precision The precision of the mean estimates. Default is 0.01
#' @param color The color to fill the histogram bars, default is "grey"
#' @param ... (various) Additional arguments passed to plot().
#' 
#' @return Histogram
#' 
#' @export


plotKS <- function(x, precision = 0.01, color = "grey", ...){
  
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
               col = color,
               yaxs ="i",
               las = 1,
               ...)
  
  plot <- box("plot", "solid")
  plot <- abline( v = minimumKS, col = "red", lwd= 2, lty = 2)
  plot <- axis(1, at = round(minimumKS, digits = 2))
               
  }

