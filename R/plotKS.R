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


plotKS <- function(x, precision = 0.01, fill_color = NULL, filename = NULL, ...){
  
  if( is.null(fill_color) ){
    fill_color <- "cadetblue4"
  }
  
  if( !(is.null(filename)) ){
    pdf(file = filename)
  }
  
  minimumESS <- minESS(precision)
  minimumKS <- ksThreshold(0.01,minimumESS)
  
  KS_values <- vector()
  
  for (i in 1:ncol(x$continuous_parameters$compare_runs)) {
    KS_values <- c(KS_values, x$continuous_parameters$compare_runs[[i]])
  }
  y_topLim <- max(hist(KS_values, plot = FALSE)$counts)
  
  par(mar = c(3.9, 2.2, 2.1, 0.1))
  plot <- hist(KS_values, 
               xlab = "KS score",
               ylab = NA,
               main = "KS histogram",
               xlim = c( (min(minimumKS, KS_values) - 0.01), (max(minimumKS, KS_values) + 0.01)  ),
               ylim = c(0,y_topLim + 1),
               col = fill_color,
               border=F,
               las = 1,
               ...)

  plot <- lines(x = c(minimumKS,minimumKS),y=c(0,y_topLim+1), col =  "antiquewhite4", lwd= 2, lty = 2)
  plot <- axis(1, at = round(minimumKS, digits = 3))
  
  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}

