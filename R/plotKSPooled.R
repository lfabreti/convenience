#' Histogram of the KS values
#' 
#' Plots the histogram of the KS values for the one-on-one comparison of runs. The MCMC must have at least 3 runs
#' 
#' @param x A list of convenience.diag type
#' @param bins The number of bins to determine the intervals of the histogram
#' @param precision The precision of the mean estimates. Default is 0.01
#' @param ... (various) Additional arguments passed to plot().
#' 
#' @return Histogram
#' 
#' @export

plotKSPooled <- function(x, bins, precision = 0.01, ...){
  
  minimumESS <- minESS(precision)
  minimumKS <- ksThreshold(0.01,minimumESS)
  
  ks_values <- list()
  for (i in 1:nrow(x$continuous_parameters$compare_runs)) {
    ks_values[[i]] <- as.numeric(x$continuous_parameters$compare_runs[i,])
  }
  
  increment <- (max( unlist(ks_values)) - min(unlist(ks_values)) ) / bins
  break_values <- seq( from = min(unlist(ks_values)), to = max( unlist(ks_values)), by = increment )
  
  colors_hist <- viridis(length(ks_values))
  
  par(mar = c(5.1, 4.1, 4.1, 8.1), xpd=F)
  plot <- hist(ks_values[[1]],
               xlab = "KS score", 
               main = "KS histogram",
               xlim = c( min(minimumKS, unlist(ks_values)), max(minimumKS, unlist(ks_values)) ),
               ylim = c( 0, length(ks_values[[1]]) ),
               breaks = break_values,
               col = scales::alpha(colors_hist[1],.5),
               border = T,
               yaxs="i",
               las = 1,
               ...)
  plot <- box("plot", "solid")
  
  for (i in 2:length(ks_values)) {
    hist(ks_values[[i]], 
         add=T, 
         col = scales::alpha(colors_hist[i],.5),
         breaks = break_values,
         border = T,
         yaxs="i",
         las = 1,
         ...)
  }
  plot <- abline( v = minimumKS, col = "red", lwd=2, lty = 2)
  plot <- axis(1, at = round(minimumKS, digits = 2))
  
  legend("topright",
         legend = rownames(x$continuous_parameters$compare_runs),
         fill = scales::alpha(colors_hist,.5),
         box.lty = 0,
         inset=c(-0.35,0),
         cex = 0.9,
         xpd=T
         )
  
}
