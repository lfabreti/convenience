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

plotKSPooled <- function(x, bins, precision = 0.01, filename = NULL, ...){
  
  if( !(is.null(filename)) ){
    pdf(file = filename, width = 6, height = 6)
  }
  
  minimumESS <- minESS(precision)
  minimumKS <- ksThreshold(0.01,minimumESS)
  
  ks_values <- list()
  for (i in 1:nrow(x$continuous_parameters$compare_runs)) {
    ks_values[[i]] <- as.numeric(x$continuous_parameters$compare_runs[i,])
  }
  
  increment <- (max( unlist(ks_values)) - min(unlist(ks_values)) ) / bins
  break_values <- seq( from = min(unlist(ks_values)), to = max( unlist(ks_values)), by = increment )
  
  colors_hist <- viridis(length(ks_values))
  
  par(mar = c(3.9, 3.9, 3.1, 5.0), xpd=F)
  plot <- plot(NA,
               xlab = "Kolmogorov-Smirnov score",
               ylab = "Counts",
               main = "Kolmogorov-Smirnov test",
               xlim = c( (min(minimumKS, unlist(ks_values)) - 0.01), (max(minimumKS, unlist(ks_values)) + 0.01) ),
               ylim = c( 0, (length(ks_values[[1]])- 2) ),
               las=1,
               bty="l")
  
  plot <- lines(x = c(minimumKS,minimumKS),y=c(0,y_topLim+1), col = col_threshold, lwd= 2, lty = 2)
  plot <- rect(xleft = min(minimumKS, unlist(ks_values)) - 0.01, ybottom = 0, xright = minimumKS, ytop = length(ks_values[[1]])- 2, border = NA, col = "gray89")
  
  
  plot <- hist(ks_values[[1]],
               xlim = c( (min(minimumKS, unlist(ks_values)) - 0.01), (max(minimumKS, unlist(ks_values)) + 0.01) ),
               ylim = c( 0, (length(ks_values[[1]])- 2) ),
               breaks = break_values,
               col = scales::alpha(colors_hist[1],.7),
               border = T,
               yaxs="i",
               las = 1,
               add=T,
               ...)
  
  for (i in 2:length(ks_values)) {
    hist(ks_values[[i]], 
         add=T, 
         col = scales::alpha(colors_hist[i],.7),
         breaks = break_values,
         border = T,
         yaxs="i",
         las = 1,
         ...)
  }
  plot <- abline( v = minimumKS, col ="antiquewhite4", lwd=2, lty = 2)
  
  legend("topright",
         legend = rownames(x$continuous_parameters$compare_runs),
         fill = scales::alpha(colors_hist,.7),
         box.lty = 0,
         inset=c(-0.18,0),
         cex = 0.9,
         xpd=T
         )
  
  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}
