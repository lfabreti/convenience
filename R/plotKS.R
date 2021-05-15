#' Histogram of the KS values
#' 
#' Plots the histogram of the KS values for the combination of all runs. The MCMC must have at least 2 runs
#' 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline hist layout legend lines par points polygon rect title
#' 
#' @param x A list of convenience.diag type
#' @param precision The precision of the mean estimates. Default is 0.01
#' @param fill_color The color to fill the histogram bars
#' @param filename The name of the file to save the plot
#' @param ... (various) Additional arguments passed to plot().
#' 
#' @return Histogram
#' 
#' @export


plotKS <- function(x, precision = 0.01, fill_color = NULL, filename = NULL, ...){
  
  col_threshold <- "gray69"
  
  if( is.null(fill_color) ){
    fill_color <- "cadetblue4"
  }
  
  if( !(is.null(filename)) ){
    pdf(file = filename, width = 4.5, height = 4.5)
  }
  
  minimumESS <- minESS(precision)
  minimumKS <- ksThreshold(0.01,minimumESS)
  
  KS_values <- vector()
  
  for (i in 1:ncol(x$continuous_parameters$compare_runs)) {
    KS_values <- c(KS_values, x$continuous_parameters$compare_runs[[i]])
  }
  y_topLim <- max(hist(KS_values, plot = FALSE)$counts)
  
  par(mar = c(4.1, 3.9, 2.1, 1.0))
  
  plot <- plot(NA,
               xlab = "Kolmogorov-Smirnov score",
               ylab = "Counts",
               main = "Kolmogorov-Smirnov test",
               cex.main = 0.9,
               xlim = c( (min(minimumKS, KS_values) - 0.01), (max(minimumKS, KS_values) + 0.01)  ),
               ylim = c(0,y_topLim + 1),
               las=1,
               bty="l")
  plot <- rect(xleft = min(minimumKS, KS_values) - 0.01, ybottom = 0, xright = minimumKS, ytop = y_topLim+1, border = NA, col = "gray89")
  plot <- lines(x = c(minimumKS,minimumKS),y=c(0,y_topLim+1), col = col_threshold, lwd= 2, lty = 2)
  
  
  plot <- hist(KS_values, 
               xlim = c( (min(minimumKS, KS_values) - 0.01), (max(minimumKS, KS_values) + 0.01)  ),
               ylim = c(0,y_topLim + 1),
               col = fill_color,
               border=F,
               las = 1,
               add=T,
               ...)

  
  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}

