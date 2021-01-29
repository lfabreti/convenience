#' Histogram of the ESS values
#' 
#' Plots the histogram of the ESS values for the splits
#' 
#' @param x A list of convenience.diag type
#' @param precision The precision of the mean estimates. Default is 0.01
#' @param color The color to fill the histogram bars, default is "grey"
#' @param ... (various) Additional arguments passed to plot().
#' 
#' @return Histogram
#' 
#' @export

plotEssSplits <- function(x, precision = 0.01, color = "grey", filename = NULL, ...){
  
  if( !(is.null(filename)) ){
    pdf(file = filename)
  }
  
  minimumESS <- minESS(precision)
  ESS_values <- vector()
  
  for (i in 1:length(x$tree_parameters$ess)) {
    ESS_values <- x$tree_parameters$ess[[i]]
    
  }
  plot <- hist(ESS_values, 
               xlab = "ESS", 
               main = "Histogram of ESS for splits", 
               xlim = c(min(minimumESS, ESS_values), (max(minimumESS, ESS_values)+1000) ),
               col = color,
               breaks = ((max(ESS_values) - min(ESS_values))/minimumESS),
               yaxs = "i",
               las = 1,
               ...)
  plot <- box("plot", "solid")
  plot <- abline(v = minimumESS, col = "red", lwd= 2, lty=2)
  plot <- axis(1, at = minimumESS)
  
  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}