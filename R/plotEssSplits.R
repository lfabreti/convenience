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

plotEssSplits <- function(x, precision = 0.01, fill_color = NULL, filename = NULL, ...){
  
  if( is.null(fill_color) ){
    fill_color <- "coral3"
  }
  
  if( !(is.null(filename)) ){
    pdf(file = filename)
  }
  
  minimumESS <- minESS(precision)
  ESS_values <- vector()
  
  for (i in 1:length(x$tree_parameters$ess)) {
    ESS_values <- x$tree_parameters$ess[[i]]
  }
  ESS_values <- ESS_values[!is.na(ESS_values)]
  y_topLim <- max(hist(ESS_values, plot = FALSE)$counts)
  
  par(mar = c(3.9, 2.2, 2.1, 0.1))
  plot <- hist(ESS_values, 
               xlab = "ESS", 
               ylab = NA,
               main = "Histogram of ESS for splits", 
               xlim = c(0, (max(minimumESS, ESS_values)+1000) ),
               ylim = c(0, y_topLim + 1),
               col = fill_color,
               las = 1,
               border=F,
               ...)
  plot <- lines(x = c(minimumESS,minimumESS),y=c(0,y_topLim+1), col =  "antiquewhite4", lwd= 2, lty=2)
  plot <- axis(1, at = minimumESS)
  
  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}