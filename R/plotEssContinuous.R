#' Histogram of the ESS values
#' 
#' Plots the histogram of the ESS values for the continuous parameters
#' 
#' @param x A list of convenience.diag type
#' @param precision The precision of the mean estimates. Default is 0.01
#' @param color The color to fill the histogram bars, default is "grey"
#' @param ... (various) Additional arguments passed to plot().
#' 
#' @return Histogram
#' 
#' @export

plotEssContinuous <- function(x, precision = 0.01, fill_color = NULL, filename = NULL, ...){
  
  if( is.null(fill_color) ){
    fill_color <- "cadetblue4"
  }
  
  if( !(is.null(filename)) ){
    pdf(file = filename, width = 4.5, height = 4.5)
  }
  
  minimumESS <- minESS(precision)
  ESS_values <- vector()
  
  for (i in 1:ncol(x$continuous_parameters$ess)) {
    ESS_values <- x$continuous_parameters$ess[,i]
  }
  
  y_topLim <- max(hist(ESS_values, breaks = 20, plot = FALSE)$counts)
  
  par(mar = c(3.9, 2.2, 2.1, 1.0))
  plot <- hist( ESS_values, 
                xlab = "ESS", 
                ylab = NA,
                main = "Histogram of ESS for continuous parameters",
                xlim = c(0, (max(minimumESS, ESS_values)+1000) ),
                ylim = c(0, y_topLim+1),
                breaks = 20,
                col = fill_color,
                las = 1,
                border=F,
                ...)
  plot <- lines(x = c(minimumESS,minimumESS),y=c(0,y_topLim+1), col = "antiquewhite4", lwd= 2, lty=2)
  #plot <- axis(1, at = minimumESS)
  
  if( !(is.null(filename)) ){
    dev.off()
  } 
}