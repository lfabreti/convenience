#' Histogram of the ESS values
#' 
#' Plots the histogram of the ESS values for the continuous parameters or the splits
#' 
#' @param x A list of convenience.diag type
#' @param trees A boolean to set if the histogram is for the splits
#' @param precision The precision of the mean estimates. Default is 0.01
#' 
#' @return Histogram
#' 
#' @export

plotESS.hist <- function(x, trees,precision = 0.01){
  
  minimumESS <- minESS(precision)
  ESS_values <- vector()
  
  if ( trees == TRUE){
    for (i in 1:length(x$tree_parameters$ess)) {
      ESS_values <- x$tree_parameters$ess[[i]]
      
    }
    plot <- hist(ESS_values, 
                 xlab = "ESS", 
                 main = "Histogram of ESS for splits", 
                 xlim = c(min(minimumESS, ESS_values), (max(minimumESS, ESS_values)+1000) ),
                 col = "grey",
                 breaks = ((max(ESS_values) - min(ESS_values))/minimumESS),
                 yaxs = "i")
    plot <- box("plot", "solid")
    plot <- abline(v = minimumESS, col = "red", lwd= 2)
    
  } else{
    for (i in 1:ncol(x$continuous_parameters$ess)) {
      ESS_values <- x$continuous_parameters$ess[,i]
      
    }
    plot <- hist( ESS_values, 
                 xlab = "ESS", 
                 main = "Histogram of ESS for continuous parameters",
                 xlim = c(0, (max(minimumESS, ESS_values)+1000) ),
                 col = "grey",
                 yaxs="i")
    plot <- box("plot", "solid")
    plot <- abline(v = minimumESS, col = "red", lwd= 2)
  }
  
}