#' Plots the difference between splits of different runs
#' 
#' @param output A list of convenience.diag type
#' @param minimumESS The threshold for the ESS, default value is 625
#' @param filename A filename to save the table, if NULL the table will be printed
#' @param ... (various) Additional arguments passed to plot()
#'  
#' @export

plotDiffSplits <- function(output, minimumESS = 625, fill_color = NULL, filename = NULL, ...){
  
  if( is.null(fill_color) ){
    fill_color <- "coral3"
  }
  
  if( !(is.null(filename)) ){
    pdf(file = filename, width = 4.5, height = 4.5)
  }
  
  if( minimumESS == 625) exp_diff_runs <- convenience::exp_diff_runs 
  else  exp_diff_runs <- expectedDiffSplits(minimumESS)
  
  frequencies <- vector()
  differences <- vector()
  for (i in 1:length(output$tree_parameters$frequencies)) {
    frequencies <- c(frequencies, as.vector(unlist(output$tree_parameters$frequencies[[i]])))
    differences <- c(differences, as.vector(unlist(output$tree_parameters$compare_runs[[i]])))
  }
  frequencies <- frequencies[!is.na(frequencies)]
  differences <- differences[!is.na(differences)]
  
  x_axis <- exp_diff_runs[1,]
  y_axis <- exp_diff_runs[2,]
  
  y_lim <- max(differences, y_axis, na.rm = T)
  y_lim <- y_lim + 0.5*y_lim
  
  par(mar = c(3.9, 4.9, 2.1, 1.3))
  plot <- plot(NA, 
               xlab = NA,
               ylab = NA,
               main = "Difference in Split Frequencies", 
               xlim = c(0.0,1.0), 
               ylim = c(0.0, y_lim),
               las = 1)
  
  plot <- lines(x_axis, y_axis, col = "antiquewhite4", lwd=3)
  title(xlab = "Split frequency", outer = T, line = -1.1)
  title(ylab = "Difference between splits", outer = T, line = -1.1)
  
  plot <- points(frequencies, differences, pch = 16, col = fill_color)
 
  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}