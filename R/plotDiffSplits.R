#' Plots the difference between splits of different runs
#' 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline hist layout legend lines par points polygon rect title
#' 
#' @param output A list of convenience.diag type
#' @param minimumESS The threshold for the ESS, default value is 625
#' @param fill_color The color for the dots on the plot
#' @param filename A filename to save the table, if NULL the table will be printed
#' @param ... (various) Additional arguments passed to plot()
#'  
#' @export

plotDiffSplits <- function(output, minimumESS = 625, fill_color = NULL, filename = NULL, ...){
  
  col_threshold <- "gray69"
  
  if( is.null(fill_color) ){
    fill_color <- "seagreen4"
  }
  
  if( !(is.null(filename)) ){
    pdf(file = filename, width = 4.5, height = 4.5)
  }
  
  if( minimumESS == 625) {
    fdir <- system.file("thresholds/expectedDiff_625.rds", package = "convenience")
    exp_diff_runs <- readRDS(fdir)
  } 
  else  exp_diff_runs <- expectedDiffSplits(minimumESS)
  
  ## Calculate the minimum ESS between runs for each split
  ess_min_between_runs <- matrix(ncol=length(output$tree_parameters$frequencies), nrow=nrow(output$tree_parameters$ess))
  n_runs <- length(output$tree_parameters$ess)
  count <-1
  for (i in 1:(n_runs-1)) {
    for (j in (i+1):n_runs) {
      for (z in 1:nrow(output$tree_parameters$ess)) {
        ess_min_between_runs[z,count] <- min(output$tree_parameters$ess[z,i], output$tree_parameters$ess[z,j], na.rm = T)
      }
      count <- count + 1
    }
  }
  row.names(ess_min_between_runs) <- row.names(output$tree_parameters$ess)
  
  frequencies <- vector()
  differences <- vector()
  freq_low_ess <- vector()
  diff_low_ess <- vector()
  for (i in 1:length(output$tree_parameters$frequencies)) {
    for (j in names(output$tree_parameters$frequencies[[i]])) {
      if ( j %in% row.names(ess_min_between_runs)){
        if ( ess_min_between_runs[j, i] < minimumESS ){
          freq_low_ess <- c( freq_low_ess, as.vector(unlist(output$tree_parameters$frequencies[[i]][j])) )
          diff_low_ess <- c( diff_low_ess, as.vector(unlist(output$tree_parameters$compare_runs[[i]][j])) )
        }
        else{
          frequencies <- c( frequencies, as.vector(unlist(output$tree_parameters$frequencies[[i]][j])) )
          differences <- c( differences, as.vector(unlist(output$tree_parameters$compare_runs[[i]][j])) )
        }
      }
      else{
        frequencies <- c( frequencies, as.vector(unlist(output$tree_parameters$frequencies[[i]][j])) )
        differences <- c( differences, as.vector(unlist(output$tree_parameters$compare_runs[[i]][j])) )
      }
    }
  }
  frequencies <- frequencies[!is.na(frequencies)]
  differences <- differences[!is.na(differences)]
  
  x_axis <- exp_diff_runs[1,]
  y_axis <- exp_diff_runs[2,]
  
  y_lim <- max(differences, y_axis, diff_low_ess, na.rm = T)
  y_lim <- y_lim + y_lim*0.5
  
  par(mar = c(4.1, 4.9, 2.1, 1.3))
  plot <- plot(NA, 
               xlab = NA,
               ylab = NA,
               main = "Difference in Split Frequencies", 
               cex.main = 0.9,
               xlim = c(0.0,1.0), 
               ylim = c(0.0, y_lim),
               las = 1)
  
  plot <- polygon(x_axis, y_axis, col = "gray89", border = NA)
  x_extra <- c(0.0, 0.01, 0.99, 1.0, 0.0)
  y_extra <- c(0.0, min(y_axis), min(y_axis), 0.0, 0.0)
  plot <- polygon(x_extra, y_extra, border = NA, col = "gray89")
  plot <- lines(x_axis, y_axis, col = col_threshold, lwd=2)
  title(xlab = "Split frequency", outer = T, line = -1.1)
  title(ylab = "Difference between splits", outer = T, line = -1.1)
  
  plot <- points(frequencies, differences, pch = 16, col = fill_color)
  if (length(freq_low_ess) > 0){
    plot <- points(freq_low_ess, diff_low_ess, pch = 16, col = "tan1")
    legend("topright",
           legend = c(paste("ESS",expression("<"),minimumESS),paste("ESS",expression(">"),minimumESS)),
           pch = c(16,16), 
           col = c("tan1",fill_color),
           box.lty = 2,
           box.col="gray7",
           cex = 0.7,
           inset = 0.01)
  }
 
  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}