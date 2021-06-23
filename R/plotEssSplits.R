#' Histogram of the ESS values
#' 
#' Plots the histogram of the ESS values for the splits
#' 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline hist layout legend lines par points polygon rect title
#' 
#' @param x A list of convenience.diag type
#' @param per_run If the plot should combine all runs or plot each run separately. Default value is FALSE
#' @param breaks The number of bins to determine the intervals of the histogram
#' @param precision The precision of the mean estimates. Default is 0.01
#' @param fill_color The color to fill the histogram bars
#' @param filename The name of the file to save the plot
#' @param ... (various) Additional arguments passed to plot().
#' 
#' @return Histogram
#' 
#' @export

plotEssSplits <- function(x, per_run = FALSE, breaks = NULL, precision = 0.01, fill_color = NULL, filename = NULL, ...){
  
  col_threshold <- "gray69"
  
  if( is.null(fill_color) ){
    fill_color <- "seagreen4"
  }
  
  if( !(is.null(filename)) ){
    pdf(file = filename, width = 4.5, height = 4.5)
  }
  
  minimumESS <- minESS(precision)
  ESS_values <- vector()
  
  if(per_run == TRUE){
    n_runs <- ncol(x$tree_parameters$ess)
    
    par(mfrow = c(n_runs/2,n_runs/2), oma = c(1.5,1.5,1.5,0) + 0.1, mar = c(2,2,2.3,2))
    layout(matrix(c(1:n_runs), nrow=n_runs/2, ncol = n_runs/2, byrow=TRUE))
    
    for (i in 1:n_runs) {
      ESS_values <- x$tree_parameters$ess[,i]
      ESS_values <- ESS_values[!is.na(ESS_values)]
      if(is.null(breaks)){
        breaks <- seq(0, (max(minimumESS, ESS_values))+50, 25)
      }
      
      y_topLim <- max(hist(ESS_values, plot = FALSE)$counts)
      x_topLim <- max(minimumESS,ESS_values) + (max(minimumESS, ESS_values))/10
      
      plot <- plot(NA,
                   xlab = NA, 
                   ylab = NA,
                   main = colnames(x$tree_parameters$ess)[i],
                   xlim = c(0, x_topLim ),
                   ylim = c(0, y_topLim+1),
                   las=1,
                   bty="l")
      plot <- rect(xleft = minimumESS, ybottom = 0, xright = x_topLim, ytop = y_topLim+1, border = NA, col = "gray89")
      plot <- lines(x = c(minimumESS,minimumESS),y=c(0,y_topLim+1), col =  col_threshold, lwd= 2, lty=2)
      
      
      
      
      plot <- hist(ESS_values, 
                   xlab = NA, 
                   ylab = NA,
                   xlim = c(0, x_topLim ),
                   ylim = c(0, y_topLim + 1),
                   breaks = breaks,
                   border = F,
                   col = fill_color,
                   add=T,
                   ...)
    }
    title(main = "Effective Sample Size for splits per run", cex.main = 0.9, xlab = "Effective Sample Size", ylab = "Counts", outer = TRUE, line = 0.5)
    
    
  } else{
    
    for (i in 1:length(x$tree_parameters$ess)) {
      ESS_values <- c(ESS_values, x$tree_parameters$ess[[i]])
    }
    ESS_values <- ESS_values[!is.na(ESS_values)]
    
    if(is.null(breaks)){
      breaks <- seq(0, (max(minimumESS, ESS_values))+50, 25)
    }
    
    x_topLim <- max(minimumESS,ESS_values) + (max(minimumESS, ESS_values))/10
    y_topLim <- max(hist(ESS_values, breaks = breaks, plot = FALSE)$counts)
    
    par(mar = c(4.1, 3.9, 2.1, 1.0))
    
    plot <- plot(NA,
                 xlab = "Effective Sample Size", 
                 ylab = "Counts",
                 main = "Effective Sample Size for splits",
                 cex.main = 0.9,
                 xlim = c(0, x_topLim ),
                 ylim = c(0, y_topLim+1),
                 las=1,
                 bty="l")
    plot <- rect(xleft = minimumESS, ybottom = 0, xright = x_topLim, ytop = y_topLim+1, border = NA, col = "gray89")
    plot <- lines(x = c(minimumESS,minimumESS),y=c(0,y_topLim+1), col = col_threshold, lwd= 2, lty=2)
    
    
    
    
    plot <- hist(ESS_values,
                 xlim = c(0, x_topLim ),
                 ylim = c(0, y_topLim + 1),
                 breaks = breaks,
                 border = F,
                 col = fill_color,
                 add=T,
                 ...)
  }

  if( !(is.null(filename)) ){
    dev.off()
  } 
  
}