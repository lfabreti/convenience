#' Remove burn-in
#'
#' Remove the initial statates from the MCMC output
#' 
#' @param output A list of rwty.chain type
#' @param burnin Fraction of generations to discard
#' 
#' @return List of type rwty.chain
#' 
#' @export


removeBurnin <- function(output, burnin) {
  
  #burnin
  if ( burnin >= length(output[[1]]$trees) | burnin >= length(output[[1]]$ptable) ) stop("Burnin larger than iterations in file")
  
  for (i in 1:length(output)) {
    
    if (burnin >= 1) {
      
      if( length(output[[i]]$trees) > 0 ) output[[i]]$trees <- output[[i]]$trees[(burnin+1):(length(output[[i]]$trees))]
      if( length(output[[i]]$ptable) > 0 ) output[[i]]$ptable <- output[[i]]$ptable[(burnin+1):(nrow(output[[i]]$ptable)),]
      
    } else if (burnin < 1 & burnin > 0) { 
      
      if( length(output[[i]]$trees) > 0 ) discard <- ceiling(burnin*(length(output[[i]]$trees)))
      else discard <- ceiling(burnin*(nrow(output[[i]]$ptable)))
      
      if( length(output[[i]]$trees) > 0 ) output[[i]]$trees <- output[[i]]$trees[(discard+1):(length(output[[i]]$trees))]
      if( length(output[[i]]$ptable) > 0 ) output[[i]]$ptable <- output[[i]]$ptable[(discard+1):(nrow(output[[i]]$ptable)),]
      
    } else if (burnin == 0) {
      
      output[[i]] <- output[[i]]
      
    } else {
      
      stop("What have you done?")
    }
    
  }
  
  return(output)
}