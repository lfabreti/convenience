#' ESS calculation for the splits
#' 
#' Calculates the Effective Sample Size for the splits in the trees
#' 
#' @param runs A list of rwty.chain type
#' @param tracer A boolean to determine if ESS should be calculated with Tracer method. If set to FALSE, ESS will be calculated with CODA
#' 
#' @return A list with the ESS for each split
#' 
#' @export

essSplitFreq <- function(runs, tracer){
  
  listESS <- list()
  names_runs <- vector()
  for (z in 1:length(runs)) {
    x <- getInfo(runs, z, trees = TRUE)
    names_runs <- c(names_runs, paste("Run_", z, sep = ""))
    clades.list <- list()
    
    for (j in 1:length(x)) {
      clades.list[[j]] <- clade.freq.tree(x[[j]])
    }
    
    clades <-  clade.freq.trees(x, start = 1, end = length(x))
    total_n_splits <- length(clades$cladenames_post)
    
    ess_splits <- vector()
    
    if( length(clades[[1]]) > 0 ){
      for (j in 1:total_n_splits) {
        is.split <- vector()
        vecNames <- vector()
        
        for (i in 1:length(clades.list)) {
          if (clades$cladenames[j] %in% clades.list[[i]]$cladenames){
            is.split <- c(is.split, 1)
          }
          else{
            is.split <- c(is.split, 0)
          }
        }
        
        vecNames <- c( vecNames, as.character(clades$cladenames_post))
        
        if( tracer == T){
          ess_splits <- c(ess_splits, essTracer(is.split))
        }else{
          ess_splits <- c(ess_splits, effectiveSize(is.split))
        }
      }
    }
    
    listESS[[z]] <- ess_splits
    if( length(ess_splits) > 0 ) names(listESS[[z]]) <- vecNames
  }
  
  names(listESS) <- names_runs
  return(listESS)
}