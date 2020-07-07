#' ESS for splits
#' 
#' Calculates the Effective Sample Size for the splits in the trees
#' 
#' @param runs The runs from loadFiles()
#' 
#' @return A list with the ESS for each split
#' 
#' @export

essSplitFreq <- function(runs){
  
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
    
    ess_splits_coda <- vector()
    
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
      
      ess_splits_coda <- c(ess_splits_coda, essTracer(is.split))
    }
    listESS[[z]] <- ess_splits_coda
    names(listESS[[z]]) <- vecNames
  }
  
  names(listESS) <- names_runs
  return(listESS)
}