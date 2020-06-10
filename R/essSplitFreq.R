
essSplitFreq <- function(runs){
  
  listESS <- list()
  names_runs <- vector()
  for (z in 1:length(runs)) {
    x <- getInfo(runs, z, trees = TRUE)
    names_runs <- c(names_runs, paste("Run", z))
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
      
      ess_splits_coda <- c(ess_splits_coda, effectiveSize(is.split))
    }
    listESS[[z]] <- ess_splits_coda
    names(listESS[[z]]) <- vecNames
  }
  
  names(listESS) <- names_runs
  return(listESS)
}


if(F){
# Load example
runs <- loadFiles("example/4_runs/")

# Get the trees for 1 run
x <- getInfo(runs, 1, trees = TRUE)

clades.list <- list()

for (i in 1:length(x)) {
  clades.list[[i]] <- clade.freq.tree(x[[i]])
}

clades <-  clade.freq.trees(x, start = start, end = end)
total_n_splits <- length(clades$cladenames)


# Let's calculate ESS
ess_splits_coda <- vector()
ess_splits_mcmcse <- vector()
ess_splits_tracer <- vector()

for (j in 1:total_n_splits) {
  is.split <- vector()
  
  for (i in 1:length(clades.list)) {
    if (clades$cladenames[j] %in% clades.list[[i]]$cladenames){
      is.split <- c(is.split, 1)
    }
    else{
      is.split <- c(is.split, 0)
    }
  }
  splits <- as.mcmc(is.split)
  traceplot(splits)
  title(main = paste(effectiveSize(is.split),ess(is.split), essTracer(is.split)))
  
  ess_splits_coda <- c(ess_splits_coda, effectiveSize(is.split))
  ess_splits_mcmcse <- c(ess_splits_mcmcse, ess(is.split))
  ess_splits_tracer <- c(ess_splits_tracer, essTracer(is.split))
  
}
  
}