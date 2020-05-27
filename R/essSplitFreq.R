runs <- loadFiles("example/4_runs/")

x <- getInfo(runs, 1, trees = TRUE)

start = 1
end = length(x)
if(class(x) == "rwty.chain"){
  x <- x$trees
}

if (length(x) == 1 && class(x[[1]]) == "multiPhylo"){
  x <- x[[1]]
}

x <- x[start:end]

clades.list <- list()

for (i in 1:length(x)) {
  clades.list[[i]] <- prop.part(x[[i]])
}

clades <-  prop.part(x)
total_n_splits <- length(clades)

ess_splits_coda <- vector()
ess_splits_mcmcse <- vector()

for (j in 1:total_n_splits) {
  is.split <- vector()
  
  for (i in 1:length(clades.list)) {
    if (clades[j] %in% clades.list[[i]]){
      is.split <- c(is.split, 1)
    }
    else{
      is.split <- c(is.split, 0)
    }
  }
  splits <- as.mcmc(is.split)
  traceplot(splits)
  title(main = paste(effectiveSize(is.split),ess(is.split)))
  
  ess_splits_coda <- c(ess_splits_coda, effectiveSize(is.split))
  ess_splits_mcmcse <- c(ess_splits_mcmcse, ess(is.split))
  
}
