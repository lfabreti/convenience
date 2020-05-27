# Function to calculate ESS according to Tracer
essTracer <- function(input,stepSize = 1){
  
  samples <- length(input)
  
  max_lag <- 2000
  max_lag <- min( (samples-1) , max_lag)
  gammaSt <- replicate(max_lag, 0)
  varStat <- 0
  
  i = 1 
  
  while (i <= max_lag) {
    
    for (j in 1:(samples-i+1)) {
      del1 <- input[j] - mean(input)
      del2 <- input[j+i-1] -mean(input)
      gammaSt[i] <- gammaSt[i] + (del1*del2)
    }
    
    gammaSt[i] <- gammaSt[i] / (samples - i + 1)
    
    if(i == 1){
      varStat <- gammaSt[1]
      
    } else if( i %% 2 == 1){
      
      if( (gammaSt[i-1]+ gammaSt[i]) > 0 ){
        varStat <- varStat + 2 * (gammaSt[i-1]+gammaSt[i])
      }
      else{
        max_lag <- i
      }
      
    }
    i = i+1
  }
  
  act <- 0
  ess <- 0
  
  if (gammaSt[1] == 0){
    act <-0
  } else {
    act <- stepSize * varStat / gammaSt[1]
  }
  
  if (act ==0){
    ess <-1
  } else {
    ess <- (stepSize*samples) / act
  }
  
  return(ess)
}

# Load example
runs <- loadFiles("example/4_runs/")

# Get the trees for 1 run
x <- getInfo(runs, 1, trees = TRUE)

# Ape functions to get info from trees
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


# Let's calculate ESS
ess_splits_coda <- vector()
ess_splits_mcmcse <- vector()
ess_splits_tracer <- vector()

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
  title(main = paste(effectiveSize(is.split),ess(is.split), essTracer(is.split)))
  
  ess_splits_coda <- c(ess_splits_coda, effectiveSize(is.split))
  ess_splits_mcmcse <- c(ess_splits_mcmcse, ess(is.split))
  ess_splits_tracer <- c(ess_splits_tracer, essTracer(is.split))
  
}
