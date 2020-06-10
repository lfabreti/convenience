runs <- loadFiles("example/4_runs/", burnin = 0.1)
  

#####   First let's check if we are sampling from the stationary distribution   #####


## KS score ##
ks_windows <- ksTest(runs, windows = TRUE)
ess_windows <- essContParam(runs, windows = T)

ks_limits_windows <-vector()
for (i in 1:length(runs)) {
  ess_runs <- ess_windows[grep(paste("Run",i), rownames(ess_windows)),]
  for (j in 1:ncol(ess_runs)) {
    ks_limits <- c(ks_limits, ksThreshold(1.95, ess_runs[1,j], ess_runs[2,j]))
  }
}
ks_limits_windows <- data.frame(matrix(unlist(ks_limits_windows), nrow = length(runs), byrow = T), stringsAsFactors = F)
colnames(ks_limits_windows) <- names(ks_windows)
for (i in 1:length(runs)) {
  rownames(ks_limits_windows)[i] <- paste("Run", i)
}

ks_limits_windows - ks_windows

## Deviation of split frequencies ##
splits_windows <- splitFreq(runs, windows = T)

exp_diff <- expectedDiffSplits(100)
  
for (i in 1:ncol(splits_windows)) {
  splits_windows[2,i]$listFrequencies <- round(splits_windows[2,i]$listFrequencies, digits = 4)
}

results_splits <- list()
for (i in 1:ncol(splits_windows)) {
  results <- vector()
  # index gets the indexes sorted
  index <- which( exp_diff[1,] %in% splits_windows[2,i]$listFrequencies )
  vecNames <- vector()
  
  for (j in index) {
    for (z in 1:length(splits_windows[2,i]$listFrequencies)) {
      if( exp_diff[1,j] == splits_windows[2,i]$listFrequencies[z]){
        results <- c(results, exp_diff[2,j] - splits_windows[1,i]$listSplits[z])
        vecNames <- c( vecNames, names(splits_windows[1,i]$listSplits[z]))
      }
    }
  }
  results_splits[[i]] <- results
  names(results_splits[[i]]) <- vecNames
}
for (i in 1:length(results_splits)) {
  names(results_splits)[i] <- paste("Run", i)
}


#####   Do we have enough samples?   #####

## ESS for continuous parameters ##
ess_runs_cont_param <- essContParam(runs)

ess_runs_cont_param - 625

## ESS for the splits ##
ess_runs_splits <- essSplitFreq(runs)



#####   Reproducible runs?   #####

## KS score ##
ks_runs <- ksTest(runs, windows = F)
ess_runs <- essContParam(runs, windows = F)

ks_limits <-vector()
count<-0
for (df1 in 1:(length(runs)-1)){
  for (df2 in (df1+1):length(runs)) {
    for (i in 1:nrow(ess_runs)) {
      ks_limits <- c(ks_limits, ksThreshold(1.95, ess_runs[i,df1], ess_runs[i,df2]))
    }
    count <- count+1
  }
}
ks_limits <- data.frame(matrix(unlist(ks_limits), nrow = count, byrow = T), stringsAsFactors = F)

colnames(ks_limits) <- names(ks_runs)

ks_limits - ks_runs


## Deviation of split frequencies ##
splits_runs <- splitFreq(runs, windows = F)

exp_diff <- expectedDiffSplits(625)

for (i in 1:ncol(splits_runs)) {
  splits_runs[2,i]$listFrequencies <- round(splits_runs[2,i]$listFrequencies, digits = 4)
}

results_splits <- list()
for (i in 1:ncol(splits_runs)) {
  results <- vector()
  index <- which( exp_diff[1,] %in% splits_runs[2,i]$listFrequencies )
  vecNames <- vector()
  
  for (j in index) {
    for (z in 1:length(splits_runs[2,i]$listFrequencies)) {
      if( exp_diff[1,j] == splits_runs[2,i]$listFrequencies[z]){
        results <- c(results, exp_diff[2,j] - splits_runs[1,i]$listSplits[z])
        vecNames <- c( vecNames, names(splits_runs[1,i]$listSplits[z]))
      }
    }
  }
  results_splits[[i]] <- results
  names(results_splits[[i]]) <- vecNames
}
