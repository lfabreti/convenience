# For the naive user

checkConvergence <- function(path, control){
  
  ##### First let's check the function arguments #####
  if ( is.null(control$precision) ){
    precision <- 0.01
  } else {
    precision <- control$precision
  }
  
  if ( is.null(control$burnin) ){
    burnin <- 0.1
  } else {
    burnin <- control$burnin
  }
  
  if ( is.null(control$namesToExclude) ){
    namesToExclude <- "br_lens|bl|Iteration|Likelihood|Posterior|Prior"
  } else {
    namesToExclude <- control$namesToExclude
  }
  #####
  
  final_output <- list()
  
  # load the mcmc output
  my_runs <- loadFiles(path, burnin, format = "revbayes")
  
  # calculating minimum ESS according to precision
  minimumESS <- minESS(precision)
  
  compar_names <- vector()
  for (r1 in 1:(length(my_runs)-1)) {
    for (r2 in (r1+1):length(my_runs)) {
      compar_names <- c(compar_names, paste("Run", r1, "_Run", r2, sep = ""))
    }
  }
  
  # list of results
  output_tree_parameters <- list()
  output_continuous_parameters <- list()
  
  #################################################
  ## IF WE HAVE TREE FILES, WE CHECK THE SPLITS  ##
  #################################################
  
  if( length(my_runs[[1]]$trees) > 0 ){ # if we have tree files, we check split freq and ESS of splits
    
    print("Analyzing tree parameters")
    
    ## ESS ##
    
    ess_runs_splits <- essSplitFreq(my_runs)
    
    for (i in 1:length(ess_runs_splits)) {
      ess_runs_splits[[i]] <- ess_runs_splits[[i]] - minimumESS
    }
    
    output_tree_parameters$ess <- ess_runs_splits
    
    ## Compare windows ##
    
    splits_windows <- splitFreq(my_runs, windows = T)
    
    exp_diff_windows <- expectedDiffSplits(100)
    
    for (i in 1:ncol(splits_windows)) {
      splits_windows[2,i]$listFrequencies <- round(splits_windows[2,i]$listFrequencies, digits = 2)
    }
    
    results_splits <- list()
    for (i in 1:ncol(splits_windows)) {
      results <- vector()
      # index gets the indexes sorted
      index <- which( exp_diff_windows[1,] %in% splits_windows[2,i]$listFrequencies )
      vecNames <- vector()
      
      for (j in index) {
        for (z in 1:length(splits_windows[2,i]$listFrequencies)) {
          if( exp_diff_windows[1,j] == splits_windows[2,i]$listFrequencies[z]){
            results <- c(results, exp_diff_windows[2,j] - splits_windows[1,i]$listSplits[z])
            vecNames <- c( vecNames, names(splits_windows[1,i]$listSplits[z]))
          }
        }
      }
      results_splits[[i]] <- results
      names(results_splits[[i]]) <- vecNames
    }
    for (i in 1:length(results_splits)) {
      names(results_splits)[i] <- paste("Run_", i, sep = "")
    }
    
    output_tree_parameters$compare_windows <- results_splits
    
    
    ## Compare runs ##
    
    splits_runs <- splitFreq(my_runs, windows = F)
    
    exp_diff_runs <- expectedDiffSplits(minimumESS)
    
    for (i in 1:ncol(splits_runs)) {
      splits_runs[2,i]$listFrequencies <- round(splits_runs[2,i]$listFrequencies, digits = 2)
    }
    
    results_splits_runs <- list()
    for (i in 1:ncol(splits_runs)) {
      results <- vector()
      index <- which( exp_diff_runs[1,] %in% splits_runs[2,i]$listFrequencies )
      vecNames <- vector()
      
      for (j in index) {
        for (z in 1:length(splits_runs[2,i]$listFrequencies)) {
          if( exp_diff_runs[1,j] == splits_runs[2,i]$listFrequencies[z]){
            results <- c(results, exp_diff_runs[2,j] - splits_runs[1,i]$listSplits[z])
            vecNames <- c( vecNames, names(splits_runs[1,i]$listSplits[z]))
          }
        }
      }
      results_splits_runs[[i]] <- results
      names(results_splits_runs[[i]]) <- vecNames
    }
    
    names(results_splits_runs) <- compar_names
    output_tree_parameters$compare_runs <- results_splits_runs
    
  }
  
  
  
  ######################################################### 
  ## IF WE HAVE LOG FILES, WE CHECK THE CONTINUOUS PARAM ##
  #########################################################
  
  if( length(my_runs[[1]]$ptable) > 0 ){ 
    
    print("Analyzing continuous parameters")
    
    ## ESS ##
    
    ess_runs_cont_param <- essContParam(my_runs)
    ess_runs_cont_param <- ess_runs_cont_param - minimumESS
    
    output_continuous_parameters$ess <- ess_runs_cont_param
    
    ## Compare windows ##
    
    ## KS score ##
    ks_windows <- ksTest(my_runs, windows = TRUE)
    ess_windows <- essContParam(my_runs, windows = T)
    
    ks_limits_windows <-vector()
    for (i in 1:length(my_runs)) {
      ess_runs <- ess_windows[grep(paste("Run_",i, sep = ""), rownames(ess_windows)),]
      for (j in 1:ncol(ess_runs)) {
        ks_limits_windows <- c(ks_limits_windows, ksThreshold(1.95, ess_runs[1,j], ess_runs[2,j]))
      }
    }
    ks_limits_windows <- data.frame(matrix(unlist(ks_limits_windows), nrow = length(my_runs), byrow = T), stringsAsFactors = F)
    colnames(ks_limits_windows) <- names(ks_windows)
    for (i in 1:length(my_runs)) {
      rownames(ks_limits_windows)[i] <- paste("Run_", i, sep = "")
    }
    
    ks_windows <- ks_limits_windows - ks_windows
    
    output_continuous_parameters$compare_windows <- ks_windows
    
    ## Compare runs ##
    
    ks_runs <- ksTest(my_runs, windows = F)
    ess_runs <- essContParam(my_runs, windows = F)
    
    ks_limits <-vector()
    count<-0
    for (df1 in 1:(length(my_runs)-1)){
      for (df2 in (df1+1):length(my_runs)) {
        for (i in 1:nrow(ess_runs)) {
          ks_limits <- c(ks_limits, ksThreshold(1.95, ess_runs[i,df1], ess_runs[i,df2]))
        }
        count <- count+1
      }
    }
    ks_limits <- data.frame(matrix(unlist(ks_limits), nrow = count, byrow = T), stringsAsFactors = F)
    
    colnames(ks_limits) <- names(ks_runs)
    
    ks_runs <- ks_limits - ks_runs
    
    rownames(ks_runs) <- compar_names
    
    output_continuous_parameters$compare_runs <- ks_runs

  }
  
  
  ##################### 
  ## DECISION MAKING ##
  #####################
  
  decision_list_trees <- list()
  decision_list_cont <- list()
  fails <- list()
  count_decision <- 0
  
  
  ##### ONLY TREE FILES #####
  
  # Check ess_runs_splits, results_splits, results_splits_runs
  if( length(my_runs[[1]]$trees) > 0 ){
    
    ess_splits <- list()
    split_freq_windows <- list()
    split_freq_runs <- list()
    
    ## ESS ##
    for (i in 1:length(output_tree_parameters[[1]])) {
      ess_splits[[i]] <- which(output_tree_parameters[[1]][[i]] < 0)
    }
    names(ess_splits) <- names(output_tree_parameters[[1]])
    
    for (i in 1:length(output_tree_parameters[[2]])) {
      split_freq_windows[[i]] <- which(output_tree_parameters[[2]][[i]] < 0)
    }
    names(split_freq_windows) <- names(output_tree_parameters[[2]])
    
    for (i in 1:length(output_tree_parameters[[3]])) {
      split_freq_runs[[i]] <- which(output_tree_parameters[[3]][[i]] < 0)
    }
    names(split_freq_runs) <- names(output_tree_parameters[[3]])
    
    decision_list_trees[[1]] <- ess_splits
    decision_list_trees[[2]] <- split_freq_windows
    decision_list_trees[[3]] <- split_freq_runs
    
    for (i in 1:length(decision_list_trees)) {
      for (j in 1:length(decision_list_trees[[i]])) {
        if( length(decision_list_trees[[i]][[j]]) > 0){
          count_decision <- count_decision + 1
          fails <- c(fails,decision_list_trees[[i]][j])
        }
      }
    }
    
  }
  
  ##### ONLY LOG FILES #####
  
  # Check ess_runs_cont_param, ks_windows, ks_runs
  if( length(my_runs[[1]]$ptable) > 0 ){
    
    decision_list_cont[[1]] <- which(output_continuous_parameters[[1]] < 0)
    
    decision_list_cont[[2]] <- which(output_continuous_parameters[[2]] < 0)
    
    decision_list_cont[[3]] <- which(output_continuous_parameters[[3]] < 0)
    
    for (i in 1:length(decision_list_cont)) {
      if( length(decision_list_cont[[i]]) > 0 ){
        count_decision <- count_decision + 1
        fails <- c(fails, decision_list_cont[[i]])
      }
    }
    
  }
  if(count_decision == 0){
    final_output$message <- "Achieved convergence"
  }else{
    final_output$message <- "Failed convergence"
    final_output$failed <- fails
  }
  
  class(output_continuous_parameters) <- "convenience.diag"
  class(output_tree_parameters) <- "convenience.diag"
  
  final_output$continuous_parameters <- output_continuous_parameters
  final_output$tree_parameters <- output_tree_parameters
  
  
  class(final_output) <- "convenience.diag"
  
  final_output
}

