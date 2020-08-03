#' Check Convergence
#' 
#' Check output from RevBayes for convergence diagnostics
#' 
#' @param path Path to directory containing all files from the same analysis
#' @param list_files List of files to check for convergence
#' @param control List of arguments to the function. Includes burn-in, precision and names of parameters to exclude from the analysis
#' 
#' @return List of type convenience.diag
#' 
#' @export


checkConvergence <- function(path = NULL, list_files = NULL, control = makeControl()){
  
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
  
  if ( !is.null(path)){
    # load the mcmc output
    my_runs <- loadFiles(path, burnin=burnin, format = "revbayes")
  }else {
    my_runs <- loadFiles(list_files = list_files, burnin = burnin, format = "revbayes")
  }

  
  # calculating minimum ESS according to precision
  minimumESS <- minESS(precision)
  minimumESS_windows <- round(minimumESS/5)
  
  if( length(my_runs) > 1 ){
    compar_names <- vector()
    for (r1 in 1:(length(my_runs)-1)) {
      for (r2 in (r1+1):length(my_runs)) {
        compar_names <- c(compar_names, paste("Run", r1, "_Run_", r2, sep = ""))
      }
    }
  }
  
  # list of results - thresholds
  output_tree_parameters <- list()
  output_continuous_parameters <- list()
  
  # list of results
  output_tree_parameters_raw <- list()
  output_continuous_parameters_raw <- list()
  
  #################################################
  ## IF WE HAVE TREE FILES, WE CHECK THE SPLITS  ##
  #################################################
  
  if( length(my_runs[[1]]$trees) > 0 ){ # if we have tree files, we check split freq and ESS of splits
    
    print("Analyzing tree parameters")
    
    ## ESS ##
    
    ess_runs_splits <- essSplitFreq(my_runs)
    output_tree_parameters_raw$ess <- ess_runs_splits
    
    for (i in 1:length(ess_runs_splits)) {
      ess_runs_splits[[i]] <- ess_runs_splits[[i]] - minimumESS
    }
    
    output_tree_parameters$ess <- ess_runs_splits
    
    ## Compare windows ##
    
    splits_windows <- splitFreq(my_runs, windows = T)
    
    exp_diff_windows <- expectedDiffSplits(minimumESS_windows)
    
    for (i in 1:ncol(splits_windows)) {
      splits_windows[2,i]$listFrequencies <- round(splits_windows[2,i]$listFrequencies, digits = 2)
    }
    
    col_names <- list()
    for (i in 1:ncol(splits_windows)) {
      col_names <- c(col_names, paste("Run_", i, sep = ""))
    }
    colnames(splits_windows) <- col_names
    output_tree_parameters_raw$compare_windows <- splits_windows
    
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
    
    if( length(my_runs) > 1 ){
      splits_runs <- splitFreq(my_runs, windows = F)
      
      exp_diff_runs <- expectedDiffSplits(minimumESS)
      
      for (i in 1:ncol(splits_runs)) {
        splits_runs[2,i]$listFrequencies <- round(splits_runs[2,i]$listFrequencies, digits = 2)
      }
      
      colnames(splits_runs) <- compar_names
      output_tree_parameters_raw$compare_runs <- splits_runs
      
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
    
  }
  
  
  
  ######################################################### 
  ## IF WE HAVE LOG FILES, WE CHECK THE CONTINUOUS PARAM ##
  #########################################################
  
  if( length(my_runs[[1]]$ptable) > 0 ){ 
    
    print("Analyzing continuous parameters")
    
    ## ESS ##
    
    ess_runs_cont_param <- essContParam(my_runs, namesToExclude = namesToExclude)
    output_continuous_parameters_raw$ess <- ess_runs_cont_param
    
    ess_runs_cont_param <- ess_runs_cont_param - minimumESS
    output_continuous_parameters$ess <- ess_runs_cont_param
    
    ## Compare windows ##
    
    ## KS score ##
    ks_windows <- ksTest(my_runs, windows = TRUE, namesToExclude = namesToExclude)
    
    ks_limits_windows <-vector()
    for (i in 1:length(my_runs)) {
      if (typeof(ks_windows) == "double"){
        ks_limits_windows <- c(ks_limits_windows, ksThreshold(1.95, minimumESS_windows, minimumESS_windows))
      }else{
        for (j in 1:ncol(ks_windows)) {
          ks_limits_windows <- c(ks_limits_windows, ksThreshold(1.95, minimumESS_windows, minimumESS_windows))
        }
      }
    }
    
    ks_limits_windows <- data.frame(matrix(unlist(ks_limits_windows), nrow = length(my_runs), byrow = T), stringsAsFactors = F)
    colnames(ks_limits_windows) <- names(ks_windows)
    for (i in 1:length(my_runs)) {
      rownames(ks_limits_windows)[i] <- paste("Run_", i, sep = "")
      rownames(ks_windows)[i] <- paste("Run_", i, sep = "")
    }
    
    output_continuous_parameters_raw$compare_windows <- ks_windows
    ks_windows <- ks_limits_windows - ks_windows
    
    output_continuous_parameters$compare_windows <- ks_windows
    
    ## Compare runs ##
    
    if( length(my_runs) > 1 ){
      ks_runs <- ksTest(my_runs, windows = F, namesToExclude = namesToExclude)
      
      ks_limits <-vector()
      count<-0
      for (df1 in 1:(length(my_runs)-1)){
        for (df2 in (df1+1):length(my_runs)) {
          for (i in 1:ncol(ks_runs)) {
            ks_limits <- c(ks_limits, ksThreshold(1.95,minimumESS, minimumESS))
          }
          count <- count+1
        }
      }
      ks_limits <- data.frame(matrix(unlist(ks_limits), nrow = count, byrow = T), stringsAsFactors = F)
      
      colnames(ks_limits) <- names(ks_runs)
      rownames(ks_runs) <- compar_names
      
      output_continuous_parameters_raw$compare_runs <- ks_runs
      ks_runs <- ks_limits - ks_runs
      
      rownames(ks_runs) <- compar_names
      
      output_continuous_parameters$compare_runs <- ks_runs
    }

  }
  
  
  ##################### 
  ## DECISION MAKING ##
  #####################
  
  decision_list_trees <- list()
  decision_list_cont <- list()
  fails <- list()
  fails_names <- list()
  count_decision <- 0
  
  
  ##### ONLY TREE FILES #####
  
  # Check ess_runs_splits, results_splits, results_splits_runs
  if( length(my_runs[[1]]$trees) > 0 ){
    
    ess_splits <- list()
    split_freq_windows <- list()
    split_freq_runs <- list()
    fails_tmp <- list()
    
    ## ESS ##
    for (i in 1:length(output_tree_parameters[[1]])) {
      ess_splits[[i]] <- names(which(output_tree_parameters[[1]][[i]] < 0))
      ess_splits[[i]] <- output_tree_parameters_raw$ess[[i]][names(output_tree_parameters_raw$ess[[i]]) %in% ess_splits[[i]]]
    }
    names(ess_splits) <- paste("ESS_of_", names(output_tree_parameters[[1]]), sep = "")
    decision_list_trees[[1]] <- ess_splits
    for (i in 1:length(decision_list_trees[[1]])) {
      if( length(decision_list_trees[[1]][[i]]) > 0 ){
        fails <- c(fails, paste(length(decision_list_trees[[1]][[i]]), "splits failed to reach", minimumESS, "for", names(decision_list_trees[[1]][i])))
      }
    }
    
    output_tree_parameters_raw$ess[[2]][which(ess_splits[[2]] == names(output_tree_parameters_raw$ess[[2]]))]
    
    ## Windows ##
    for (i in 1:length(output_tree_parameters[[2]])) {
      split_freq_windows[[i]] <- names(which(output_tree_parameters[[2]][[i]] < 0))
    }
    names(split_freq_windows) <- paste("Windows_of_", names(output_tree_parameters[[2]]), sep = "")
    decision_list_trees[[2]] <- split_freq_windows
    for (i in 1:length(decision_list_trees[[2]])) {
      if ( length(decision_list_trees[[2]][[i]]) > 0 ){
        fails <- c(fails, paste(length(decision_list_trees[[2]][[i]]), "splits failed the split difference test between windows for", names(decision_list_trees[[2]][i])))
      }
    }
    
    ## Runs ##
    if (length(my_runs)>1){
      for (i in 1:length(output_tree_parameters[[3]])) {
        split_freq_runs[[i]] <- names(which(output_tree_parameters[[3]][[i]] < 0))
      }
      names(split_freq_runs) <- paste("Between_", names(output_tree_parameters[[3]]), sep = "")
      decision_list_trees[[3]] <- split_freq_runs
      for (i in 1:length(decision_list_trees[[3]])) {
        if ( length(decision_list_trees[[3]][[i]]) > 0 ){
          fails <- c(fails, paste(length(decision_list_trees[[3]][[i]]), "splits failed the split difference test between runs for", names(decision_list_trees[[3]][i])))
        }
      }
    }

    
    for (i in 1:length(decision_list_trees)) {
      for (j in 1:length(decision_list_trees[[i]])) {
        if( length(decision_list_trees[[i]][[j]]) > 0){
          count_decision <- count_decision + 1
          fails_tmp <- c(fails_tmp,decision_list_trees[[i]][j])
        }
      }
    }
    if (length(fails_tmp) > 0){
      fails_names$tree_parameters <- fails_tmp
    }
    
  }
  
  ##### ONLY LOG FILES #####
  
  # Check ess_runs_cont_param, ks_windows, ks_runs
  if( length(my_runs[[1]]$ptable) > 0 ){
    
    ess_cont <- list()
    ess_value <- list()
    ks_windows <- list()
    ks_runs <- list()
    fails_tmp <- list()
    
    for (i in 1:length(output_continuous_parameters[[1]])) {
      ess_cont[[i]] <- row.names(output_continuous_parameters[[1]])[which(output_continuous_parameters[[1]][i] < 0 )]
      ess_value[[i]] <- output_continuous_parameters_raw$ess[ess_cont[[i]], i]
      names(ess_value[[i]]) <- ess_cont[[i]]
    }
    names(ess_value) <- names(output_continuous_parameters[[1]])
    decision_list_cont[[1]] <- ess_value
    for (i in 1:length(decision_list_cont[[1]])) {
      if( length(decision_list_cont[[1]][[i]]) > 0 ){
        fails <- c(fails, paste(length(decision_list_cont[[1]][[i]]), "parameters failed to reach", minimumESS, "for", names(decision_list_cont[[1]][i])))
      }
    }
    
    for (i in 1:nrow(output_continuous_parameters[[2]])) {
      ks_windows[[i]] <- colnames(output_continuous_parameters[[2]])[which(output_continuous_parameters[[2]][i,] < 0)]
    }
    names(ks_windows) <- rownames(output_continuous_parameters[[2]])
    decision_list_cont[[2]] <- ks_windows
    for (i in 1:length(decision_list_cont[[2]])) {
      if ( length(decision_list_cont[[2]][[i]]) > 0 ){
        fails <- c(fails, paste(length(decision_list_cont[[2]][[i]]), "parameters failed KS test between windows for", names(decision_list_cont[[2]][i])))
      }
    }
    
    
    if( length(my_runs) > 1){
      for (i in 1:nrow(output_continuous_parameters[[3]])) {
        ks_runs[[i]] <- colnames(output_continuous_parameters[[3]])[which(output_continuous_parameters[[3]][i,] < 0)]
      }
      names(ks_runs) <- rownames(output_continuous_parameters[[3]])
      decision_list_cont[[3]] <- ks_runs
      
      for (i in 1:length(decision_list_cont[[3]])) {
        if ( length(decision_list_cont[[3]][[i]]) > 0 ){
          fails <- c(fails, paste(length(decision_list_cont[[3]][[i]]), "parameters failed KS test between runs for", names(decision_list_cont[[3]][i])))
        }
      }
    }


    for (i in 1:length(decision_list_cont)) {
      for (j in 1:length(decision_list_cont[[i]])) {
        if( length(decision_list_cont[[i]][[j]]) > 0){
          count_decision <- count_decision + 1
          fails_tmp <- c(fails_tmp,decision_list_cont[[i]][j])
        }
      }
    }
    
    if (length(fails_tmp) > 0){
      fails_names$continuous_parameters <- fails_tmp
    }
    
  }
  
  if( length(fails) > 0 ){
    tmp <- list()
    for (i in 1:(length(fails))) {
      #tmp <- paste(tmp,fails[i], "\n", fails[i+1], "\n")
      tmp <- paste(tmp,fails[i], "\n")
    }
    fails <- tmp
  }
  
  if(count_decision == 0){
    final_output$message <- "Achieved convergence"
    final_output$converged <- TRUE
  }else{
    final_output$message <- "Failed convergence"
    final_output$converged <- FALSE
    final_output$failed <- fails
    final_output$failed_names <- fails_names

    class(final_output$failed) <- "list.fails"
   
  }
  
  class(output_continuous_parameters_raw) <- "convenience.diag"
  class(output_tree_parameters_raw) <- "convenience.diag"
  
  final_output$continuous_parameters <- output_continuous_parameters_raw
  final_output$tree_parameters <- output_tree_parameters_raw
  
  
  class(final_output) <- "convenience.diag"
  
  final_output
}