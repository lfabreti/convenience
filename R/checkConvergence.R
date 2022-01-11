#' Check Convergence full assessment
#' 
#' Check output from phylogenetic MCMC analysis for convergence diagnostics
#' 
#' @importFrom plyr ldply
#' @importFrom stats var setNames
#' 
#' @param path Path to directory containing all files from the same analysis
#' @param list_files List of files to check for convergence
#' @param format The format of the phylogenetic output. Current supported formats are: "revbayes", "mb", "beast", "*beast", "phylobayes"
#' @param control List of arguments to the function. Includes tracer (to calculate the ESS with Tracer method), burn-in, precision and names of parameters to exclude from the analysis
#' 
#' @return List of type convenience.diag
#' 
#' @examples 
#' \dontrun{
#' checkConvergence( path ,
#'  format = "revbayes" )
#' }
#' 
#' @export


checkConvergence <- function(path = NULL, list_files = NULL, format = "revbayes", control = makeControl()){
  
  ##### First let's check the function arguments #####
  if( is.null(control$tracer) ){
    tracer <- TRUE
  } else {
    tracer <- control$tracer
  }
  
  if ( is.null(control$precision) ){
    precision <- 0.01
  } else {
    precision <- control$precision
  }
  
  if ( is.null(control$burnin) ){
    burnin <- 0.0
  } else {
    burnin <- control$burnin
  }
  
  if ( is.null(control$namesToExclude) ){
    namesToExclude <- "br_lens|bl|Iteration|Likelihood|Posterior|Prior|Gen|LnL|LnPr|state|joint|prior|likelihood|time|loglik|iter|topo|Replicate_ID|Sample|posterior|it"
  } else {
    namesToExclude <- control$namesToExclude
  }
  #####
  
  final_output <- list()
  
  # calculating minimum ESS according to precision
  minimumESS <- minESS(precision)
  minimumESS_windows <- round(minimumESS/5)
  
  # list of results - thresholds
  output_tree_parameters <- list()
  output_continuous_parameters <- list()
  
  # list of results
  output_tree_parameters_raw <- list()
  output_continuous_parameters_raw <- list()
  
  # Load MCMC output
  if ( !is.null(path)){
    # load the mcmc output
    my_runs <- loadFiles(path, format = format)
  }else {
    my_runs <- loadFiles(list_files = list_files, format = format)
  }
  if(burnin > 0){
    my_runs <- removeBurnin(my_runs, burnin)
  }
  
  #####   BURN-IN   #####
  print("Calculating burn-in")
  my_runs_aux <- my_runs
  while (burnin <= 0.5) {
    
    list_control <- 0 
    
    # tree files #
    if( length(my_runs[[1]]$trees) > 0 ){
      
      splits_windows <- splitFreq(my_runs, windows = T)
      if( minimumESS_windows == 125 ) {
        fdir <- system.file("thresholds/expectedDiff_125.rds", package = "convenience")
        exp_diff_windows <- readRDS(fdir)
      }
      else exp_diff_windows <- expectedDiffSplits(minimumESS_windows)
      
      for (i in 1:ncol(splits_windows)) {
        splits_windows[2,i]$listFrequencies <- round(splits_windows[2,i]$listFrequencies, digits = 2)
      }
      
      col_names <- list()
      for (i in 1:ncol(splits_windows)) {
        col_names <- c(col_names, paste("Run_", i, sep = ""))
      }
      colnames(splits_windows) <- col_names
      
      
      results_splits <- list()
      for (i in 1:ncol(splits_windows)) {
        results <- vector()
        # index gets the indexes sorted
        index <- which( exp_diff_windows[1,] %in% splits_windows[2,i]$listFrequencies )
        vecNames <- vector()
        
        for (j in index) {
          for (z in 1:length(splits_windows[2,i]$listFrequencies)) {
            if( exp_diff_windows[1,j] == splits_windows[2,i]$listFrequencies[z]){
              results <- c(results, exp_diff_windows[2,j] - splits_windows[1,i]$listDiffSplits[z])
              vecNames <- c( vecNames, names(splits_windows[1,i]$listDiffSplits[z]))
            }
          }
        }
        results_splits[[i]] <- results
        names(results_splits[[i]]) <- vecNames
      }
      for (i in 1:length(results_splits)) {
        names(results_splits)[i] <- paste("Run_", i, sep = "")
      }
      
    }
    
    # log files #
    if( length(my_runs[[1]]$ptable) > 0 ){
      
      ## KS score ##
      ks_windows <- ksTest(my_runs, windows = TRUE, namesToExclude = namesToExclude)
      
      ks_limits_windows <-vector()
      for (i in 1:length(my_runs)) {
        if (typeof(ks_windows) == "double"){
          ks_limits_windows <- c(ks_limits_windows, ksThreshold(0.01, minimumESS_windows))
        }else{
          for (j in 1:ncol(ks_windows)) {
            ks_limits_windows <- c(ks_limits_windows, ksThreshold(0.01, minimumESS_windows))
          }
        }
      }
      
      ks_limits_windows <- data.frame(matrix(unlist(ks_limits_windows), nrow = length(my_runs), byrow = T), stringsAsFactors = F)
      colnames(ks_limits_windows) <- names(ks_windows)
      for (i in 1:length(my_runs)) {
        rownames(ks_limits_windows)[i] <- paste("Run_", i, sep = "")
        rownames(ks_windows)[i] <- paste("Run_", i, sep = "")
      }
      
      ks_windows <- ks_limits_windows - ks_windows
    }
    
    if( length(my_runs[[1]]$trees) > 0 & length(my_runs[[1]]$ptable) > 0 ){
      for (i in 1:length(results_splits)) {
        if( length(results_splits[[i]]) > 0 | any(ks_windows[i,] < 0) ) list_control <- list_control+1
      }
    } else if ( length(my_runs[[1]]$ptable) == 0){
      for (i in 1:length(results_splits)) {
        if( length(results_splits[[i]]) > 0 ) list_control <- list_control+1
      }
    }
    
    if(list_control > 0){
      burnin <- burnin + 0.1
      my_runs <- removeBurnin(my_runs_aux, burnin)
    } 
    else break
    
  }
  if( burnin > 0.5 & list_control > 0) stop("Burn-in too large")
  
  # Name runs
  if( length(my_runs) > 1 ){
    compar_names <- vector()
    for (r1 in 1:(length(my_runs)-1)) {
      for (r2 in (r1+1):length(my_runs)) {
        compar_names <- c(compar_names, paste("Run_", r1, "_Run_", r2, sep = ""))
      }
    }
  }
  
  #################################################
  ## IF WE HAVE TREE FILES, WE CHECK THE SPLITS  ##
  #################################################
  
  if( length(my_runs[[1]]$trees) > 0 ){ # if we have tree files, we check split freq and ESS of splits
    
    print("Analyzing tree parameters")
    
    ## Check splits with freq above 0.975 or below 0.025 ##
    output_tree_parameters_raw$exclude_high <- check.clades.freq(my_runs, 0.975)
    output_tree_parameters_raw$exclude_low <- check.clades.freq(my_runs, 0.025)
    
    ## Split frequency ##
    
    output_tree_parameters_raw$frequencies <- splitFreq(my_runs)
    
    ## ESS ##
    
    ess_runs_splits <- essSplitFreq(my_runs, tracer)
    output_tree_parameters_raw$ess <- ess_runs_splits
    ess_run_splits_aux <- ess_runs_splits
    
    for (i in 1:length(ess_runs_splits)) {
      ess_runs_splits[[i]] <- ess_runs_splits[[i]] - minimumESS
    }
    
    output_tree_parameters$ess <- ess_runs_splits
    
    ## Compare runs ##
    
    if( length(my_runs) > 1 ){
      splits_runs <- output_tree_parameters_raw$frequencies
      
      if( minimumESS == 625) {
        fdir <- system.file("thresholds/expectedDiff_625.rds", package = "convenience")
        exp_diff_runs <- readRDS(fdir)
      }
      else  exp_diff_runs <- expectedDiffSplits(minimumESS)
      
      results_splits_runs <- list()
      if( length(splits_runs) > 0 ){
        for (i in 1:ncol(splits_runs)) {
          splits_runs[2,][[i]] <- round(splits_runs[2,][[i]], digits = 2)
        }
        
        colnames(splits_runs) <- compar_names
        output_tree_parameters_raw$compare_runs <- splits_runs
        
        for (i in 1:ncol(splits_runs)) {
          results <- vector()
          index <- which( exp_diff_runs[1,] %in% splits_runs[2,][[i]] )
          vecNames <- vector()
          
          for (j in index) {
            for (z in 1:length(splits_runs[2,][[i]])) {
              if( exp_diff_runs[1,j] == splits_runs[2,][[i]][z]){
                results <- c(results, (exp_diff_runs[2,j] - splits_runs[1,][[i]][z]))
                vecNames <- c( vecNames, names(splits_runs[1,][[i]][z]))
              }
            }
          }
          results_splits_runs[[i]] <- results
          names(results_splits_runs[[i]]) <- vecNames
        }
        
        names(results_splits_runs) <- compar_names
      }
      
      output_tree_parameters$compare_runs <- results_splits_runs
    }
    
    # Split frequency per run #
    if(length(my_runs) > 1 ){
      all_df <- vector("list", length = 0)
      for (i in 1:length(my_runs)){
        
        x <- getInfo(my_runs, i, trees = TRUE)
        cladefreqs <- clade.freq.named(x, start = 1, end = length(x))
        
        all_df[[i]] <- as.numeric(cladefreqs$cladefreqs_post)
        names(all_df[[i]]) <- as.character(cladefreqs$cladenames_post)
      }
      df_freqs <- plyr::ldply(all_df, rbind)
      df_freqs <- t(df_freqs)
      colnames(df_freqs) <- names(my_runs)
      output_tree_parameters_raw$freq_per_run <- df_freqs
    }
    
  }
  
  
  ######################################################### 
  ## IF WE HAVE LOG FILES, WE CHECK THE CONTINUOUS PARAM ##
  #########################################################
  
  if( length(my_runs[[1]]$ptable) > 0 ){ 
    
    print("Analyzing continuous parameters")
    
    ## Check parameters with var = 0 and exclude from convergence assessment ##
    
    param_names <- vector("list", length = length(my_runs))
    for (i in 1:length(my_runs)) {
      aux <- vector()
      for ( j in 1:length(my_runs[[i]]$ptable) ) {
        if( var(my_runs[[i]]$ptable[j]) == 0 ){
          aux <- c(aux, j)
        }
      }
      if (length(aux) > 0 ){
        param_names[[i]] <- names(my_runs[[i]]$ptable[aux])
        my_runs[[i]]$ptable <- my_runs[[i]]$ptable[,-aux]
      }
    }
    output_continuous_parameters_raw$exclude <- param_names
    
    ## Means of continuous parameters ##
    
    output_continuous_parameters_raw$means <- meanContParam(my_runs, namesToExclude = namesToExclude)
    
    ## ESS ##
    
    ess_runs_cont_param <- essContParam(my_runs, namesToExclude = namesToExclude, tracer = tracer)
    output_continuous_parameters_raw$ess <- ess_runs_cont_param
    
    ess_runs_cont_param <- ess_runs_cont_param - minimumESS
    output_continuous_parameters$ess <- ess_runs_cont_param
    
    ## Compare runs ##
    
    if( length(my_runs) > 1 ){
      ks_runs <- ksTest(my_runs, windows = F, namesToExclude = namesToExclude)
      
      ks_limits <-vector()
      count<-0
      for (df1 in 1:(length(my_runs)-1)){
        for (df2 in (df1+1):length(my_runs)) {
          for (i in 1:ncol(ks_runs)) {
            ks_limits <- c(ks_limits, ksThreshold(0.01,minimumESS))
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
      if(length(output_tree_parameters[[1]][[i]])>0) ess_splits[[i]] <- names(which(output_tree_parameters[[1]][[i]] < 0))
      if(length(output_tree_parameters[[1]][[i]])>0) ess_splits[[i]] <- output_tree_parameters_raw$ess[[i]][names(output_tree_parameters_raw$ess[[i]]) %in% ess_splits[[i]]]
    }
    if(length(ess_splits) > 0) names(ess_splits) <- paste("ESS_of_", names(output_tree_parameters[[1]]), sep = "")
    decision_list_trees[[1]] <- ess_splits
    if( length(decision_list_trees[[1]]) > 0 ) for (i in 1:length(decision_list_trees[[1]])) {
      if( length(decision_list_trees[[1]][[i]]) > 0 ){
        fails <- c(fails, paste(length(decision_list_trees[[1]][[i]]), "splits failed to reach", minimumESS, "for", names(decision_list_trees[[1]][i])))
      }
    }
    
    ## Runs ##
    if (length(my_runs) > 1){
      if ( length(output_tree_parameters_raw) > 3 & length(output_tree_parameters[[2]]) > 0 ){
        for (i in 1:length(output_tree_parameters[[2]])) {
          split_freq_runs[[i]] <- names(which(output_tree_parameters[[2]][[i]] <= 0))
        }
        names(split_freq_runs) <- paste("Between_", names(output_tree_parameters[[2]]), sep = "")
        decision_list_trees[[2]] <- split_freq_runs
        for (i in 1:length(decision_list_trees[[2]])) {
          if ( length(decision_list_trees[[2]][[i]]) > 0 ){
            fails <- c(fails, paste(length(decision_list_trees[[2]][[i]]), "splits failed the split difference test between runs for", names(decision_list_trees[[2]][i])))
          }
        }
      }
    }
    
    for (i in 1:length(decision_list_trees)) {
      if(length(decision_list_trees[[i]]) > 0 ) for (j in 1:length(decision_list_trees[[i]])) {
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
    
    
    if( length(my_runs) > 1){
      for (i in 1:nrow(output_continuous_parameters[[2]])) {
        ks_runs[[i]] <- colnames(output_continuous_parameters[[2]])[which(output_continuous_parameters[[2]][i,] < 0)]
      }
      names(ks_runs) <- rownames(output_continuous_parameters[[2]])
      decision_list_cont[[2]] <- ks_runs
      
      for (i in 1:length(decision_list_cont[[2]])) {
        if ( length(decision_list_cont[[2]][[i]]) > 0 ){
          fails <- c(fails, paste(length(decision_list_cont[[2]][[i]]), "parameters failed KS test between runs for", names(decision_list_cont[[2]][i])))
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
      tmp <- paste(tmp,fails[i], " \n ", sep = "")
    }
    fails <- tmp
  }
  
  
  ##### SUMMARIZING RESULTS #####
  
  final_output$burnin <- burnin
  
  message_list <- list()
  message_complete <- list()
  if(count_decision == 0){
    message_list <- paste(message_list, "ACHIEVED CONVERGENCE", "\n")
    message_list <- paste(message_list, " \n")
  }else{
    message_list <- paste(message_list, "FAILED CONVERGENCE", "\n")
    message_list <- paste(message_list, " \n")
    message_list <- paste(message_list, fails)
    message_list <- paste(message_list, " \n")
  } 
  
  message_list <- paste(message_list, "BURN-IN SET AT", burnin, "\n")
  message_list <- paste(message_list, " \n")
  
  message_complete <- paste(message_list)
  message_complete <- paste(message_complete, " \n")
  # reporting split with too high or too low frequency
  if( length(my_runs[[1]]$trees) > 0 ){
    if ( length(output_tree_parameters_raw$exclude_high[[1]]) > 0 ) message_complete <- paste(message_complete, "SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT \n")
    for (i in 1:length(my_runs)) {
      if ( length(output_tree_parameters_raw$exclude_high[[i]]) > 0 ) message_complete <- paste(message_complete, "FREQUENCY HIGHER THAN 0.975 FOR RUN", i, "\n")
      if ( length(output_tree_parameters_raw$exclude_high[[i]]) > 0 ){
        for (j in 1:length(output_tree_parameters_raw$exclude_high[[i]])) {
          message_complete <- paste(message_complete, "    ", output_tree_parameters_raw$exclude_high[[i]][j], "\n")
        }
        message_complete <- paste(message_complete, " \n")
      }
      if ( length(output_tree_parameters_raw$exclude_low[[i]]) > 0 ) message_complete <- paste(message_complete, "FREQUENCY LOWER THAN 0.025 FOR RUN", i, "\n")
      if ( length(output_tree_parameters_raw$exclude_low[[i]]) > 0 ){
        for (j in 1:length(output_tree_parameters_raw$exclude_low[[i]])) {
          message_complete <- paste(message_complete, "    ", output_tree_parameters_raw$exclude_low[[i]][j], "\n")
        }
        message_complete <- paste(message_complete, " \n")
      }
    }
  }
  
  # reportin cont param with no variance
  if( length(my_runs[[1]]$ptable) > 0 ){
    if( length(output_continuous_parameters_raw$exclude[[1]]) > 0 ) message_complete <- paste(message_complete, "CONTINUOUS PARAMETERS WITH NO VARIANTION AND EXCLUDED FROM CONVERGENCE ASSESSMENT \n")
    for (i in 1:length(my_runs)) {
      if( length(output_continuous_parameters_raw$exclude[[i]]) > 0 ) message_complete <- paste(message_complete, "RUN", i, "\n")
      if( length(output_continuous_parameters_raw$exclude[[i]]) > 0 ){
        for (j in 1:length(output_continuous_parameters_raw$exclude[[i]])) {
          message_complete <- paste(message_complete, "    ", output_continuous_parameters_raw$exclude[[i]][j], "\n")
        }
        message_complete <- paste(message_complete, " \n")
      }
    }
  }
  
  # reporting split with lowest ESS
  if( length(my_runs[[1]]$trees) > 0 ){
    message_list <- paste(message_list, "LOWEST SPLIT ESS \n")
    message_complete <- paste(message_complete, "LOWEST SPLIT ESS \n")
    for (i in 1:length(my_runs)) {
      message_list <- paste(message_list, "     RUN", i, "->", names(which.min(output_tree_parameters_raw$ess[[i]])), round(min(output_tree_parameters_raw$ess[[i]], na.rm = T), digits = 2), "\n")
      message_complete <- paste(message_complete, "     RUN", i, "->", names(which.min(output_tree_parameters_raw$ess[[i]])), round(min(output_tree_parameters_raw$ess[[i]], na.rm = T), digits = 2), "\n")
    }
    message_list <- paste(message_list, " \n")
    message_complete <- paste(message_complete, " \n")
  }
  
  # reporting cont param with lowest ESS
  if( length(my_runs[[1]]$ptable) > 0 ){
    message_list <- paste(message_list, "LOWEST CONTINUOUS PARAMETER ESS \n")
    message_complete <- paste(message_complete, "LOWEST CONTINUOUS PARAMETER ESS \n")
    for (i in 1:length(my_runs)){
      message_list <- paste(message_list, "     RUN", i, "->", row.names(output_continuous_parameters_raw$ess)[which.min(output_continuous_parameters_raw$ess[[i]])],  round(min(output_continuous_parameters_raw$ess[[i]]), digits = 2), "\n")
      message_complete <- paste(message_complete, "     RUN", i, "->", row.names(output_continuous_parameters_raw$ess)[which.min(output_continuous_parameters_raw$ess[[i]])],  round(min(output_continuous_parameters_raw$ess[[i]]), digits = 2), "\n")
    }
    message_list <- paste(message_list, " \n")
    message_complete <- paste(message_complete, " \n")
  }
  
  # reporting commands that the user can use to visualize the output
  if( length(my_runs[[1]]$ptable) > 0 ){
    message_list <- paste(message_list, "To check the calculated parameters for the continuous parameters type: \n")
    message_list <- paste(message_list, "     Means: output$continuous_parameters$means \n")
    message_list <- paste(message_list, "     ESS: output$continuous_parameters$ess \n")
    if( length(my_runs) > 0 ) message_list <- paste(message_list, "     KS score: output$continuous_parameters$compare_runs \n")
    message_list <- paste(message_list, " \n")
  }
  if( length(my_runs[[1]]$trees) > 0 ){
    message_list <- paste(message_list, "To check the calculated parameters for the splits type: \n")
    message_list <- paste(message_list, "     Frequencies of splits: output$tree_parameters$frequencies \n")
    message_list <- paste(message_list, "     ESS: output$tree_parameters$ess \n")
    if( length(my_runs) > 0 ) message_list <- paste(message_list, "     Difference in frequencies: output$tree_parameters$compare_runs \n")
    message_list <- paste(message_list, " \n")
  }
  message_list <- paste(message_list, "To check the full summary message with splits and parameters excluded from the analysis type: \n")
  message_list <- paste(message_list, "     output$message_complete \n")
  
  class(message_list) <- "list.fails"
  class(message_complete) <- "list.fails"
  
  if(count_decision == 0){
    final_output$message <- message_list
    final_output$message_complete <- message_complete
    final_output$converged <- TRUE
  }else{
    final_output$message <- message_list
    final_output$message_complete <- message_complete
    final_output$converged <- FALSE
    final_output$failed <- fails
    final_output$failed_names <- fails_names
    
    class(final_output$failed) <- "list.fails"
    
  }
  
  ##### Changing formats for better looking output #####
  
  # format of diff splits comparison between runs
  tmp_freq <- list()
  tmp_diff <- list()
  if( length(my_runs[[1]]$trees) > 0 ){
    if( length(my_runs) > 1){
      if( length(output_tree_parameters_raw$compare_runs) > 0 ){
        for (i in 1:ncol(output_tree_parameters_raw$compare_runs)) {
          tmp_diff <- c(tmp_diff, output_tree_parameters_raw$compare_runs[1,i])
          tmp_freq <- c(tmp_freq, output_tree_parameters_raw$compare_runs[2,i])
        }
        names(tmp_diff) <- colnames(output_tree_parameters_raw$compare_runs)
        names(tmp_freq) <- colnames(output_tree_parameters_raw$compare_runs)
        output_tree_parameters_raw$compare_runs <- tmp_diff
        output_tree_parameters_raw$frequencies <- tmp_freq
        
      }
    }
  }
  #######
  
  
  # format of ESS for splits
  if(length(my_runs[[1]]$trees) > 0){

    df_2 <- plyr::ldply(ess_run_splits_aux, rbind)
    output_tree_parameters_raw$ess <- setNames(data.frame(t(df_2[,-1]), row.names = colnames(df_2)[-1]), df_2[,1])
  }
  
  
  ####
  
  final_output$continuous_parameters <- output_continuous_parameters_raw[-1]
  final_output$tree_parameters <- output_tree_parameters_raw[-(1:2)]
  
  class(final_output$continuous_parameters) <- "convenience.table"
  class(final_output$tree_parameters) <- "convenience.table"
  class(final_output) <- "convenience.diag"
  
  final_output
}