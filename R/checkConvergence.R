# For the naive user

checkConvergence <- function(runs, burnin = 0.1, min_split= 0.05, percent =0.01, min_stats = 0.5, max_psrf=1.05){
  
  # load the mcmc output
  my_runs <- loadFiles(runs, burnin)
  
  #initializing a counter
  count = 0
  
  name <- vector("list", length = 0)
  
  ## If we have tree files, we check split freq ##
  
  if( length(my_runs[[1]]$trees) > 0 ){ # if we have tree files, we check split freq
    
    ## First we check if the runs are ok, by comparing windows of the same run ##
    
    vec_splits_windows <- splitFreq(my_runs, windows = T)
    
    for (i in 1:length(vec_splits_windows)) {
      if ( vec_splits_windows[i] > min_split ){
        stop("Try running your MCMC for more iterations")
      }
    }
    
    ## If the runs are ok, we compare the split freq between them ##
    
    vec_split_runs <- splitFreq(my_runs)
    
    for (i in 1:length(vec_split_runs)) {
      if(vec_split_runs[1] < min_split){
        count <- count + 1
      }
    }
    print(paste((count/length(vec_split_runs))*100,"% of the splits are above the minimum value"))
  }
  
  ## If we have log files, we check the cont parameters ##
  
  if( length(my_runs[[1]]$ptable) > 0 ){ # if we have log files, we check the cont parameters
    
    # check if min ESS is achieved for all cont parameters
    min_ess <- minESS(percent)
    ess_runs <- essContParam(runs)
    
    for (i in 1:length(ess_runs)) {
      for (j in 1:length(ess_runs[[1]])) {
        if (ess_runs[[i]][j] < min_ess){
          stop("Try running your MCMC for more iterations")
        }
      }
    }
    
    # check if std error is within 1% of the 95% interval
    std_error <- stderrContParam(runs)
    min_std_error <- stderrMin(runs)
    
    comp <- mapply("-", std_error,min_std_error)
    
    for (i in 1:length(comp)) {
      for (j in 1:length(comp[[1]])) {
        if (comp[[i]][j] < 0){
          stop("Try running your MCMC for more iterations")
        }
      }
    }
    
    #check if stats within runs are ok
    stats_values <- statsContParam(my_runs, windows = TRUE)
    
    for (i in 1:length(stats_values)) {
      for (j in 1:length(stats_values[[1]])){
        for (k in 1:nrow(stats_values[[1]][1])) {
          if (stats_values[[i]][[j]][k] > min_stats){
            stop("Try running your MCMC for more iterations")
          }
        }
      }
    }
    
    #check if psrf within runs is ok
    psrf_windows <- psrfContParams(my_runs, windows = TRUE)
    
    for (i in 1:length(psrf_windows)) {
      for (j in 1:(length(psrf[[1]]$psrf)/2)){
        if( psrf[[i]]$psrf[j] > max_psrf ){
          stop("Try running your MCMC for more iterations")
        }
      }
    }
    
    # check stats between runs
    stats_runs <- statsContParam(my_runs)
    count <- 0
    
    for (i in 1:length(stats_runs)) {
      for (j in 1:length(stats_runs[[1]])){
        for (k in 1:nrow(stats_runs[[1]][1])) {
          if (stats_runs[[i]][[j]][k] > min_stats){
            append(name, rownames(stats_runs[[i]][j][k]))
            print("a")
          }
        }
      }
    }
    
    
    
        
    #check psrf between runs
    psrf_runs <- psrfContParams(my_runs)
    
    
    
  }
  

  
}

