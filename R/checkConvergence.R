# For the naive user

checkConvergence <- function(runs, burnin = 0.1, min_split = 0.05, percent = 0.01, min_stats = 0.5, max_psrf = 1.05){
  
  # load the mcmc output
  my_runs <- loadFiles(runs, burnin, format = "revbayes")
  
  # initializing a counter
  count = 0
  
  # initializing strings for name of parameters
  name_psrf <- character(0)
  names_ess <- character(0)
  ks_result <- vector("list", length = 0)
  
  
  #################################################
  ## IF WE HAVE TREE FILES, WE CHECK THE SPLITS  ##
  #################################################
  
  if( length(my_runs[[1]]$trees) > 0 ){ # if we have tree files, we check split freq
    
    ## First we check if the runs are ok, by comparing windows of the same run ##
    
    vec_splits_windows <- splitFreq(my_runs, windows = T)
    
    for (i in 1:length(vec_splits_windows)) {
      if ( vec_splits_windows[i] > min_split ){
        stop("Try running your MCMC for more iterations, splits between windows failed")
      }
    }
    
    ## If the runs are ok, we compare the split freq between them ##
    
    vec_split_runs <- splitFreq(my_runs)
    
    for (i in 1:length(vec_split_runs)) {
      if(vec_split_runs[1] < min_split){
        count <- count + 1
      }
    }
    print(paste((count/length(vec_split_runs))*100,"% of the splits are below the minimum value"))
  }
  
  
  #########################################################
  ## IF WE HAVE LOG FILES, WE CHECK THE CONTINUOUS PARAM ##
  #########################################################
  
  if( length(my_runs[[1]]$ptable) > 0 ){ 
    
    # check if min ESS is achieved for all cont parameters
    min_ess <- minESS(percent)
    ess_runs <- essContParam(my_runs)
    
    
    for (i in 1:length(ess_runs)) {
      for (j in 1:length(ess_runs[[1]])) {
        if (ess_runs[[i]][j] < min_ess){
          names_ess <- c(names_ess, (rownames(ess_runs)[j]))
        }
      }
    }
    
    if (length(names_ess)>0){
      print(paste("This parameters have ESS below", min_ess, ":"))
      print(unique(names_ess))
    }
    
    
    # check psrf between runs
    psrf_runs <- psrfContParams(my_runs)
    
    for (i in 1:(length(psrf_runs[[1]])/2)) {
      if( psrf_runs[[1]][i] > max_psrf){
        print(psrf_runs[[1]][i])
        name_psrf <- c(name_psrf, (rownames(psrf_runs[[1]])[i]))
        print( rownames(psrf_runs[[1]])[i] )
      }
    }
    
    if (length(name_psrf)>0){
      print(paste("The following parameters failed in PSRF:", name_psrf))
    }
    
    else{
      print(paste("PSRF values are below", max_psrf))
    }
    
  }
  
  all_df <- vector("list", length = 0)
  
  for (i in 1:length(my_runs)) {
    
    #get the cont_param for each run
    all_df[[i]] <- getInfo(runs, i, namesToExclude)
  }
  
  for (df1 in 1:(length(all_df)-1)){
    for (i in 1:length(all_df[[df1]])) {
      ks_result <-c( ks_result, (ks.test(all_df[[df1]][[i]], all_df[[df1+1]][[i]])))
    }
  }

  
}

