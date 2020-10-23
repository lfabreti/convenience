#' Kolmogorov-Smirnov test for the continuous parameters
#' 
#' Calculates the KS test for the comparison between runs or windows of the same run
#' 
#' @param runs A list of rwty.chain type
#' @param windows A boolean to set if the calculation is within runs or between runs
#' @param namesToExclude Column names to exclude from the calculations, default = bl, Iteration, Likelihood, Posterior, Prior
#' 
#' @return 
#' 
#' @export


ksTest <- function(runs, windows = FALSE, namesToExclude){
  
  if(!windows){
    all_df <- vector("list", length = 0)
    
    for (i in 1:length(runs)) {
      #get the cont_param for each run
      all_df[[i]] <- getInfo(runs, i, namesToExclude)
    }
    names_parameters <- names(all_df[[1]])
    
    vecAll <- vector("list", length = 0)
    count <- 0
    
    # In case there is only one continuous parameter
    if( typeof(all_df[[1]]) == "double" ){
      for (df1 in 1:(length(all_df)-1)) {
        for (df2 in (df1+1):length(all_df)) {
          
          ksTest <- ks.test( all_df[[df1]], all_df[[df2]])
          vecAll <- c(vecAll, ksTest$statistic)
          count <- count + 1
        }
      }
    } else {
      for (df1 in 1:(length(all_df)-1)){
        for (df2 in (df1+1):length(all_df)) {
          vecKS <- vector()
          
          for (i in 1:length(all_df[[1]])) {
            
            ksTest <- ks.test( all_df[[df1]][[i]] , all_df[[df2]][[i]] )
            vecKS <- c(vecKS, ksTest$statistic)
            
          }
          vecAll <- c(vecAll, vecKS)
          count <- count+1
        }
      }
    }
    
    
    df_ks <- data.frame(matrix(unlist(vecAll), nrow = count, byrow=T), stringsAsFactors = F)
    colnames(df_ks) <- names_parameters
    
  }
  
  else{
    
    vecAll <- vector()
    
    for (i in 1:length(runs)) {
      
      cont_param <- getInfo(runs, i, splitWindows = TRUE, namesToExclude)
      names_parameters <- names(cont_param[[1]])
      vecKS <- vector()
      
      # In case there is only one continuous parameter
      if ( typeof(cont_param[[1]]) == "double"){
        ks_test <- ks.test( cont_param[[1]], cont_param[[2]])
        vecAll <- c(vecAll, ks_test$statistic)
      }
      
      else{
        for (j in 1:length(cont_param[[1]])) {
          ks_test <- ks.test( cont_param[[1]][[j]] , cont_param[[2]][[j]] )
          vecKS <- c(vecKS, ks_test$statistic)
        }
        vecAll <- c(vecAll, vecKS)
      }
    }
    
    df_ks <- data.frame(matrix(unlist(vecAll), nrow = length(runs), byrow=T), stringsAsFactors = F)
    colnames(df_ks) <- names_parameters
    
  }
  
  return(df_ks)
  
}
