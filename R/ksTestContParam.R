#' Kolmogorov-Smirnov test for the Continuous Parameters
#' 
#' 
#' @param runs A list of rwty.chain type
#' @param namesToExclude Column names to exclude from the calculations, default = bl, Iteration, Likelihood, Posterior, Prior
#' 
#' @return 
#' 
#' @example 
#' # statsContParam(runs, windows = FALSE, namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior")
#' 
#' @export


ksTest <- function(runs, namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior"){
  
  all_df <- vector("list", length = 0)
  
  for (i in 1:length(runs)) {
    #get the cont_param for each run
    all_df[[i]] <- getInfo(runs, i, namesToExclude)
  }
  names_parameters <- names(all_df[[1]])
  
  vecAll <- vector("list", length = 0)
  
  for (df1 in 1:(length(all_df)-1)){
    vecKS <- vector()
    
    for (i in 1:length(all_df[[1]])) {
      
      ksTest <- ks.test( all_df[[df1]][[i]] , all_df[[df1+1]][[i]] )
      vecKS <- c(vecKS, ksTest$statistic)
      
    }
    vecAll <- c(vecAll, vecKS)
  }
  
  df_ks <- data.frame(matrix(unlist(vecAll), nrow = length(runs)-1, byrow=T), stringsAsFactors = F)
  colnames(df_ks) <- names_parameters
  
  return(df_ks)
  
}
