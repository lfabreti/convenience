#' Measures of Statistical Dispersion for the Continuous Parameters
#' 
#' Calculates mean, median, standard deviation, 2.5% quantile and 97.5% quantile for all continuous parameters
#' 
#' @param runs A list of rwty.chain type
#' @param windows A boolean to set if the calculations are in between windows of the same run, default = FALSE
#' @param namesToExclude Column names to exclude from the calculations, default = bl, Iteration, Likelihood, Posterior, Prior
#' 
#' @return a list of data frames with the calculated values
#' 
#' @example 
#' # statsContParam(runs, windows = FALSE, namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior")
#' 
#' @export

statsContParam <- function(runs, windows = FALSE, namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior"){
  
  list_df <- list()
  
  # comparison between runs
  if(!windows){
    
    all_df <- vector("list", length = 0)
    
    for (i in 1:length(runs)) {
      
      #get the cont_param for each run
      all_df[[i]] <- getInfo(runs, i, namesToExclude)
    }
    
    vecAll <- vector("list", length = 0)
    
    for (df1 in 1:(length(all_df)-1)){
      
      vecStats <- vector(length = 0)
      
      vecStats <- c(vecStats,calcRelativeDiff(all_df[[df1]],all_df[[df1+1]],mean))
      
      vecStats <- c(vecStats,calcRelativeDiff(all_df[[df1]],all_df[[df1+1]],median))
      
      vecStats <- c(vecStats,calcRelativeDiff(all_df[[df1]],all_df[[df1+1]],sd))
      
      vecStats <- c(vecStats,calcRelativeDiff(all_df[[df1]],all_df[[df1+1]],quants))
      
      vecAll <- c(vecAll, vecStats)
      
    }
    
    size_df = (length(vecAll))/(length(runs)-1) # get the size of each data frame
    
    n_row <- size_df/5 # get the number of continuous parameters
    
    names_param <- names(vecAll[1:n_row]) # get the names of the continuous parameters
    
    name_stats <- c("Mean", "Median", "Std deviation", "2.5% Quantile", "97.5% Quantile") # name of the stats calculated
    
    # get all the information into a list of data frames
    for (i in 1:(length(runs)-1)){
      
      start <- ((i-1)*size_df)+1
      end <- size_df*i
      
      df_stats <- data.frame(matrix(unlist(vecAll[start:end]), nrow = n_row, byrow=T), stringsAsFactors = F)
      rownames(df_stats) <- names_param
      colnames(df_stats) <- name_stats
      
      list_df[[i]] <- df_stats
      
    }
    
  }
  
  # comparison between windows of the same run
  else{
    
    for (i in 1:length(runs)) {
      
      cont_param <- getInfo(runs, i, splitWindows = TRUE,  namesToExclude)
      
      vecStats <- vector(length = 0)
      
      vecStats <- c(vecStats,calcRelativeDiff(cont_param[[1]],cont_param[[2]],mean))
      
      vecStats <- c(vecStats,calcRelativeDiff(cont_param[[1]],cont_param[[2]],median))
      
      vecStats <- c(vecStats,calcRelativeDiff(cont_param[[1]],cont_param[[2]],sd))
      
      vecStats <- c(vecStats,calcRelativeDiff(cont_param[[1]],cont_param[[2]],quants))
      
      size_df <- length(vecStats) # get the size of each data frame
      
      n_row <- size_df/5 # get the number of continuous parameters
      
      names_param <- names(vecStats[1:n_row]) # get the names of the continuous parameters
      
      name_stats <- c("Mean", "Median", "Std deviation", "2.5% Quantile", "97.5% Quantile") # name of the stats calculated
      
      df_stats <- data.frame(matrix(unlist(vecStats), nrow = n_row, byrow=T), stringsAsFactors = F)
      rownames(df_stats) <- names_param
      colnames(df_stats) <- name_stats
      
      list_df[[i]] <- df_stats
      
    }
  }
  
  return(list_df)
}
