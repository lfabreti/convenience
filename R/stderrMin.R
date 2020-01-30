#' Minimum Standard Error
#' 
#' Calculates the minimum standard error of the continuous parameters, based on the 95% interval of the distribution of the continuous parameter
#' 
#' @param runs A list of rwty.chain type
#' @param namesToExlude Column names to exclude from the calculations, default = bl, Iteration, Likelihood, Posterior, Prior
#' 
#' @return 
#' 
#' @example 
#' 
#' @export

stderrMin <- function(runs, percent=0.01 ,namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior"){
  
  vecstdDev <- vector("double", length = 0)
  
  list_df <- vector("double", length = 0)
  
  for (i in 1:length(runs)) {
    
    cont_param <- getInfo(runs, i, namesToExclude)
    #vecstdDev <- sapply(cont_param, abc)# continue
    vecstdDev <- mapply(abc, cont_param, percent = percent)
    
    names_param <- names(vecstdDev) # get the names of the continuous parameters
    
    df <- data.frame(matrix(unlist(vecstdDev), nrow = length(vecstdDev), byrow = T), stringsAsFactors = F)
    rownames(df) <- names_param
    colnames(df) <- "Min std error"
    
    list_df[[i]] <- df
  }
  
  return(list_df)
}
