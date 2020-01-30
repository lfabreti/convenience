#' Standard error of the contiuoun parameters
#' 
#' Calculates the standard error of the continuous parameters
#' 
#' @param runs A list of rwty.chain type
#' @param namesToExlude Column names to exclude from the calculations, default = bl, Iteration, Likelihood, Posterior, Prior
#' 
#' @return 
#' 
#' @example 
#' 
#' @export

stderrContParam <- function(runs,namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior"){

  vecStderr <- vector("double", length = 0)
  
  list_df <- vector("double", length = 0)
  
  for (i in 1:length(runs)) {
    
    cont_param <- getInfo(runs, i, namesToExclude)
    vecStderr <- (sapply(cont_param, se))
    vecstdDev <- sapply(cont_param, sd)# continue
    
    names_param <- names(vecStderr) # get the names of the continuous parameters
    
    df <- data.frame(matrix(unlist(vecStderr), nrow = length(vecStderr), byrow = T), stringsAsFactors = F)
    rownames(df) <- names_param
    colnames(df) <- "Std error"
    
    list_df[[i]] <- df
  }
  
  return(list_df)
}