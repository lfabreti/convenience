# Calculate ess for the continuous parameters

#' ESS for the continuous parameters
#' 
#' Calculates the Effective Sample Size for the continuous parameters
#' 
#' @param runs A list of rwty.chain type
#' @param namesToExclude Column names to exclude from the calculations, default = bl, Iteration, Likelihood, Posterior, Prior
#' 
#' @return
#' 
#' @example 
#' 
#' @export


essContParam <- function(runs, windows=FALSE, namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior") {

  if(!windows){
    vecEss <- vector("list", length = 0)
    for (i in 1:length(runs)) {
      cont_param <- getInfo(runs, i, namesToExclude)
      
      for (rows in 1:length(cont_param)) {
        
        if ( is.nan( ess(cont_param[[rows]]) ) ) {
          print(paste(names(cont_param[rows]), " ess is not a number!"))
        }
      }
      
      vecEss <- c(vecEss, ess(cont_param))
    }
    
    n_param <- length(vecEss)/length(runs)
    name_param <- names(vecEss[1:n_param])
    
    df_ess <- data.frame(matrix(unlist(vecEss), nrow = n_param, byrow = T), stringsAsFactors = F)
    rownames(df_ess) <- name_param
    
    for (i in 1:ncol(df_ess)) {
      colnames(df_ess)[i] <- paste("ESS run", i)
      
    }
  }
  
  else{
    vecESS <- vector()
    
    for (i in 1:length(runs)) {
      cont_param <- getInfo(runs, i, splitWindows = TRUE, namesToExclude)
      vecESS <- c(vecESS, effectiveSize(cont_param[[1]]))
      vecESS <- c(vecESS, effectiveSize(cont_param[[2]]))
    }
    
    n_param <- length(vecESS)/(length(runs)*2)
    name_param <- names(vecESS[1:n_param])
    
    df_ess <- data.frame(matrix(unlist(vecESS), nrow = length(runs)*2, byrow = T), stringsAsFactors = F)
    colnames(df_ess) <- name_param
    name_runs <- vector()
    for (i in 1:length(runs)) {
      name_runs <- c(name_runs, paste("Run", i, " window 1"), paste("Run", i, " window 2") )
    }
    rownames(df_ess) <- name_runs
  }
  
  return(df_ess)
}
