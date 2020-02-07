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


essContParam <- function(runs, namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior") {

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

  return(df_ess)
}
