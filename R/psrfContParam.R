#' PSRF for the Continuous Parameters
#' 
#' Calculate the Potential Scale Reduction Factor for the continuous parameters
#' 
#' @param runs A list of rwty.chain type
#' @param windows
#' @param conf_interval
#' @param nameToExlude Column names to exclude from the calculations, default = bl, Iteration, Likelihood, Posterior, Prior
#' 
#' @return 
#' 
#' @example 
#' 
#' @export

psrfContParams <- function(runs, windows=FALSE, conf_interval = 0.95, namesToExclude = "bl|Iteration|Likelihood|Posterior|Prior"){

  if (!windows){
    
    list_runs <- vector("list", length = 0)
    
    for (i in 1:length(runs)) {

      cont_param <- getInfo(runs, i, namesToExclude)
      list_runs[[i]] <- as.mcmc(cont_param)
    }
    
    psrf_diag_runs <- gelman.diag(list_runs, confidence = conf_interval, transform = FALSE, autoburnin = FALSE, multivariate = FALSE)

    
    return(psrf_diag_runs)
    
  }

  else{
    
    list_windows <- vector("list", length = 0)
    
    list_gelman <- vector("list", length = 0)
    
    for (i in 1:length(runs)) {

      cont_param <- getInfo(runs, i, splitWindows = TRUE, namesToExclude)
      list_windows[[1]] <- as.mcmc(cont_param[[1]])
      list_windows[[2]] <- as.mcmc(cont_param[[2]])

      psrf_windows <- gelman.diag(list_windows,confidence = conf_interval, transform = FALSE, autoburnin = FALSE, multivariate = FALSE)
      
      list_gelman[[i]] <- psrf_windows
      
    }
    
    return(list_gelman)
  }
  
}
