#' ESS calculation for the continuous parameters
#' 
#' Calculates the Effective Sample Size for the continuous parameters
#' 
#' @param runs A list of rwty.chain type
#' @param windows A boolean to set if the calculation is within runs or between runs
#' @param namesToExclude Column names to exclude from the calculations
#' @param tracer A boolean to determine if ESS should be calculated with Tracer method. If set to FALSE, ESS will be calculated with CODA
#' 
#' @return A data-frame with ESS values for the continuous parameters
#' 
#' @export


essContParam <- function(runs, windows=FALSE, namesToExclude, tracer) {

  if(!windows){
    vecEss <- vector("list", length = 0)
    for (i in 1:length(runs)) {
      cont_param <- getInfo(runs, i, namesToExclude)
      
      if ( typeof(cont_param) == "list" ){
        for (rows in 1:length(cont_param)) {
          
          if ( is.nan( effectiveSize(cont_param[[rows]]) ) ) {
            print(paste(names(cont_param[rows]), " ess is not a number!"))
          }
        }
      }
      
      if( tracer == T){
        vecEss <- c(vecEss, sapply(cont_param, essTracerC))
      } else{
        vecEss <- c(vecEss, effectiveSize(cont_param))
      } 
    }
    
    n_param <- length(vecEss)/length(runs)
    name_param <- names(vecEss[1:n_param])
    
    df_ess <- data.frame(matrix(unlist(vecEss), nrow = n_param, byrow = F), stringsAsFactors = F)
    rownames(df_ess) <- name_param
    
    for (i in 1:ncol(df_ess)) {
      colnames(df_ess)[i] <- paste("ESS_run_", i, sep = "")
      
    }
  }
  
  else{
    vecESS <- vector()
    
    for (i in 1:length(runs)) {
      cont_param <- getInfo(runs, i, splitWindows = TRUE, namesToExclude)
      
      if( tracer == T ){
        vecESS <- c(vecESS, sapply(cont_param[[1]], essTracerC))
        vecESS <- c(vecESS, sapply(cont_param[[2]], essTracerC))
      }else{
        vecESS <- c(vecESS, effectiveSize(cont_param[[1]]))
        vecESS <- c(vecESS, effectiveSize(cont_param[[2]]))
      }
    }
    
    n_param <- length(vecESS)/(length(runs)*2)
    name_param <- names(vecESS[1:n_param])
    
    df_ess <- data.frame(matrix(unlist(vecESS), nrow = length(runs)*2, byrow = T), stringsAsFactors = F)
    colnames(df_ess) <- name_param
    name_runs <- vector()
    for (i in 1:length(runs)) {
      name_runs <- c(name_runs, paste("Run_", i, "_window_1",sep = ""), paste("Run_", i, "window_2", sep = "") )
    }
    rownames(df_ess) <- name_runs
  }
  
  return(df_ess)
}
