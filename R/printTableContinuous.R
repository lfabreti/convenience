#' Prints the means and the ESS of the continuous parameters
#' 
#' @param output A list of convenience.diag type
#' @param filename A filename to save the table, if NULL the table will be printed
#'  
#' @export

printTableContinuous <- function(output, filename = NULL){
  
  df_continuous <- as.data.frame(output$continuous_parameters$means)
  colnames(df_continuous) <- "means"
  df_continuous$ESS <- rowSums(output$continuous_parameters$ess)
  
  if( is.null(filename) ){
    return(df_continuous)
  } else{
    write.csv(df_continuous, file = filename)
  }
}