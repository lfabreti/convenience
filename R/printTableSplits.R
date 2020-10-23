#' Prints the frequencies and ESS of the splits
#' 
#' @param output A list of convenience.diag type
#' @param filename A filename to save the table, if NULL the table will be printed
#'  
#' @export

printTableSplits <- function(output, filename = NULL){
  
  if( ncol(output$tree_parameters$ess) == 1 ){
    names_spits <- output$tree_parameters$frequencies[1,]
    freq_splits <- as.numeric(output$tree_parameters$frequencies[2,])
    df_splits <- data.frame(row.names = names_spits, "frequencies" = freq_splits)
  }else{
    df_splits <- as.data.frame(output$tree_parameters$frequencies)
  }
  colnames(df_splits) <- "frequencies"
  df_splits$ESS <- rowSums(output$tree_parameters$ess)
  
  if( is.null(filename) ){
    return(df_splits)
  } else{
    write.csv(df_splits, file = filename)
  }
}