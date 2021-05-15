#' Prints the frequencies and ESS of the splits
#' 
#' @importFrom utils write.csv
#' 
#' @param output A list of convenience.diag type
#' @param splits_per_run When set to TRUE, shows the information for the splits for each run
#' @param filename A filename to save the table, if NULL the table will be printed
#'  
#' @export

printTableSplits <- function(output, splits_per_run = FALSE,filename = NULL){
  
  if(splits_per_run == TRUE){
    df_splits <- output$tree_parameters$freq_per_run
    
  } else{
    if( length(output$tree_parameters$frequencies) == 1 ){
      names_spits <- names(output$tree_parameters$frequencies[[1]])
      freq_splits <- as.vector(unlist(output$tree_parameters$frequencies[[1]]))
      df_splits <- data.frame(row.names = names_spits, "frequencies" = freq_splits)
      df_splits <- as.data.frame(df_splits[row.names(output$tree_parameters$ess),])
      row.names(df_splits) <- row.names(output$tree_parameters$ess)
    }else{
      df_3 <- plyr::ldply(output$tree_parameters$frequencies, rbind)
      df_3 <- setNames(data.frame(t(df_3[,-1]), row.names = colnames(df_3)[-1]), df_3[,1])
      
      df_3 <- (rowSums(as.data.frame(df_3))*2)/(length(output$tree_parameters$ess)-1)
      df_3 <- (df_3)/4
      df_splits <- as.data.frame(df_3)
      df_splits <- as.data.frame(df_splits[row.names(output$tree_parameters$ess),])
      row.names(df_splits) <- row.names(output$tree_parameters$ess)
    }
    colnames(df_splits) <- "frequencies"
    df_splits$ESS <- rowSums(output$tree_parameters$ess)
    df_splits <- df_splits[order(df_splits$frequencies),]
  }
  
  if( is.null(filename) ){
    return(df_splits)
  } else{
    write.csv(df_splits, file = filename)
  }
}