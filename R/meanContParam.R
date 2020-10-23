#' Calculates the means of the continuous parameters#' 
#' 
#' @param runs A list of rwty.chain type
#' @param namesToExclude Column names to exclude from the calculations
#' 
#' @return A vector with the means
#' 
#' @export

meanContParam <- function(runs, namesToExclude){
  all_df <- vector("list", length = 0)
  
  for (i in 1:length(runs)) {
    #get the cont_param for each run
    all_df[[i]] <- getInfo(runs, i, namesToExclude)
  }
  names_parameters <- names(all_df[[1]])
  
  
  if ( length(runs) >  1 ){
    df_aux <- data.frame()
    for (i in 1:length(runs)) {
      df_aux <- rbind(df_aux, all_df[[i]])
    }
    all_df <- df_aux
  }
  
  all_df <- as.data.frame(all_df)
  vecMeans <- sapply(all_df, mean)
  
  return(vecMeans)
}