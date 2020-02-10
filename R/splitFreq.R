# Calculate split frequencies between windows of the same run or between runs

#' Split Frequencies
#' 
#' Calculate the split frequencies
#' 
#' @param runs A list of rwty.chain type
#' @param windows
#' 
#' @return 
#' 
#' @example 
#' 
#' @export

splitFreq <- function(runs, windows=FALSE){
  
  vecSplits <- vector(length = 0)
  
  if(!windows){
    
    all_df <- vector("list", length = 0)
    for (i in 1:length(runs)){
      
      x <- getInfo(runs, i, trees = TRUE)
      cladefreqs <- clade.freq.named(x, start = 1, end = length(x))
      all_df[[i]] <- cladefreqs
    }
    
    for (r1 in 1:(length(all_df)-1)){
      #print(r1)
      for (z in 1:length(all_df[[r1]]$cladenames)){
        for (j in 1:length(all_df[[r1+1]]$cladenames)){
          if( as.character(all_df[[r1]]$cladenames[z]) == as.character(all_df[[r1+1]]$cladenames[j]) ){
            vecSplits <- c ( vecSplits, ( abs (all_df[[r1]]$cladefreqs[z] - all_df[[r1+1]]$cladefreqs[j] ) ) )
          }
        }
      }
    }
    
  }
  
  else{
    
    for (i in 1:length(runs)){
      
      x <- getInfo(runs, i, trees=TRUE, splitWindows = TRUE)
      
      compar_1 <- clade.freq.named(x[[1]], start = 1, end = length(x[[1]]))
      compar_2 <- clade.freq.named(x[[2]], start = 1, end = length(x[[2]]))
      
      for (z in 1:length(compar_1[[1]])) {
        for (j in 1:length(compar_2[[1]])) {
          if( as.character(compar_1$cladenames[z]) == as.character(compar_2$cladenames[j]) ){
            vecSplits <- c(vecSplits, ( abs(compar_1$cladefreqs[z] - compar_2$cladefreqs[j])) )
          }
        }
      }
    }
  }
  
  return(vecSplits)
  
}

