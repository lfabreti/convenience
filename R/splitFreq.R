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
  
  if(!windows){
    
    all_df <- vector("list", length = 0)
    
    for (i in 1:length(runs)){
      
      x <- getInfo(runs, i, trees = TRUE)
      cladefreqs <- clade.freq.named(x, start = 1, end = length(x))
      all_df[[i]] <- cladefreqs
    }
    
    listSplits <- list()
    listFrequencies <- list()
    count <- 0
    
    for (r1 in 1:(length(all_df)-1)){
      for (r2 in (r1+1):length(all_df)) {
        
        vecSplits <- vector()
        vecNames <- vector()
        vecFreqs <- vector()
        
        for (z in 1:length(all_df[[r1]]$cladenames_post)){
          
          for (j in 1:length(all_df[[r2]]$cladenames_post)){
            
            if( as.character(all_df[[r1]]$cladenames_post[z]) == as.character(all_df[[r2]]$cladenames_post[j]) ){
              
              exp_diff <- 0
              
              vecSplits <- c ( vecSplits, ( abs (all_df[[r1]]$cladefreqs_post[z] - all_df[[r2]]$cladefreqs_post[j] ) ) )
              
              vecNames <- c( vecNames, as.character(all_df[[r1]]$cladenames_post[z]) )
              
              vecFreqs <- c( vecFreqs, ((all_df[[r1]]$cladefreqs_post[z] + all_df[[r2]]$cladefreqs_post[j])/2))
            }
          }
        }
        
        count <- count+1
        listSplits[[count]] <- vecSplits
        names(listSplits[[count]]) <- vecNames
        listFrequencies[[count]] <- vecFreqs
        names(listFrequencies[[count]]) <- vecNames
      }
      
    }
    listAll <- rbind(listSplits,listFrequencies)
  }
  
  else{
    
    listSplits <- list()
    listFrequencies <- list()
    count <- 0
    
    for (i in 1:length(runs)){
      
      vecSplits <- vector()
      vecFreqs <- vector()
      vecNames <- vector()
      
      x <- getInfo(runs, i, trees=TRUE, splitWindows = TRUE)
      
      compar_1 <- clade.freq.named(x[[1]], start = 1, end = length(x[[1]]))
      compar_2 <- clade.freq.named(x[[2]], start = 1, end = length(x[[2]]))
      
      for (z in 1:length(compar_1[[1]])) {
        for (j in 1:length(compar_2[[1]])) {
          
          if( as.character(compar_1$cladenames[z]) == as.character(compar_2$cladenames[j]) ){
            
            exp_diff <- 0
            
            vecSplits <- c ( vecSplits, ( abs (compar_1$cladefreqs[z] - compar_2$cladefreqs[j] ) ) )
            
            vecFreqs <- c( vecFreqs, ((compar_1$cladefreqs[z] + compar_2$cladefreqs[j])/2))
            
            vecNames <- c( vecNames, as.character(compar_1$cladenames[z]))
          }
        }
      }
      
      count <- count+1
      listSplits[[count]] <- vecSplits
      names(listSplits[[count]]) <- vecNames
      listFrequencies[[count]] <- vecFreqs
      names(listFrequencies[[count]]) <- vecNames
      
    }
    listAll <- rbind(listSplits,listFrequencies)
  }
  
  return(listAll)
  
}

