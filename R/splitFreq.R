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

splitFreq <- function(runs, windows=FALSE, ESS = 100){
  
  if(!windows){
    
    all_df <- vector("list", length = 0)
    
    for (i in 1:length(runs)){
      
      x <- getInfo(runs, i, trees = TRUE)
      cladefreqs <- clade.freq.named(x, start = 1, end = length(x))
      all_df[[i]] <- cladefreqs
    }
    
    listSplits <- list()
    listExpectedDiff <- list()
    count <- 0
    
    for (r1 in 1:(length(all_df)-1)){
      for (r2 in (r1+1):length(all_df)) {
        
        vecSplits <- vector()
        vecExpectedDiff <- vector()
        
        for (z in 1:length(all_df[[r1]]$cladenames)){
          
          for (j in 1:length(all_df[[r2]]$cladenames)){
            
            if( as.character(all_df[[r1]]$cladenames[z]) == as.character(all_df[[r2]]$cladenames[j]) ){
              
              exp_diff <- 0
              
              vecSplits <- c ( vecSplits, ( abs (all_df[[r1]]$cladefreqs[z] - all_df[[r2]]$cladefreqs[j] ) ) )
              
              for (f1 in 0:ESS) {
                for (f2 in 0:ESS) {
                  
                  exp_diff <- exp_diff + abs( (f1/ESS) - (f2/ESS) ) * dbinom(f1,size=ESS,prob=all_df[[r1]]$cladefreqs[z]) * dbinom(f2,size=ESS,prob=all_df[[r2]]$cladefreqs[j])
                }
              }
              vecExpectedDiff <- c(vecExpectedDiff, exp_diff)
              
            }
          }
        }
        count <- count+1
        listSplits[[count]] <- vecSplits
        listExpectedDiff[[count]] <- vecExpectedDiff
      }
      
    }
    
  }
  
  else{
    
    listSplits <- list()
    listExpectedDiff <- list()
    count <- 0
    
    for (i in 1:length(runs)){
      
      vecSplits <- vector()
      vecExpectedDiff <- vector()
      
      x <- getInfo(runs, i, trees=TRUE, splitWindows = TRUE)
      
      compar_1 <- clade.freq.named(x[[1]], start = 1, end = length(x[[1]]))
      compar_2 <- clade.freq.named(x[[2]], start = 1, end = length(x[[2]]))
      
      for (z in 1:length(compar_1[[1]])) {
        for (j in 1:length(compar_2[[1]])) {

          if( as.character(compar_1$cladenames[z]) == as.character(compar_2$cladenames[j]) ){
            
            exp_diff <- 0
            
            vecSplits <- c ( vecSplits, ( abs (compar_1$cladefreqs[z] - compar_2$cladefreqs[j] ) ) )
            
            for (f1 in 0:ESS) {
              for (f2 in 0:ESS) {
                
                exp_diff <- exp_diff + abs( (f1/ESS) - (f2/ESS) ) * dbinom(f1,size=ESS,prob=compar_1$cladefreqs[z]) * dbinom(f2,size=ESS,prob=compar_2$cladefreqs[j])
              }
            }
            vecExpectedDiff <- c(vecExpectedDiff, exp_diff)
            
          }
        }
      }
      count <- count+1
      listSplits[[count]] <- vecSplits
      listExpectedDiff[[count]] <- vecExpectedDiff
      
    }
  }
  
  listAll <- rbind(listSplits,listExpectedDiff)
  return(listAll)
  
}

