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
  
  vecSplits <- vector(length = 0)
  
  if(!windows){
    
    all_df <- vector("list", length = 0)
    for (i in 1:length(runs)){
      
      x <- getInfo(runs, i, trees = TRUE)
      cladefreqs <- clade.freq.named(x, start = 1, end = length(x))
      all_df[[i]] <- cladefreqs
    }
    
    # initiate some vectors to plot splits against runs
    vec1 <- vector()
    vec2 <- vector()
    
    for (r1 in 1:(length(all_df)-1)){
      #print(r1)
      for (z in 1:length(all_df[[r1]]$cladenames)){
        for (j in 1:length(all_df[[r1+1]]$cladenames)){
          if( as.character(all_df[[r1]]$cladenames[z]) == as.character(all_df[[r1+1]]$cladenames[j]) ){
            #vec1 <- c(vec1, all_df[[r1]]$cladefreqs[z])
            #vec2 <- c(vec2, all_df[[r1+1]]$cladefreqs[j])
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
          sum_freq <- 0
          exp_diff <- 0 
          obs_diff <-0
          if( as.character(compar_1$cladenames[z]) == as.character(compar_2$cladenames[j]) ){

            obs_diff <- ( abs(compar_1$cladefreqs[z] - compar_2$cladefreqs[j]))
            
            sum_freq <- ( (compar_1$cladefreqs[z]*length(x[[1]])) + (compar_2$cladefreqs[j]*length(x[[2]])) ) / ( length(x[[1]])+length(x[[2]]) )
            
            for (i in 1:ESS){
	            exp_diff <- exp_diff + abs(i/ESS-sum_freq)* dbinom(i, size = ESS, prob = sum_freq)
            }

            vecSplits <- c(vecSplits, (exp_diff-obs_diff))
          }
        }
      }
    }
  }
  
  return(vecSplits)
  
}

