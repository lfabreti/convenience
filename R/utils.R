#' Utilities
#' 
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @export makeControl
#' @export ksThreshold

#' @useDynLib convenience
#' @importFrom Rcpp sourceCpp
NULL
#> NULL

# Calculates the relative difference of some stats for 2 dataframes
calcRelativeDiff <- function(dataframe1, dataframe2, stats){
  return( (abs (sapply(dataframe1,stats) - sapply(dataframe2, stats)) ) / ( (abs (sapply(dataframe1,stats) + sapply(dataframe2,stats)) ) / 2) )
}

# subtract one data frame from other
check <- function(df1, df2){
  
  df <- vector("double", length = 0)
  
  for (i in 1:length(df1)){
    df[[i]] <- df1[[i]] - df2[[i]]
  }
  
  return(df)
}

check.clades.freq <- function(runs, freq){
  
  cladenames_post <- vector("list", length = length(runs))
  for (z in 1:length(runs)) {
    
    x <- getInfo(runs, z, trees = TRUE)
    start <- 1
    end <- length(x)
    
    if(class(x) == "rwty.chain"){
      x <- x$trees
    }
    
    if (length(x) == 1 && class(x[[1]]) == "multiPhylo"){
      x <- x[[1]]
    }
    
    x <- x[start:end]
    
    clades <-  prop.part(x)
    
    cladefreqs <- as.numeric(as.character(attr(clades, which="number")[1:length(clades)] ))
    
    cladefreqs <- cladefreqs/length(x)
    
    tiplabels <- as.character(x[[1]]$tip.label)
    
    cladenames <- rep(NA, length(clades))
    
    for(i in 1:length(clades)){
      taxa_indices <- clades[[i]]
      taxa <- tiplabels[taxa_indices]
      sorted_taxa <- sort(taxa)
      cladenames[i] <- paste(sorted_taxa, collapse=" ")
    }
    aux <- vector()
    for (i in 1:length(cladenames)) {
      if(freq > 0.5) {
        if( cladefreqs[i] > freq ) aux <- c(aux, cladenames[i])
      } else{
        if( cladefreqs[i] < freq ) aux <- c(aux, cladenames[i])
      }
    }
    cladenames_post[[z]] <- aux
  }
  
  return(cladenames_post)
}

# Adjust the names of tips for split frequencies
clade.freq.named <- function (x, start, end, rooted=FALSE, ...) {
  
  if(class(x) == "rwty.chain"){
    x <- x$trees
  }
  
  if (length(x) == 1 && class(x[[1]]) == "multiPhylo"){
    x <- x[[1]]
  }
  
  x <- x[start:end]
  
  clades <-  prop.part(x)
  
  cladefreqs <- as.numeric(as.character(attr(clades, which="number")[1:length(clades)] ))
  
  cladefreqs <- cladefreqs/length(x)
  
  tiplabels <- as.character(x[[1]]$tip.label)
  
  cladenames <- rep(NA, length(clades))
  
  for(i in 1:length(clades)){
    taxa_indices <- clades[[i]]
    taxa <- tiplabels[taxa_indices]
    sorted_taxa <- sort(taxa)
    cladenames[i] <- paste(sorted_taxa, collapse=" ")
  }
  
  cladenames_post <- vector()
  cladefreqs_post <- vector()
  for (i in 1:length(cladenames)) {
    cladenames_post <- c(cladenames_post, cladenames[i])
    cladefreqs_post <- c(cladefreqs_post, cladefreqs[i])
  }
  
  clade.df <- data.frame(cladenames_post, cladefreqs_post)
  
  return(clade.df)
}

clade.freq.tree <- function (x, rooted=FALSE, ...) {
  
  clades <-  prop.part(x)
  
  cladefreqs <- as.numeric(as.character(attr(clades, which="number")[1:length(clades)] ))
  
  tiplabels <- as.character(x$tip.label)
  
  cladenames <- rep(NA, length(clades))
  
  for(i in 1:length(clades)){
    taxa_indices <- clades[[i]]
    taxa <- tiplabels[taxa_indices]
    sorted_taxa <- sort(taxa)
    cladenames[i] <- paste(sorted_taxa, collapse=" ")
  }
  
  clade.df <- data.frame(cladenames, cladefreqs)
  
  return(clade.df)
}

clade.freq.trees <- function (x, start, end, rooted=FALSE, ...) {
  
  if(class(x) == "rwty.chain"){
    x <- x$trees
  }
  
  if (length(x) == 1 && class(x[[1]]) == "multiPhylo"){
    x <- x[[1]]
  }
  
  x <- x[start:end]
  
  clades <-  prop.part(x)
  
  cladefreqs <- as.numeric(as.character(attr(clades, which="number")[1:length(clades)] ))
  
  clade_frequencies <- cladefreqs/length(x)
  
  tiplabels <- as.character(x[[1]]$tip.label)
  
  cladenames <- rep(NA, length(clades))
  
  for(i in 1:length(clades)){
    taxa_indices <- clades[[i]]
    taxa <- tiplabels[taxa_indices]
    sorted_taxa <- sort(taxa)
    cladenames[i] <- paste(sorted_taxa, collapse=" ")
  }
  
  cladenames_post <- vector()
  cladefreqs_post <- vector()
  for (i in 1:length(cladenames)) {
    if(  clade_frequencies[i] <= 0.975 & clade_frequencies[i] >= 0.025 ){
      cladenames_post <- c(cladenames_post, cladenames[i])
      cladefreqs_post <- c(cladefreqs_post, cladefreqs[i])
    }
  }
  
  clade.df <- data.frame(cladenames_post, cladefreqs_post)
  
  return(clade.df)
}

# Function to calculate ESS according to Tracer
essTracer <- function(input,stepSize = 1){
  
  samples <- length(input)
  
  max_lag <- 2000
  max_lag <- min( (samples-1) , max_lag)
  gammaSt <- replicate(max_lag, 0)
  varStat <- 0
  
  i = 1 
  
  while (i <= max_lag) {
    
    for (j in 1:(samples-i+1)) {
      del1 <- input[j] - mean(input)
      del2 <- input[j+i-1] -mean(input)
      gammaSt[i] <- gammaSt[i] + (del1*del2)
    }
    
    gammaSt[i] <- gammaSt[i] / (samples - i + 1)
    
    if(i == 1){
      varStat <- gammaSt[1]
      
    } else if( i %% 2 == 1){
      
      if( (gammaSt[i-1]+ gammaSt[i]) > 0 ){
        varStat <- varStat + 2 * (gammaSt[i-1]+gammaSt[i])
      }
      else{
        max_lag <- i
      }
      
    }
    i = i+1
  }
  
  act <- 0
  ess <- 0
  
  if (gammaSt[1] == 0){
    act <-0
  } else {
    act <- stepSize * varStat / gammaSt[1]
  }
  
  if (act ==0){
    ess <-1
  } else {
    ess <- (stepSize*samples) / act
  }
  
  return(ess)
}

# Calculates expected difference for splits frequencies
expectedDiffSplits <- function(ess){
  
  prob <- vector()
  threshold <- vector()
  for (p in 1:99/100) {
    exp_diff <- 0
    probs <- array(0,ess+1)
    
    for (f1 in 0:ess) {
      for (f2 in 0:ess) {
        exp_diff <- exp_diff + abs( (f1/ess) - (f2/ess) ) * dbinom(f1, size = ess, prob = p) * dbinom(f2, size = ess, prob = p)
        diff <- abs(f1-f2)
        probs[diff+1] <- probs[diff+1] + dbinom(f1, size = ess, prob = p) * dbinom(f2, size = ess, prob = p)
      }
    }
    prob <- c(prob, p)
    cdf_diff <- cumsum(probs)
    thresh <- (min(which( cdf_diff >= 0.95 ))-1) / ess
    threshold <- c(threshold, thresh)
  }
  return(rbind(prob,threshold))
}

# This function takes the name of a format and returns a list containing important info about file suffixes and whatnot
get.format <- function(format){
  
  # Default behavior for MrBayes files
  if(format == "mb"){
    return(list(
      trees.suffix = ".t",
      log.suffix = ".p",
      type = "nexus",
      skip = 1
    ))
  }
  
  # Default behavior for *BEAST files
  if(format == "*beast"){
    return(list(
      trees.suffix = ".species.trees",
      log.suffix = ".log",
      type = "nexus",
      skip = 2
    ))
  }
  
  # Default behavior for BEAST files
  if(format == "beast"){
    return(list(
      trees.suffix = ".trees",
      log.suffix = ".log",
      type = "nexus",
      skip = 2
    ))
  }
  
  # Default behavior for revbayes files
  if(format == "revbayes"){
    return(list(
      trees.suffix = ".trees",
      log.suffix = ".log",
      type = "revbayes",
      skip = 0
    ))
  }
  
  # Default behavior for pyhlobayes files
  if(format == "phylobayes"){
    return(list(
      trees.suffix = ".treelist",
      log.suffix = ".trace",
      type = "newick",
      skip = 0
    ))
  }
  
}

# Get the continuous parameters or the trees of a list of rwty.trees, returns a dataframe
getInfo <- function(all_runs, run, namesToExclude, trees = FALSE, splitWindows = FALSE){
  
  if (!trees){
    cont_param <- all_runs[[run]]["ptable"]
    cont_param <- as.data.frame(all_runs[[run]]["ptable"])
    names(cont_param) <- gsub("ptable.","",names(cont_param),fixed=TRUE)
    column <- grep(pattern = namesToExclude, names(cont_param), value = T)
    cont_param <- cont_param[-match(column, names(cont_param))]
    if( !splitWindows ){
      return(cont_param)
    }
    else{
      if( typeof(cont_param) == "list"){
        all_wind <- vector("list", length = 0)
        len_run <- length(cont_param[[1]])
        second <- (1*(0.2*len_run))
        third <- (2*(0.2*len_run))
        fourth <- (3*(0.2*len_run))
        fifth <- (4*(0.2*len_run))
        
        # gets the third window of the run
        all_wind[[1]] <- as.data.frame(cont_param[1:second,]) # check first and last  window
        #gets the fifth window of the run
        all_wind[[2]] <- as.data.frame(cont_param[fifth:len_run,])
        names(all_wind[[2]]) <- names(cont_param)
        return(all_wind)
      } 
      
      else{
        all_wind <- vector("list", length = 0)
        len_run <- length(cont_param)
        second <- (1*(0.2*len_run))
        third <- (2*(0.2*len_run))
        fourth <- (3*(0.2*len_run))
        fifth <- (4*(0.2*len_run))
        
        # gets the third window of the run
        all_wind[[1]] <- cont_param[1:second]
        #gets the fifth window of the run
        all_wind[[2]] <- cont_param[fifth:len_run]
        return(all_wind)
      }
    }
  }
  
  else{
    x <- all_runs[[run]]$trees
    if( !splitWindows ){
      return(x)
    }
    else{
      all_wind <- vector("list", length = 0)
      len_run <- length(x)
      second <- (1*(0.2*len_run))
      third <- (2*(0.2*len_run))
      fourth <- (3*(0.2*len_run))
      fifth <- (4*(0.2*len_run))
      
      #gets the third window of the run
      all_wind[[1]] <- x[1:second]
      #gets the fifth window of the run
      all_wind[[2]] <- x[fifth:len_run]
      return(all_wind)
    }
  }
}

# Calculates the threshold for the KS test
ksThreshold <- function(alpha, ess){
  c_alpha <- sqrt( -(log(alpha/2)) * 0.5) 
  return( (c_alpha*sqrt((2*ess)/(ess*ess))) )
}

# Function to create control argument for checkConvergence
makeControl <- function( tracer = NULL, burnin = NULL, precision = NULL, namesToExclude = NULL ){
  control <- vector(mode = "list", length = 4)
  names(control) <- c("tracer", "burnin", "precision", "namesToExclude")
  control$tracer <- tracer
  control$burnin <- burnin
  control$precision <- precision
  control$namesToExclude <- namesToExclude
  #names(control) <- c("tracer", "burnin", "precision", "namesToExclude")
  return(control)
}

# Calculates min ESS according to the std error of the mean
minESS <- function(per){
  return((1/(per*4))^2)
}

# Calculates the 2.5 and 97.5 quantiles of a dataframe
quants <- function(x){
  return(quantile(x,probs = c(0.025, 0.975)))
}

# From RWTY
isTree<-function(x) {
  !is.null(try(read.tree(text=x)))
}
  
# From RWTY
read.revbayestrees<-function(file) {
  filelines<-readLines(file)
  column.names<-strsplit(filelines[1], split="\t")[[1]]
  data<-strsplit(filelines[-1], split="\t")
  samplerow<-data[[1]]
  treecheck<-unlist(lapply(samplerow, FUN=isTree))
  param<-matrix(ncol=sum(!treecheck), nrow=length(data))
  colnames(param)<-column.names[!treecheck]
  for(i in 1:length(data))
    param[i,]<-as.numeric(data[[i]][!treecheck])
  tree<-list()
  for(i in 1:length(data)) {
    tree[[i]]<-read.tree(text=data[[i]][treecheck])
  }
  class(tree)<-"multiPhylo"
  return(list(tree=tree, param=param))
}

# Calculate standard error
se <- function(x){
  effSamSize <- ess(x)
  return(sd(x)/sqrt(effSamSize))
}
