### Utils ###

abc <- function(x, percent){
  return(sd(x)*4*percent)
}

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
  
  if(!rooted){
    clades <- postprocess.prop.part(clades)
  }
  
  # recover()
  
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
  
  clade.df <- data.frame(cladenames, cladefreqs)
  
  return(clade.df)
}

# Get the continuous parameters or the trees of a list of rwty.trees, returns a dataframe
getInfo <- function(all_runs, run, namesToExclude, trees = FALSE, splitWindows = FALSE){
  
  if (!trees){
    cont_param <- all_runs[[run]]["ptable"]
    cont_param <- as.data.frame(all_runs[[run]]["ptable"])
    names(cont_param) <- gsub("ptable.","",names(cont_param),fixed=TRUE)
    column <- grepl(pattern = namesToExclude, names(cont_param))
    cont_param <- cont_param[,!column]
    if( !splitWindows ){
      return(cont_param)
    }
    else{
      all_wind <- vector("list", length = 0)
      len_run <- length(cont_param[[1]])
      fourth <- (3*(0.2*len_run))
      fifth <- (4*(0.2*len_run))
      all_wind[[1]] <- cont_param[fourth:fifth,]
      all_wind[[2]] <- cont_param[fifth:len_run,]
      return(all_wind)
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
      fourth <- (3*(0.2*len_run))
      fifth <- (4*(0.2*len_run))
      all_wind[[1]] <- x[fourth:fifth]
      all_wind[[2]] <- x[fifth:len_run]
      return(all_wind)
    }
  }
}

# Calculates the 2.5 and 97.5 quantiles of a dataframe
quants <- function(x){
  return(quantile(x,probs = c(0.025, 0.975)))
}

# Calculate standard error
se <- function(x){
  effSamSize <- ess(x)
  return(sd(x)/sqrt(effSamSize))
}