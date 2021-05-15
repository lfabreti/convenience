#' Modified from RWTY
#'
#' Loads trees, looks for a log file of tree likelihoods and parameter values, returns an rwty.chain object containing both
#' 
#' @param file The tree file 
#' @param type The type of tree file
#' @param format The format of the phylogenetic output. Current supported formats are: "revbayes", "mb", "beast", "*beast", "phylobayes"
#' @param gens.per.tree Number of generations between sampled trees
#' @param trim Increment for the sequence of trees to read
#' @param logfile The log files from the MCMC
#' @param skip Number of lines to skip before reading the data from the log file
#'
#' @return output A list of rwty.chain objects containing the multiPhylos and the tables of values from the log files if available
#'
#' @export

loadTrees <- function(file, type=NA, format = "mb", gens.per.tree=NA, trim=1, logfile=NA, skip=NA){
  
  format <- tolower(format)
  format_choices <- c("mb", "beast", "*beast", "revbayes", "mrbayes", "phylobayes")
  format <- match.arg(format, format_choices)
  if(format=="mrbayes") format="mb"
  if(!is.na(type)){
    type <- tolower(type)
    type_choices <- c("nexus", "newick")
    type <- match.arg(type, type_choices)
  }
  
  file.format <- get.format(format)
  
  if(is.na(type)){
    type <- file.format$type
  }
  
  if(is.na(skip)){
    skip <- file.format$skip
  }
  
  
  # Read in trees
  print("Reading trees...")
  if(type == "nexus") {
    treelist <- read.nexus(file=file)
  } else if(type=="revbayes") {
    tmp <- read.revbayestrees(file=file)
    treelist <- tmp$tree
    rb_ptable <- tmp$param
  } else {
    treelist <- read.tree(file=file)
  }
  
  treelist <- treelist[seq(from=1, to=length(treelist), by=trim)]
  
  if(is.na(gens.per.tree)){
    if(type=="revbayes") {
      gens.per.tree <- rb_ptable[2,"Iteration"] - rb_ptable[1,"Iteration"]
    } else {
      #   "beast" | "*beast" | "mb" | "phylobayes" 
      if(!is.null(names(treelist))){
        gens.per.tree <- as.numeric(utils::tail(strsplit(x=names(treelist)[3], split="[[:punct:]]")[[1]], 1)) -
          as.numeric(utils::tail(strsplit(x=names(treelist)[2], split="[[:punct:]]")[[1]], 1))
      }
      else gens.per.tree <- 1
    }
  }
  
  print(paste(gens.per.tree, "generations per tree..."))
  
  # Unroot all trees.  Can't use lapply because it was
  # deleting tip labels.
  if(is.rooted(treelist[[1]])){
    print("Unrooting, this may take a while...")
    treelist <- unroot(treelist)
    #for(i in 1:length(treelist)){
    #  treelist[[i]] <- unroot(treelist[[i]])
    #}
  }
  else{print("Trees are unrooted...")}
  
  
  # Reset class
  class(treelist) <- "multiPhylo"
  treelist <- ape::.compressTipLabel(treelist) # saves memory
  
  
  ptable <- NULL
  
  # Check whether log file path has been supplied and doesn't exist
  if(!is.na(logfile) && !file.exists(logfile)){
    stop(paste("Logfile not found at", logfile))
  }
  
  # logfile path has been supplied and file exists
  if(!is.na(logfile) && file.exists(logfile)){
    print(paste("Reading parameter values from", basename(logfile)))
    ptable <- read.table(logfile, skip=skip, header=TRUE, comment.char = "")
    ptable <- ptable[seq(from=1, to=length(ptable[,1]), by=trim),]
  }
  
  # If logfile hasn't been supplied try to find it by searching
  if(is.na(logfile)){
    
    if(grepl(paste(file.format$trees.suffix, "$"), file)){
      logfile <- sub(pattern = paste0(file.format$trees.suffix, "$"), file.format$log.suffix, file)
    }
    
    if(!is.na(logfile)){
      if(file.exists(logfile)){
        print(paste("Reading parameter values from", basename(logfile)))
        ptable <- read.table(logfile, skip=skip, header=TRUE)
        ptable <- ptable[seq(from=1, to=length(ptable[,1]), by=trim),]
      } else {
        print(paste("Couldn't find", basename(logfile)))
      }
    }
  }
  
  # add any columns from rb_ptable (from treefile) that are not already in ptable (from log)
  if(format=="revbayes") {
    to_add<-!(colnames(rb_ptable) %in% colnames(ptable))
    if(sum(to_add)>0)
      ptable<-cbind(rb_ptable[,to_add], ptable)
    print("rerooting trees...")
    outgroup_taxon = sort(treelist[[1]]$tip.label)[1]
    print(paste("Outgroup",outgroup_taxon))
    treelist <- root(treelist, outgroup=outgroup_taxon)
    treelist <- unroot(treelist)
    #for(i in 1:length(treelist)){
    #  tmp = root(treelist[[i]],outgroup=outgroup_taxon)
    #  treelist[[i]] <- unroot( tmp )
    #}
  }
  
  output <- list(
    "trees" = treelist,
    "ptable" = ptable,
    "gens.per.tree" = gens.per.tree)
  
  class(output) <- "rwty.chain"
  
  output
}