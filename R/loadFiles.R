#' Load Files
#'
#' Load files from your output, they can be log files, tree files or both
#' 
#' @param path The path to the folder with all the output
#' @param format The format of the output, default = "revbayes"
#' @param burnin Fraction of generations to discard, default = 0.1
#' @param tree_name The name of the column containing the trees, default = "psi"
#' @param log_ext The extension of your log files, default = "*.log"
#' @param tree_ext The extension of your tree files, default = "*.trees"
#' 
#' @return List of type rwty.chain
#' 
#' @example 
#' # loadFiles(path, format="revbayes", burnin = 0.1, tree_name =  "psi", log_ext = "*.log", tree_ext="*.trees")
#' 
#' @export


loadFiles <- function(path,format="revbayes", burnin = 0.1, tree_name =  "psi", log_ext = "*.log", tree_ext="*.trees"){
  
  output <- list()
  files <- list.files(path, recursive=F)
  files <- files[ grepl("*run*", files) ]
  
  logFiles <- files[ grepl(log_ext, files) ]
  treeFiles <- files[ grepl(tree_ext, files) ]
  
  if ( length(logFiles) == 0 & length(treeFiles) == 0 ){
    
    stop("No files to read")
  }
  
  else if ( length(logFiles) > 0 & length(treeFiles) == 0 ){
    
    output <- list()
    for (i in 1:length(logFiles)){
      
      output[[i]] <- readTrace(logFiles[i], burnin = burnin)
    }
    
  }
  
  else{
    
    all_vecs <- vector("list", length = 0)
    files <- list.files(path, pattern = "*trees")
    
    for(i in 1:length(files)){
      all_vecs[[i]] <- paste("run_", i)
    }
    vec <- unlist(all_vecs)
    
    output <- load.multi(path , format = format, labels=vec) # add burnin
    
    count_sizes <- vector("double", length=0)
    
    for (i in 1:length(output)){
      
      count_sizes <- c(count_sizes, length(output[[i]]$trees))
    }
    
    mean <- mean(count_sizes)
    
    for (i in 1:length(output)) {
      if( (length(output[[i]]$trees)) > mean ){
        output <- output[-i]
      }
      
    }
    
    #burnin
    if (burnin >= length(output[[1]]$trees)) stop("Burnin larger than iterations in file")
    
    for (i in 1:length(output)) {
      
      if (burnin >= 1) {
        
        output[[i]]$trees <- output[[i]]$trees[(burnin+1):(length(output[[i]]$trees))]
        output[[i]]$ptable <- output[[i]]$ptable[(burnin+1):(nrow(output[[i]]$ptable)),]
        
      } else if (burnin < 1 & burnin > 0) { #not working
        
        discard <- ceiling(burnin*(length(output[[i]]$trees)))
        output[[i]]$trees <- output[[i]]$trees[(discard+1):(length(output[[i]]$trees))]
        output[[i]]$ptable <- output[[i]]$ptable[(discard+1):(nrow(output[[i]]$ptable)),]
        
      } else if (burnin == 0) {
        
        output[[i]] <- output[[1]]
        
      } else {
        
        stop("What have you done?")
      }
    }
  }
  return(output)
}

