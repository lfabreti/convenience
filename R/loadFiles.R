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


loadFiles <- function( path = NULL, list_files = NULL, format = "revbayes", burnin, tree_name =  "psi", log_ext = "*.log", tree_ext="*.trees") {
  
  options(warn = -1)
  
  output <- list()
  
  if ( !is.null(path) ){
    files <- list.files(path, recursive=F)
    files <- files[ grepl("*run*|*joint*|*stone*", files) ]
    
    logFiles <- files[ grepl(log_ext, files) ]
    treeFiles <- files[ grepl(tree_ext, files) ]
  }else {
    
    # enforce argument matching
    character_paths_are_strings <- is.character(list_files)
    if ( any(character_paths_are_strings == FALSE) == TRUE ) {
      # print out the ones that are not character strings
      cat( "Some paths are not character strings:",
           paste0("\t",paths[character_paths_are_strings == FALSE]), sep="\n")
      stop()
    }
    
    do_files_exist <- file.exists(list_files)
    if ( any(do_files_exist == FALSE) == TRUE ) {
      # print out paths to files that don't exist
      cat( "Some files do not exist:",
           paste0("\t",list_files[do_files_exist == FALSE]), sep="\n")
      stop()
    }
    
    logFiles <- list_files[ grepl(log_ext, list_files) ]
    treeFiles <- list_files[ grepl(tree_ext, list_files) ]
    
  }

  # No files
  if ( length(logFiles) == 0 & length(treeFiles) == 0 ){
    
    stop("No files to read")
  }
  
  # Only log files
  else if ( length(logFiles) > 0 & length(treeFiles) == 0 ){
    
    # path is provided
    if(!is.null(path)){
      setwd(path)
      #output <- list()
      
      for (i in 1:length(logFiles)){
        
        output[[i]] <- readTrace(logFiles[i], burnin = burnin)
      }
      setwd("..")
    }
    # list of files is provided
    else{
      output <- list()
      
      for (i in 1:length(logFiles)){
        
        output[[i]] <- readTrace(logFiles[i], burnin = burnin)
      }
    }
    
  }
  
  # Log files and tree files
  else{

    all_vecs <- vector("list", length = 0)
    
    for(i in 1:length(treeFiles)){
      all_vecs[[i]] <- paste("run_", i, sep = "")
    }
    vec <- unlist(all_vecs)
    
    if ( !is.null(path)){
      output <- loadMulti(path , format = format, labels=vec)
    }else {
      output <- loadMulti( tree_files = treeFiles, log_files = logFiles, format = format, labels = vec)
    }
    
    
    # exclude combined files
    count_sizes <- vector("double", length=0)
    
    for (i in 1:length(output)){
      
      count_sizes <- c(count_sizes, length(output[[i]]$trees))
    }
    
    mean <- mean(count_sizes)
    
    output_tmp <- output
    
    for (i in 1:length(output)) {
      if( (length(output[[i]]$trees)) > mean ){
        output_tmp <- output[-i]
      }
    }
    output <- output_tmp
    
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
    
    # exclude logs for when there is only tree files
    if (length(logFiles) == 0){
      
      for(i in 1:length(output)){
        
        output[[i]]$ptable <- vector("list", length = 0)
        
      }
    }
  }
  
  
  if ( length(output[[1]]$ptable) > 0 ){
    output_tmp <- output
    # Exclude continuous parameters that are fixed
    for (i in 1:length(output)) {
      for (j in 1:length(output[[i]]$ptable)) {
        if( var(output[[i]]$ptable[j]) == 0 ){
          #output[[i]]$ptable[j] <- NULL
          output_tmp[[i]]$ptable[j] <- NULL
        }
      }
    }
    output <- output_tmp
  }

  
  return(output)
  
}
