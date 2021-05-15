#' Load Files
#'
#' Load files from the output of a MCMC. The input can be log files, tree files or both
#' 
#' @param path The path to the folder with all the output
#' @param list_files List of files with the output
#' @param format The format of the phylogenetic output. Current supported formats are: "revbayes", "mb", "beast", "*beast", "phylobayes"
#' @param tree_name The name of the column containing the trees, default = "psi"
#' 
#' @return List of type rwty.chain
#' 
#' @export


loadFiles <- function( path = NULL, list_files = NULL, format, tree_name =  "psi") {
  
  skip = 0
  
  if(format == "revbayes"){
    log_ext = "\\.log$"
    tree_ext="\\.trees$"
    
  }else if(format == "mb"){
    skip = 1
    log_ext = "\\.p$"
    tree_ext="\\.t$"
    
  }else if(format == "beast"){
    skip = 2
    log_ext = "\\.log$"
    tree_ext="\\.trees$"
    
  }else if(format == "*beast"){
    skip = 2
    log_ext = "\\.log$"
    tree_ext="\\.species.trees$"
    
  }else if(format == "phylobayes"){
    log_ext = "\\.trace$"
    tree_ext="\\.treelist$"
    
  }else{
    stop("Provide format!")
  }
  
  options(warn = -1)
  
  output <- list()
  
  if ( !is.null(path) ){
    files <- list.files(path, recursive=F)
    #files <- files[ grepl("*run*|*joint*|*stone*", files) ]
    
    logFiles <- files[ grepl(log_ext, files) ]
    treeFiles <- files[ grepl(tree_ext, files) ]
  }else {
    
    # enforce argument matching
    character_paths_are_strings <- is.character(list_files)
    if ( any(character_paths_are_strings == FALSE) == TRUE ) {
      # print out the ones that are not character strings
      cat( "Some paths are not character strings:",
           paste0("\t",path[character_paths_are_strings == FALSE]), sep="\n")
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
        
        output[[i]] <- readTrace(logFiles[i], burnin = 0, skip = skip)
      }
      setwd("..")
    }
    # list of files is provided
    else{
      output <- vector("list", length(logFiles))
      
      for (i in 1:length(logFiles)){
        
        output[[i]] <- readTrace(logFiles[i], burnin = 0, skip = skip)
      }
    }
    
  }
  
  # Log files and tree files
  else{

    all_vecs <- vector("list", length(treeFiles))
    
    for(i in 1:length(treeFiles)){
      all_vecs[[i]] <- paste("run_", i, sep = "")
    }
    vec <- unlist(all_vecs)
    
    if ( !is.null(path)){
      output <- loadMulti(path , format = format, labels=vec, skip = skip)
    }else {
      output <- loadMulti( tree_files = treeFiles, log_files = logFiles, format = format, skip = skip, labels = vec)
    }
    
    
    # exclude combined files
    count_sizes <- vector("double", length=0)
    
    for (i in 1:length(output)){
      
      count_sizes <- c(count_sizes, length(output[[i]]$trees))
    }
    
    mean <- mean(count_sizes)
    
    output_tmp <- output
    
    for (i in 1:length(output)) {
      if( (length(output[[i]]$trees)) > mean*2 ){
        output_tmp <- output[-i]
      }
    }
    output <- output_tmp
    
    # exclude logs for when there is only tree files
    if (length(logFiles) == 0){
      
      for(i in 1:length(output)){
        
        output[[i]]$ptable <- vector("list", length = 0)
        
      }
    }
  }
 
  return(output)
  
}
