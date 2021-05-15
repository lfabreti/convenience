#' Make control
#'
#' Function to create control argument for checkConvergence
#' 
#' @param tracer To choose the ESS calculation using CODA, turn this parameter off with tracer = FALSE
#' @param burnin Sets a burn-in value to exclude the initial samples
#' @param precision Changes the precision for the threshold calculations
#' @param namesToExclude The names of parameters to exclude from the convergence assessment 
#' 
#' @return List
#' 
#' @export

makeControl <- function( tracer = NULL, burnin = NULL, precision = NULL, namesToExclude = NULL ){
  control <- vector(mode = "list", length = 4)
  names(control) <- c("tracer", "burnin", "precision", "namesToExclude")
  control$tracer <- tracer
  control$burnin <- burnin
  control$precision <- precision
  control$namesToExclude <- namesToExclude

  return(control)
}