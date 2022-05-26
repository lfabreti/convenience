#' Modified from RevGadgets
#' 
#' Read trace
#'
#' Reads in MCMC log files
#'
#' Reads in one or multiple MCMC log files from the same analysis
#' and discards a user-specified burn-in, compatible with multiple monitor types.
#' If the trace contains vectors of vectors and the user does not specify
#' format = "complex", readTrace() will read in those columns as factors
#' rather than as numeric vectors.
#' 
#' @importFrom utils read.table
#'
#' @param paths (vector of character strings; no default) File path(s) to trace file.
#' @param format (single character string; default = simple) Indicates type of
#' MCMC trace, complex indicates cases where trace contains vectors of vectors/
#' matrices - mnStochasticVariable monitor will sometimes be of this type.
#' @param delim (single character string) Delimiter of file.
#' @param burnin (single numeric value; default = 0.1) Fraction of generations to
#' discard (if value provided is between 0 and 1) or number of generations (if
#' value provided is greater than 1).
#' @param check.names (logical; default = FALSE) Passed to utils::read.table(); indicates
#' if utils::read.table() should check column names and replace syntactically invalid
#' characters.
#' @param skip Number of lines to skip before reading the data from the log file.
#' @param ... (various) Additional arguments passed to utils::read.table().
#'
#' @return List of dataframes (of length 1 if only 1 log file provided).
#'
#' @examples
#'
#' \dontrun{
#' file <- system.file("extdata",
#'     "sub_models/primates_cytb_GTR_mini.p", package="RevGadgets")
#' one_trace <- readTrace(paths = file, delim = "\t", skip = 0)
#' multi_trace <- readTrace(paths = c(file, file))
#' }
#' @export

readTrace <- function(paths, format = "simple",
                      delim, burnin = 0, check.names = FALSE, ...){

  # read output using RevGadgets  
  output <- RevGadgets::readTrace(paths, format, delim, burnin, check.names, ...)

  # convert to rwty.chain object  
  treelist <- NULL
  gens.per.tree <- NULL
  ptable <- output[[1]]
  output <- list("trees" = treelist, "ptable" = ptable, "gens.per.tree" = gens.per.tree)
  class(output) <- "rwty.chain"
  
  return(output)
}
