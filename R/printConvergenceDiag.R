#' Print method for convenience.diag
#' 
#' @param x A list of convenience.diag type
#' @param ... (various) Additional arguments
#' 
#' @export

print.convenience.diag <- function(x, ...){
  print(x$message)
}
