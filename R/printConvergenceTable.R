#' Print method for convenience.table
#' 
#' @param x A list of convenience.diag type
#' 
#' @export

print.convenience.table <- function(x){
  print(summary(x))
}