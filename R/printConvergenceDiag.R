#' ESS for split frequencies
#' 
#' explain
#' 
#' @param 
#' 
#' @return 
#' 
#' @example 
#' 
#' @export

print.convenience.diag <- function(x){
  if ( length(x)>3 ){
    print(x$message)
  }
  print(summary(x))
}
