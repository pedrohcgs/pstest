#' @title Summary
#'
#' @description Summary of a pstest object
#'
#' @param object A pstest object
#' @param ... Other params (required as generic function, but not used)
#'
#' @export
#' @noRd
# Define new summary function
summary.pstest <- function(object, ...){
  pstest.obj <- object
  print(pstest.obj)

}
