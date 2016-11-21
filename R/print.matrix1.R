#' @title print.matrix1
#'
#' @description Helper function to print a matrix; used by the print methods
#'
#' @param m Some matrix
#'
#' @noRd
#' @importFrom utils write.table
#'

# Got this from Brant QTE package
print.matrix1 <- function(m){
  utils::write.table(format(m, justify="right", digits=2, nsmall=2),
                     row.names=F, col.names=F, quote=F, sep="\t")
}
