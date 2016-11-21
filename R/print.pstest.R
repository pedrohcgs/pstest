#' @title Print
#'
#' @description Prints a pstest Object
#'
#' @param x A pstest object
#' @param ... Other params (required as generic function, but not used)
#'
# Define new print function
print.pstest <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  ks.mat <- data.frame(x$kstest, x$pvks)
  colnames(ks.mat) <- c("Test statistic", "Bootstrapped P-value")
  rownames(ks.mat) <- c("")

  cvm.mat <- data.frame(x$cvmtest, x$pvcvm)
  colnames(cvm.mat) <- c("Test statistic", "Bootstrapped P-value")
  rownames(cvm.mat) <- c("")

  cat("\n Kolmogorov-Smirnov test:\n")
  print(ks.mat)

  cat("\n Cramer-von Mises test:\n")
  print(cvm.mat)
}
