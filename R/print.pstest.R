#' @title Print
#'
#' @description Prints a pstest Object
#'
#' @param x A pstest object
#' @param ... Other params (required as generic function, but not used)
#'
#' @export
#' @noRd
# Define new print function
print.pstest <- function(x, ...){
    #-----------------------------------------------------------------------------
    # Preliminaries
    # Put test and pvalues in matrix
    ks.mat <- data.frame(x$kstest, x$pvks)
    colnames(ks.mat) <- c("Test statistic", "Bootstrapped P-value")
    rownames(ks.mat) <- c("")

    cvm.mat <- data.frame(x$cvmtest, x$pvcvm)
    colnames(cvm.mat) <- c("Test statistic", "Bootstrapped P-value")
    rownames(cvm.mat) <- c("")
    #-----------------------------------------------------------------------------
    # Weight function used
    if(x$argu$w == "ind")       ww <- "Indicator function, w(q,u) = 1(q<=u)"
    if(x$argu$w == "exp")       ww <- "Exponential function, w(q,u) = exp(qu)"
    if(x$argu$w == "logistic")  ww <- "Logistic function, w(q,u) = 1/[1+exp(-qu)]"
    if(x$argu$w == "sin")       ww <- "Sine function, w(q,u) = sin(qu)"
    if(x$argu$w == "sincos")    ww <- "Sine and cosine function, w(q,u) = sin(qu)+ cos(qu)"

    cat(" Call:\n")
    print(x$call)

    cat("\n Sant'Anna and Song (2016) specification test for the propensity score.\n")

    cat("\n Weight function:", ww, "\n")
    cat("\n Number of Boostrap draws:", x$argu$nboot, "\n")

    header <- c("", "Test statistic", "Bootstrapped P-value")
    body <- cbind(c("Kolmogorov-Smirnov", "Cramer-von Mises"),
                  c(x$kstest, x$cvmtest), c(x$pvks, x$pvcvm))
    colnames(body) <- header
    print.matrix1(rbind(header, body))

    #cat("\n Kolmogorov-Smirnov test:\n")
    #print(ks.mat)

    #cat("\n Cramer-von Mises test:\n")
    #print(cvm.mat)
}


#' @title print.matrix1
#'
#' @description Helper function to print a matrix; used by the print methods
#'
#' @param m Some matrix
#'
#' @noRd
#'@importFrom utils write.table
#'
print.matrix1 <- function(m){
  utils::write.table(format(m, justify="right", digits=2, nsmall=2),
              row.names=F, col.names=F, quote=F, sep="\t")
  ##print(m, print.gap=3, right=T)
}
