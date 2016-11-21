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
    # Weight function used
    if(x$argu$w == "ind")       ww <- "Indicator function, w(q,u) = 1(q<=u)"
    if(x$argu$w == "exp")       ww <- "Exponential function, w(q,u) = exp(qu)"
    if(x$argu$w == "logistic")  ww <- "Logistic function, w(q,u) = 1/[1+exp(-qu)]"
    if(x$argu$w == "sin")       ww <- "Sine function, w(q,u) = sin(qu)"
    if(x$argu$w == "sincos")    ww <- "Sine and cosine function, w(q,u) = sin(qu)+ cos(qu)"

    #Creat parameters for the Table
    header <- c("", "Test statistic", "Bootstrapped P-value")
    body <- cbind(c("Kolmogorov-Smirnov", "Cramer-von Mises"),
                  c(round(x$kstest, digits = 4), round(x$cvmtest, digits = 4)),
                  c(round(x$pvks, digits = 4), round(x$pvcvm, digits = 4)))

    colnames(body) <- header
    #-----------------------------------------------------------------------------
    #Output

    cat(" Call:\n")
    print(x$call)

    cat("\n Sant'Anna and Song (2016) specification test for the propensity score.\n")

    cat("\n Weight function:", ww, "\n")
    cat("Number of Boostrap draws:", x$argu$nboot, "\n")
    cat("\n \n")

    print.matrix1(rbind(header, body))
}


#' @title print.matrix1
#'
#' @description Helper function to print a matrix; used by the print methods
#'
#' @param m Some matrix
#'
#' @noRd
#' @importFrom utils write.table
#'
print.matrix1 <- function(m){
  utils::write.table(format(m, justify="right", digits=2, nsmall=2),
              row.names=F, col.names=F, quote=F, sep="\t")
  ##print(m, print.gap=3, right=T)
}
