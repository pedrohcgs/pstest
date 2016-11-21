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
    if(x$argu$w == "ind")       ww <- "Indicator function: w(q,u) = 1(q<=u)"
    if(x$argu$w == "exp")       ww <- "Exponential function: w(q,u) = exp(qu)"
    if(x$argu$w == "logistic")  ww <- "Logistic function: w(q,u) = 1/[1+exp(-qu)]"
    if(x$argu$w == "sin")       ww <- "Sine function: w(q,u) = sin(qu)"
    if(x$argu$w == "sincos")    ww <- "Sine and cosine function: w(q,u) = sin(qu)+ cos(qu)"

    #Creat parameters for the Table
    header <- c("", "Test statistic", "Bootstrapped P-value")
    body <- cbind(c("Kolmogorov-Smirnov", "Cramer-von Mises"),
                  c(round(x$kstest, digits = 4), round(x$cvmtest, digits = 4)),
                  c(round(x$pvks, digits = 4), round(x$pvcvm, digits = 4)))

    colnames(body) <- header
    #-----------------------------------------------------------------------------
    #Output
    cat(" Call:\n")
    cat(" "); print(x$call)
    cat("-------------------------------------------------------------------------")
    cat("\n Sant'Anna and Song (2016) specification tests for the propensity score:\n")
    print.matrix1(rbind(header, body))
    cat("-------------------------------------------------------------------------")
    cat("\n Weight function:", "\t \t", ww)
    cat("\n Number of Boostrap draws:", "\t", x$argu$nboot)
    cat("\n Boostrap draws from", x$argu$dist, 'distribution')

}
