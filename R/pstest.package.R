#' pstest: An R Package for assessing the goodness-of-fit of parametric propensity score models.
#'
#'@description
#' The propensity score is one of the most widely used tools in studying the causal effect
#' of a treatment, intervention, or policy. Given that the propensity score is usually unknown,
#' it has to be estimated, implying that the reliability of many treatment effect estimators depends
#' on the correct specification of the (parametric) propensity score. This package provides
#' data-driven nonparametric diagnostic tools for detecting propensity score misspecification.
#'
#'
#'@details
#' This R package implements the class of specification test for the propensity score
#' proposed in Sant'Anna and Song (2016), 'Specification Tests for the Propensity Score',
#' available at Pedro H.C. Sant'Anna webpage, \url{http://sites.google.com/site/pedrohcsantanna/}.
#'
#' In short, this package implements Kolmogorov-Smirnov and Cramer-von Mises type tests
#' for parametric propensity score models with either logistic ('logit'), or
#' normal ('probit') link function. Critical values are computed with the assistance of a
#' multiplier bootstrap.
#'
#' The tests are based on the integrated conditional moment approach, where the weight function
#' used is based on an orthogonal projection onto the tangent space of nuisance parameters.
#' As a result, the tests (a) enjoy improved power properties, (b) do not suffer from the
#' 'curse of dimensionality' when the vector of covariates is of high-dimensionality,
#' (c) are fully data-driven, (e) do not require tuning parameters such as bandwidths, and
#' (e) are able to detect a broad class of local alternatives converging to the null at the
#' parametric rate. These appealing features highlight that the tests can be of great use
#' in practice.
#'
#' It is worth stressing that this package implements in a unified manner a large class of
#' specification tests, depending on the chosen weight function \eqn{w(q,u)}:
#' \itemize{
#' \item{'ind' - the indicator weight function \eqn{w(q,u)=1(q \le u)}. This is the default.}
#' \item{'exp' - the exponential weight function \eqn{w(q,u)=exp(qu)}.}
#' \item{'logistic' - the logistic weight function \eqn{w(q,u)=1/[1+exp(-qu)]}.}
#' \item{'sin' - the sine weight function \eqn{w(q,u)=sin(qu)}.}
#' \item{'sincos' - the sine and cosine weight function \eqn{w(q,u)=sin(qu)+cos(qu)}.}
#' }
#'
#'Different weight functions \eqn{w(q,u)} have different power properties, and therefore,
#'being able to choose different \eqn{w(q,u)} gives us flexibility to direct power in desired
#'directions.
#'
#' @references
#'       Sant'Anna, Pedro H. C, and Song, Xiaojun (2016), \emph{Specification Tests for the Propensity Score},
#'       available at \url{http://sites.google.com/site/pedrohcsantanna/}.
#'
#'@docType package
#'@name pstest-package
NULL
