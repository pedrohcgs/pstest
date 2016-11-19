#' pstest: An R Package for assessing the goodness-of-fit of parametric propensity score models.
#'
#'
#'An R package for implementing all the tests proposed in Sant'Anna and
#'Song (2016), 'Specification Tests for the Propensity Score', available at
#'Pedro H.C. Sant'Anna webpage, \url{http://sites.google.com/site/pedrohcsantanna/}.
#'
#'In short, this package implements Kolmogorov-Smirnov
#'and Cramer-von Mises type tests for parametric propensity score models with
#'either logistic ('logit'), or standard normal ('probit') link function. Critical values are
#'computed with the assistance of a multiplier bootstrap.
#'
#'The tests are based on the integrated conditional moment approach, where the weight function
#'used is based on an orthogonal projection onto the tangent space of nuisance parameters.
#'As a result, the tests (a) enjoy improved power properties, (b) do not suffer from the
#''curse of dimensionality' when the vector of covariates is of high-dimensionality,
#'(c) are fully data-driven, (e) do not require tuning parameters such as bandwidths, and
#'(e) are able to detect a broad class of local alternatives converging to the null at the
#' parametric rate. These appealing features highlight that the tests can be of great use
#' in practice.
#'
#' @section Authors:
#'       The \emph{pstest} package was written by Pedro H. C. Sant'Anna (Vanderbilt University),
#'       \email{pedro.h.santanna@@vanderbilt.edu}, and Xiaojun Song (Peking University),
#'       \email{sxj@@gsm.pku.edu.cn}.
#'
#'       In case you have any question, please do not hesitate to contact us.
#' @references
#'       Sant'Anna, Pedro H. C, and Song, Xiaojun (2016), \emph{Specification Tests for the Propensity Score},
#'       available at \url{http://sites.google.com/site/pedrohcsantanna/}.
#'
#'
"_PACKAGE"
