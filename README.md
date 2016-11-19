# pstest: An R Package to assess the goodness-of-fit of parametric propensity score models.

## Description 
The `pstest` R package implements the specification tests proposed in Sant'Anna and Song (2016), 'Specification Tests for the Propensity Score', available at Pedro H.C. Sant'Anna webpage, http://sites.google.com/site/pedrohcsantanna/ .

In short, this package implement bootstrapped based Kolmogorov-Smirnov and Cramer-von Mises type tests for parametric propensity score models with either logistic ('logit'), or standard normal ('probit') link function.

The tests are based on the integrated conditional moment approach, where the weight function used is based on an orthogonal projection onto the tangent space of nuisance parameters. As a result, the tests (a) enjoy improved power properties, (b) do not suffer from the "curse of dimensionality" when the vector of covariates is of high-dimensionality, (c) are fully data-driven, (d) do not require tuning parameters such as bandwidths, and (e) are able to detect a broad class of local alternatives converging to the null at the parametric rate. These appealing features highlight that the tests can be of great use in practice.

As of now, the `pstest` function accommodates the following weight functions w:
* `ind` - use the indicator weight function w(q,u)=1(q<=u). This is the default.
* `exp` - use the exponential weight function w(q,u)=exp(qu).
* `logistic` - use the logistic weight function w(q,u)=1/[1+exp(-qu)].
* `sin` - use the sine weight function w(q,u)=sin(qu).
* `sincos` - use the sine and cosine weight function w(q,u)=sin(qu)+cos(qu).

For further details, please see the paper. In case you have doubts, don't hesitate to contact us.

## Installing pstest
This github website hosts the source code. To install the `pstest` package, you simply need to run the following two lines of code in R:

        library(devtools)
        install_github("pedrohcgs/pstest")

A help file with examples will be made available in the near future. 

This package was written by Pedro H. C. Sant'Anna (Vanderbilt University), and Xiaojun Song (Peking University).
