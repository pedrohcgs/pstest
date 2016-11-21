# pstest: An R Package to assess the goodness-of-fit of parametric propensity score models.

## Overview 
The propensity score is one of the most widely used tools in studying the causal effect
of a treatment, intervention, or policy. Given that the propensity score is usually unknown,
it has to be estimated, implying that the reliability of many treatment effect estimators depends
on the correct specification of the (parametric) propensity score. This package provides
data-driven nonparametric diagnostic tools for detecting propensity score misspecification.

In short, this package implements the class of specification test for the propensity score
proposed in Sant'Anna and Song (2016), [Specification Tests for the Propensity Score](https://papers.ssrn.com/abstract=2872084 " target="_blank). The package accommodates Logit and Probit propensity score specifications. Critical values are computed with the assistance of a multiplier bootstrap.

The tests are based on the integrated conditional moment approach, where the weight function
used is based on an orthogonal projection onto the tangent space of nuisance parameters.
As a result, the tests (a) enjoy improved power properties, (b) do not suffer from the
'curse of dimensionality' when the vector of covariates is of high-dimensionality,
(c) are fully data-driven, (e) do not require tuning parameters such as bandwidths, and
(e) are able to detect a broad class of local alternatives converging to the null at the
parametric rate. These appealing features highlight that the tests can be of great use
in practice.

It is worth stressing that the `pstest` package implements in a unified manner a large class of specification tests, depending on the chosen weight function w(q,u). Current choices are:

* `ind` - the indicator weight function w(q,u)=1(q<=u). This is the default.
* `exp` - the exponential weight function w(q,u)=exp(qu).
* `logistic` - the logistic weight function w(q,u)=1/[1+exp(-qu)].
* `sin` - the sine weight function w(q,u)=sin(qu).
* `sincos` - the sine and cosine weight function w(q,u)=sin(qu)+cos(qu).

Different weight functions w(q,u) have different power properties. Thus, being able to choose between different w(q,u) gives one the flexibility to direct power in alternative directions.

For further details, please see the paper Sant'Anna and Song (2016), [Specification Tests for the Propensity Score](https://papers.ssrn.com/abstract=2872084 " target="_blank).

## Installing pstest
This github website hosts the source code. The difference between what is here what is in CRAN is that here we always have the most updated version of the package.

To install the `pstest` package, you have two options: (a) install the CRAN version, or (b) instal the GitHub (most updated) version. Both procedures are very easy!

You can install the package from CRAN as follows:

        install.packages("pstest")

Alternatively, you can install the most updated version of the `pstest` package from GitHub:

        devtools::install_github("pedrohcgs/pstest")
For this, the package devtools is required.          
        
## Authors 

Pedro H. C. Sant'Anna, Vanderbilt University, Nashville, TN. E-mail: pedro.h.santanna [at] vanderbilt [dot] edu.

Xiaojun Song, Peking University, Beijing, China. E-mail: sxj [at] gsm [dot] pku [dot] edu [dot] cn.

In case you have questions, please do not hesitate to contact us.

