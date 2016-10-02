#' pstest: Tests for the Propensity Score
#'
#' \emph{pstest} computes Kolmogorov-Smirnov and Cramer-von Mises type tests
#' for the null hypothesis that a parametric model for the propensity score is
#' is correctly specified. For details of the testing procedure, see
#' Sant'Anna and Song (2016),'Specification Tests for the Propensity Score'.
#'
#'@param d vector containing the binary treatment indicator
#'@param pscore vector containing the fitted propensity score
#'@param xpscore matrix (or data frame) containing the covariates (and their
#'               transformations) included in the propensity score
#'               estimation. You should always include a constant term.
#'@param model functional form (link function) employed while estimating the propensity score.
#'              The alternatives are: 'logit' (default), and 'probit'.
#'@param nboot number of bootstrap draws. Default is 1,000.
#'@param cores number of cores to use during the bootstrap (default is 1).
#'        If cores>1, the bootstrap is conducted using parLapply, instead
#'        of lapply type call.
#'
#'@return a list containing the Kolmogorov-Smirnov and Cramer-von Mises test
#'        statistics for the null hypothesis of correctly specified propensity
#'        score model (kstest and cvmtest, respectively), and their associate bootstrapped
#'        p-values, pvks and pvcvm, respectively.
#'
#'@export
#'@importFrom stats binomial glm rbinom runif
#'@importFrom MASS ginv
#'@importFrom parallel makeCluster parLapply stopCluster
#'@importFrom harvestr gather


###############################################################
pstest = function(d, pscore, xpscore, model = c("logit", "probit"), nboot = 1000, cores = 1) {
    n <- length(d)
    xx <- as.matrix(xpscore)
    pscore.fit <- pscore
    uhat <- d - pscore.fit
    if (model == "logit") {
        g <- pscore.fit * (1 - pscore.fit) * xx
    }
    if (model == "probit") {
        beta.x <- stats::qnorm(pscore.fit)
        g <- stats::dnorm(beta.x) * xx
    }
    gg <- (t(g) %*% g)
    gginv <- solve(gg)
    w <- (outer(pscore.fit, unique(pscore.fit), "<="))
    Gw <- t(g) %*% w
    beta <- gginv %*% Gw
    w1 <- w - g %*% beta

    Rw <- colSums(uhat * w1)/n
    cvmtest <- sum(Rw^2)
    kstest <- sqrt(n) * max(abs(Rw))


    # Use the Mammen(1993) binary V's
    k1 <- 0.5 * (1 - 5^0.5)
    k2 <- 0.5 * (1 + 5^0.5)
    pkappa <- 0.5 * (1 + 5^0.5)/(5^0.5)

    ## Define seeds
    ss <- floor(stats::runif(1) * 10000)
    seed.temp <- harvestr::gather(nboot, seed = ss)

    Seed <- matrix(nrow = nboot, ncol = 6)
    for (i in 1:nboot) {
        Seed[i, ] <- seed.temp[[i]][2:7]
    }

    bootapply <- function(nn, n, pkappa, k1, k2, uhat, Seed) {
        # to make each run fully reproducible, we set the seed
        seed.run <- Seed[nn, ]
        set.seed(seed.run, "L'Ecuyer-CMRG")
        v <- stats::rbinom(n, 1, pkappa)
        v <- ifelse(v == 1, k1, k2)
        # Bootstrapped emprirical process
        Rwb <- colSums(uhat * v * w1)/n
        # KS test
        ksb <- sqrt(n) * max(abs(Rwb))
        # Cramer-von Mises test
        cvmb <- sum(Rwb^2)
        # Return both tests
        return(cbind(ksb, cvmb))
    }

    if (cores==1){
      boottests <- lapply(1:nboot, bootapply, n, pkappa, k1, k2,
                          uhat, Seed)
    }

    if (cores>1){
      cl <- parallel::makeCluster(cores)
      boottests <- parallel::parLapply(cl, 1:nboot, bootapply,
                                     n, pkappa, k1, k2,
                                     uhat, Seed)
    parallel::stopCluster(cl)
    }

    # Put the Bootstrap resuls in a matrix
    boottest <- t(matrix(unlist(boottests), 2, nboot))

    # Name the Columns
    colnames(boottest) <- c("ksb", "cvmb")

    # compute the Bootstrap P-value
    pvksb <- sum((boottest[, 1] > kstest))/nboot
    pvcvmb <- sum((boottest[, 2] > cvmtest))/nboot

    list(kstest = kstest, cvmtest = cvmtest,
         pvks = pvksb, pvcvm = pvcvmb)
}
