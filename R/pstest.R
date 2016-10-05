#' pstest: Tests for the Propensity Score
#'
#' \emph{pstest} computes Kolmogorov-Smirnov and Cramer-von Mises type tests
#' for the null hypothesis that a parametric model for the propensity score is
#' is correctly specified. For details of the testing procedure, see
#' Sant'Anna and Song (2016),'Specification Tests for the Propensity Score'.
#'
#'@param d a vector containing the binary treatment indicator.
#'@param pscore a vector containing the estimated propensity scores.
#'@param xpscore a matrix (or data frame) containing the covariates (and their
#'               transformations) included in the propensity score
#'               estimation. It should also include the constant term.
#'@param model  a description of the functional form (link function) used
#'              to estimated propensity score. The alternatives are:
#'              'logit' (default), and 'probit'.
#'@param nboot number of bootstrap replicates to perform. Default is 1,000.
#'@param cores number of cores to use during the bootstrap. Default is 1.
#'              If cores is greater than 1, the bootstrap is conducted using
#'              parLapply, instead of lapply type call.
#'@param chunk a value that determine the size of each 'tile'. Such argument is used
#'              to splitthe original data into chunks, saving memory.
#'              Default value is 1,000. If the \emph{pstest} function throw a
#'              memory error, you should choose a smaller value for \emph{chunk}.
#'
#'@return a list containing the Kolmogorov-Smirnov and Cramer-von Mises test
#'        statistics for the null hypothesis of correctly specified propensity
#'        score model (kstest and cvmtest, respectively), and their associate
#'        bootstrapped p-values, pvks and pvcvm, respectively.
#'
#'@export
#'@importFrom stats binomial rbinom runif
#'@importFrom parallel makeCluster parLapply stopCluster
#'@importFrom harvestr gather
#-------------------------------------------------------------------------------
pstest = function(d, pscore, xpscore, model = c("logit", "probit"),
                  nboot = 1000, cores = 1, chunk = 1000) {
    #-----------------------------------------------------------------------------
    # Define some underlying variables
    n <- length(d)
    xx <- as.matrix(xpscore)
    pscore.fit <- pscore
    uhat <- d - pscore.fit
    #-----------------------------------------------------------------------------
    # #Define the score variables for the projection
    if (model == "logit") {
        g <- pscore.fit * (1 - pscore.fit) * xx
    }
    if (model == "probit") {
        beta.x <- stats::qnorm(pscore.fit)
        g <- stats::dnorm(beta.x) * xx
        rm(beta.x)
    }
    gg <- crossprod(g)
    #-----------------------------------------------------------------------------
    # Define variables to be used in the loop
    # Number of covariates
    k.dim = dim(xx)[2]

    # unique pscores
    un.pscores <- unique(pscore.fit)
    n.unique <- length(un.pscores)

    # Initialize `beta` matrix (K coefficients for each of n.unique responses)
    beta <- matrix(0, k.dim, n.unique)

    # Initialize `Rw` row vector (n.unique dimension)
    Rw <- matrix(0, 1, n.unique)

    # We split n columns into l tiles, each with chunk columns
    l <- floor(n.unique/chunk) + 1

    # Initialize the bootststrap vector
    #ksb1 <- matrix(0, nboot, l)
    #cvmb1 <- matrix(0, nboot, l)
    ksb1 <- matrix(0, nboot, 1)
    cvmb1 <- matrix(0, nboot, 1)
    #-----------------------------------------------------------------------------
    # Let's define some parameters for the bootstrap
    # Better to define these outside the loop that will follow.

    # Use the Mammen(1993) binary V's
    k1 <- 0.5 * (1 - 5^0.5)
    k2 <- 0.5 * (1 + 5^0.5)
    pkappa <- 0.5 * (1 + 5^0.5)/(5^0.5)

    # function for the bootstrap
    bootapply <- function(nn, n, pkappa, k1, k2, uhat, w1.temp, Seed) {
        # to make each run fully reproducible, we set the seed
        seed.run <- Seed[nn, ]
        set.seed(seed.run, "L'Ecuyer-CMRG")
        v <- stats::rbinom(n, 1, pkappa)
        v <- ifelse(v == 1, k1, k2)
        # Bootstrapped emprirical process
        Rwb <- colSums(uhat * v * w1.temp)/n
        # KS test
        ksb <- sqrt(n) * max(abs(Rwb))
        # Cramer-von Mises test
        cvmb <- sum(Rwb^2)
        # Return both tests
        return(cbind(ksb, cvmb))
    }
    #-----------------------------------------------------------------------------
    # Define seeds: Guarantee reproducibility
    ss <- floor(stats::runif(1) * 10000)
    seed.temp <- harvestr::gather(nboot, seed = ss)

    Seed <- matrix(nrow = nboot, ncol = 6)
    for (i in 1:nboot) {
      Seed[i, ] <- seed.temp[[i]][2:7]
    }
    #-----------------------------------------------------------------------------
    # If we are going to use paralell coding, initialize the cores
    if (cores > 1) {
        cl <- parallel::makeCluster(cores)
    }
    #-----------------------------------------------------------------------------
    # Start the loop to compute the tests (this is more memory efficient)
    for (i in 1:l) {
        start <- min(chunk * (i - 1) + 1, n.unique)
        end <- min(chunk * i, n.unique)
        w.temp <- outer(pscore.fit, un.pscores[start:end], "<=")
        Gw <- crossprod(g, w.temp)
        beta[, start:end] <- solve(gg, Gw)
        w1.temp <- (w.temp - g %*% beta[, start:end])
        Rw[start:end] <- colSums(uhat * w1.temp)/n
        # Now the bootstrapped test in the chunk
        if (cores == 1) {
            boot.chunk <- lapply(1:nboot, bootapply, n, pkappa, k1, k2,
                                 uhat, w1.temp, Seed)
        }
        if (cores > 1) {
            boot.chunk <- parallel::parLapply(cl, 1:nboot, bootapply, n,
                                              pkappa, k1, k2, uhat, w1.temp, Seed)
        }
        # Put the Bootstrap resuls in a matrix
        boot.chunk <- t(matrix(unlist(boot.chunk), 2, nboot))
        # Compute the KSb and CvMb over chunks
        if (1000 * (i - 1) + 1 <= n.unique) {
          ksb1 <- pmax(ksb1,boot.chunk[, 1])
          cvmb1 <- cvmb1+boot.chunk[, 2]
            #ksb1[, i] <- boot.chunk[, 1]
            #cvmb1[, i] <- boot.chunk[, 2]
        }
    }
    #-----------------------------------------------------------------------------
    # close the clusters, if we used paralell
    if (cores > 1) {
        parallel::stopCluster(cl)
    }
    #-----------------------------------------------------------------------------
    # Compute our test statistics
    cvmtest1 <- sum(Rw^2)
    kstest1 <- sqrt(n) * max(abs(Rw))
    #-----------------------------------------------------------------------------
    # Put the bootstrap tests in a matrix
    boottest <- matrix(0, nboot, 2)
    #boottest[, 1] <- apply(ksb1, 1, max)
    #boottest[, 2] <- apply(cvmb1, 1, sum)
    boottest[, 1] <- max(ksb1)
    boottest[, 2] <- sum(cvmb1)
    #-----------------------------------------------------------------------------
    # Name the Columns
    colnames(boottest) <- c("ksb", "cvmb")
    #-----------------------------------------------------------------------------
    # compute the Bootstrap P-value
    pvksb <- sum((boottest[, 1] > kstest1))/nboot
    pvcvmb <- sum((boottest[, 2] > cvmtest1))/nboot
    #---------------------------------------------------------------------
    # Return these variables
    list(kstest = kstest1, cvmtest = cvmtest1, pvks = pvksb, pvcvm = pvcvmb)
}
