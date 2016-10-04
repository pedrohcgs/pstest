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
#'@param model functional form (link function) employed while estimating the
#'              propensity score. The alternatives are: 'logit' (default),
#'              and 'probit'.
#'@param nboot number of bootstrap draws. Default is 1,000.
#'@param cores number of cores to use during the bootstrap (default is 1).
#'              If cores>1, the bootstrap is conducted using parLapply, instead
#'              of lapply type call.
#'@param big logical value indication if you have "big data". Default is
#'        FALSE. It is recommended to first try with the default. If function
#'        return a memory error, try to set big = TRUE to see if solves the
#'        issue.
#'
#'@return a list containing the Kolmogorov-Smirnov and Cramer-von Mises test
#'        statistics for the null hypothesis of correctly specified propensity
#'        score model (kstest and cvmtest, respectively), and their associate
#'        bootstrapped p-values, pvks and pvcvm, respectively.
#'
#'@export
#'@importFrom stats binomial glm rbinom runif
#'@importFrom MASS ginv
#'@importFrom parallel makeCluster parLapply stopCluster
#'@importFrom harvestr gather


#-------------------------------------------------------------------------------
pstest = function(d, pscore, xpscore, model = c("logit", "probit"),
                  nboot = 1000, cores = 1, big = FALSE) {

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
        rm(beta.x)
    }
    gg <- crossprod(g)

    #----------------------------------------------------------------------------
    # Handle first the case in which I have enough memory to solve
    if (big == FALSE) {
        w <- (outer(pscore.fit, unique(pscore.fit), "<="))
        Gw <- crossprod(g, w)
        beta <- solve(gg, Gw)
        w1 <- w - g %*% beta
        # Get the functions we need
        Rw <- colSums(uhat * w1)/n
        cvmtest1 <- sum(Rw^2)
        kstest1 <- sqrt(n) * max(abs(Rw))

        # Use the Mammen(1993) binary V's
        k1 <- 0.5 * (1 - 5^0.5)
        k2 <- 0.5 * (1 + 5^0.5)
        pkappa <- 0.5 * (1 + 5^0.5)/(5^0.5)

        # Define seeds
        ss <- floor(stats::runif(1) * 10000)
        seed.temp <- harvestr::gather(nboot, seed = ss)

        Seed <- matrix(nrow = nboot, ncol = 6)
        for (i in 1:nboot) {
            Seed[i, ] <- seed.temp[[i]][2:7]
        }

        bootapply <- function(nn, n, pkappa, k1, k2, uhat, w1, Seed) {
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

        if (cores == 1) {
            boottests <- lapply(1:nboot, bootapply, n, pkappa, k1, k2,
                                uhat, w1, Seed)
        }

        if (cores > 1) {
            cl <- parallel::makeCluster(cores)
            boottests <- parallel::parLapply(cl, 1:nboot, bootapply, n,
                                             pkappa, k1, k2,
                uhat, w1, Seed)
            parallel::stopCluster(cl)
        }

        # Put the Bootstrap resuls in a matrix
        boottest <- t(matrix(unlist(boottests), 2, nboot))

        # Name the Columns
        colnames(boottest) <- c("ksb", "cvmb")

        # compute the Bootstrap P-value
        pvksb <- sum((boottest[, 1] > kstest1))/nboot
        pvcvmb <- sum((boottest[, 2] > cvmtest1))/nboot
    }
    #------------------------------------------------------------------------

    # Now, let us try to solve the case with 'big data'
    if (big == TRUE) {
        # Number of covariates
        k.dim = dim(xx)[2]

        # unique pscores
        un.pscores <- unique(pscore.fit)
        n.unique <- length(un.pscores)
        # Initialize `beta` matrix (K coefficients for each of
        #n.unique responses)
        beta <- matrix(0, k.dim, n.unique)
        # Initialize `Rw` row vector (n.unique dimension)
        Rw <- matrix(0, 1, n.unique)

        # Initialize the bootststrap vector
        boottest <- matrix(0, nboot, 2)

        # We split n columns into l tiles, each with 1000 columns
        l <- floor(n.unique/1000) + 1

        # Let's define some parameters for the bootstrap
        # Better to define these outside the
        # loop, so we define only once

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

        bootapply <- function(nn, n, pkappa, k1, k2, uhat, w1, Seed) {
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

        for (i in 1:l) {
            start <- min(1000 * (i - 1) + 1, n.unique)  ## chunk start
            end <- min(1000 * i, n.unique)  ## chunk end
            w.temp <- outer(pscore.fit, un.pscores[start:end], "<=")
            Gw <- crossprod(g, w.temp)
            beta[, start:end] <- solve(gg, Gw)
            w1.temp <- uhat * (w.temp - g %*% beta[, start:end])
            Rw[start:end] <- colSums(w1.temp)/n
            # Now the bootstrapped test in the chunk
            if (cores == 1) {
                boot.chunk <- lapply(1:nboot, bootapply, n, pkappa, k1,
                                     k2, uhat, w1.temp, Seed)
            }

            if (cores > 1) {
                cl <- parallel::makeCluster(cores)
                boot.chunk <- parallel::parLapply(cl, 1:nboot, bootapply,
                                                  n, pkappa, k1,
                                                  k2, uhat, w1.temp, Seed)
            }

            # Put the Bootstrap resuls in a matrix
            boot.chunk1 <- t(matrix(unlist(boot.chunk), 2, nboot))

            # Compute the KSb and CvMb over chunks
            if (1000 * (i - 1) + 1 <= n.unique) {
                # First KsB (maximum between the KSB in this chunk and
                # the previous maximum)
                boottest[, 1] <- pmax(boottest[, 1], boot.chunk1[, 1])
                # Now the Cvmb, which i just need to sum
                boottest[, 2] <- sum(boottest[, 2], boot.chunk1[, 2])
            }

        }
        # close the clusters
        if (cores > 1) {
          parallel::stopCluster(cl)
        }


        cvmtest1 <- sum(Rw^2)
        kstest1 <- sqrt(n) * max(abs(Rw))
        # Name the Columns
        colnames(boottest) <- c("ksb", "cvmb")

        # compute the Bootstrap P-value
        pvksb <- sum((boottest[, 1] > kstest1))/nboot
        pvcvmb <- sum((boottest[, 2] > cvmtest1))/nboot

    }

    #---------------------------------------------------------------------
    # Return these

    list(kstest = kstest1, cvmtest = cvmtest1, pvks = pvksb, pvcvm = pvcvmb)
}
