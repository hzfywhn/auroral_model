# a simplified rewrite of the original LatticeKrig package on Cartesian geometry
# no standard error estimate yet, implementing

library(package = "spam")

LKfit <- function(lambda, wy, wU, wX, Gchol) {
    # calculate coefficients
    A <- t(wU) %*% (wU - wX %*% backsolve(Gchol, forwardsolve(Gchol, t(wX) %*% wU))) / lambda
    b <- t(wU) %*% (wy - wX %*% backsolve(Gchol, forwardsolve(Gchol, t(wX) %*% wy))) / lambda
    d <- solve(A) %*% b
    r <- wy - wU %*% d
    c <- forwardsolve(Gchol, t(wX) %*% r)
    rhoMLE <- (sum(r^2) - sum(c^2)) / lambda / length(wy)
    c <- backsolve(Gchol, c)
    return (list(d = d, c = c, rhoMLE = rhoMLE))
}

LK <- function(lambda, wy, wU, wX, Q, Qchol, rho, weights, reps) {
    # LK for fixed lambda

    n <- length(wy)
    m <- nrow(Q)

    G <- t(wX) %*% wX + lambda * Q

    # pivoting can be turned on depending on performance
    Gchol <- chol(G, pivot = FALSE)

    LKfit_res <- LKfit(lambda, wy, wU, wX, Gchol)
    residual <- (wy - (wU %*% LKfit_res$d + wX %*% LKfit_res$c)) / sqrt(weights)

    # likelihood
    lnDetCov <- 2*sum(log(diag(Gchol))) - 2*sum(log(diag(Qchol))) + (n-m)*log(lambda) - sum(log(weights))
    lnProfileLike <- -n/2 - n/2*log(2*pi) - n/2*log(LKfit_res$rhoMLE) - 1/2*lnDetCov
    lnLike <- -LKfit_res$rhoMLE*n/(2*rho) - n/2*log(2*pi) - n/2*log(rho) - 1/2*lnDetCov

    est <- NA
    SE <- NA
    GCV <- NA
    if (reps >= 1) {
        wEy <- matrix(data = rnorm(reps * n), nrow = n, ncol = reps) * sqrt(weights)
        trAest <- LKfit(lambda, wEy, wU, wX, Gchol)
        wEyhat <- wU %*% trAest$d + wX %*% trAest$c
        info <- colSums(wEy * wEyhat / weights)
        est <- mean(info)
        SE <- sqrt(var(info) / length(info))
        GCV <- sum(weights * residual^2) / n / (1 - est/n)^2
    }

    return (c(LKfit_res, list(lnProfileLike = lnProfileLike, lnLike = lnLike, est = est, SE = SE, GCV = GCV)))
}

LKfixed <- function(loc, Z) {
    # include possible covariates
    return (cbind(1, loc, Z))
}

LKbasis <- function(obsloc, basis, rho, normalize, Qchol) {
    # pairwise distance and basis evaluation (phi matrix)

    n <- nrow(obsloc)
    nlev <- length(basis)

    mall <- 0
    for (ilev in 1: nlev)
        mall <- mall + nrow(basis[[ilev]]$loc)
    Xall <- spam(x = 0, nrow = n, ncol = mall)

    mall <- 0
    for (ilev in 1: nlev) {
        m <- nrow(basis[[ilev]]$loc)

        distance <- spam(x = 0, nrow = n, ncol = m)
        for (i in 1: n) {
            for (j in 1: m) {
                dist2 <- sum((obsloc[i, ] - basis[[ilev]]$loc[j, ])^2)
                if (dist2 <= basis[[ilev]]$delta^2) distance[i, j] <- sqrt(dist2) / basis[[ilev]]$delta
            }
        }
        phi <- (1 - distance)^6 * (35*distance^2 + 18*distance + 3) / 3

        X <- sqrt(basis[[ilev]]$alpha) * phi

        if (normalize) {
            normweight <- colSums(forwardsolve(Qchol[[ilev]], t(phi))^2)
            normweight[normweight == 0] <- 1
            X <- diag.spam(1/sqrt(normweight)) %*% X
        }

        Xall[, (mall+1): (mall+m)] <- X
        mall <- mall + m
    }

    return (sqrt(rho) * Xall)
}

LatticeKriging <- function(obs, basis, Z, rho, normalize, reps) {
    # entrance of LK, optimize with lambda

    n <- nrow(obs$loc)
    nlev <- length(basis)

    # need to verify
    weights <- 1 / obs$err^2

    wy <- sqrt(weights) * obs$val

    # fixed part
    wU <- sqrt(weights) * LKfixed(obs$loc, Z)

    # precision matrix
    mall <- 0
    for (ilev in 1: nlev)
        mall <- mall + nrow(basis[[ilev]]$loc)
    Qall <- spam(x = 0, nrow = mall, ncol = mall)
    Qchol <- vector(mode = "list", length = nlev)
    Qcholall <- spam(x = 0, nrow = mall, ncol = mall)

    mall <- 0
    for (ilev in 1: nlev) {
        m <- nrow(basis[[ilev]]$loc)

        B <- spam(x = 0, nrow = m, ncol = m)
        for (j in 1: m) {
            B[j, j] <- basis[[ilev]]$centerweight
            con <- basis[[ilev]]$connect[j, ]
            B[j, con[!is.na(con)]] <- -1
        }
        Q <- t(B) %*% B

        # disable pivoting for matrix assignment
        Qchol[[ilev]] <- chol(Q, pivot = FALSE)

        Qall[(mall+1): (mall+m), (mall+1): (mall+m)] <- Q
        Qcholall[(mall+1): (mall+m), (mall+1): (mall+m)] <- as.spam(Qchol[[ilev]])
        mall <- mall + m
    }

    # random part
    X <- LKbasis(obs$loc, basis, rho, normalize, Qchol)
    wX <- diag.spam(sqrt(weights)) %*% X

    # main loop for optimization
    opt <- optimize(
        function(loglambda) {
            LK_res <- LK(exp(loglambda), wy, wU, wX, Qall, Qcholall, rho, weights, 0)
            return (LK_res$lnProfileLike)
        }, interval = c(-9, 4), maximum = TRUE, tol = 5e-3)

    lambda <- exp(opt$maximum)
    LK_res <- LK(lambda, wy, wU, wX, Qall, Qcholall, rho, weights, reps)

    return (c(LK_res, list(Qchol = Qchol, lambda = lambda)))
}

predictMean <- function(loc, basis, Z, rho, normalize, Qchol, d, c) {
    return (LKfixed(loc, Z) %*% d + LKbasis(loc, basis, rho, normalize, Qchol) %*% c)
}

predictSD <- function() {
    # to be completed
    return ()
}