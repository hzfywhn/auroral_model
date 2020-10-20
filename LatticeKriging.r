# a simplified rewrite of the original LatticeKrig package on Cartesian geometry

library(package = "spam")

SAR <- function(basis) {
    # spatial autogression
    m <- nrow(basis$loc)
    B <- spam(x = 0, nrow = m, ncol = m)
    for (j in 1: m) {
        B[j, j] <- basis$centerweight
        con <- basis$connect[j, ]
        B[j, con[!is.na(con)]] <- -1
    }
    Bs <- B / sqrt(basis$alpha)
    return (Bs)
}

precision <- function(basis, rho) {
    B <- SAR(basis)
    Q <- t(B) %*% B / rho
    return (Q)
}

regression <- function(loc, basis) {
    n <- nrow(loc)
    m <- nrow(basis$loc)
    delta2 <- basis$delta^2
    dist <- spam(x = 0, nrow = n, ncol = m)
    # pairwise distance
    for (i in 1: n) {
        for (j in 1: m) {
            dist2 <- sum((loc[i, ] - basis$loc[j, ])^2)
            if (dist2 <= delta2) dist[i, j] <- sqrt(dist2)
        }
    }
    d <- dist / basis$delta
    phi <- (1 - d)^6 * (35*d^2 + 18*d + 3) / 3
    return (phi)
}

combineMR <- function(loc, basis, normalize, rho) {
    n <- nrow(loc)
    nlev <- length(basis)

    m0 <- array(dim = nlev)
    for (ilev in 1: nlev)
        m0[ilev] <- nrow(basis[[ilev]]$loc)
    m <- cumsum(c(0, m0))

    phi <- spam(x = 0, nrow = n, ncol = m[nlev+1])
    Qc <- spam(x = 0, nrow = m[nlev+1], ncol = m[nlev+1])
    for (ilev in 1: nlev) {
        start_index <- m[ilev] + 1
        end_index <- m[ilev+1]

        Q <- precision(basis[[ilev]], rho)
        Qc0 <- chol(Q, pivot = FALSE)
        phi0 <- regression(loc, basis[[ilev]])

        if (normalize) {
            # still don't understand how it is related to Sec 2.6
            normweight <- rho * colSums(forwardsolve(Qc0, t(phi0))^2)
            stopifnot(normweight != 0)
            phi0 <- diag.spam(1/sqrt(normweight)) %*% phi0
        }

        phi[, start_index: end_index] <- phi0
        Qc[start_index: end_index, start_index: end_index] <- as.spam(Qc0)
    }

    return (list(phi = phi, Qc = Qc))
}

covariate <- function(loc) {
    Z <- cbind(1, loc)
    return (Z)
}

errorCov <- function(err) {
    W <- diag.spam(x = 1/err^2)
    return (W)
}

LatticeKriging <- function(obs, basis, normalize, lambda, rho) {
    n <- nrow(obs$loc)

    MR <- combineMR(obs$loc, basis, normalize, rho)
    Qinv <- chol2inv(MR$Qc)

    Z <- covariate(obs$loc)
    y <- obs$val
    W <- errorCov(obs$err)

    M <- MR$phi %*% Qinv %*% t(MR$phi) + lambda * chol2inv(chol(W, pivot = FALSE))

    # here is the most expansive computation taken place
    Mc <- chol(M, pivot = FALSE)
    Minv <- chol2inv(Mc)
    tZMinvZinv <- chol2inv(chol(t(Z) %*% Minv %*% Z, pivot = FALSE))

    d <- tZMinvZinv %*% t(Z) %*% Minv %*% y
    r <- y - Z %*% d
    c <- Qinv %*% t(MR$phi) %*% Minv %*% r
    rhoMLE <- t(r) %*% Minv %*% r / n
    likelihood <- -n/2 - n/2*log(rhoMLE) - sum(log(diag(Mc)))
    # discrepancy happens, the last term is +n/2*log(pi) in Eq 7
    likelihood <- likelihood - n/2*log(2*pi)

    # collapse 1x1 matrix to scalar
    rhoMLE <- rhoMLE[1, 1]
    return (c(MR, list(Qinv = Qinv, Z = Z, W = W, d = d, c = c, rhoMLE = rhoMLE,
        likelihood = likelihood, Minv = Minv, tZMinvZinv = tZMinvZinv)))
}

predictLK <- function(obs, basis, loc, lambda, normalize, rho) {
    fit <- LatticeKriging(obs, basis, normalize, lambda, rho)

    Z1 <- covariate(loc)
    phi1 <- combineMR(loc, basis, normalize, rho)$phi

    mean <- Z1 %*% fit$d + phi1 %*% fit$c

    # standard deviation prediction is confusing with no reference found, need verification
    y1 <- fit$phi %*% fit$Qinv %*% t(phi1)
    d1 <- fit$tZMinvZinv %*% t(fit$Z) %*% fit$Minv %*% y1
    r1 <- y1 - fit$Z %*% d1
    c1 <- fit$Qinv %*% t(fit$phi) %*% fit$Minv %*% r1
    resid <- y1 - fit$Z %*% d1 - fit$phi %*% c1
    joint <- colSums(t(Z1) * (fit$tZMinvZinv %*% t(Z1))) - 2 * colSums(t(Z1) * d1)
    # seems like only diagonals of error covariance matrix is used in marginal variance calculation
    marginal <- rho * colSums(forwardsolve(fit$Qc, t(phi1))^2) - colSums(y1*resid) / lambda
    sd <- sqrt(fit$rhoMLE * abs(joint + marginal))

    return (list(mean = mean, sd = sd))
}


# codes below are for testing

# testing domain
xylim <- c(-2, 2)

# fake observation
xy <- seq(from = xylim[1], to = xylim[2], by = 0.5)
coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]
z <- x^2 + y^2
obs <- list(loc = cbind(x, y), val = z, err = rep(1, times = length(z)))

# distance between basis functions
delta <- 1
xy <- seq(from = xylim[1], to = xylim[2], by = delta)
nc <- length(xy)

coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]

# constructing connectivity matrix is simple when the whole domain is filled
connect <- array(dim = c(nc^2, 4))
for (i in 1: nc) {
    for (j in 1: nc) {
        idx <- (i-1)*nc + j
        if (i >= 2) connect[idx, 1] <- (i-2)*nc + j
        if (j >= 2) connect[idx, 2] <- (i-1)*nc + j-1
        if (i <= nc-1) connect[idx, 3] <- i*nc + j
        if (j <= nc-1) connect[idx, 4] <- (i-1)*nc + j+1
    }
}
basis1 <- list(loc = cbind(x, y), connect = connect, centerweight = 4.01, delta = delta*2.5, alpha = 1)
basis <- list(basis1)

normalize <- TRUE
rho <- 1
lambda <- exp(optimize(
    function(loglambda) {
        fit <- LatticeKriging(obs, basis, normalize, exp(loglambda), rho)
        return (fit$likelihood)
    }, interval = c(-9, 4), maximum = TRUE, tol = 5e-3)$maximum)
fit <- LatticeKriging(obs, basis, normalize, lambda, rho)

xy <- seq(from = xylim[1], to = xylim[2], by = 0.2)
coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]
LK <- predictLK(obs, basis, cbind(x, y), lambda, normalize, rho)