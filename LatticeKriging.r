# a simplified rewrite of the original LatticeKrig package on Cartesian geometry

library(package = "spam")

SAR <- function(basis) {
    # spatial autogression
    m <- nrow(basis$loc)
    m0 <- ncol(basis$connect)

    ia <- array(dim = m*m0)
    ja <- array(dim = m*m0)
    a <- array(dim = m*m0)

    k <- 1
    for (j in 1: m) {
        ia[k] <- j
        ja[k] <- j
        a[k] <- basis$centerweight
        k <- k + 1

        for (j0 in 1: m0) {
            js <- basis$connect[j, j0]
            if (!is.na(js)) {
                ia[k] <- j
                ja[k] <- js
                a[k] <- -1
                k <- k + 1
            }
        }
    }

    ia <- ia[1: (k-1)]
    ja <- ja[1: (k-1)]
    a <- a[1: (k-1)] / sqrt(basis$alpha)

    return (list(i = ia, j = ja, values = a))
}

precision <- function(basis, rho) {
    B <- spam(x = SAR(basis))
    Q <- t(B) %*% B / rho
    return (triplet(Q, tri=TRUE))
}

regression <- function(loc, basis, derivative) {
    n <- nrow(loc)
    m <- nrow(basis$loc)
    ndim <- ncol(loc)

    nr <- n
    if (derivative) nr <- ndim * n

    ia <- array(dim = nr*m)
    ja <- array(dim = nr*m)
    a <- array(dim = nr*m)

    k <- 1
    for (i in 1: n) {
        for (j in 1: m) {
            displace <- loc[i, ] - basis$loc[j, ]
            d <- norm(x = displace, type = "2") / basis$delta
            if (d <= 1) {
                if (derivative) {
                    r <- -(1 - d)^5 * (5*d + 1) * 56/3 / basis$delta^2
                    for (idim in 1: ndim) {
                        ia[k] <- i + (idim-1)*n
                        ja[k] <- j
                        a[k] <- r * displace[idim]
                        k <- k + 1
                    }
                } else {
                    ia[k] <- i
                    ja[k] <- j
                    a[k] <- (1 - d)^6 * (35*d^2 + 18*d + 3) / 3
                    k <- k + 1
                }
            }
        }
    }

    ia <- ia[1: (k-1)]
    ja <- ja[1: (k-1)]
    a <- a[1: (k-1)]

    return (list(i = ia, j = ja, values = a))
}

combineMR <- function(loc, basis, derivative) {
    n <- nrow(loc)
    ndim <- ncol(loc)
    nlev <- length(basis)

    nr <- n
    if (derivative) nr <- ndim * n

    m <- 0
    for (ilev in 1: nlev) m <- m + nrow(basis[[ilev]]$loc)

    iphi <- array(dim = nr*m)
    jphi <- array(dim = nr*m)
    phi <- array(dim = nr*m)
    kphi <- 1

    iQ <- array(dim = m^2)
    jQ <- array(dim = m^2)
    Q <- array(dim = m^2)
    kQ <- 1

    for (ilev in 1: nlev) {
        phi0 <- regression(loc, basis[[ilev]], derivative)
        k0 <- length(phi0$values)
        iphi[kphi: (kphi+k0-1)] <- phi0$i
        jphi[kphi: (kphi+k0-1)] <- phi0$j
        phi[kphi: (kphi+k0-1)] <- phi0$values
        kphi <- kphi + k0

        Q0 <- precision(basis[[ilev]], rho)
        k0 <- length(Q0$values)
        iQ[kQ: (kQ+k0-1)] <- Q0$i
        jQ[kQ: (kQ+k0-1)] <- Q0$j
        Q[kQ: (kQ+k0-1)] <- Q0$values
        kQ <- kQ + k0
    }

    iphi <- iphi[1: (kphi-1)]
    jphi <- jphi[1: (kphi-1)]
    phi <- phi[1: (kphi-1)]

    iQ <- iQ[1: (kQ-1)]
    jQ <- jQ[1: (kQ-1)]
    Q <- Q[1: (kQ-1)]

    phi <- spam(x = list(i = iphi, j = jphi, values = phi))
    Q <- spam(x = list(i = iQ, j = jQ, values = Q))

    return (list(phi = phi, Q = Q))
}

covariate <- function(loc) cbind(1, loc)

errorCov <- function(err) diag.spam(x = 1/err^2)

constants <- function(obs, basis, normalization, rho, derivative) {
    y <- c(obs$val)
    W <- errorCov(c(obs$err))

    # concatenate vector observation into an array (x1 ... xn y1 ... yn ...)
    Z <- covariate(obs$loc)
    if (derivative) {
        # extend covariate vector to length of obs
        Z <- covariate(do.call(what = rbind, args = replicate(n = ncol(obs$loc), expr = obs$loc, simplify = FALSE)))
    }

    MR <- combineMR(obs$loc, basis, derivative)

    if (normalization) {
        # still don't understand how it is related to Sec 2.6
        normweight <- rho * colSums(forwardsolve(chol(MR$Q), t(MR$phi))^2)
        stopifnot(normweight != 0)
        MR$phi <- diag.spam(x = 1/sqrt(normweight)) %*% MR$phi
    }

    return (c(list(y = y, W = W, Z = Z), MR))
}

kriging <- function(lambda, y, W, Z, phi, Q) {
    n <- length(y)

    # auxiliary matrix M_lambda
    M <- phi %*% solve(Q, t(phi)) + lambda * chol2inv(chol(W))

    # here are the most expansive computations taken place
    d <- solve(t(Z) %*% solve(M, Z), t(Z)) %*% solve(M, y)
    r <- y - Z %*% d
    c <- solve(Q, t(phi)) %*% solve(M, r)
    rhoMLE <- t(r) %*% solve(M, r) / n
    likelihood <- -n/2 - n/2*log(rhoMLE) - sum(log(diag(chol(M))))
    # the last term +n/2*log(pi) in Eq 7 was wrong, corrected here
    likelihood <- likelihood - n/2*log(2*pi)

    return (list(d = drop(d), c = drop(c), rhoMLE = drop(rhoMLE), likelihood = drop(likelihood), M = M))
}

prediction <- function(loc, basis, normalization, rho, lambda, Z, phi, Q, M, d, c, rhoMLE) {
    # no vector simulation for now
    Z1 <- covariate(loc)
    phi1 <- combineMR(loc, basis, FALSE)$phi

    normweight <- rho * colSums(forwardsolve(chol(Q), t(phi1))^2)
    stopifnot(normweight != 0)
    if (normalization) phi1 <- diag.spam(x = 1/sqrt(normweight)) %*% phi1

    m <- Z1 %*% d + phi1 %*% c

    # standard deviation prediction is confusing with no reference found, need verification
    y1 <- phi %*% solve(Q, t(phi1))
    ZMZ <- t(Z) %*% solve(M, Z)
    d1 <- solve(ZMZ, t(Z)) %*% solve(M, y1)
    r1 <- y1 - Z %*% d1
    c1 <- solve(Q, t(phi)) %*% solve(M, r1)
    residual <- y1 - Z %*% d1 - phi %*% c1
    joint <- colSums(t(Z1) * solve(ZMZ, t(Z1))) - 2 * colSums(t(Z1) * d1)
    # seems like only diagonals of error covariance matrix is used in marginal variance calculation
    # TODO: the second term left out error covariance
    marginal <- normweight - colSums(y1 * residual) / lambda
    sd <- sqrt(rhoMLE * abs(joint + marginal))

    return (list(m = drop(m), sd = sd))
}

# codes below are for testing
if (FALSE) {
    # testing domain
    x1 <- -2
    x2 <- 2
    y1 <- -1
    y2 <- 1

    # fake observation
    delta <- 0.2
    x0 <- seq(from = x1, to = x2, by = delta)
    y0 <- seq(from = y1, to = y2, by = delta)
    loc <- data.matrix(frame = expand.grid(x0, y0))
    x <- loc[, 1]
    y <- loc[, 2]
    z <- 1 / (1 + (x-1)^2 + y^2) - 1 / (1 + (x+1)^2 + y^2)
    vx <- -2*(x-1) / (1 + (x-1)^2 + y^2)^2 + 2*(x+1) / (1 + (x+1)^2 + y^2)^2
    vy <- -2*y / (1 + (x-1)^2 + y^2)^2 + 2*y / (1 + (x+1)^2 + y^2)^2
    val <- z
    obs <- list(loc = loc, val = val, err = array(data = 1, dim = length(val)))

    # basis function
    delta <- 0.5
    x0 <- seq(from = x1, to = x2, by = delta)
    y0 <- seq(from = y1, to = y2, by = delta)
    nx <- length(x0)
    ny <- length(y0)
    # constructing connectivity matrix is simple when the whole domain is filled
    connect <- array(dim = c(nx*ny, 4))
    for (i in 1: nx) {
        for (j in 1: ny) {
            idx <- (i-1)*ny + j
            if (i >= 2) connect[idx, 1] <- (i-2)*ny + j
            if (j >= 2) connect[idx, 2] <- (i-1)*ny + j-1
            if (i <= nx-1) connect[idx, 3] <- i*ny + j
            if (j <= ny-1) connect[idx, 4] <- (i-1)*ny + j+1
        }
    }
    basis1 <- list(loc = data.matrix(frame = expand.grid(x0, y0)), connect = connect,
        centerweight = 4.01, delta = delta*2.5, alpha = 1)
    # combine different levels
    basis <- list(basis1)

    normalization <- TRUE
    rho <- 1
    derivative <- FALSE
    # interval needs to be adjusted for specific purposes
    con <- constants(obs, basis, normalization, rho, derivative)
    lambda <- exp(optimize(
        function(l) kriging(exp(l), con$y, con$W, con$Z, con$phi, con$Q)$likelihood,
        interval = c(-12, -8), maximum = TRUE)$maximum)
    fit <- kriging(lambda, con$y, con$W, con$Z, con$phi, con$Q)
    pred <- prediction(loc, basis, normalization, rho, lambda, con$Z, con$phi, con$Q, fit$M, fit$d, fit$c, fit$rhoMLE)
    cat(x, "\n", y, "\n", pred$m, "\n", z, file = "test.txt")
}