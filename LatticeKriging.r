# a simplified rewrite of LatticeKrig package on Cartesian geometry

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
    a <- a[1: (k-1)]

    return (list(i = ia, j = ja, values = a))
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

combineMR <- function(loc, basis, normalization, rho, derivative) {
    n <- nrow(loc)
    ndim <- ncol(loc)
    nlev <- length(basis)

    nr <- n
    if (derivative) nr <- ndim * n

    m0 <- array(dim = nlev)
    for (ilev in 1: nlev) m0[ilev] <- nrow(basis[[ilev]]$loc)
    m1 <- cumsum(c(0, m0))
    m <- m1[nlev+1]

    iB <- array(dim = m^2)
    jB <- array(dim = m^2)
    B <- array(dim = m^2)
    kB <- 1

    iphi <- array(dim = nr*m)
    jphi <- array(dim = nr*m)
    phi <- array(dim = nr*m)
    kphi <- 1

    for (ilev in 1: nlev) {
        B0 <- SAR(basis[[ilev]])
        phi0 <- regression(loc, basis[[ilev]], derivative)

        if (normalization) {
            # don't understand how it is related to Sec 2.6
            B1 <- spam(x = B0, nrow = m0[ilev], ncol = m0[ilev])
            Q1 <- t(B1) %*% B1
            phi1 <- spam(x = phi0, nrow = nr, ncol = m0[ilev])
            normweight <- colSums(forwardsolve(chol(Q1), t(phi1))^2)
            stopifnot(normweight != 0)
            phi1 <- diag.spam(x = 1/sqrt(normweight)) %*% phi1
            phi0 <- triplet(phi1, tri = TRUE)
        }

        k0 <- length(B0$values)
        iB[kB: (kB+k0-1)] <- B0$i + m1[ilev]
        jB[kB: (kB+k0-1)] <- B0$j + m1[ilev]
        # this is different from Eq 11, there is a scaling factor in front
        B[kB: (kB+k0-1)] <- B0$values
        kB <- kB + k0

        k0 <- length(phi0$values)
        iphi[kphi: (kphi+k0-1)] <- phi0$i
        jphi[kphi: (kphi+k0-1)] <- phi0$j + m1[ilev]
        # scaling factor of precision matrix moved here, not clear why
        phi[kphi: (kphi+k0-1)] <- phi0$values * sqrt(basis[[ilev]]$alpha)
        kphi <- kphi + k0
    }

    iB <- iB[1: (kB-1)]
    jB <- jB[1: (kB-1)]
    B <- B[1: (kB-1)]

    iphi <- iphi[1: (kphi-1)]
    jphi <- jphi[1: (kphi-1)]
    # scaling factor of precision matrix moved here, not clear why
    phi <- phi[1: (kphi-1)] * sqrt(rho)

    B <- spam(x = list(i = iB, j = jB, values = B), nrow = m, ncol = m)
    Q <- t(B) %*% B
    phi <- spam(x = list(i = iphi, j = jphi, values = phi), nrow = nr, ncol = m)

    return (list(Q = Q, phi = phi))
}

constants <- function(obs, basis, normalization, rho, derivative) {
    # concatenate vector observation into an array (x1 ... xn y1 ... yn ...)
    y <- c(obs$val)
    W <- diag.spam(x = 1/c(obs$err)^2)

    Z <- cbind(1, obs$loc)
    # extend covariate vector to length of obs
    if (derivative) Z <- cbind(1, do.call(what = rbind,
        args = replicate(n = ncol(obs$loc), expr = obs$loc, simplify = FALSE)))

    MR <- combineMR(obs$loc, basis, normalization, rho, derivative)

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
    Z1 <- cbind(1, loc)
    phi1 <- combineMR(loc, basis, normalization, rho, FALSE)$phi

    m <- Z1 %*% d + phi1 %*% c

    # standard deviation prediction is confusing with no reference found
    # don't trust for now, verification needed
    y1 <- phi %*% solve(Q, t(phi1))
    ZMZ <- t(Z) %*% solve(M, Z)
    d1 <- solve(ZMZ, t(Z)) %*% solve(M, y1)
    r1 <- y1 - Z %*% d1
    c1 <- solve(Q, t(phi)) %*% solve(M, r1)
    residual <- y1 - Z %*% d1 - phi %*% c1
    joint <- colSums(t(Z1) * solve(ZMZ, t(Z1))) - 2 * colSums(t(Z1) * d1)
    marginal <- colSums(forwardsolve(chol(Q), t(phi1))^2) - colSums(y1 * residual) / lambda
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
    delta <- 0.1
    x0 <- seq(from = x1, to = x2, by = delta)
    y0 <- seq(from = y1, to = y2, by = delta)
    loc <- expand.grid(x0, y0)
    x <- loc[, 1]
    y <- loc[, 2]
    z <- 1 / (1 + (x-1)^2 + y^2) - 1 / (1 + (x+1)^2 + y^2)
    vx <- -2*(x-1) / (1 + (x-1)^2 + y^2)^2 + 2*(x+1) / (1 + (x+1)^2 + y^2)^2
    vy <- -2*y / (1 + (x-1)^2 + y^2)^2 + 2*y / (1 + (x+1)^2 + y^2)^2
    obs <- list(loc = cbind(x, y), val = z, err = array(data = 1, dim = length(z)))

    # basis set, column major order for comparison with original LatticeKrig package

    delta <- 0.2
    x0 <- seq(from = x1, to = x2, by = delta)
    y0 <- seq(from = y1, to = y2, by = delta)
    nx <- length(x0)
    ny <- length(y0)
    NC <- max(c(nx, ny))
    loc <- array(dim = c(nx*ny, 2))
    con <- array(dim = c(nx*ny, 4))
    for (j in 1: ny) {
        for (i in 1: nx) {
            idx <- (j-1)*nx + i
            loc[idx, 2] <- y0[j]
            loc[idx, 1] <- x0[i]
            if (j >= 2) con[idx, 1] <- (j-2)*nx + i
            if (i >= 2) con[idx, 2] <- (j-1)*nx + i-1
            if (j <= ny-1) con[idx, 3] <- j*nx + i
            if (i <= nx-1) con[idx, 4] <- (j-1)*nx + i+1
        }
    }
    basis1 <- list(loc = loc, connect = con, centerweight = 4.01, delta = delta*2.5, alpha = 0.8)

    delta <- 0.1
    x0 <- seq(from = x1, to = x2, by = delta)
    y0 <- seq(from = y1, to = y2, by = delta)
    nx <- length(x0)
    ny <- length(y0)
    loc <- array(dim = c(nx*ny, 2))
    con <- array(dim = c(nx*ny, 4))
    for (j in 1: ny) {
        for (i in 1: nx) {
            idx <- (j-1)*nx + i
            loc[idx, 2] <- y0[j]
            loc[idx, 1] <- x0[i]
            if (j >= 2) con[idx, 1] <- (j-2)*nx + i
            if (i >= 2) con[idx, 2] <- (j-1)*nx + i-1
            if (j <= ny-1) con[idx, 3] <- j*nx + i
            if (i <= nx-1) con[idx, 4] <- (j-1)*nx + i+1
        }
    }
    basis2 <- list(loc = loc, connect = con, centerweight = 4.01, delta = delta*2.5, alpha = 0.2)

    # combine different levels
    basis <- list(basis1, basis2)

    normalization <- TRUE
    rho <- 1
    derivative <- FALSE
    cons <- constants(obs, basis, normalization, rho, derivative)
    # interval needs to be adjusted for specific purposes
    lambda <- exp(optimize(
        function(l) kriging(exp(l), cons$y, cons$W, cons$Z, cons$phi, cons$Q)$likelihood,
        interval = c(-9, 5), maximum = TRUE, tol = 5e-3)$maximum)
    fit <- kriging(lambda, cons$y, cons$W, cons$Z, cons$phi, cons$Q)
    # lk <- LatticeKrig::LatticeKrig(x = obs$loc, y = obs$val, nlevel = 2,
    #     NC = NC, NC.buffer = 0, normalize = normalization)

    delta <- 0.1
    x0 <- seq(from = x1, to = x2, by = delta)
    y0 <- seq(from = y1, to = y2, by = delta)
    loc <- expand.grid(x0, y0)
    x <- loc[, 1]
    y <- loc[, 2]
    z <- 1 / (1 + (x-1)^2 + y^2) - 1 / (1 + (x+1)^2 + y^2)
    pred <- prediction(cbind(x, y), basis, normalization, rho, lambda,
        cons$Z, cons$phi, cons$Q, fit$M, fit$d, fit$c, fit$rhoMLE)
    # z and pred$m should be quite close
}