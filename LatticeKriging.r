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

regression <- function(loc, basis, derivative) {
    n <- nrow(loc)
    m <- nrow(basis$loc)
    ndim <- ncol(loc)

    nr <- n
    if (derivative) nr <- ndim * n

    phi <- spam(x = 0, nrow = nr, ncol = m)
    for (i in 1: n) {
        for (j in 1: m) {
            displace <- loc[i, ] - basis$loc[j, ]
            d <- norm(x = displace, type = "2") / basis$delta
            if (d <= 1) {
                if (derivative) {
                    r <- -(1 - d)^5 * (5*d + 1) * 56/3 / basis$delta^2
                    for (idim in 1: ndim) phi[i+(idim-1)*n, j] <- r * displace[idim]
                } else
                    phi[i, j] <- (1 - d)^6 * (35*d^2 + 18*d + 3) / 3
            }
        }
    }

    return (phi)
}

combineMR <- function(loc, basis, derivative) {
    n <- nrow(loc)
    ndim <- ncol(loc)
    nlev <- length(basis)

    nr <- n
    if (derivative) nr <- ndim * n

    m0 <- array(dim = nlev)
    for (ilev in 1: nlev) m0[ilev] <- nrow(basis[[ilev]]$loc)
    m <- cumsum(c(0, m0))

    phi <- spam(x = 0, nrow = nr, ncol = m[nlev+1])
    Q <- spam(x = 0, nrow = m[nlev+1], ncol = m[nlev+1])
    for (ilev in 1: nlev) {
        start_index <- m[ilev] + 1
        end_index <- m[ilev+1]

        phi[, start_index: end_index] <- regression(loc, basis[[ilev]], derivative)
        Q[start_index: end_index, start_index: end_index] <- precision(basis[[ilev]], rho)
    }

    return (list(phi = phi, Q = Q))
}

covariate <- function(loc) cbind(1, loc)

errorCov <- function(err) diag.spam(x = 1/err^2)

LatticeKriging <- function(obs, basis, normalization, lambda, rho, derivative) {
    n <- nrow(obs$loc)
    ndim <- ncol(obs$loc)

    nr <- n
    z <- obs$loc
    if (derivative) {
        nr <- ndim * n
        # extend covariate vector to length of obs
        z <- do.call(what = rbind, args = replicate(n = ndim, expr = obs$loc, simplify = FALSE))
    }

    MR <- combineMR(obs$loc, basis, derivative)
    Qc <- chol(MR$Q)

    if (normalization) {
        # still don't understand how it is related to Sec 2.6
        normweight <- rho * colSums(forwardsolve(Qc, t(MR$phi))^2)
        stopifnot(normweight != 0)
        MR$phi <- diag.spam(1/sqrt(normweight)) %*% MR$phi
    }

    # concatenate vector observation into an array (x1 ... xn y1 ... yn ...)
    Z <- covariate(z)
    y <- c(obs$val)
    W <- errorCov(c(obs$err))

    # auxiliary matrix M_lambda
    M <- MR$phi %*% solve(MR$Q, t(MR$phi)) + lambda * chol2inv(chol(W))

    # here are the most expansive computations taken place
    d <- solve(t(Z) %*% solve(M, Z), t(Z)) %*% solve(M, y)
    r <- y - Z %*% d
    c <- solve(MR$Q, t(MR$phi)) %*% solve(M, r)
    rhoMLE <- t(r) %*% solve(M, r) / nr
    likelihood <- -nr/2 - nr/2*log(rhoMLE) - 1/2*determinant(M)$modulus
    # the last term +nr/2*log(pi) in Eq 7 was wrong, corrected here
    likelihood <- likelihood - nr/2*log(2*pi)

    return (c(list(d = drop(d), c = drop(c), rhoMLE = drop(rhoMLE),
        likelihood = drop(likelihood), Z = Z, W = W, M = M, Qc = Qc), MR))
}

predictLK <- function(obs, basis, loc, lambda, normalization, rho, derivative) {
    fit <- LatticeKriging(obs, basis, normalization, lambda, rho, derivative)

    # no vector simulation for now
    Z1 <- covariate(loc)
    phi1 <- combineMR(loc, basis, FALSE)$phi

    normweight <- rho * colSums(forwardsolve(fit$Qc, t(phi1))^2)
    stopifnot(normweight != 0)
    if (normalization) phi1 <- diag.spam(1/sqrt(normweight)) %*% phi1

    m <- Z1 %*% fit$d + phi1 %*% fit$c

    # standard deviation prediction is confusing with no reference found, need verification
    y1 <- fit$phi %*% solve(fit$Q, t(phi1))
    ZMZ <- t(fit$Z) %*% solve(fit$M, fit$Z)
    d1 <- solve(ZMZ, t(fit$Z)) %*% solve(fit$M, y1)
    r1 <- y1 - fit$Z %*% d1
    c1 <- solve(fit$Q, t(fit$phi)) %*% solve(fit$M, r1)
    residual <- y1 - fit$Z %*% d1 - fit$phi %*% c1
    joint <- colSums(t(Z1) * solve(ZMZ, t(Z1))) - 2 * colSums(t(Z1) * d1)
    # seems like only diagonals of error covariance matrix is used in marginal variance calculation
    # TODO: the second term left error covariance
    marginal <- normweight - colSums(y1 * residual) / lambda
    sd <- sqrt(fit$rhoMLE * abs(joint + marginal))

    return (list(m = drop(m), sd = sd))
}

# codes below are for testing
if (FALSE) {
    # testing domain
    x1 <- -2
    x2 <- 2
    y1 <- -1
    y2 <- 1

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

    # fake observation
    delta <- 0.5
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

    normalization <- TRUE
    rho <- 1
    derivative <- FALSE
    # interval needs to be adjusted for specific purposes
    lambda <- exp(optimize(
        function(l) LatticeKriging(obs, basis, normalization, exp(l), rho, derivative)$likelihood,
        interval = c(-10, 10), maximum = TRUE)$maximum)
    fit <- LatticeKriging(obs, basis, normalization, lambda, rho, derivative)
    LK <- predictLK(obs, basis, loc, lambda, normalization, rho, derivative)
    cat(x, "\n", y, "\n", LK$m, "\n", z, file = "test.txt")
}