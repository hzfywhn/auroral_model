source(file = "LatticeKriging.r")

xylim <- c(-2, 2)

# small test set for 1-1 comparison with LatticeKrig package
xy <- seq(from = xylim[1], to = xylim[2], by = 0.5)
coor <- expand.grid(xy, xy)
loc <- cbind(coor$Var1, coor$Var2)
obs <- list(loc = loc, val = loc[, 1]^2 + loc[, 2]^2, err = rep(1, times = nrow(loc)))

# build up basis structure
delta <- 1
xy <- seq(from = xylim[1], to = xylim[2], by = 1)
nc <- length(xy)
m <- nc^2
loc <- array(dim = c(m, 2))
connect <- array(dim = c(m, 4))
for (i in 1: nc) {
    for (j in 1: nc) {
        idx <- (i-1)*nc + j

        # column-major order for 1-1 comparison with LatticeKrig package
        loc[idx, 1] <- xy[j]
        loc[idx, 2] <- xy[i]

        # connectivity matrix, at most 2/4/6 elements in 1/2/3 dimensions
        if (i >= 2) connect[idx, 1] <- (i-2)*nc + j
        if (j >= 2) connect[idx, 2] <- (i-1)*nc + j-1
        if (i <= nc-1) connect[idx, 3] <- i*nc + j
        if (j <= nc-1) connect[idx, 4] <- (i-1)*nc + j+1
    }
}
basis1 <- list(loc = loc, connect = connect, centerweight = 4.01, delta = delta*2.5, alpha = 1)

# combine basis sets into one large basis structure
basis <- list(basis1)

LatticeKriging_res <- LatticeKriging(obs, basis, Z = NULL, rho = 1, normalize = TRUE, reps = 20)

xy <- seq(from = xylim[1], to = xylim[2], by = 0.25)
coor <- expand.grid(xy, xy)
loc <- cbind(coor$Var1, coor$Var2)
predictMean_res <- predictMean(loc, basis, Z = NULL, normalize = TRUE, Qchol = LatticeKriging_res$Qchol, rho = 1, LatticeKriging_res$d, LatticeKriging_res$c)
dif <- loc[, 1]^2 + loc[, 2]^2 - predictMean_res
