source(file = "LatticeKriging.r")

xylim <- c(-2, 2)

# small test set for 1-1 comparison with LatticeKrig package
xy <- seq(from = xylim[1], to = xylim[2], by = 0.5)
coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]
z <- x^2 + y^2
obs <- list(loc = cbind(x, y), val = z, err = rep(1, times = length(z)))

write.csv(data.frame(x, y, z), file = "obs.csv")

# build up basis structure
delta <- 1
xy <- seq(from = xylim[1], to = xylim[2], by = delta)
nc <- length(xy)

# column-major order for 1-1 comparison with LatticeKrig package
coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]
write.csv(data.frame(x, y), file = "basis.csv")

# connectivity matrix, at most 2/4/6 elements in 1/2/3 dimensions
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

# combine basis sets into one large basis structure
basis <- list(basis1)

LatticeKriging_res <- LatticeKriging(obs, basis, NULL, 1, TRUE, 20)

xy <- seq(from = xylim[1], to = xylim[2], by = 0.25)
coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]
z <- predictMean(cbind(x, y), basis, NULL, 1, TRUE, LatticeKriging_res$Qchol, LatticeKriging_res$d, LatticeKriging_res$c)

write.csv(data.frame(x, y, z), file = "predict.csv")

dif <- x^2 + y^2 - z
