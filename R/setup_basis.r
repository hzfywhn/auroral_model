setup_basis <- function(bndry, delta, centerweight, overlap, weight) {
    nt <- length(bndry$rmax)
    basis <- vector(mode = "list", length = nt)

    for (it in 1: nt) {
        # fill up full domain
        xy <- seq(from = -bndry$rmax[it], to = bndry$rmax[it], by = delta)
        nxy <- length(xy)
        coor <- expand.grid(xy, xy)

        # select valid basis functions (within auroral boundaries)
        radius2 <- array(data = coor[, 1]^2 + coor[, 2]^2, dim = c(nxy, nxy))
        valid <- radius2 >= bndry$rmin[it]^2 & radius2 <= bndry$rmax[it]^2

        x <- array(dim = nxy^2)
        y <- array(dim = nxy^2)
        index <- array(dim = c(nxy, nxy))
        # setup index matrix: mark index on valid basis function
        idx <- 1
        for (i in 1: nxy) {
            for (j in 1: nxy) {
                if (valid[i, j]) {
                    x[idx] <- xy[i]
                    y[idx] <- xy[j]
                    index[i, j] <- idx
                    idx <- idx + 1
                }
            }
        }

        m <- idx - 1
        x <- x[1: m]
        y <- y[1: m]
        con <- array(dim = c(m, 4))
        # construct connectivity matrix using index matrix
        idx <- 1
        for (i in 1: nxy) {
            for (j in 1: nxy) {
                if (valid[i, j]) {
                    if (i >= 2 && valid[i-1, j]) con[idx, 1] <- index[i-1, j]
                    if (j >= 2 && valid[i, j-1]) con[idx, 2] <- index[i, j-1]
                    if (i <= nxy-1 && valid[i+1, j]) con[idx, 3] <- index[i+1, j]
                    if (j <= nxy-1 && valid[i, j+1]) con[idx, 4] <- index[i, j+1]
                    idx <- idx + 1
                }
            }
        }

        basis1 <- list(loc = cbind(x + bndry$am[it], y + bndry$bm[it]), connect = con, centerweight = centerweight, delta = overlap, alpha = weight)

        # combine different levels, single level for now
        basis[[it]] <- list(basis1)
    }

    return (basis)
}