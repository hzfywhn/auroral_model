auroral_boundary <- function(x, y, flux, lb, a, b, r) {
    nx <- nrow(x)
    ny <- ncol(y)

    na <- length(a)
    nb <- length(b)
    nr <- length(r)

    r2 <- r^2
    rlb2 <- (r[1] - (r[2] - r[1]) / 2)^2
    rub2 <- (r[nr] + (r[nr] - r[nr-1]) / 2)^2

    nt <- length(flux[1, 1, ])
    am <- array(dim = nt)
    bm <- array(dim = nt)
    rmin2 <- array(dim = nt)
    rmax2 <- array(dim = nt)

    for (it in 1: nt) {
        aurora <- flux[, , it] >= lb
        aurora[is.na(aurora)] <- FALSE
        vote <- array(data = 0, dim = c(na, nb, nr))

        for (ix in 1: nx) {
            for (iy in 1: ny) {
                if (aurora[ix, iy]) {
                    for (ia in 1: na) {
                        for (ib in 1: nb) {
                            d2 <- (x[ix, iy] - a[ia])^2 + (y[ix, iy] - b[ib])^2
                            if (d2 >= rlb2 && d2 <= rub2) {
                                ir <- which.min(abs(r2 - d2))
                                vote[ia, ib, ir] <- vote[ia, ib, ir] + 1
                            }
                        }
                    }
                }
            }
        }

        idx <- which(vote == max(vote), arr.ind = TRUE)
        am[it] <- a[idx[1, 1]]
        bm[it] <- b[idx[1, 2]]

        aur <- flux[, , it] > 0
        aur[is.na(aur)] <- FALSE
        rm2 <- (x[aur] - am[it])**2 + (y[aur] - bm[it])**2
        rmin2[it] <- min(rm2)
        rmax2[it] <- max(rm2)
    }

    return (data.frame(am, bm, rmin = sqrt(rmin2), rmax = sqrt(rmax2)))
}