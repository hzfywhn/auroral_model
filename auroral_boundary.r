auroral_boundary <- function(x, y, ut, flux, time, ub, lb, a, b, r) {

    nx <- nrow(x)
    ny <- ncol(y)

    na <- length(a)
    nb <- length(b)
    nr <- length(r)

    r2 <- r^2
    rlb2 <- (r[1] - (r[2] - r[1]) / 2)^2
    rub2 <- (r[nr] + (r[nr] - r[nr-1]) / 2)^2

    nt <- length(time)
    am <- array(dim = c(nt))
    bm <- array(dim = c(nt))
    rmin2 <- array(dim = c(nt))
    rmax2 <- array(dim = c(nt))

    for (i in 1: nt) {
        if (i == 1) t1 <- time[1] - (time[2] - time[1]) / 2
        else t1 <- (time[i-1] + time[i]) / 2

        if (i == nt) t2 <- time[nt] + (time[nt] - time[nt-1]) / 2
        else t2 <- (time[i] + time[i+1]) / 2

        f <- flux
        f[ut < t1 | ut > t2] <- NA
        f <- rowMeans(f, na.rm = TRUE, dims = 2)

        aurora <- f >= ub
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
        am[i] <- a[idx[1, 1]]
        bm[i] <- b[idx[1, 2]]

        aur <- f >= lb
        aur[is.na(aur)] <- FALSE
        rm2 <- (x[aur] - am[i])**2 + (y[aur] - bm[i])**2
        rmin2[i] <- min(rm2)
        rmax2[i] <- max(rm2)
    }

    return (data.frame(am, bm, rmin = sqrt(rmin2), rmax = sqrt(rmax2)))
}