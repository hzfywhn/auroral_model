ground <- function(timein, mlat, mlt, flux, energy, timeout) {
    ntin <- length(timein)
    ntout <- length(timeout)
    nmlat <- length(mlat[1, , 1])
    nmlt <- length(mlt[, 1, 1])
    mlatout <- array(dim = c(nmlt, nmlat, ntout))
    mltout <- array(dim = c(nmlt, nmlat, ntout))
    fluxout <- array(dim = c(nmlt, nmlat, ntout))
    energyout <- array(dim = c(nmlt, nmlat, ntout))

    s <- 1
    while (s <= ntout & timeout[s] < timein[1]) s <- s + 1

    j <- 1
    for (i in s: ntout) {
        while (j <= ntin - 1 & timein[j + 1] <= timeout[i]) j <- j + 1
        if (j == ntin) break
        mlatout[, , i] <- mlat[, , j]
        mltout[, , i] <- mlt[, , j]
        fluxout[, , i] <- flux[, , j]
        energyout[, , i] <- energy[, , j]
    }

    return (list(mlat = mlatout, mlt = mltout, flux = fluxout, energy = energyout))
}