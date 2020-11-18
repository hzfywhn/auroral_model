interp_satellite <- function(ut, flux, energy, time) {
    nmlt <- length(ut[, 1, 1])
    nmlat <- length(ut[1, , 1])
    nt <- length(time)
    flux_interp <- array(dim = c(nmlt, nmlat, nt))
    energy_interp <- array(dim = c(nmlt, nmlat, nt))

    for (imlt in 1: nmlt) {
        for (imlat in 1: nmlat) {
            x <- ut[imlt, imlat, ]
            y1 <- flux[imlt, imlat, ]
            y2 <- energy[imlt, imlat, ]
            valid <- !(is.na(x) | is.na(y1) | is.na(y2))
            x <- x[valid]
            y1 <- y1[valid]
            y2 <- y2[valid]
            if (length(x) >= 2) {
                ix <- sort(x, index.return = TRUE)$ix
                x <- x[ix]

                # TODO: control of time interval
                flux_interp[imlt, imlat, ] <- approx(x = x, y = y1[ix], xout = time)$y
                energy_interp[imlt, imlat, ] <- approx(x = x, y = y2[ix], xout = time)$y
            }
        }
    }

    return (list(flux = flux_interp, energy = energy_interp))
}