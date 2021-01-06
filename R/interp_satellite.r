interp_satellite <- function(ut, flux, energy, time, max_interval) {
    nmlat <- nrow(ut)
    nmlt <- ncol(ut)
    nt <- length(time)
    flux_interp <- array(dim = c(nmlat, nmlt, nt))
    energy_interp <- array(dim = c(nmlat, nmlt, nt))

    for (imlat in 1: nmlat) {
        for (imlt in 1: nmlt) {
            x <- ut[imlat, imlt, ]
            y1 <- flux[imlat, imlt, ]
            y2 <- energy[imlat, imlt, ]
            valid <- !(is.na(x) | is.na(y1) | is.na(y2))
            x <- x[valid]
            y1 <- y1[valid]
            y2 <- y2[valid]
            if (length(x) >= 2) {
                ix <- sort(x, index.return = TRUE)$ix
                x <- x[ix]

                if (min(diff(x)) < max_interval) {
                    flux_interp[imlat, imlt, ] <- approx(x = x, y = y1[ix], xout = time)$y
                    energy_interp[imlat, imlt, ] <- approx(x = x, y = y2[ix], xout = time)$y
                }
            }
        }
    }

    return (list(flux = flux_interp, energy = energy_interp))
}