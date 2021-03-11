interp_satellite <- function(ut, flux, energy, time) {
    nmlat <- nrow(ut)
    nmlt <- ncol(ut)
    nt <- length(time)
    flux_interp <- array(dim = c(nmlat, nmlt, nt))
    energy_interp <- array(dim = c(nmlat, nmlt, nt))

    for (imlat in 1: nmlat) {
        for (imlt in 1: nmlt) {
            t <- ut[imlat, imlt, ]

            f <- flux[imlat, imlt, ]
            valid <- !(is.na(t) | is.na(f))
            x <- t[valid]
            if (length(unique(x)) >= 2) {
                flux_interp[imlat, imlt, ] <- approx(x = x, y = f[valid], xout = time)$y
            }

            e <- energy[imlat, imlt, ]
            valid <- !(is.na(t) | is.na(e))
            x <- t[valid]
            if (length(unique(x)) >= 2) {
                energy_interp[imlat, imlt, ] <- approx(x = x, y = e[valid], xout = time)$y
            }
        }
    }

    return (list(flux = flux_interp, energy = energy_interp))
}