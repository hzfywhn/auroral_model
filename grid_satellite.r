grid_satellite <- function(mlat, mlt, ut, flux, energy, time, coverage) {
    reps <- length(ut[1, 1, ])
    mlat <- replicate(n = reps, expr = mlat)
    mlt <- replicate(n = reps, expr = mlt)

    nt <- length(time)
    gridded <- vector(mode = "list", length = nt)

    for (i in 1: nt) {
        mask <- ut > time[i] - coverage & ut < time[i] + coverage
        mask[is.na(mask)] = FALSE
        gridded[[i]] <- data.frame(mlat = mlat[mask], mlt = mlt[mask], flux = flux[mask], energy = energy[mask])
    }

    return (gridded)
}