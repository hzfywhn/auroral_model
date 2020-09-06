grid_satellite <- function(mlat, mlt, ut, flux, energy, time) {
    reps <- length(ut[1, 1, ])
    mlat <- replicate(n = reps, expr = mlat)
    mlt <- replicate(n = reps, expr = mlt)

    nt <- length(time)
    gridded <- list()

    # Low efficiency due to unfixed length of each data frame
    for (i in 1: nt) {
        if (i == 1) t1 <- time[1] - (time[2] - time[1]) / 2
        else t1 <- (time[i-1] + time[i]) / 2

        if (i == nt) t2 <- time[nt] + (time[nt] - time[nt-1]) / 2
        else t2 <- (time[i] + time[i+1]) / 2

        mask <- ut > t1 & ut < t2
        mask[is.na(mask)] = FALSE
        gridded[[i]] <- data.frame(mlat = mlat[mask], mlt = mlt[mask], flux = flux[mask], energy = energy[mask])
    }

    return (gridded)
}