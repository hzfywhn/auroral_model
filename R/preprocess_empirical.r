preprocess_empirical <- function(emp, aurtype, mlat_emp, mlt_emp, interp, grnd, alphaRadius) {
    # keep only one type of aurora
    emp$flux <- emp$flux[, aurtype, , ]
    emp$energy <- emp$energy[, aurtype, , ]

    # add ghost points in mlt
    emp$mlt <- c(emp$mlt, 24)

    nt <- nrow(emp$flux)
    nmlt <- length(mlt_emp)
    nmlat <- length(mlat_emp)
    nmlt0 <- length(emp$mlt)
    nmlat0 <- length(emp$mlat)

    # add ghost points in flux and energy
    flux <- array(dim = c(nt, nmlt0, nmlat0))
    flux[, 1: (nmlt0-1), ] <- emp$flux
    flux[, nmlt0, ] <- emp$flux[, 1, ]
    emp$flux <- flux
    energy <- array(dim = c(nt, nmlt0, nmlat0))
    energy[, 1: (nmlt0-1), ] <- emp$energy
    energy[, nmlt0, ] <- emp$energy[, 1, ]
    emp$energy <- energy

    # interpolate empirical data to new grids
    aux1 <- array(dim = c(nmlt0, nmlat))
    aux2 <- array(dim = c(nmlt0, nmlat))
    flux <- array(dim = c(nt, nmlt, nmlat))
    energy <- array(dim = c(nt, nmlt, nmlat))
    for (it in 1: nt) {
        for (imlt in 1: nmlt0) {
            aux1[imlt, ] <- approx(x = emp$mlat, y = emp$flux[it, imlt, ], xout = mlat_emp)$y
            aux2[imlt, ] <- approx(x = emp$mlat, y = emp$energy[it, imlt, ], xout = mlat_emp)$y
        }
        for (imlat in 1: nmlat) {
            flux[it, , imlat] <- approx(x = emp$mlt, y = aux1[, imlat], xout = mlt_emp)$y
            energy[it, , imlat] <- approx(x = emp$mlt, y = aux2[, imlat], xout = mlt_emp)$y
        }
    }
    emp$flux <- flux
    emp$energy <- energy

    # reformat mlt, mlat into 2D array
    coor <- expand.grid(mlt_emp, mlat_emp)
    emp$mlt <- array(data = coor[, 1], dim = c(nmlt, nmlat))
    emp$mlat <- array(data = coor[, 2], dim = c(nmlt, nmlat))

    # empirical data is valid only if it is out of the observation region
    r <- pi/2 - emp$mlat*pi/180
    t <- emp$mlt * pi/12
    loc0 <- cbind(c(r * cos(t)), c(r * sin(t)))

    emp_new <- vector(mode = "list", length = nt)
    for (it in 1: nt) {
        valid_interp <- interp$flux[, , it] > 0
        valid_interp[is.na(valid_interp)] <- FALSE
        valid_grnd <- grnd$flux[, , it] > 0
        valid_grnd[is.na(valid_grnd)] <- FALSE
        mlt_obs <- c(interp$mlt[valid_interp], grnd$mlt[, , it][valid_grnd])
        mlat_obs <- c(interp$mlat[valid_interp], grnd$mlat[, , it][valid_grnd])
        loc_obs <- unique(cbind(mlt_obs, mlat_obs))
        r <- pi/2 - loc_obs[, 2]*pi/180
        t <- loc_obs[, 1] * pi/12
        valid <- !inahull(ahull(x = r * cos(t), y = r * sin(t), alpha = alphaRadius), loc0)
        mlt <- atan2(loc0[valid, 2], loc0[valid, 1]) * 12/pi
        mlat <- 90 - sqrt(loc0[valid, 1]^2 + loc0[valid, 2]^2)*180/pi
        flux <- c(emp$flux[it, , ])[valid]
        energy <- c(emp$energy[it, , ])[valid]
        emp_new[[it]] <- data.frame(mlt, mlat, flux, energy)
    }

    return (emp_new)
}