get_HP <- function(kp) {
    A <- c(16.8244, 0.323365, -4.86128)
    B <- c(1.82336, 0.613192, 26.1798)
    if (kp <= 6.9) {
        HP <- A[1] * exp(A[2] * kp) + A[3]
    } else {
        HP <- B[1] * exp(B[2] * kp) + B[3]
    }
    return (HP)
}

get_fitting_value <- function(t, r, AA, BB) {
    ang <- (0: (nrow(AA)-1)) * t
    cosang <- cos(ang)
    sinang <- sin(ang)

    ncoeff <- ncol(AA)
    coeff <- array(dim = ncoeff)
    for (i in 1: ncoeff) coeff[i] <- sum(AA[, i] * cosang + BB[, i] * sinang)

    F1 <- exp((r - coeff[2]) / coeff[3])
    F2 <- 1 + exp((r - coeff[2]) / coeff[4])
    value <- coeff[1] * F1 / F2^2

    return (value)
}

get_guvi_kp_model <- function(kp, mlat, mlt) {
    nc <- ncdf4::nc_open(filename = 'guvi_auroral_model.nc')
    kp_list <- ncdf4::ncvar_get(nc = nc, varid = "kp")
    AA_flux <- ncdf4::ncvar_get(nc = nc, varid = "AA_flux")
    BB_flux <- ncdf4::ncvar_get(nc = nc, varid = "BB_flux")
    AA_energy <- ncdf4::ncvar_get(nc = nc, varid = "AA_energy")
    BB_energy <- ncdf4::ncvar_get(nc = nc, varid = "BB_energy")
    ncdf4::nc_close(nc = nc)

    ida <- which.min(abs(kp_list - kp))
    if (kp_list[ida] > kp) ida <- ida - 1
    idb <- ida + 1

    x0 <- (kp_list[idb] - kp) / (kp_list[idb] - kp_list[ida])
    x1 <- (kp - kp_list[ida]) / (kp_list[idb] - kp_list[ida])

    HPa <- get_HP(kp_list[ida])
    HPb <- get_HP(kp_list[idb])
    HP <- get_HP(kp)

    xa <- (HPb - HP) / (HPb - HPa)
    xb <- (HP - HPa) / (HPb - HPa)

    n <- length(mlat)
    flux <- array(dim = n)
    energy <- array(dim = n)

    for (i in 1: n) {
        r <- 90 - mlat[i]
        t <- mlt[i] * pi/12
        Va <- get_fitting_value(t, r, AA_flux[, ida, ], BB_flux[, ida, ])
        Vb <- get_fitting_value(t, r, AA_flux[, idb, ], BB_flux[, idb, ])
        flux[i] <- Va*xa + Vb*xb
        Va <- get_fitting_value(t, r, AA_energy[, ida, ], BB_energy[, ida, ])
        Vb <- get_fitting_value(t, r, AA_energy[, idb, ], BB_energy[, idb, ])
        energy[i] <- Va*x0 + Vb*x1
    }

    return (list(flux = flux, energy = energy))
}