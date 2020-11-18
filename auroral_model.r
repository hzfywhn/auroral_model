library(package = "ncdf4")

source(file = "grid_satellite.r")
source(file = "interp_satellite.r")
source(file = "ground.r")
source(file = "empirical.r")
source(file = "auroral_boundary.r")
source(file = "LatticeKriging.r")

# adjustable parameters
doy <- 51
time <- seq(from = 0, to = 24, by = 1/6)
coverage <- 1/2
ratio <- c(1, 1/6, 1/6, 1/36)
delta <- pi/180
xy <- seq(from = -pi*2/9, to = pi*2/9, by = delta*2)
normalization <- TRUE
rho <- 1
derivative <- FALSE

nc <- nc_open(filename = "ssusi/ssusi_combine_20140219_20140221.nc")
mlat_sat <- ncvar_get(nc = nc, varid = "mlat")
mlt_sat <- ncvar_get(nc = nc, varid = "mlt")
ut_sat <- ncvar_get(nc = nc, varid = "ut_n") - 24 # satellite data time offset
flux_sat <- ncvar_get(nc = nc, varid = "flux_n")
energy_sat <- ncvar_get(nc = nc, varid = "energy_n")
nc_close(nc = nc)

nc <- nc_open(filename = "themis/themis2014022006_mag.nc")
time_grnd <- ncvar_get(nc = nc, varid = "time") / 60 + 6 # ground data time offset
mlat_grnd <- ncvar_get(nc = nc, varid = "mlat")
mlt_grnd <- ncvar_get(nc = nc, varid = "mlt")
flux_grnd <- ncvar_get(nc = nc, varid = "flux")
energy_grnd <- ncvar_get(nc = nc, varid = "energy") * 1e-3
nc_close(nc = nc)

obs <- grid_satellite(mlat = mlat_sat, mlt = mlt_sat, ut = ut_sat, flux = flux_sat, energy = energy_sat, time = time, coverage = coverage)

interp <- interp_satellite(ut = ut_sat, flux = flux_sat, energy = energy_sat, time = time)
lon2 <- mlt_sat * 15
lat2 <- mlat_sat

grnd <- ground(timein = time_grnd, mlat = mlat_grnd, mlt = mlt_grnd, flux = flux_grnd, energy = energy_grnd, timeout = time)

# simplified dFdt calculation
f <- read.table(file = "ovation/omni2.lst")
hour <- (f[, 2] - doy) * 24 + f[, 3]
dFdt <- f[, 6]^(4/3) * sqrt(f[, 4]^2 + f[, 5]^2)^(2/3) * abs(sin(atan2(f[, 4], f[, 5])/2))^(8/3)
dFdt <- approx(x = hour, y = dFdt, xout = time)$y
emp <- empirical(doy = doy, dFdt, hemi = "N", premodel = "ovation/ovation-prime_2.3.0/premodel")

# interpolate empirical data in mlt to be comparable with mlat resolution
dmlt <- 1/30
mlt_emp_0 <- c(emp$mlt, 24)
mlt_emp_1 <- seq(from = 0, to = 24-dmlt, by = dmlt)
nt <- length(time)
ntype <- length(emp$type)
nmlt_0 <- length(mlt_emp_0)
nmlt_1 <- length(mlt_emp_1)
nmlat <- length(emp$mlat)
flux_emp_0 <- array(dim = c(nt, ntype, nmlt_0, nmlat))
flux_emp_0[, , 1: (nmlt_0-1), ] <- emp$flux
flux_emp_0[, , nmlt_0, ] <- emp$flux[, , 1, ]
energy_emp_0 <- array(dim = c(nt, ntype, nmlt_0, nmlat))
energy_emp_0[, , 1: (nmlt_0-1), ] <- emp$energy
energy_emp_0[, , nmlt_0, ] <- emp$energy[, , 1, ]
flux_emp_1 <- array(dim = c(nt, ntype, nmlt_1, nmlat))
energy_emp_1 <- array(dim = c(nt, ntype, nmlt_1, nmlat))
for (it in 1: nt) {
    for (itype in 1: ntype) {
        for (imlat in 1: nmlat) {
            flux_emp_1[it, itype, , imlat] <- approx(x = mlt_emp_0, y = flux_emp_0[it, itype, , imlat], xout = mlt_emp_1)
            energy_emp_1[it, itype, , imlat] <- approx(x = mlt_emp_0, y = energy_emp_0[it, itype, , imlat], xout = mlt_emp_1)
        }
    }
}
emp$mlt <- mlt_emp_1
emp$flux <- flux_emp_1
emp$energy <- energy_emp_1

coor <- expand.grid(emp$mlt * 15, emp$mlat)
lon4 <- coor[, 1]
lat4 <- coor[, 2]

r2 <- pi/2 - lat2 * pi/180
t2 <- lon2 * pi/180
x2 <- r2 * cos(t2)
y2 <- r2 * sin(t2)
bndry <- auroral_boundary(x = x2, y = y2, ut = ut_sat, flux = flux_sat, time = time, coverage = coverage, ub = 5, lb = 1,
    a = seq(from = -pi/18, to = pi/18, length.out = 21),
    b = seq(from = -pi/18, to = pi/18, length.out = 21),
    r = seq(from = pi/18, to = 4*pi/18, length.out = 31))
# center: bndry$am, bndry$bm; boundary: bndry$rmin, bndry$rmax

nt <- length(time)
nxy <- length(xy)
fluxsim <- array(dim = c(nt, nxy, nxy))
energysim <- array(dim = c(nt, nxy, nxy))
locsim <- expand.grid(xy, xy)

for (it in 1: nt) {
    lon1 <- obs[[it]]$mlt * 15
    lat1 <- obs[[it]]$mlat
    flux1 <- obs[[it]]$flux
    energy1 <- obs[[it]]$energy

    flux2 <- interp$flux[, , it]
    energy2 <- interp$energy[, , it]

    lon3 <- grnd$mlt[, , it] * 15
    lat3 <- grnd$mlat[, , it]
    flux3 <- grnd$flux[, , it]
    energy3 <- grnd$energy[, , it]

    flux4 <- c(emp$flux[it, 1, , ])
    energy4 <- c(emp$energy[it, 1, , ])

    lon <- list(lon1, lon2, lon3, lon4)
    lat <- list(lat1, lat2, lat3, lat4)
    flux <- list(flux1, flux2, flux3, flux4)
    energy <- list(energy1, energy2, energy3, energy4)

    # remove NA and 0
    for (j in 1: 4) {
        valid <- flux[[j]] > 0 & energy[[j]] > 0
        valid[is.na(valid)] <- FALSE
        lon[[j]] <- lon[[j]][valid]
        lat[[j]] <- lat[[j]][valid]
        flux[[j]] <- flux[[j]][valid]
        energy[[j]] <- energy[[j]][valid]
    }

    # downsampling
    n1 <- length(lon[[1]])
    for (j in 2: 4) {
        nj <- length(lon[[j]])
        len <- as.integer(n1 * ratio[j])

        # two ways of downsampling make no significant difference
        # idx <- as.integer(seq(from = 1, to = nj, length.out = len))
        idx <- sample(1: nj, size = min(len, nj))

        lon[[j]] <- lon[[j]][idx]
        lat[[j]] <- lat[[j]][idx]
        flux[[j]] <- flux[[j]][idx]
        energy[[j]] <- energy[[j]][idx]
    }

    loc <- cbind(unlist(lon), unlist(lat))
    flux <- unlist(flux)
    energy <- unlist(energy)

    # remove duplicate
    valid <- !duplicated(loc)
    loc <- loc[valid, ]
    flux <- flux[valid]
    energy <- energy[valid]

    r <- pi/2 - loc[, 2] * pi/180
    t <- loc[, 1] * pi/180
    # translate observation to auroral center
    loc <- cbind(r*cos(t) - bndry$am[it], r*sin(t) - bndry$bm[it])
    obsflux <- list(loc = loc, val = log(flux), err = rep(1, times = length(flux)))
    obsenergy <- list(loc = loc, val = energy, err = rep(1, times = length(energy)))

    # translate simulation locations to auroral center
    locsim_trans <- cbind(locsim[, 1] - bndry$am[it], locsim[, 2] - bndry$bm[it])


    # basis function setup

    # fill up full domain
    basisxy <- seq(from = -bndry$rmax[it], to = bndry$rmax[it], by = delta)
    nxy <- length(basisxy)
    coor <- expand.grid(basisxy, basisxy)

    # select valid basis functions (within auroral boundaries)
    radius <- array(data = sqrt(coor[, 1]^2 + coor[, 2]^2), dim = c(nxy, nxy))
    valid <- radius >= bndry$rmin[it] & radius <= bndry$rmax[it]

    x <- array(dim = nxy^2)
    y <- array(dim = nxy^2)
    index <- array(dim = c(nxy, nxy))
    # setup index matrix: mark index on valid basis function
    idx <- 1
    for (i in 1: nxy) {
        for (j in 1: nxy) {
            if (valid[i, j]) {
                x[idx] <- basisxy[i]
                y[idx] <- basisxy[j]
                index[i, j] <- idx
                idx <- idx + 1
            }
        }
    }

    m <- idx - 1
    x <- x[1: m]
    y <- y[1: m]
    con <- array(dim = c(m, 4))
    # construct connectivity matrix using index matrix
    idx <- 1
    for (i in 1: nxy) {
        for (j in 1: nxy) {
            if (valid[i, j]) {
                if (i >= 2 && valid[i-1, j]) con[idx, 1] <- index[i-1, j]
                if (j >= 2 && valid[i, j-1]) con[idx, 2] <- index[i, j-1]
                if (i <= nxy-1 && valid[i+1, j]) con[idx, 3] <- index[i+1, j]
                if (j <= nxy-1 && valid[i, j+1]) con[idx, 4] <- index[i, j+1]
                idx <- idx + 1
            }
        }
    }

    basis1 <- list(loc = cbind(x, y), connect = con, centerweight = 4.01, delta = delta*2.5, alpha = 1)
    # combine different levels, single level for now
    basis <- list(basis1)

    cons <- constants(obsflux, basis, normalization, rho, derivative)
    lambda <- exp(optimize(
        function(l) kriging(exp(l), cons$y, cons$W, cons$Z, cons$Q, cons$phi)$likelihood,
        interval = c(-9, 5), maximum = TRUE, tol = 5e-3)$maximum)
    fit <- kriging(lambda, cons$y, cons$W, cons$Z, cons$Q, cons$phi)
    pred <- prediction(locsim_trans, basis, normalization, rho, lambda, cons$Z, cons$Q, cons$phi, fit$M, fit$d, fit$c, fit$rhoMLE)
    fluxsim[it, , ] <- exp(pred$m)

    cons <- constants(obsenergy, basis, normalization, rho, derivative)
    lambda <- exp(optimize(
        function(l) kriging(exp(l), cons$y, cons$W, cons$Z, cons$Q, cons$phi)$likelihood,
        interval = c(-9, 5), maximum = TRUE, tol = 5e-3)$maximum)
    fit <- kriging(lambda, cons$y, cons$W, cons$Z, cons$Q, cons$phi)
    pred <- prediction(locsim_trans, basis, normalization, rho, lambda, cons$Z, cons$Q, cons$phi, fit$M, fit$d, fit$c, fit$rhoMLE)
    energysim[it, , ] <- pred$m
}

dim_time <- ncdim_def(name = "time", units = "hours since 2014-02-20 00:00:00", vals = time)
dim_x <- ncdim_def(name = "x", units = "", vals = xy)
dim_y <- ncdim_def(name = "y", units = "", vals = xy)
var_flux <- ncvar_def(name = "flux", units = "erg/g/cm2", dim = list(dim_time, dim_x, dim_y))
var_energy <- ncvar_def(name = "energy", units = "keV", dim = list(dim_time, dim_x, dim_y))
nc <- nc_create(filename = "auroral_model.nc", vars = list(var_flux, var_energy), force_v4 = TRUE)
ncvar_put(nc = nc, varid = var_flux, vals = fluxsim)
ncvar_put(nc = nc, varid = var_energy, vals = energysim)
nc_close(nc = nc)