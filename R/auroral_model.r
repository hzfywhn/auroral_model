library(package = "ncdf4")
library(package = "abind")
library(package = "class")
library(package = "alphahull")
# library(package = "ggplot2")

source(file = "grid_satellite.r")
source(file = "interp_satellite.r")
source(file = "ground.r")
source(file = "empirical.r")
source(file = "preprocess_empirical.r")
source(file = "auroral_boundary.r")
source(file = "setup_basis.r")
source(file = "LatticeKriging.r")

# adjustable parameters

doy <- 51
hemi <- "N"

# ut to simulate (hour)
time <- seq(from = 6, to = 12, by = 1/3)

# satellite data files
filename_sat <- c("../f16_20140220.nc", "../f17_20140220.nc", "../f18_20140220.nc")

# time interval to gather satellite data to one time point (hour)
coverage <- 1/4

# ground data files
filename_grnd <- "../themis20140220.nc"

# select only diffuse aurora in empirical data
aurtype <- 1

# simplified dFdt calculation
f <- read.table(file = "../omni2.lst")
dFdt <- approx(x = (f[, 2] - doy) * 24 + f[, 3], y = f[, 6]^(4/3) * sqrt(f[, 4]^2 + f[, 5]^2)^(2/3) * abs(sin(atan2(f[, 4], f[, 5])/2))^(8/3), xout = time)$y

# directory to hold coefficients
premodel <- "../premodel"

# interpolate empirical data to approximate equidistance in mlt and mlat
mlt_emp <- seq(from = 0, to = 23.9, by = 0.1)
mlat_emp <- seq(from = 50, to = 89, by = 1)

# ratio of samples among satellite, interpolated satellite, ground and empirical data
ratio <- c(1, 1/3, 1/3, 1/9)

# downsampling method
downsampling <- "random"

# parameter for determining the region with observation
alphaRadius <- 0.2

# parameters for auroral boundaries
lb <- 1
a <- seq(from = -pi/18, to = pi/18, by = pi/180)
b <- seq(from = -pi/18, to = pi/18, by = pi/180)
r <- seq(from = pi/18, to = 4*pi/18, by = pi/180)
lev_flux <- c(0.5, 1)
lev_energy <- c(1, 5)
KN <- 10

# modeling domain
mlt_sim <- seq(from = 0, to = 23.9, by = 0.1)
mlat_sim <- seq(from = 50, to = 89, by = 1)
coor_sim <- expand.grid(mlt_sim * pi/12, pi/2 - mlat_sim*pi/180)
loc_sim <- cbind(coor_sim[, 2] * cos(coor_sim[, 1]), coor_sim[, 2] * sin(coor_sim[, 1]))

min_data_points <- 50

# parameters for Lattice Kriging
sponge <- pi/18
delta <- pi/180
centerweight <- 4.01
overlap <- delta * 2.5
weight <- 1
normalization <- FALSE
rho <- 1
derivative <- FALSE

filename_output <- "auroral_model.nc"

# read satellite data
nc <- nc_open(filename = filename_sat[1])
mlat_sat <- ncvar_get(nc = nc, varid = "mlat")
mlt_sat <- ncvar_get(nc = nc, varid = "mlt")
ut_sat <- ncvar_get(nc = nc, varid = "ut_n")
flux_sat <- ncvar_get(nc = nc, varid = "flux_n")
energy_sat <- ncvar_get(nc = nc, varid = "energy_n")
nc_close(nc = nc)
for (ifile in 2: length(filename_sat)) {
    nc <- nc_open(filename = filename_sat[ifile])
    ut_sat <- abind(ut_sat, ncvar_get(nc = nc, varid = "ut_n"))
    flux_sat <- abind(flux_sat, ncvar_get(nc = nc, varid = "flux_n"))
    energy_sat <- abind(energy_sat, ncvar_get(nc = nc, varid = "energy_n"))
    nc_close(nc = nc)
}

# gather satellite data within certain time period
sat <- grid_satellite(mlat_sat, mlt_sat, ut_sat, flux_sat, energy_sat, time, coverage)

# interpolate satellite data to given time
interp <- interp_satellite(ut_sat, flux_sat, energy_sat, time)
interp$mlt <- mlt_sat
interp$mlat <- mlat_sat

# read ground data
nc <- nc_open(filename = filename_grnd)
time_grnd <- ncvar_get(nc = nc, varid = "time") / 3600
mlat_grnd <- ncvar_get(nc = nc, varid = "mlat")
mlt_grnd <- ncvar_get(nc = nc, varid = "mlt")
flux_grnd <- ncvar_get(nc = nc, varid = "flux")
energy_grnd <- ncvar_get(nc = nc, varid = "energy")
nc_close(nc = nc)

# select ground data at given time
grnd <- ground(time_grnd, mlat_grnd, mlt_grnd, flux_grnd, energy_grnd, time)

# calculate empirical data with given dFdt
emp <- empirical(doy, dFdt, hemi, premodel)

# regrid empirical data to new grids and remove values in the common region with observations
emp <- preprocess_empirical(emp, aurtype, mlat_emp, mlt_emp, interp, grnd, alphaRadius)

# coordinates for boundary calculation
r_sat <- pi/2 - mlat_sat*pi/180
t_sat <- mlt_sat * pi/12
x_sat <- r_sat * cos(t_sat)
y_sat <- r_sat * sin(t_sat)

nt <- length(time)
nmlt_sim <- length(mlt_sim)
nmlat_sim <- length(mlat_sim)

# probabilistic boundary
loc_sat <- cbind(c(x_sat), c(y_sat))
nlev_flux <- length(lev_flux)
nlev_energy <- length(lev_energy)
prob_flux <- array(dim = c(nt, nmlt_sim*nmlat_sim))
prob_energy <- array(dim = c(nt, nmlt_sim*nmlat_sim))
for (it in 1: nt) {
    flux_interp <- c(interp$flux[, , it])
    valid <- flux_interp > 0
    valid[is.na(valid)] <- FALSE
    flux_interp <- flux_interp[valid]
    label <- array(dim = length(flux_interp))
    label[flux_interp < lev_flux[1]] <- 0
    for (ilev in 1: (nlev_flux-1)) {
        label[flux_interp >= lev_flux[ilev] & flux_interp < lev_flux[ilev+1]] <- ilev / nlev_flux
    }
    label[flux_interp >= lev_flux[nlev_flux]] <- 1
    prob <- knn(train = loc_sat[valid, ], test = loc_sim, cl = label, k = KN, prob = TRUE)
    prob_flux[it, ] <- attr(prob, "prob") #* as.numeric(levels(prob))[prob]

    energy_interp <- c(interp$energy[, , it])
    valid <- energy_interp > 0
    valid[is.na(valid)] <- FALSE
    energy_interp <- energy_interp[valid]
    label <- array(dim = length(energy_interp))
    label[energy_interp < lev_energy[1]] <- 0
    for (ilev in 1: (nlev_energy-1)) {
        label[energy_interp >= lev_energy[ilev] & energy_interp < lev_energy[ilev+1]] <- ilev / nlev_energy
    }
    label[energy_interp >= lev_energy[nlev_energy]] <- 1
    prob <- knn(train = loc_sat[valid, ], test = loc_sim, cl = label, k = KN, prob = TRUE)
    prob_energy[it, ] <- attr(prob, "prob") #* as.numeric(levels(prob))[prob]
}

# deterministic boundary: center (bndry$am, bndry$bm), boundary (bndry$rmin, bndry$rmax)
bndry <- auroral_boundary(x_sat, y_sat, interp$flux, lb, a, b, r)

# setup basis based on deterministic boundaries
basis <- setup_basis(bndry, sponge, delta, centerweight, overlap, weight)

nsim <- nmlt_sim * nmlat_sim
flux_sim <- array(dim = c(nt, nmlt_sim, nmlat_sim))
energy_sim <- array(dim = c(nt, nmlt_sim, nmlat_sim))

for (it in 1: nt) {
    print(it)
    mlt <- list(sat[[it]]$mlt, interp$mlt, grnd$mlt[, , it], emp[[it]]$mlt)
    mlat <- list(sat[[it]]$mlat, interp$mlat, grnd$mlat[, , it], emp[[it]]$mlat)
    flux <- list(sat[[it]]$flux, interp$flux[, , it], grnd$flux[, , it], emp[[it]]$flux)
    energy <- list(sat[[it]]$energy, interp$energy[, , it], grnd$energy[, , it], emp[[it]]$energy)

    # downsampling
    n1 <- length(mlt[[1]])
    for (j in 2: 4) {
        nj <- length(mlt[[j]])
        len <- as.integer(n1 * ratio[j])

        # downsampling
        idx <- switch (downsampling, sequential = as.integer(seq(from = 1, to = nj, length.out = len)), random = sample(1: nj, size = min(len, nj)))
        mlt[[j]] <- mlt[[j]][idx]
        mlat[[j]] <- mlat[[j]][idx]
        flux[[j]] <- flux[[j]][idx]
        energy[[j]] <- energy[[j]][idx]
    }

    # remove NA and 0
    for (j in 1: 4) {
        valid <- flux[[j]] > 0 & energy[[j]] > 0
        valid[is.na(valid)] <- FALSE
        mlt[[j]] <- mlt[[j]][valid]
        mlat[[j]] <- mlat[[j]][valid]
        flux[[j]] <- flux[[j]][valid]
        energy[[j]] <- energy[[j]][valid]
    }

    loc <- cbind(unlist(mlt), unlist(mlat))
    flux <- unlist(flux)
    energy <- unlist(energy)

    # remove duplicate
    valid <- !duplicated(loc)
    loc <- loc[valid, ]
    flux <- flux[valid]
    energy <- energy[valid]

    if (nrow(loc) <= min_data_points) next

    # reformat for simulation
    r <- pi/2 - loc[, 2]*pi/180
    t <- loc[, 1] * pi/12
    loc <- cbind(r*cos(t), r*sin(t))
    Z1 <- rep(1, times = nsim)
    MR1 <- combineMR(list(loc = loc_sim), basis[[it]], normalization, rho, derivative)

    n <- length(flux)
    yfit <- log(flux)
    W <- diag.spam(x = 1, nrow = n, ncol = n)
    Z <- rep(1, times = n)
    MR <- combineMR(list(loc = loc), basis[[it]], normalization, rho, derivative)
    lambda <- exp(optimize(
        function(l) kriging(exp(l), yfit, W, Z, MR$Q, MR$phi)$likelihood,
        interval = c(-9, 5), maximum = TRUE, tol = 5e-3)$maximum)
    fit <- kriging(lambda, yfit, W, Z, MR$Q, MR$phi)
    pred <- prediction(Z1, MR1$phi1, lambda, Z, MR$Q, MR$phi, fit$M, fit$d, fit$c, fit$rhoMLE)

    # ggplot(data = data.frame(x = loc[, 1], y = loc[, 2], c = flux$val)) + geom_point(mapping = aes(x, y, colour = c))
    # ggplot(data = data.frame(x = loc_sim[, 1], y = loc_sim[, 2], c = pred$m)) + geom_point(mapping = aes(x, y, colour = c))
    flux_sim[it, , ] <- array(data = exp(pred$m) * prob_flux[it, ], dim = c(nmlt_sim, nmlat_sim))

    n <- length(energy)
    yfit <- log(energy)
    W <- diag.spam(x = 1, nrow = n, ncol = n)
    Z <- rep(1, times = n)
    MR <- combineMR(list(loc = loc), basis[[it]], normalization, rho, derivative)
    lambda <- exp(optimize(
        function(l) kriging(exp(l), yfit, W, Z, MR$Q, MR$phi)$likelihood,
        interval = c(-9, 5), maximum = TRUE, tol = 5e-3)$maximum)
    fit <- kriging(lambda, yfit, W, Z, MR$Q, MR$phi)
    pred <- prediction(Z1, MR1$phi1, lambda, Z, MR$Q, MR$phi, fit$M, fit$d, fit$c, fit$rhoMLE)

    # ggplot(data = data.frame(x = loc[, 1], y = loc[, 2], c = energy$val)) + geom_point(mapping = aes(x, y, colour = c))
    # ggplot(data = data.frame(x = loc_sim[, 1], y = loc_sim[, 2], c = pred$m)) + geom_point(mapping = aes(x, y, colour = c))
    energy_sim[it, , ] <- array(data = exp(pred$m) * prob_energy[it, ], dim = c(nmlt_sim, nmlat_sim))
}

dim_time <- ncdim_def(name = "time", units = "", vals = time)
dim_mlt <- ncdim_def(name = "mlt", units = "", vals = mlt_sim)
dim_mlat <- ncdim_def(name = "mlat", units = "", vals = mlat_sim)
var_flux <- ncvar_def(name = "flux", units = "", dim = list(dim_time, dim_mlt, dim_mlat))
var_energy <- ncvar_def(name = "energy", units = "", dim = list(dim_time, dim_mlt, dim_mlat))
nc <- nc_create(filename = filename_output, vars = list(var_flux, var_energy), force_v4 = TRUE)
ncvar_put(nc = nc, varid = var_flux, vals = flux_sim)
ncvar_put(nc = nc, varid = var_energy, vals = energy_sim)
nc_close(nc = nc)