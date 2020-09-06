library(package = "ncdf4")
library(package = "LatticeKrig")

source(file = "grid_satellite.r")
source(file = "interp_satellite.r")
source(file = "ground.r")
source(file = "empirical.r")

# adjustable parameters
time <- seq(from = 0, to = 24, by = 4)
ratio <- c(1, 1/6, 1/6, 1/36)
nxy <- 101
xy <- seq(from = -pi*2/9, to = pi*2/9, length.out = nxy)

nc <- nc_open(filename = "ssusi/ssusi_combine_20140219_20140221.nc")
mlat_sat <- ncvar_get(nc = nc, varid = "mlat")
mlt_sat <- ncvar_get(nc = nc, varid = "mlt")
ut_sat <- ncvar_get(nc = nc, varid = "ut_n") - 24
flux_sat <- ncvar_get(nc = nc, varid = "flux_n")
energy_sat <- ncvar_get(nc = nc, varid = "energy_n")
nc_close(nc = nc)

nc <- nc_open(filename = "themis/themis2014022006_mag.nc")
time_grnd <- ncvar_get(nc = nc, varid = "time") / 60 + 6
mlat_grnd <- ncvar_get(nc = nc, varid = "mlat")
mlt_grnd <- ncvar_get(nc = nc, varid = "mlt")
flux_grnd <- ncvar_get(nc = nc, varid = "flux")
energy_grnd <- ncvar_get(nc = nc, varid = "energy") * 1e-3
nc_close(nc = nc)

obs <- grid_satellite(mlat = mlat_sat, mlt = mlt_sat, ut = ut_sat, flux = flux_sat, energy = energy_sat, time = time)

interp <- interp_satellite(ut = ut_sat, flux = flux_sat, energy = energy_sat, time = time)
r2 <- pi/2 - interp$mlat * pi/180
t2 <- interp$mlt * pi/12
x2 <- r2 * cos(t2)
y2 <- r2 * sin(t2)

grnd <- ground(timein = time_grnd, mlat = mlat_grnd, mlt = mlt_grnd, flux = flux_grnd, energy = energy_grnd, timeout = time)

# simplified dFdt calculation
f <- read.table(file = "ovation/omni2.lst")
hour <- (f$V2 - 51) * 24 + f$V3
dFdt <- f$V6^(4/3) * sqrt(f$V4^2 + f$V5^2)^(2/3) * abs(sin(atan2(f$V4, f$V5)/2))^(8/3)
dFdt <- approx(x = hour, y = dFdt, xout = time)$y
emp <- empirical(doy = 51, dFdt, hemi = "N", premodel = "ovation/ovation-prime_2.3.0/premodel")
coor <- expand.grid(emp$mlt * pi/12, pi/2 - emp$mlat * pi/180)
t4 <- coor$Var1
r4 <- coor$Var2
x4 <- r4 * cos(t4)
y4 <- r4 * sin(t4)

# deleted probabilistic boundary generation due to unnecessity and high cost
# probably add back later in a low cost way

flux <- array(dim = c(nt, nxy, nxy))
energy <- array(dim = c(nt, nxy, nxy))
for (i in 1: length(time)) {
    r1 <- pi/2 - obs[[i]]$mlat * pi/180
    t1 <- obs[[i]]$mlt * pi/12
    x1 <- r1 * cos(t1)
    y1 <- r1 * sin(t1)
    f1 <- obs[[i]]$flux
    e1 <- obs[[i]]$energy

    f2 <- interp$flux[, , i]
    e2 <- interp$energy[, , i]

    r3 <- pi/2 - grnd$mlat[, , i] * pi/180
    t3 <- grnd$mlt[, , i] * pi/12
    x3 <- r3 * cos(t3)
    y3 <- r3 * sin(t3)
    f3 <- grnd$flux[, , i]
    e3 <- grnd$energy[, , i]

    f4 <- c(emp$flux[i, 1, , ])
    e4 <- c(emp$energy[i, 1, , ])

    x <- list(x1, x2, x3, x4)
    y <- list(y1, y2, y3, x4)
    f <- list(f1, f2, f3, x4)
    e <- list(e1, e2, e3, x4)

    # remove NA and 0
    for (j in 1: 4) {
        valid <- f[[j]] > 0 & e[[j]] > 0
        valid[is.na(valid)] <- FALSE
        x[[j]] <- x[[j]][valid]
        y[[j]] <- y[[j]][valid]
        f[[j]] <- f[[j]][valid]
        e[[j]] <- e[[j]][valid]
    }

    # downsampling
    n1 <- length(x[[1]])
    for (j in 2: 4) {
        nj <- length(x[[j]])
        len <- as.integer(n1 * ratio[j])

        # two ways of downsampling make no significant difference
        # idx <- as.integer(seq(from = 1, to = nj, length.out = len))
        idx <- sample(1: nj, size = min(len, nj))

        x[[j]] <- x[[j]][idx]
        y[[j]] <- y[[j]][idx]
        f[[j]] <- f[[j]][idx]
        e[[j]] <- e[[j]][idx]
    }

    xyfe <- data.frame(x = unlist(x), y = unlist(y), f = unlist(f), e = unlist(e))

    # remove duplicate
    xyfe <- xyfe[!duplicated(xyfe[c("x", "y")]), ]

    modelf <- LatticeKrig(x = xyfe[c("x", "y")], y = log(xyfe["f"]), nlevel = 1, NC = 50)
    predf <- predictSurface(object = modelf, grid.list = list(xy, xy))
    flux[i, , ] <- exp(predf$z)

    modele <- LatticeKrig(x = xyfe[c("x", "y")], y = xyfe["e"], nlevel = 1, NC = 50)
    prede <- predictSurface(object = modele, grid.list = list(xy, xy))
    energy[i, , ] <- prede$z
}

dim_t <- ncdim_def(name = "time", units = "hours since 2014-02-20 00:00:00", vals = time)
dim_x <- ncdim_def(name = "x", units = "", vals = xy)
dim_y <- ncdim_def(name = "y", units = "", vals = xy)
var_f <- ncvar_def(name = "flux", units = "erg/g/cm2", dim = list(dim_t, dim_x, dim_y))
var_e <- ncvar_def(name = "energy", units = "keV", dim = list(dim_t, dim_x, dim_y))
nc <- nc_create(filename = "auroral_model.nc", vars = c(var_f, var_e), force_v4 = TRUE)
ncvar_put(nc = nc, varid = var_f, vals = flux)
ncvar_put(nc = nc, varid = var_e, vals = energy)
nc_close(nc = nc)
