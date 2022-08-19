source(file = "empirical.r")
library(package = "spam64")

# variable to be modeled
varname <- "flux"

# flag for SE calculation
SE <- TRUE

# error percentage
error_ratio <- c(0.15, 0.2, 0.3, 1)

# error threshold
thres_e <- 1e-2

# resolution tag
res <- "high"

# parameters for LatticeKrig
nlevel <- 3
NC <- 25
reps <- 100

# flux threshold for significant activity
flux_label <- 2

# parameter for KNN
K <- 10

# modeling domain
nxy <- 401
xy <- seq(from = -pi*2/9, to = pi*2/9, length.out = nxy)
coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]
loc_new <- cbind(x, y)
mlat_new <- 90 - sqrt(x^2 + y^2) * 180/pi
mlt_new <- atan2(y, x) * 12/pi

src <- c("sat", "grnd", "interp", "emp")
nsrc <- length(src)
nrec <- vector(mode = "list", length = nsrc)
mlat <- vector(mode = "list", length = nsrc)
mlt <- vector(mode = "list", length = nsrc)
variable <- vector(mode = "list", length = nsrc)
nc <- ncdf4::nc_open(filename = paste(varname, "in.nc", sep = "_"))
time <- ncdf4::ncvar_get(nc = nc, varid = "time")
for (i in 1: nsrc) {
    nrec[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], "nrec", sep = "_"))
    mlat[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], "mlat", sep = "_"))
    mlt[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], "mlt", sep = "_"))
    variable[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], varname, sep = "_"))
}
kp <- ncdf4::ncvar_get(nc = nc, varid = "kp")
ncdf4::nc_close(nc = nc)

ntime <- length(time)
weight <- array(dim = c(nxy, nxy, ntime))
var_m <- array(dim = c(nxy, nxy, ntime))
var_se <- array(dim = c(nxy, nxy, ntime))
for (itime in 1: ntime) {
    mlat_i <- vector(mode = "list", length = nsrc)
    mlt_i <- vector(mode = "list", length = nsrc)
    var_i <- vector(mode = "list", length = nsrc)
    err_i <- vector(mode = "list", length = nsrc)
    for (i in 1: nsrc) {
        n <- nrec[[i]][itime]
        if (n == 0) next
        mlat_i[[i]] <- mlat[[i]][1: n, itime]
        mlt_i[[i]] <- mlt[[i]][1: n, itime]
        var_i[[i]] <- variable[[i]][1: n, itime]
        if (i == 2) {
            emp <- get_guvi_kp_model(kp[itime], mlat_i[[2]], mlt_i[[2]])
            valid <- var_i[[2]] > emp[[varname]]
            mlat_i[[2]] <- mlat_i[[2]][valid]
            mlt_i[[2]] <- mlt_i[[2]][valid]
            var_i[[2]] <- var_i[[2]][valid]
        }
        err_i[[i]] <- var_i[[i]] * error_ratio[i]
    }
    err_i[[nsrc]] <- pmax(err_i[[nsrc]], unlist(err_i[1: (nsrc-1)]))
    mlat_i <- unlist(mlat_i)
    mlt_i <- unlist(mlt_i)
    var_i <- unlist(var_i)
    err_i <- unlist(err_i)

    # remove 0
    valid <- var_i > 0
    coor <- cbind(mlat_i[valid], mlt_i[valid])
    var_i <- var_i[valid]
    err_i <- err_i[valid]

    # remove duplicate
    valid <- !duplicated(coor)
    mlat_i <- coor[valid, 1]
    mlt_i <- coor[valid, 2]
    r <- pi/2 - mlat_i*pi/180
    t <- mlt_i * pi/12
    loc <- cbind(r*cos(t), r*sin(t))
    var_i <- var_i[valid]
    err_i <- err_i[valid]

    # empirical model as background model
    emp <- get_guvi_kp_model(kp[itime], mlat_i, mlt_i)
    emp_new <- get_guvi_kp_model(kp[itime], mlat_new, mlt_new)

    y <- log(var_i)
    ye <- pmax(log(err_i), thres_e)
    # ye <- 1
    Z <- log(emp[[varname]])
    Znew <- log(emp_new[[varname]])

    fit <- LatticeKrig::LatticeKrig(x = loc, y = y, weights = ye, Z = Z, nlevel = nlevel, NC = NC)

    m <- LatticeKrig::predict.LKrig(object = fit, xnew = loc_new, Znew = Znew)
    var_m[, , itime] <- exp(m)

    if (SE) {
        # this method breaks down with large amount of input data
        # se <- LatticeKrig::predictSE.LKrig(object = fit, xnew = loc_new, Znew = Znew)
        # var_se[, , itime] <- exp(se)

        sim <- LatticeKrig::LKrig.sim.conditional(LKrigObj = fit, M = reps, x.grid = loc_new, Z.grid = Znew)
        # var_m[, , itime] <- exp(sim$ghat)
        var_se[, , itime] <- exp(sim$SE)
    }

    if (varname == "flux") {
        # post-processing weighting
        cl <- class::knn(train = loc, test = loc_new, cl = var_i > flux_label, k = K, prob = TRUE)
        prob <- attributes(cl)$prob
        factor <- as.numeric(cl) - 1
        weight[, , itime] <- array(data = factor*prob+(1-factor)*(1-prob), dim = c(nxy, nxy))
    }
}

dim_time <- ncdf4::ncdim_def(name = "time", units = "", vals = time)
dim_x <- ncdf4::ncdim_def(name = "x", units = "", vals = xy)
dim_y <- ncdf4::ncdim_def(name = "y", units = "", vals = xy)

m_out <- ncdf4::ncvar_def(name = varname, units = "", dim = list(dim_x, dim_y, dim_time))
se_out <- ncdf4::ncvar_def(name = paste(varname, "se", sep = "_"), units = "", dim = list(dim_x, dim_y, dim_time))
weight_out <- ncdf4::ncvar_def(name = "weight", units = "", dim = list(dim_x, dim_y, dim_time))

if (varname == "flux") {
    if (SE) {
        vars_out <- list(m_out, se_out, weight_out)
    } else {
        vars_out <- list(m_out, weight_out)
    }
} else {
    if (SE) {
        vars_out <- list(m_out, se_out)
    } else {
        vars_out <- m_out
    }
}
nc <- ncdf4::nc_create(filename = paste(varname, "_out_", res, ".nc", sep = ""), vars = vars_out, force_v4 = TRUE)
ncdf4::ncvar_put(nc = nc, varid = m_out, vals = var_m)
if (SE) ncdf4::ncvar_put(nc = nc, varid = se_out, vals = var_se)
if (varname == "flux") ncdf4::ncvar_put(nc = nc, varid = weight_out, vals = weight)
ncdf4::nc_close(nc = nc)
