source(file = "empirical.r")
source(file = "LatticeKriging.r")

# predict auroral energy flux and mean energy based on satellite and ground observations

hemi <- "north"

# variable to be modeled
varname <- "flux"

# resolution in degree, defined as the spacing between basis functions at each level
res <- 5

# input and output files
input <- paste(hemi, varname, "nz.nc", sep = "_")
output <- paste(hemi, "_", varname, "_out_", res, "d.nc", sep = "")

# whether mean will be calculated
pred_m <- TRUE

# since the amount of auroral data is huge, the calculation of variance using the standard method (current package) will fail
# only Monte Carlo simulation (implemented in the original LatticeKrig package) will successfully produce results
# current package hasn't implemented Monte Carlo simulation, therefore set it to FALSE
pred_se <- FALSE

# error percentage with respect to the value
# this proxy is used since errors are usually not provided in the data
error_ratio <- c(0.15, 0.2, 0.3, 1)

# minimum error to be used in the model
thres_e <- 0.2

# modeling domain, (-range, range)
range <- pi*2/9

# flux threshold for significant activity, only used for calculating post-processing weights
flux_label <- 2

# parameter for KNN, only used for calculating post-processing weights
K <- 10

# read input data, definitions of variables are given in prepare_input.py
src <- c("sat", "grnd", "interp", "emp")
nsrc <- length(src)
nrec <- vector(mode = "list", length = nsrc)
mlat <- vector(mode = "list", length = nsrc)
mlt <- vector(mode = "list", length = nsrc)
variable <- vector(mode = "list", length = nsrc)
nc <- ncdf4::nc_open(filename = input)
time <- ncdf4::ncvar_get(nc = nc, varid = "time")
for (i in 1: nsrc) {
    nrec[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], "nrec", sep = "_"))
    mlat[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], "mlat", sep = "_"))
    mlt[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], "mlt", sep = "_"))
    variable[[i]] <- ncdf4::ncvar_get(nc = nc, varid = paste(src[i], varname, sep = "_"))
}
kp <- ncdf4::ncvar_get(nc = nc, varid = "kp")
ncdf4::nc_close(nc = nc)

# set up the basis functions
res <- res * pi/180
nres <- length(res)
x0 <- vector(mode = "list", length = nres)
connect <- vector(mode = "list", length = nres)
width <- array(dim = nres)
alpha <- array(dim = nres)

for (ires in 1: nres) {
# within each level, the basis functions are equally spaced
    xy <- seq(from = -range[ires], to = range[ires], by = res[ires])
    nxy <- length(xy)
    loc <- array(dim = c(nxy^2, 2))
    con <- array(dim = c(nxy^2, 4))
    for (j in 1: nxy) {
        for (i in 1: nxy) {
            idx <- (j-1)*nxy + i
            loc[idx, 2] <- xy[j]
            loc[idx, 1] <- xy[i]
            if (j >= 2) con[idx, 1] <- (j-2)*nxy + i
            if (i >= 2) con[idx, 2] <- (j-1)*nxy + i-1
            if (j <= nxy-1) con[idx, 3] <- j*nxy + i
            if (i <= nxy-1) con[idx, 4] <- (j-1)*nxy + i+1
        }
    }
    x0[[ires]] <- loc
    connect[[ires]] <- con
# increase the width of basis functions to allow overlapping
    width[ires] <- res[ires] * 2.5
# weighting of each level can be adjusted based on actual problems
    alpha[ires] <- 1 / nres
}

# center weight is usually 2^d+0.01 where d is dimension
centerweight <- 4.01

# initial guesss of the marginal variance, will be re-estimated later
rho <- 1

# number of points to be estimated along each side in the modeling domain
nxy <- 401

# create the grid to be predicted
xy <- seq(from = -range, to = range, length.out = nxy)
coor <- expand.grid(xy, xy)
x <- coor[, 1]
y <- coor[, 2]
xnew <- cbind(x, y)
mlat_new <- 90 - sqrt(x^2 + y^2) * 180/pi
mlt_new <- atan2(y, x) * 12/pi

# save essential modeling parameters for diagnostics
ntime <- length(time)
lambda <- array(dim = ntime)
d <- array(dim = ntime)
rhoMLE <- array(dim = ntime)
like <- array(dim = ntime)
var_m <- array(dim = c(nxy, nxy, ntime))
var_se <- array(dim = c(nxy, nxy, ntime))
weight <- array(dim = c(nxy, nxy, ntime))
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
        err_i[[i]] <- var_i[[i]] * error_ratio[i]
    }

    mlat_i <- unlist(mlat_i)
    mlt_i <- unlist(mlt_i)
    var_i <- unlist(var_i)
    err_i <- unlist(err_i)

# from geomagnetic coordinates to model coordinates
    r <- pi/2 - mlat_i*pi/180
    t <- mlt_i * pi/12
    loc <- cbind(r*cos(t), r*sin(t))

# empirical model as background model
    emp <- get_guvi_kp_model(kp[itime], mlat_i, mlt_i)
    emp_new <- get_guvi_kp_model(kp[itime], mlat_new, mlt_new)

    y <- var_i
    w <- pmax(err_i, thres_e)
    Z <- emp[[varname]]
    Znew <- emp_new[[varname]]

    lk <- LatticeKriging(loc, NULL, y, w, Z,
        x0, width, alpha, connect, centerweight, rho,
        xnew, Znew, pred_m, pred_se)

    lambda[itime] <- lk$lambda
    rhoMLE[itime] <- lk$rhoMLE
    like[itime] <- lk$like
    d[itime] <- lk$d
    if (pred_m) var_m[, , itime] <- array(data = lk$m, dim = c(nxy, nxy))
    if (pred_se) var_se[, , itime] <- array(data = lk$se, dim = c(nxy, nxy))

    if (varname == "flux") {
# post-processing weighting is calculated based on flux
        cl <- class::knn(train = loc, test = xnew, cl = var_i > flux_label, k = K, prob = TRUE)
        prob <- attributes(cl)$prob

# KNN predicts labels for every point, and the label starts at 1. label-1 gives 0/1 labels.
        factor <- as.numeric(cl) - 1

# label*probability is the actual weight
        weight[, , itime] <- array(data = factor*prob+(1-factor)*(1-prob), dim = c(nxy, nxy))
    }
}

# note that the model might give negative values, set to 0
var_m[var_m < 0] <- 0

# save outputs
dim_time <- ncdf4::ncdim_def(name = "time", units = "", vals = time)
dim_x <- ncdf4::ncdim_def(name = "x", units = "", vals = xy)
dim_y <- ncdf4::ncdim_def(name = "y", units = "", vals = xy)

lambda_out <- ncdf4::ncvar_def(name = "lambda", units = "", dim = dim_time)
rho_out <- ncdf4::ncvar_def(name = "rho", units = "", dim = dim_time)
like_out <- ncdf4::ncvar_def(name = "like", units = "", dim = dim_time)
d_out <- ncdf4::ncvar_def(name = "d", units = "", dim = dim_time)

# these are mandatory outputs
vars_out <- list(lambda_out, rho_out, like_out, d_out)

# output mean, variance and weight if flagged
if (pred_m) {
    m_out <- ncdf4::ncvar_def(name = varname, units = "", dim = list(dim_x, dim_y, dim_time))
    vars_out <- c(vars_out, list(m_out))
}
if (pred_se) {
    se_out <- ncdf4::ncvar_def(name = paste(varname, "se", sep = "_"), units = "", dim = list(dim_x, dim_y, dim_time))
    vars_out <- c(vars_out, list(se_out))
}
if (varname == "flux") {
    weight_out <- ncdf4::ncvar_def(name = "weight", units = "", dim = list(dim_x, dim_y, dim_time))
    vars_out <- c(vars_out, list(weight_out))
}

nc <- ncdf4::nc_create(filename = output, vars = vars_out, force_v4 = TRUE)
ncdf4::ncvar_put(nc = nc, varid = m_out, vals = var_m)
if (pred_m) ncdf4::ncvar_put(nc = nc, varid = m_out, vals = var_m)
if (pred_se) ncdf4::ncvar_put(nc = nc, varid = se_out, vals = var_se)
if (varname == "flux") ncdf4::ncvar_put(nc = nc, varid = weight_out, vals = weight)
ncdf4::nc_close(nc = nc)