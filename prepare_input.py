from numpy import arange, loadtxt, concatenate, ndarray, logical_not, logical_and, nan, nanmean, full, \
    isnan, count_nonzero, argsort, argmin, abs, interp, pi, meshgrid, linspace, rad2deg, sqrt, arctan2
from netCDF4 import Dataset
from guvi_auroral_model import get_guvi_kp_model


# combine auroral data from satellite and ground, interpolate satellite data, and create model data
# input: SSUSI data with orbits combined, THEMIS data, and 3 hourly kp indices
# output: auroral data to be used in downsampling


hemi = 'north'
caseid = 20140220

# time to simulate (hour)
step = 1/6
time = arange(start=0, stop=24+step, step=step)
ntime = len(time)

# satellite code and data file patterns
code = [16, 17, 18]
file_sat = 'f{:d}_{:d}.nc'

# time interval (hour) to gather satellite data
coverage = 1/6

# max time interval (hour) to perform interpolation
max_interval = 2

# ground data files
file_grnd = 'themis{:d}.nc'.format(caseid)

# Kp for empirical model
data = loadtxt(fname='omni2.lst')
dt2 = (data[::3, 1]-data[0, 1])*24 + data[::3, 2]
kp2 = data[::3, 3]/10

# read satellite data, note that data from one single satellite has been combined
data_in = Dataset(filename=file_sat.format(code[0], caseid))
mlat_sat = data_in['mlat'][:].filled()
mlt_sat = data_in['mlt'][:].filled()
ut_sat = data_in['ut_'+hemi[0]][:].filled()
flux_sat = data_in['flux_'+hemi[0]][:].filled()
energy_sat = data_in['energy_'+hemi[0]][:].filled()
data_in.close()
for c in code[1:]:
    data_in = Dataset(filename=file_sat.format(c, caseid))
    ut_sat = concatenate((ut_sat, data_in['ut_'+hemi[0]][:].filled()))
    flux_sat = concatenate((flux_sat, data_in['flux_'+hemi[0]][:].filled()))
    energy_sat = concatenate((energy_sat, data_in['energy_'+hemi[0]][:].filled()))
    data_in.close()

# read ground data, only northern hemisphere has THEMIS data
if hemi == 'north':
    data_in = Dataset(filename=file_grnd)
    time_grnd = data_in['time'][:].filled() / 3600
    mlat_grnd = data_in['mlat'][:].filled()
    mlt_grnd = data_in['mlt'][:].filled()
    flux_grnd = data_in['flux'][:].filled()
    energy_grnd = data_in['energy'][:].filled()
    data_in.close()

# gather satellite data within certain time period
nmlt_sat, nmlat_sat = mlat_sat.shape
flux_grid = ndarray(shape=(ntime, nmlt_sat, nmlat_sat))
energy_grid = ndarray(shape=(ntime, nmlt_sat, nmlat_sat))
for itime in range(ntime):
    invalid = logical_not(logical_and(ut_sat > time[itime]-coverage, ut_sat < time[itime]+coverage))

# this format with nan kept is easy for plotting purposes
    flux = flux_sat.copy()
    flux[invalid] = nan
    flux_grid[itime, :, :] = nanmean(flux, axis=0)

    energy = energy_sat.copy()
    energy[invalid] = nan
    energy_grid[itime, :, :] = nanmean(energy, axis=0)

# interpolate satellite data to given time
flux_interp = full(shape=(ntime, nmlt_sat, nmlat_sat), fill_value=nan)
energy_interp = full(shape=(ntime, nmlt_sat, nmlat_sat), fill_value=nan)
for imlt in range(nmlt_sat):
    for imlat in range(nmlat_sat):
# remove filling values
        valid = logical_not(isnan(ut_sat[:, imlt, imlat]))
        if count_nonzero(valid) <= 1:
            continue
        ut = ut_sat[valid, imlt, imlat]
        flux = flux_sat[valid, imlt, imlat]
        energy = energy_sat[valid, imlt, imlat]

# ut is not in order by now, sort it for interpolation
        ind = argsort(ut)
        ut = ut[ind]
        flux = flux[ind]
        energy = energy[ind]

        nt = len(ut)
        for itime in range(ntime):
# no extrapolation
            if time[itime] <= ut[0] or time[itime] >= ut[nt-1]:
                continue

            idx = argmin(abs(ut - time[itime]))
            if ut[idx] > time[itime]:
                idx -= 1

# if the consecutive orbits are seperated too far, skip the interpolation
            if ut[idx+1] - ut[idx] > max_interval:
                continue

            flux_interp[itime, imlt, imlat] = interp(x=time[itime], xp=ut[idx: idx+1], fp=flux[idx: idx+1])
            energy_interp[itime, imlt, imlat] = interp(x=time[itime], xp=ut[idx: idx+1], fp=energy[idx: idx+1])

# select ground data at given time
# ground data has different mlat/mlt at different time
if hemi == 'north':
    _, nmlat_grnd, nmlt_grnd = mlat_grnd.shape
    dt_grnd = time_grnd[1] - time_grnd[0]
    mlat_grnd_s = full(shape=(ntime, nmlat_grnd, nmlt_grnd), fill_value=nan)
    mlt_grnd_s = full(shape=(ntime, nmlat_grnd, nmlt_grnd), fill_value=nan)
    flux_grnd_s = full(shape=(ntime, nmlat_grnd, nmlt_grnd), fill_value=nan)
    energy_grnd_s = full(shape=(ntime, nmlat_grnd, nmlt_grnd), fill_value=nan)
    for itime in range(ntime):
        idx = argmin(abs(time_grnd - time[itime]))
        if abs(time_grnd[idx] - time[itime]) > dt_grnd:
            continue
        mlat_grnd_s[itime, :, :] = mlat_grnd[idx, :, :]
        mlt_grnd_s[itime, :, :] = mlt_grnd[idx, :, :]
        flux_grnd_s[itime, :, :] = flux_grnd[idx, :, :]
        energy_grnd_s[itime, :, :] = energy_grnd[idx, :, :]
else:
    nmlat_grnd = 1
    nmlt_grnd = 1
    mlat_grnd_s = 90
    mlt_grnd_s = 0
    flux_grnd_s = 0
    energy_grnd_s = 0

# create empirical data at a similar resolution
nxy_emp = mlat_sat.shape[0]
xy_emp = linspace(start=-pi*2/9, stop=pi*2/9, num=nxy_emp)
x_emp, y_emp = meshgrid(xy_emp, xy_emp)
mlat_emp = 90 - rad2deg(sqrt(x_emp**2 + y_emp**2))
mlt_emp = arctan2(y_emp, x_emp) * 12/pi

# the empirical model executes slowly, so calculate model data at a coarse time interval
ndt2 = len(dt2)
flux2 = ndarray(shape=(ndt2, nxy_emp, nxy_emp))
energy2 = ndarray(shape=(ndt2, nxy_emp, nxy_emp))
for itime in range(ndt2):
    flux, energy = get_guvi_kp_model(kp2[itime], mlat_emp.flatten(), mlt_emp.flatten())
    flux2[itime, :, :] = flux.reshape((nxy_emp, nxy_emp))
    energy2[itime, :, :] = energy.reshape((nxy_emp, nxy_emp))

# interpolate empirical data to the selected time
kp = interp(x=time, xp=dt2, fp=kp2)
flux_emp = ndarray(shape=(ntime, nxy_emp, nxy_emp))
energy_emp = ndarray(shape=(ntime, nxy_emp, nxy_emp))
for iy in range(nxy_emp):
    for ix in range(nxy_emp):
        flux_emp[:, iy, ix] = interp(x=time, xp=dt2, fp=flux2[:, iy, ix])
        energy_emp[:, iy, ix] = interp(x=time, xp=dt2, fp=energy2[:, iy, ix])

data_out = Dataset(filename=hemi+'_in.nc', mode='w')
data_out.createDimension(dimname='time', size=ntime)
data_out.createDimension(dimname='sat_nmlat', size=nmlat_sat)
data_out.createDimension(dimname='sat_nmlt', size=nmlt_sat)
data_out.createDimension(dimname='grnd_nmlat', size=nmlat_grnd)
data_out.createDimension(dimname='grnd_nmlt', size=nmlt_grnd)
data_out.createDimension(dimname='emp_x', size=nxy_emp)
data_out.createDimension(dimname='emp_y', size=nxy_emp)
time_out = data_out.createVariable(varname='time', datatype='f4', dimensions='time')
kp_out = data_out.createVariable(varname='kp', datatype='f4', dimensions='time')
sat_mlat_out = data_out.createVariable(varname='sat_mlat', datatype='f4', dimensions=('sat_nmlt', 'sat_nmlat'))
sat_mlt_out = data_out.createVariable(varname='sat_mlt', datatype='f4', dimensions=('sat_nmlt', 'sat_nmlat'))
sat_flux_out = data_out.createVariable(varname='sat_flux', datatype='f4', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
sat_energy_out = data_out.createVariable(varname='sat_energy', datatype='f4', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
grnd_mlat_out = data_out.createVariable(varname='grnd_mlat', datatype='f4', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
grnd_mlt_out = data_out.createVariable(varname='grnd_mlt', datatype='f4', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
grnd_flux_out = data_out.createVariable(varname='grnd_flux', datatype='f4', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
grnd_energy_out = data_out.createVariable(varname='grnd_energy', datatype='f4', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
interp_mlat_out = data_out.createVariable(varname='interp_mlat', datatype='f4', dimensions=('sat_nmlt', 'sat_nmlat'))
interp_mlt_out = data_out.createVariable(varname='interp_mlt', datatype='f4', dimensions=('sat_nmlt', 'sat_nmlat'))
interp_flux_out = data_out.createVariable(varname='interp_flux', datatype='f4', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
interp_energy_out = data_out.createVariable(varname='interp_energy', datatype='f4', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
emp_mlat_out = data_out.createVariable(varname='emp_mlat', datatype='f4', dimensions=('emp_y', 'emp_x'))
emp_mlt_out = data_out.createVariable(varname='emp_mlt', datatype='f4', dimensions=('emp_y', 'emp_x'))
emp_flux_out = data_out.createVariable(varname='emp_flux', datatype='f4', dimensions=('time', 'emp_y', 'emp_x'))
emp_energy_out = data_out.createVariable(varname='emp_energy', datatype='f4', dimensions=('time', 'emp_y', 'emp_x'))
time_out[:] = time
kp_out[:] = kp
sat_mlat_out[:] = mlat_sat
sat_mlt_out[:] = mlt_sat
sat_flux_out[:] = flux_grid
sat_energy_out[:] = energy_grid
grnd_mlat_out[:] = mlat_grnd_s
grnd_mlt_out[:] = mlt_grnd_s
grnd_flux_out[:] = flux_grnd_s
grnd_energy_out[:] = energy_grnd_s
interp_mlat_out[:] = mlat_sat
interp_mlt_out[:] = mlt_sat
interp_flux_out[:] = flux_interp
interp_energy_out[:] = energy_interp
emp_mlat_out[:] = mlat_emp
emp_mlt_out[:] = mlt_emp
emp_flux_out[:] = flux_emp
emp_energy_out[:] = energy_emp
data_out.close()