from numpy import arange, loadtxt, concatenate, ndarray, logical_not, logical_and, nan, nanmean, full, \
    isnan, count_nonzero, argsort, argmin, abs, interp, pi, meshgrid, linspace, rad2deg, sqrt, arctan2
from netCDF4 import Dataset
from guvi_auroral_model import get_guvi_kp_model


caseid = 20140220

# ut to simulate (hour)
step = 1/6
time = arange(start=0, stop=24+step, step=step)
ntime = len(time)

# satellite code
code = [16, 17, 18]
# satellite data files
file_sat = '../ssusi/f{:d}_{:d}.nc'

# time interval to gather satellite data
coverage = 1/6

# max time interval to perform interpolation
max_interval = 2

# ground data files
file_grnd = '../themis/themis{:d}.nc'.format(caseid)

# Kp for empirical model
data = loadtxt(fname='../omni2.lst')
dt2 = (data[::3, 1]-data[0, 1])*24 + data[::3, 2]
kp2 = data[::3, 3]/10
kp = interp(x=time, xp=dt2, fp=kp2)

# read satellite data
data_in = Dataset(filename=file_sat.format(code[0], caseid))
mlat_sat = data_in['mlat'][:].filled()
mlt_sat = data_in['mlt'][:].filled()
ut_sat = data_in['ut_n'][:].filled()
flux_sat = data_in['flux_n'][:].filled()
energy_sat = data_in['energy_n'][:].filled()
data_in.close()
for c in code[1:]:
    data_in = Dataset(filename=file_sat.format(c, caseid))
    ut_sat = concatenate((ut_sat, data_in['ut_n'][:].filled()))
    flux_sat = concatenate((flux_sat, data_in['flux_n'][:].filled()))
    energy_sat = concatenate((energy_sat, data_in['energy_n'][:].filled()))
    data_in.close()

# read ground data
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
        valid = logical_not(isnan(ut_sat[:, imlt, imlat]))
        if count_nonzero(valid) <= 1:
            continue
        ut = ut_sat[valid, imlt, imlat]
        flux = flux_sat[valid, imlt, imlat]
        energy = energy_sat[valid, imlt, imlat]
        ind = argsort(ut)
        ut = ut[ind]
        flux = flux[ind]
        energy = energy[ind]
        nt = len(ut)
        for itime in range(ntime):
            if time[itime] <= ut[0] or time[itime] >= ut[nt-1]:
                continue
            idx = argmin(abs(ut - time[itime]))
            if ut[idx] > time[itime]:
                idx -= 1
            if ut[idx+1] - ut[idx] > max_interval:
                continue
            flux_interp[itime, imlt, imlat] = interp(x=time[itime], xp=ut[idx: idx+1], fp=flux[idx: idx+1])
            energy_interp[itime, imlt, imlat] = interp(x=time[itime], xp=ut[idx: idx+1], fp=energy[idx: idx+1])

# select ground data at given time
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

# create empirical data at a similar resolution
nxy_emp = mlat_sat.shape[0]
xy_emp = linspace(start=-pi*2/9, stop=pi*2/9, num=nxy_emp)
x_emp, y_emp = meshgrid(xy_emp, xy_emp)
mlat_emp = 90 - rad2deg(sqrt(x_emp**2 + y_emp**2))
mlt_emp = arctan2(y_emp, x_emp) * 12/pi
flux_emp = ndarray(shape=(ntime, nxy_emp, nxy_emp))
energy_emp = ndarray(shape=(ntime, nxy_emp, nxy_emp))
for itime in range(ntime):
    flux, energy = get_guvi_kp_model(kp[itime], mlat_emp.flatten(), mlt_emp.flatten())
    flux_emp[itime, :, :] = flux.reshape((nxy_emp, nxy_emp))
    energy_emp[itime, :, :] = energy.reshape((nxy_emp, nxy_emp))

data_out = Dataset(filename='input.nc', mode='w')
data_out.createDimension(dimname='time', size=ntime)
data_out.createDimension(dimname='sat_nmlat', size=nmlat_sat)
data_out.createDimension(dimname='sat_nmlt', size=nmlt_sat)
data_out.createDimension(dimname='grnd_nmlat', size=nmlat_grnd)
data_out.createDimension(dimname='grnd_nmlt', size=nmlt_grnd)
data_out.createDimension(dimname='emp_x', size=nxy_emp)
data_out.createDimension(dimname='emp_y', size=nxy_emp)
time_out = data_out.createVariable(varname='time', datatype='d', dimensions='time')
kp_out = data_out.createVariable(varname='kp', datatype='d', dimensions='time')
sat_mlat_out = data_out.createVariable(varname='sat_mlat', datatype='d', dimensions=('sat_nmlt', 'sat_nmlat'))
sat_mlt_out = data_out.createVariable(varname='sat_mlt', datatype='d', dimensions=('sat_nmlt', 'sat_nmlat'))
sat_flux_out = data_out.createVariable(varname='sat_flux', datatype='d', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
sat_energy_out = data_out.createVariable(varname='sat_energy', datatype='d', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
grnd_mlat_out = data_out.createVariable(varname='grnd_mlat', datatype='d', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
grnd_mlt_out = data_out.createVariable(varname='grnd_mlt', datatype='d', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
grnd_flux_out = data_out.createVariable(varname='grnd_flux', datatype='d', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
grnd_energy_out = data_out.createVariable(varname='grnd_energy', datatype='d', dimensions=('time', 'grnd_nmlat', 'grnd_nmlt'))
interp_mlat_out = data_out.createVariable(varname='interp_mlat', datatype='d', dimensions=('sat_nmlt', 'sat_nmlat'))
interp_mlt_out = data_out.createVariable(varname='interp_mlt', datatype='d', dimensions=('sat_nmlt', 'sat_nmlat'))
interp_flux_out = data_out.createVariable(varname='interp_flux', datatype='d', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
interp_energy_out = data_out.createVariable(varname='interp_energy', datatype='d', dimensions=('time', 'sat_nmlt', 'sat_nmlat'))
emp_mlat_out = data_out.createVariable(varname='emp_mlat', datatype='d', dimensions=('emp_y', 'emp_x'))
emp_mlt_out = data_out.createVariable(varname='emp_mlt', datatype='d', dimensions=('emp_y', 'emp_x'))
emp_flux_out = data_out.createVariable(varname='emp_flux', datatype='d', dimensions=('time', 'emp_y', 'emp_x'))
emp_energy_out = data_out.createVariable(varname='emp_energy', datatype='d', dimensions=('time', 'emp_y', 'emp_x'))
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