from glob import glob
from netCDF4 import Dataset
from numpy import ndarray, logical_or, nan
from datetime import date


# combine orbits from single satellite for easy comparison of time evolutions
# input: SSUSI auroral EDR data grouped by satellites
# output: combined SSUSI data from one single satellite


Year = 2014
Month = 2
Day = 20
code = 16

files = glob(pathname='f{0:d}/PS.APL_V0105S027C?00??_SC.U_DI.A_GP.F{0:d}-SSUSI_PA.APL-EDR-AURORA_DD.????????_SN.?????-0?_DF.NC'.format(code))
files.sort()
norbit = len(files)

reference_doy = date(year=Year, month=Month, day=Day).timetuple().tm_yday

ssusi_in = Dataset(filename=files[0])
nmlt = ssusi_in.dimensions['N_GEOMAGNETIC_LONGITUDE'].size
nmlat = ssusi_in.dimensions['N_GEOMAGNETIC_LATITUDE'].size
mlt = ssusi_in['MLT_GRID_MAP'][:].filled()
mlat = ssusi_in['LATITUDE_GEOMAGNETIC_GRID_MAP'][:].filled()
flux_thres = ssusi_in['ELECTRON_ENERGY_FLUX_THRESHOLDS'][0].filled()
ssusi_in.close()

ut_n = ndarray(shape=(norbit, nmlt, nmlat))
ut_s = ndarray(shape=(norbit, nmlt, nmlat))
flux_n = ndarray(shape=(norbit, nmlt, nmlat))
flux_s = ndarray(shape=(norbit, nmlt, nmlat))
energy_n = ndarray(shape=(norbit, nmlt, nmlat))
energy_s = ndarray(shape=(norbit, nmlt, nmlat))

for iorbit in range(norbit):
    filename = files[iorbit]
    print(filename)

    ssusi_in = Dataset(filename=filename)
    starting_time = ssusi_in.STARTING_TIME
    stopping_time = ssusi_in.STOPPING_TIME
    ut_n_0 = ssusi_in['UT_N'][:].filled()
    ut_s_0 = ssusi_in['UT_S'][:].filled()
    flux_n_0 = ssusi_in['ENERGY_FLUX_NORTH_MAP'][:].filled()
    flux_s_0 = ssusi_in['ENERGY_FLUX_SOUTH_MAP'][:].filled()
    energy_n_0 = ssusi_in['ELECTRON_MEAN_NORTH_ENERGY_MAP'][:].filled()
    energy_s_0 = ssusi_in['ELECTRON_MEAN_SOUTH_ENERGY_MAP'][:].filled()
    ssusi_in.close()

# stop time doesn't pass the current day
    doy = int(stopping_time[4:7])
    day_inc = doy - reference_doy

# valid UT is non-zero, reset to nan
    invalid = ut_n_0 == 0
    ut_n_0[invalid] = nan
    flux_n_0[invalid] = nan
    energy_n_0[invalid] = nan

    invalid = ut_s_0 == 0
    ut_s_0[invalid] = nan
    flux_s_0[invalid] = nan
    energy_s_0[invalid] = nan

# if the current orbit starts at the previous day, subtract 24 hours to be consistent with other orbits
    if int(starting_time[4:7]) == doy-1:
        ut_n_0[ut_n_0 > 22] -= 24
        ut_s_0[ut_s_0 > 22] -= 24

    ut_n[iorbit, :, :] = ut_n_0 + day_inc*24
    ut_s[iorbit, :, :] = ut_s_0 + day_inc*24
    flux_n[iorbit, :, :] = flux_n_0
    flux_s[iorbit, :, :] = flux_s_0
    energy_n[iorbit, :, :] = energy_n_0
    energy_s[iorbit, :, :] = energy_s_0

ssusi_out = Dataset(filename='f{:d}_{:4d}{:02d}{:02d}.nc'.format(code, Year, Month, Day), mode='w')
ssusi_out.createDimension(dimname='orbit', size=norbit)
ssusi_out.createDimension(dimname='nmlt', size=nmlt)
ssusi_out.createDimension(dimname='nmlat', size=nmlat)
flux_thres_out = ssusi_out.createVariable(varname='flux_thres', datatype='f4')
mlt_out = ssusi_out.createVariable(varname='mlt', datatype='f4', dimensions=('nmlt', 'nmlat'))
mlat_out = ssusi_out.createVariable(varname='mlat', datatype='f4', dimensions=('nmlt', 'nmlat'))
ut_n_out = ssusi_out.createVariable(varname='ut_n', datatype='f4', dimensions=('orbit', 'nmlt', 'nmlat'))
ut_s_out = ssusi_out.createVariable(varname='ut_s', datatype='f4', dimensions=('orbit', 'nmlt', 'nmlat'))
flux_n_out = ssusi_out.createVariable(varname='flux_n', datatype='f4', dimensions=('orbit', 'nmlt', 'nmlat'))
flux_s_out = ssusi_out.createVariable(varname='flux_s', datatype='f4', dimensions=('orbit', 'nmlt', 'nmlat'))
energy_n_out = ssusi_out.createVariable(varname='energy_n', datatype='f4', dimensions=('orbit', 'nmlt', 'nmlat'))
energy_s_out = ssusi_out.createVariable(varname='energy_s', datatype='f4', dimensions=('orbit', 'nmlt', 'nmlat'))

timeunits = 'hours since {:4d}-{:02d}-{:02d} 00:00:00'.format(Year, Month, Day)
ut_n_out.units = timeunits
ut_s_out.units = timeunits
flux_thres_out[:] = flux_thres
mlt_out[:] = mlt
mlat_out[:] = mlat
ut_n_out[:] = ut_n
ut_s_out[:] = ut_s
flux_n_out[:] = flux_n
flux_s_out[:] = flux_s
energy_n_out[:] = energy_n
energy_s_out[:] = energy_s

ssusi_out.close()