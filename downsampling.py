from netCDF4 import Dataset
from numpy import size, ndarray, logical_not, isnan, arange, max
from random import sample


varname = 'flux'
ratio = [1, 1, 1/3, 1/20]

src = ['sat', 'grnd', 'interp', 'emp']
nsrc = len(src)

mlat = [None] * nsrc
mlt = [None] * nsrc
var = [None] * nsrc

data_in = Dataset(filename='input.nc')
time = data_in['time'][:].filled()
kp = data_in['kp'][:].filled()
for isrc in range(nsrc):
    mlat[isrc] = data_in[src[isrc]+'_mlat'][:].filled()
    mlt[isrc] = data_in[src[isrc]+'_mlt'][:].filled()
    var[isrc] = data_in[src[isrc]+'_'+varname][:].filled()
data_in.close()

ntime = len(time)
max_nrec = size(mlat[1][0, :, :])
for isrc in range(nsrc):
    if isrc == 1:
        continue
    if max_nrec < size(mlat[isrc]):
        max_nrec = size(mlat[isrc])

nrec_ds = ndarray(shape=(nsrc, ntime), dtype=int)
mlat_ds = ndarray(shape=(nsrc, ntime, max_nrec))
mlt_ds = ndarray(shape=(nsrc, ntime, max_nrec))
var_ds = ndarray(shape=(nsrc, ntime, max_nrec))

for itime in range(ntime):
    for isrc in range(nsrc):
        if isrc == 1:
            mlat_i = mlat[isrc][itime, :, :]
            mlt_i = mlt[isrc][itime, :, :]
        else:
            mlat_i = mlat[isrc]
            mlt_i = mlt[isrc]
        var_i = var[isrc][itime, :, :]

        valid = logical_not(isnan(var_i))
        mlat_i = mlat_i[valid]
        mlt_i = mlt_i[valid]
        var_i = var_i[valid]

        cnt = len(var_i)
        nsample = int(cnt * ratio[isrc])
        nrec_ds[isrc, itime] = nsample
        if nsample == 0:
            continue

        indices = sample(population=list(arange(start=0, stop=cnt, dtype=int)), k=nsample)
        mlat_ds[isrc, itime, 0: nsample] = mlat_i[indices]
        mlt_ds[isrc, itime, 0: nsample] = mlt_i[indices]
        var_ds[isrc, itime, 0: nsample] = var_i[indices]
max_nrec_ds = max(nrec_ds)

data_out = Dataset(filename=varname+'_in.nc', mode='w')
data_out.createDimension(dimname='time', size=ntime)
data_out.createDimension(dimname='max_nrec', size=max_nrec_ds)
time_out = data_out.createVariable(varname='time', datatype='d', dimensions='time')
kp_out = data_out.createVariable(varname='kp', datatype='d', dimensions='time')
for isrc in range(nsrc):
    data_out.createVariable(varname=src[isrc]+'_nrec', datatype='i', dimensions='time')
    data_out.createVariable(varname=src[isrc]+'_mlat', datatype='d', dimensions=('time', 'max_nrec'))
    data_out.createVariable(varname=src[isrc]+'_mlt', datatype='d', dimensions=('time', 'max_nrec'))
    data_out.createVariable(varname=src[isrc]+'_'+varname, datatype='d', dimensions=('time', 'max_nrec'))
time_out[:] = time
kp_out[:] = kp
for isrc in range(nsrc):
    data_out[src[isrc]+'_nrec'][:] = nrec_ds[isrc, :]
    data_out[src[isrc]+'_mlat'][:] = mlat_ds[isrc, :, 0: max_nrec_ds]
    data_out[src[isrc]+'_mlt'][:] = mlt_ds[isrc, :, 0: max_nrec_ds]
    data_out[src[isrc]+'_'+varname][:] = var_ds[isrc, :, 0: max_nrec_ds]
data_out.close()