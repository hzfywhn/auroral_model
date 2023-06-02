from netCDF4 import Dataset
from numpy import size, ndarray, logical_not, isnan, arange
from random import sample


# random downsampling of auroral data from different sources
# input: full input data containing flux and energy
# output: downsampled flux or energy


hemi = 'north'

# flux or energy
varname = 'flux'

# designed downsampling ratio for each data source
ratio = [1, 1, 1/3, 1/20]

src = ['sat', 'grnd', 'interp', 'emp']
nsrc = len(src)

mlat = [None] * nsrc
mlt = [None] * nsrc
var = [None] * nsrc

data_in = Dataset(filename=hemi+'_in.nc')
time = data_in['time'][:].filled()
kp = data_in['kp'][:].filled()
for isrc in range(nsrc):
    mlat[isrc] = data_in[src[isrc]+'_mlat'][:].filled()
    mlt[isrc] = data_in[src[isrc]+'_mlt'][:].filled()
    var[isrc] = data_in[src[isrc]+'_'+varname][:].filled()
data_in.close()

ntime = len(time)
maxnrec = size(mlat[1][0, :, :])
for isrc in range(nsrc):
    if isrc == 1:
        continue
    if maxnrec < size(mlat[isrc]):
        maxnrec = size(mlat[isrc])

nrec_ds = ndarray(shape=(nsrc, ntime), dtype=int)
mlat_ds = ndarray(shape=(nsrc, ntime, maxnrec))
mlt_ds = ndarray(shape=(nsrc, ntime, maxnrec))
var_ds = ndarray(shape=(nsrc, ntime, maxnrec))

for itime in range(ntime):
    for isrc in range(nsrc):
        if isrc == 1:
# note that ground data has a different format (varying mlat/mlt with time), dealt separately
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

# randomly choose nsample from cnt
        indices = sample(population=list(arange(start=0, stop=cnt, dtype=int)), k=nsample)
        mlat_ds[isrc, itime, 0: nsample] = mlat_i[indices]
        mlt_ds[isrc, itime, 0: nsample] = mlt_i[indices]
        var_ds[isrc, itime, 0: nsample] = var_i[indices]
maxnrec_ds = nrec_ds.max()

data_out = Dataset(filename=hemi+'_'+varname+'_ds.nc', mode='w')
data_out.createDimension(dimname='time', size=ntime)
data_out.createDimension(dimname='nrec', size=maxnrec_ds)
time_out = data_out.createVariable(varname='time', datatype='f4', dimensions='time')
kp_out = data_out.createVariable(varname='kp', datatype='f4', dimensions='time')
for isrc in range(nsrc):
    data_out.createVariable(varname=src[isrc]+'_nrec', datatype='i4', dimensions='time')
    data_out.createVariable(varname=src[isrc]+'_mlat', datatype='f4', dimensions=('time', 'nrec'))
    data_out.createVariable(varname=src[isrc]+'_mlt', datatype='f4', dimensions=('time', 'nrec'))
    data_out.createVariable(varname=src[isrc]+'_'+varname, datatype='f4', dimensions=('time', 'nrec'))
time_out[:] = time
kp_out[:] = kp
for isrc in range(nsrc):
    data_out[src[isrc]+'_nrec'][:] = nrec_ds[isrc, :]
    data_out[src[isrc]+'_mlat'][:] = mlat_ds[isrc, :, 0: maxnrec_ds]
    data_out[src[isrc]+'_mlt'][:] = mlt_ds[isrc, :, 0: maxnrec_ds]
    data_out[src[isrc]+'_'+varname][:] = var_ds[isrc, :, 0: maxnrec_ds]
data_out.close()