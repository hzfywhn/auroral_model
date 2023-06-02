from netCDF4 import Dataset
from numpy import ndarray


# remove zeros in the original data to better assist spatial modeling
# input: downsampled input data of flux or energy
# output: flux or energy with zeros removed


hemi = 'north'

# flux or energy
varname = 'flux'

src = ['sat', 'grnd', 'interp', 'emp']
nsrc = len(src)

data_in = Dataset(filename=hemi+'_'+varname+'_ds.nc')
maxnrec = data_in.dimensions['nrec'].size
time = data_in['time'][:].filled()

ntime = len(time)
nrec = ndarray(shape=(nsrc, ntime), dtype=int)
mlat = ndarray(shape=(nsrc, ntime, maxnrec))
mlt = ndarray(shape=(nsrc, ntime, maxnrec))
var = ndarray(shape=(nsrc, ntime, maxnrec))

kp = data_in['kp'][:].filled()
for isrc in range(nsrc):
    nrec[isrc, :] = data_in[src[isrc]+'_nrec'][:]
    mlat[isrc, :, :] = data_in[src[isrc]+'_mlat'][:]
    mlt[isrc, :, :] = data_in[src[isrc]+'_mlt'][:]
    var[isrc, :, :] = data_in[src[isrc]+'_'+varname][:]
data_in.close()

nrec_nz = ndarray(shape=(nsrc, ntime), dtype=int)
mlat_nz = ndarray(shape=(nsrc, ntime, maxnrec))
mlt_nz = ndarray(shape=(nsrc, ntime, maxnrec))
var_nz = ndarray(shape=(nsrc, ntime, maxnrec))
for isrc in range(nsrc):
    for itime in range(ntime):
        n = nrec[isrc, itime]
        mlat_i = mlat[isrc, itime, 0: n]
        mlt_i = mlt[isrc, itime, 0: n]
        var_i = var[isrc, itime, 0: n]

        valid = var_i > 0
        mlat_i = mlat_i[valid]
        mlt_i = mlt_i[valid]
        var_i = var_i[valid]

        cnt = len(var_i)
        nrec_nz[isrc, itime] = cnt
        mlat_nz[isrc, itime, 0: cnt] = mlat_i
        mlt_nz[isrc, itime, 0: cnt] = mlt_i
        var_nz[isrc, itime, 0: cnt] = var_i
maxnrec_nz = nrec_nz.max()

data_out = Dataset(filename=hemi+'_'+varname+'_nz.nc', mode='w')
data_out.createDimension(dimname='time', size=ntime)
data_out.createDimension(dimname='nrec', size=maxnrec_nz)
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
    data_out[src[isrc]+'_nrec'][:] = nrec_nz[isrc, :]
    data_out[src[isrc]+'_mlat'][:] = mlat_nz[isrc, :, 0: maxnrec_nz]
    data_out[src[isrc]+'_mlt'][:] = mlt_nz[isrc, :, 0: maxnrec_nz]
    data_out[src[isrc]+'_'+varname][:] = var_nz[isrc, :, 0: maxnrec_nz]
data_out.close()