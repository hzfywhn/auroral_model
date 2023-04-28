from numpy import exp, linspace, cos, sin, ndarray, sum, argmin, abs, pi
from netCDF4 import Dataset


# a rewrite of Zhang and Paxton (2008) auroral model in python


def get_HP(kp):
    A = [16.8244, 0.323365, -4.86128]
    B = [1.82336, 0.613192, 26.1798]
    if kp <= 6.9:
        HP = A[0] * exp(A[1] * kp) + A[2]
    else:
        HP = B[0] * exp(B[1] * kp) + B[2]
    return HP


def get_fitting_value(t, r, AA, BB):
    ncoeff, nkp = AA.shape
    ang = linspace(start=0, stop=nkp-1, num=nkp) * t
    cosang = cos(ang)
    sinang = sin(ang)

    coeff = ndarray(shape=ncoeff)
    for i in range(ncoeff):
        coeff[i] = sum(AA[i, :]*cosang + BB[i, :]*sinang)

    F1 = exp((r - coeff[1]) / coeff[2])
    F2 = 1 + exp((r - coeff[1]) / coeff[3])
    value = coeff[0] * F1 / F2**2

    return value


def get_guvi_kp_model(kp, mlat, mlt):
# input: kp, mlat, mlt (of same length)
# output: energy flux and mean energy

# guvi_auroral_model.nc contains paramters of the auroral model
    data = Dataset(filename='guvi_auroral_model.nc')
    kp_list = data['kp'][:].filled()
    AA_energy = data['AA_energy'][:].filled()
    BB_energy = data['BB_energy'][:].filled()
    AA_flux = data['AA_flux'][:].filled()
    BB_flux = data['BB_flux'][:].filled()
    data.close()

# simplified kp index calculation
    ida = argmin(abs(kp_list - kp))
    if kp_list[ida] > kp:
        ida -= 1
    idb = ida + 1

    x0 = (kp_list[idb] - kp) / (kp_list[idb] - kp_list[ida])
    x1 = (kp - kp_list[ida]) / (kp_list[idb] - kp_list[ida])

    HPa = get_HP(kp_list[ida])
    HPb = get_HP(kp_list[idb])
    HP = get_HP(kp)

    xa = (HPb - HP) / (HPb - HPa)
    xb = (HP - HPa) / (HPb - HPa)

    n = len(mlat)
    flux = ndarray(shape=n)
    energy = ndarray(shape=n)

    for i in range(n):
        r = 90 - mlat[i]
        t = mlt[i] * pi/12
        Va = get_fitting_value(t, r, AA_flux[:, ida, :], BB_flux[:, ida, :])
        Vb = get_fitting_value(t, r, AA_flux[:, idb, :], BB_flux[:, idb, :])
        flux[i] = Va*xa + Vb*xb
        Va = get_fitting_value(t, r, AA_energy[:, ida, :], BB_energy[:, ida, :])
        Vb = get_fitting_value(t, r, AA_energy[:, idb, :], BB_energy[:, idb, :])
        energy[i] = Va*x0 + Vb*x1

    return flux, energy
