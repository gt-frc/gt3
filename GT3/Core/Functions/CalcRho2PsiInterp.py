#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata,  interp1d
from shapely.geometry import Point, LineString


def calc_rho2psi_interp(pts, psi_data, sep_val):
    rho_vals = np.linspace(0, sep_val, 100)

    ptsRZ = np.zeros((len(rho_vals), 2))

    obmp_line = LineString([Point(pts.axis.mag), Point(pts.obmp)])
    for i, rho in enumerate(rho_vals):
        ptsRZ[i] = np.asarray(obmp_line.interpolate(rho, normalized=True).coords)

    psi_vals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                        psi_data.psi.flatten(),
                        ptsRZ,
                        method='cubic')

    psinorm_vals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                            psi_data.psi_norm.flatten(),
                            ptsRZ,
                            method='cubic')

    psinorm_vals[0] = 0
    psinorm_vals[-1] = sep_val

    # For some reason, psi_vals will sometimes not be monotonically increasing, especially near the magnetic axis.
    # This prevents UnivariateSpline from working. To prevent this, we're just going to delete any points that are
    # lower than the previous point, prior to doing the fit. As long as the psi data isn't too wonky, this should be
    # fine.
    psivals_mi = []  # psivals that are monotonically increasing
    psinormvals_mi = []  # psinormvals corresponding to monotonically increasing psi
    rhovals_mi = []  # rhovals corresponding to monotonically increasing psi
    for i, psi_val in enumerate(psi_vals):
        if i == 0:
            rhovals_mi.append(0)
            psivals_mi.append(psi_vals[0])  # this is probably supposed to be zero as well, but we'll leave it for now
            psinormvals_mi.append(0)
        elif psi_val > psivals_mi[-1]:
            rhovals_mi.append(rho_vals[i])
            psivals_mi.append(psi_vals[i])
            psinormvals_mi.append(psinorm_vals[i])

    rhovals_mi = np.asarray(rhovals_mi)
    psivals_mi = np.asarray(psivals_mi)
    psinormvals_mi = np.asarray(psinormvals_mi)

    rho2psi = interp1d(rhovals_mi, psivals_mi, fill_value='extrapolate')
    rho2psinorm = interp1d(rhovals_mi, psinormvals_mi, fill_value='extrapolate')
    psi2rho = interp1d(psivals_mi, rhovals_mi, fill_value='extrapolate')
    psi2psinorm = interp1d(psivals_mi, psinormvals_mi, fill_value='extrapolate')
    psinorm2rho = interp1d(psinormvals_mi, rhovals_mi, fill_value='extrapolate')
    psinorm2psi = interp1d(psinormvals_mi, psivals_mi, fill_value='extrapolate')

    # rho2psi = interp1d(rho_vals, psi_vals, fill_value='extrapolate')
    # psi2rho = interp1d(psi_vals, rho_vals, fill_value='extrapolate')
    return rho2psi, rho2psinorm, psi2rho, psi2psinorm, psinorm2rho, psinorm2psi
