#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata


def calc_grad(quant, psi, R, Z):
    psi_vals = psi.psi_norm[:,0]
    vals = quant[:,0]


    # For some reason, psi_vals will sometimes not be monotonically increasing, especially near the magnetic axis.
    # This prevents UnivariateSpline from working. To prevent this, we're just going to delete any points that are
    # lower than the previous point, prior to doing the fit. As long as the psi data isn't too wonky, this should be
    # fine.
    psivals_mi = []  # psivals that are monotonically increasing
    vals_mi = []  # surf_area corresponding to monotonically increasing psi
    for i, psi_val in enumerate(psi_vals):
        if i == 0:
            psivals_mi.append(0)
            vals_mi.append(0)
        elif psi_val > psivals_mi[-1]:
            psivals_mi.append(psi_vals[i])
            vals_mi.append(vals[i])


    # calculate values as function of psi and get psi derivative function
    #psi_fit = UnivariateSpline(psivals_mi, vals_mi, k=3, s=0)
    psi_fit = UnivariateSpline(psi_vals, vals, k=3, s=0)
    d_dpsi_fit = psi_fit.derivative()

    # get value of dval_dpsi on the main computational grid
    dval_dpsi = d_dpsi_fit(psi.psi)

    # calculate dpsi_norm_dr everywhere on the main computational grid
    dpsi_dr = griddata(np.column_stack((psi.psi_data.R.flatten(),
                                        psi.psi_data.Z.flatten())),
                       psi.psi_data.dpsidr.flatten(),
                       (R,Z),
                       method='linear')

    dval_dr = dval_dpsi * dpsi_dr

    return dval_dr