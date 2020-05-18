#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
from scipy import constants
from scipy.interpolate import griddata

e = constants.elementary_charge

def calc_svcx_st(T):
    tint = np.array([-1, 0, 1, 2, 3])
    tnnt = np.array([0, 1, 2])

    cx = np.array([[-1.4097E+01, -1.3921E+01, -1.3553E+01, -1.4097E+01, -1.3921E+01],
                   [-1.3553E+01, -1.3921E+01, -1.3824E+01, -1.3538E+01, -1.3553E+01],
                   [-1.3538E+01, -1.3432E+01, -1.3553E+01, -1.3538E+01, -1.3432E+01]])

    interp1 = interp2d(tint, tnnt, cx)

    Ti_exps = np.linspace(-1, 3, 100)
    Tn_exps = np.linspace(0, 2, 100)
    svcx_vals = 10.0 ** (interp1(Ti_exps, Tn_exps))  # in m^3/s

    Ti_vals = np.logspace(-1, 3, 100) * e  # in joules
    Tn_vals = np.logspace(0, 2, 100) * e  # in joules

    dsvcx_dTi_vals = np.gradient(svcx_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * e, T.i.ev * e)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * e

    sv_cx = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                          svcx_vals.flatten(),
                          (Ti_mod, Tn_mod),
                          method='linear', rescale=False)

    dsv_cx_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                              dsvcx_dTi_vals.flatten(),
                              (Ti_mod, Tn_mod),
                              method='linear', rescale=False)

    return sv_cx, dsv_cx_dT