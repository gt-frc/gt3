#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy import constants
from scipy.interpolate import interp2d
from scipy.interpolate import griddata


e = constants.elementary_charge


def calc_svel_st(T):
    tint = np.array([-1, 0, 1, 2, 3])
    tnnt = np.array([0, 1, 2])

    elast = np.array([[-1.3569E+01, -1.3337E+01, -1.3036E+01, -1.3569E+01, -1.3337E+01],
                      [-1.3036E+01, -1.3337E+01, -1.3167E+01, -1.3046E+01, -1.3036E+01],
                      [-1.3046E+01, -1.2796E+01, -1.3036E+01, -1.3046E+01, -1.2796E+01]])

    interp1 = interp2d(tint, tnnt, elast)

    Ti_exps = np.linspace(-1, 3, 100)
    Tn_exps = np.linspace(0, 2, 100)
    svel_vals = 10.0 ** (interp1(Ti_exps, Tn_exps))  # in m^3/s

    Ti_vals = np.logspace(-1, 3, 100) * e  # in joules
    Tn_vals = np.logspace(0, 2, 100) * e  # in joules

    dsvel_dTi_vals = np.gradient(svel_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * e, T.i.ev * e)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * e

    sv_el = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                          svel_vals.flatten(),
                          (Ti_mod, Tn_mod),
                          method='linear', rescale=False)
    dsv_el_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                              dsvel_dTi_vals.flatten(),
                              (Ti_mod, Tn_mod),
                              method='linear', rescale=False)

    return sv_el, dsv_el_dT