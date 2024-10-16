#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from GT3.IOL.Functions.CalcVSep import calc_vsep

def calc_iol_mono_en(z, m, param, thetapts, v_mono, coslist, numcos, *args, **kwargs):
    """calculates IOL for a monoenergetic species with isotropic launch angle (i.e. fast alphas)"""

    if kwargs.get("iol_override"):
        if kwargs.get("iol_override").get("vflip"):
            v_sep, v_sep_min, _, _, _ = calc_vsep(z, m, param, alternate=True)
        else:
            v_sep, v_sep_min, _, _, _ = calc_vsep(z, m, param)
    else:
        v_sep, v_sep_min, _, _, _ = calc_vsep(z, m, param)

    # Create the launch angle matrix (used in the calculation of M_orb)
    zeta_matrix = np.zeros(v_sep_min.shape)
    for indx, column in enumerate(zeta_matrix.T):
        zeta_matrix[:, indx] = coslist[indx]

    # F_orb calculation
    #integrand_f = np.where(v_sep_min <= v_mono, 1.0, 0)
    integrand_f = np.heaviside(v_mono - np.abs(v_sep_min), 1)
    F_orb_1D = np.sum(integrand_f, axis=1) * (2 / numcos) / 2
    F_orb_1D = np.nan_to_num(F_orb_1D)
    F_orb = np.repeat(F_orb_1D.reshape(-1, 1), thetapts, axis=1)

    # M_orb calculation
    #integrand_m = zeta_matrix * np.where(v_sep_min <= v_mono, 1, 0)
    integrand_m = zeta_matrix * np.heaviside(v_mono - np.abs(v_sep_min), 1)
    M_orb_1D = np.sum(integrand_m, axis=1) * (2 / numcos) / 2
    M_orb_1D = np.nan_to_num(M_orb_1D)
    M_orb = np.repeat(M_orb_1D.reshape(-1, 1), thetapts, axis=1)

    # E_orb calculation
    # E_orb is mathematically identical to F_orb for monoenergetic, isotropic launch angle species
    E_orb = F_orb

    return F_orb, M_orb, E_orb
