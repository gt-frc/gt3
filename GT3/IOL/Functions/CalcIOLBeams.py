#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp1d
from GT3.IOL.Functions.CalcVSep import calc_vsep


def calc_iol_beams(z, m, param, thetapts, v_mono, zeta_beam, coslist):
    """calculates IOL for a monoenergetic species with a single known launch angle (i.e. beam ions)
    :param z: The Z value of the beam ions
    :type z: float
    :param m: The mass of the beam ions in kg
    :type m: float
    :param param: The IOL parameters
    :param thetapts: The number of theta points
    :type thetapts: int
    :param v_mono: The particle velocity in m/s
    :type v_mono: float
    :param zeta_beam: The zeta value
    :param coslist:
    :return:
    """

    v_sep, v_sep_min = calc_vsep(z, m, param)

    # Obtain v_sep_min(zeta_beam) for each rho value
    v_sep_min_zeta = np.zeros(len(v_sep_min))

    for i,v in enumerate(v_sep_min):
        # create interpolation function of v_sep_min(rho) vs zeta
        v_sep_min_zeta_interp = interp1d(coslist, v, fill_value='extrapolate')
        # get v_sep_min(zeta_beam)
        v_sep_min_zeta[i] = v_sep_min_zeta_interp(zeta_beam)

    # F_orb calculation
    F_orb_1D = np.heaviside(v_mono - v_sep_min_zeta, 1)
    F_orb_1D = np.nan_to_num(F_orb_1D)
    F_orb = np.repeat(F_orb_1D.reshape(-1, 1), thetapts, axis=1)

    # M_orb calculation
    M_orb_1D = zeta_beam * np.heaviside(v_mono - v_sep_min_zeta, 1)
    M_orb_1D = np.nan_to_num(M_orb_1D)
    M_orb = np.repeat(M_orb_1D.reshape(-1, 1), thetapts, axis=1)

    # E_orb calculation
    # E_orb is mathematically identical to F_orb for monoenergetic, monodirectional species
    E_orb = F_orb

    return F_orb, M_orb, E_orb
