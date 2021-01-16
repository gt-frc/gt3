#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import e
from GT3.IOL.Functions.CalcVSep import calc_vsep
from scipy.special import gammaincc



def calc_iol_maxwellian(z, m, param, thetapts, Tprofile, coslist, numcos):
    """Calculates eps_min for IOL of species treated with a truncated maxwellian."""

    v_sep, v_sep_min = calc_vsep(z, m, param)

    # Define the Maxwellian for every point in the plasma based on its temperature

    T_matrix = np.zeros(v_sep_min.shape)
    for indx, row in enumerate(T_matrix):
        T_matrix[indx, :] = Tprofile[indx]

    # Create the launch angle matrix (used in the calculation of M_orb)
    zeta_matrix = np.zeros(v_sep_min.shape)
    for indx, column in enumerate(zeta_matrix.T):
        zeta_matrix[:, indx] = coslist[indx]

    eps_min = m * v_sep_min**2 / (2.*T_matrix*1E3*e)

    # F_orb calculation
    # note: the use of gammaincc renders the denominator in Dr. Stacey's equations obsolete.
    integrand_f = gammaincc(3./2., eps_min)
    F_orb_1D = np.sum(integrand_f, axis=1)*(2./numcos)/2.
    F_orb_1D = np.nan_to_num(F_orb_1D)
    F_orb = np.repeat(F_orb_1D.reshape(-1, 1), thetapts, axis=1)

    # M_orb calculation
    integrand_m = zeta_matrix*gammaincc(2, eps_min)
    M_orb_1D = np.sum(integrand_m, axis=1)*(2./numcos)/2.
    M_orb_1D = np.nan_to_num(M_orb_1D)
    M_orb = np.repeat(M_orb_1D.reshape(-1, 1), thetapts, axis=1)

    # E_orb calculation
    integrand_e = gammaincc(5./2, eps_min)
    E_orb_1D = np.sum(integrand_e, axis=1)*(2./numcos)/2.
    E_orb_1D = np.nan_to_num(E_orb_1D)
    E_orb = np.repeat(E_orb_1D.reshape(-1, 1), thetapts, axis=1)

    return F_orb, M_orb, E_orb
