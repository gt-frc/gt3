#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy import constants
from scipy.interpolate import UnivariateSpline
from GT3.Core.Functions.ProfileClasses import PressureProfiles, TemperatureProfiles

e = constants.elementary_charge


def calc_chi_jet(T: TemperatureProfiles, p: PressureProfiles, a, q, B_T, m_i, rho):
    def calc_bohm(T: TemperatureProfiles, p: PressureProfiles, a, q, B_T, m_i, rho):
        """"""
        cs = np.sqrt(T.i.J / m_i)  # sound speed
        rho_s = cs * m_i / (e * B_T)  # gyroradius

        Te_interp = UnivariateSpline(rho[:, 0], T.e.J.to1D(), k=3, s=0)
        delta_Te = (Te_interp(0.8) - Te_interp(1.0)) / Te_interp(1.0)
        chi_bohm = rho_s * cs * q ** 2 * a * p.e.L * delta_Te
        return chi_bohm

    def calc_gyro_bohm(T, B_T, m_i):
        cs = np.sqrt(T.i.J / m_i)  # sound speed
        rho_s = cs * m_i / (e * B_T)  # gyroradius
        chi_gyro_bohm = rho_s ** 2 * cs * T.e.J.L
        return chi_gyro_bohm

    chi_bohm = calc_bohm(T, p, a, q, B_T, m_i, rho)
    chi_gyro_bohm = calc_gyro_bohm(T, B_T, m_i)
    chi_i_jet = 1.6E-4 * chi_bohm + 1.75E-2 * chi_gyro_bohm
    chi_e_jet = 8E-5 * chi_bohm + 3.5E-2 * chi_gyro_bohm
    return chi_i_jet, chi_e_jet
