#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline
from scipy import constants
import numpy as np

e = constants.elementary_charge

def calc_svfus(T, mode='dd'):
    def sigv(T, mode):  # function takes T in kev
        if mode == 'dt':
            B_G = 34.3827
            m_rc2 = 1124656

            C1 = 1.17302E-9
            C2 = 1.51361E-2
            C3 = 7.51886E-2
            C4 = 4.60643E-3
            C5 = 1.35000E-2
            C6 = -1.06750E-4
            C7 = 1.36600E-5

            theta = T / (1.0 - (T * (C2 + T * (C4 + T * C6))) / (1.0 + T * (C3 + T * (C5 + T * C7))))
            xi = (B_G ** 2.0 / (4.0 * theta)) ** (1.0 / 3.0)
            sigv = C1 * theta * np.sqrt(xi / (m_rc2 * T ** 3.0)) * np.exp(-3.0 * xi)
            sigv = sigv / 1.0E6  # convert from cm^3/s to m^3/s

        elif mode == 'dd':

            B_G = 31.3970
            m_rc2 = 937814

            # first for the D(d, p)T reaction
            C1_1 = 5.65718E-12
            C2_1 = 3.41267E-3
            C3_1 = 1.99167E-3
            C4_1 = 0.0
            C5_1 = 1.05060E-5
            C6_1 = 0.0
            C7_1 = 0.0

            theta_1 = T / (
                        1.0 - (T * (C2_1 + T * (C4_1 + T * C6_1))) / (1.0 + T * (C3_1 + T * (C5_1 + T * C7_1))))
            xi_1 = (B_G ** 2.0 / (4.0 * theta_1)) ** (1.0 / 3.0)
            sigv_1 = C1_1 * theta_1 * np.sqrt(xi_1 / (m_rc2 * T ** 3.0)) * np.exp(-3.0 * xi_1)

            # then for the D(d, n)He3 reaction

            C1_2 = 5.43360E-12
            C2_2 = 5.85778E-3
            C3_2 = 7.68222E-3
            C4_2 = 0.0
            C5_2 = -2.96400E-6
            C6_2 = 0.0
            C7_2 = 0.0

            theta_2 = T / (
                        1.0 - (T * (C2_2 + T * (C4_2 + T * C6_2))) / (1.0 + T * (C3_2 + T * (C5_2 + T * C7_2))))
            xi_2 = (B_G ** 2.0 / (4.0 * theta_2)) ** (1.0 / 3.0)
            sigv_2 = C1_2 * theta_2 * np.sqrt(xi_2 / (m_rc2 * T ** 3.0)) * np.exp(-3.0 * xi_2)

            sigv = (0.5 * sigv_1 + 0.5 * sigv_2) / 1.0E6  # convert from cm^3/s to m^3/s
        elif mode == 'dHe3':
            B_G = 68.7508
            m_rc2 = 1124572

            C1 = 5.51036E-10
            C2 = 6.41918E-3
            C3 = -2.02896E-3
            C4 = -1.91080E-5
            C5 = 1.35776E-4
            C6 = 0
            C7 = 0

            theta = T / (1.0 - (T * (C2 + T * (C4 + T * C6))) / (1.0 + T * (C3 + T * (C5 + T * C7))))
            xi = (B_G ** 2.0 / (4.0 * theta)) ** (1.0 / 3.0)
            sigv = C1 * theta * np.sqrt(xi / (m_rc2 * T ** 3.0)) * np.exp(-3.0 * xi)
            sigv = sigv / 1.0E6  # convert from cm^3/s to m^3/s
        return sigv

    # create logspace over the relevant temperature range
    # (bosch hale technically only valid over 0.2 - 100 kev)
    Ti_range = np.logspace(-1, 2, 1000)  # values in kev
    sigv_fus_range = sigv(Ti_range, mode=mode)  # in m^3/s
    sigv_fus_interp = UnivariateSpline(Ti_range * 1.0E3 * e, sigv_fus_range, s=0)  # converted to Joules
    sv_fus = sigv_fus_interp(T.i.J)
    dsv_fus_dT = sigv_fus_interp.derivative()(T.i.J)

    return sv_fus, dsv_fus_dT
