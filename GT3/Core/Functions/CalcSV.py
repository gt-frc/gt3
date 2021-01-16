#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d, UnivariateSpline, griddata
from scipy import constants

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

    Ti_mod = np.where(T.i.ev.val > 1E3, 1E3 * e, T.i.ev * e)
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

    Ti_mod = np.where(T.i.ev.val > 1E3, 1E3 * e, T.i.ev * e)
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
    sv_fus = sigv_fus_interp(T.i.J.val)
    dsv_fus_dT = sigv_fus_interp.derivative()(T.i.J.val)

    return sv_fus, dsv_fus_dT

def calc_svion_st(T):
    # TODO: configure so it can use any of the cross section libraries
    # currently using the Stacey-Thomas cross sections
    T_exps_fit = np.array([-1, 0, 1, 2, 3, 4, 5])
    sigv_exps_fit = np.array([-2.8523E+01, -1.7745E+01, -1.3620E+01,
                              -1.3097E+01, -1.3301E+01, -1.3301E+01, -1.3301E+01])
    interp1 = UnivariateSpline(T_exps_fit, sigv_exps_fit, s=0)

    T_exps_range = np.linspace(-1, 5, 1000)
    sigv_vals_range = 10.0 ** interp1(T_exps_range)  # in m^3/s

    T_vals_range = np.logspace(-1, 5, 1000) * e  # in joules
    interp2 = UnivariateSpline(T_vals_range, sigv_vals_range, s=0)

    sv_ion = interp2(T.i.J.val)
    dsv_ion_dT = interp2.derivative()(T.i.J.val)

    return sv_ion, dsv_ion_dT

def calc_svrec_st(n, T):
    # # TODO: check this calculation. -MH
    # znint = np.array([16, 18, 20, 21, 22])
    # Tint = np.array([-1, 0, 1, 2, 3])
    #
    # rec = np.array([[-1.7523E+01, -1.6745E+01, -1.5155E+01, -1.4222E+01, -1.3301E+01],
    #                 [-1.8409E+01, -1.8398E+01, -1.8398E+01, -1.7886E+01, -1.7000E+01],
    #                 [-1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01],
    #                 [-2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01],
    #                 [-2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01]])
    #
    # interp1 = interp2d(znint, Tint, rec)
    #
    # zni_exps = np.linspace(16, 22, 100)
    # Ti_exps = np.linspace(-1, 3, 100)
    # svrec_vals = 10.0 ** (interp1(zni_exps, Ti_exps))  # in m^3/s
    #
    # zni_vals = np.logspace(16, 22, 100)
    # Ti_vals = np.logspace(-1, 3, 100) * e  # in joules
    #
    # dsvrec_dTi_vals = np.gradient(svrec_vals, Ti_vals, axis=0)
    #
    # zni_vals2d, Ti_vals2d = np.meshgrid(zni_vals, Ti_vals)
    #
    # zni_mod = np.where(n.i > 1E22, 1E22, n.i)
    # zni_mod = np.where(n.i < 1E16, 1E16, zni_mod)
    # Ti_mod = np.where(T.i.ev > 1E3, 1E3 * e, T.i.ev * e)
    # Ti_mod = np.where(T.i.ev < 1E-1, 1E-1 * e, Ti_mod)
    #
    # plt.semilogx(zni_vals2d.flatten(), Ti_vals2d.flatten())
    # plt.show()
    # print np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten()))
    # sys.exit()
    # sv_rec = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())),
    #                   svrec_vals.flatten(),
    #                   (zni_mod, Ti_mod),
    #                   method='linear', rescale=False)
    #
    # dsv_rec_dT = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())),
    #                            dsvrec_dTi_vals.flatten(),
    #                            (zni_mod, Ti_mod),
    #                            method='linear', rescale=False)

    # return sv_rec, dsv_rec_dT
    return 0, 0