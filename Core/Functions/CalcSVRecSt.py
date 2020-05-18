#!/usr/bin/env python2
# -*- coding: utf-8 -*-

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