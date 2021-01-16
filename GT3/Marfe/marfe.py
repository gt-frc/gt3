#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:21:23 2018

@author: max
"""

import numpy as np
from GT3.ImpRadiation import ImpRad
from collections import namedtuple
import sys


def calc_z_0(n):
    z_0 = n.C * 6.0 ** 2 / n.i
    return z_0


def calc_z_eff(n):
    z_eff = (n.i * 1.0 ** 2 + n.C * 6.0 ** 2) / n.e
    return z_eff


# TODO: Generalize this to non carbon impurities. Refer to equations 14.39 in Stacey's book.
def calc_Ci2(n):
    z_0 = calc_z_0(n)
    Ci2 = 1.56 * (1.0 + np.sqrt(2) * z_0) * (1.0 + 0.52 * z_0) / \
          ((1 + 2.65 * z_0) * (1 + 0.285 * z_0) * (z_0 + np.sqrt(0.5 * (1 + (1 / 6)))))
    return Ci2


def calc_fz(n):
    fz = n.C / n.i  # ne?
    return fz


def calc_f0(n):
    f0 = (n.n.s + n.n.t) / n.i  # ne?
    return f0


def calc_f0c(n):
    f0c = n.n.s / n.i  # ne?
    return f0c


def calc_Ce2(n):
    z_eff = calc_z_eff(n)
    Ce2 = 1.5 * (1 - 0.6934 / 1.3167 ** z_eff)
    return Ce2


def calc_C2(n):
    z_0 = calc_z_0(n)

    Ce2 = calc_Ce2(n)
    Ci2 = calc_Ci2(n)
    C2 = Ce2 - z_0 * Ci2
    return C2


def calc_E_ion():
    return 15.466 * 1.6021E-19  # ionization energy of deuterium in Joules


def calc_Lz(n, T):

    # calculate a weighted average neutral temperature
    Tn = (n.n.s*T.n.s + n.n.t*T.n.t) / (n.n.s + n.n.t)

    # calculate a total neutral fraction
    nf = 1.0E-7  #(n.n.s + n.n.t) / n.e

    # obtain interpolations from ImpRad
    imp_C = ImpRad(z=6)
    Lz_interp_C = imp_C.Lz
    dLzdT_interp_C = imp_C.dLzdT

    # obtain values from the interpolations
    Lz_C = Lz_interp_C(np.log10(Tn), np.log10(nf), np.log10(T.e.kev))
    dLzdT_C = dLzdT_interp_C(np.log10(Tn), np.log10(nf), np.log10(T.e.kev))

    return Lz_C, dLzdT_C


def calc_n_marfe(n, sv, T, L, Lz, chi_r):
    """
    Calculates the density limit for MARFE onset

    :param Lz:
    :param n: class instance or namedtuple with density attributes
    :param sv: class instance or namedtuple with cross section attributes
    :param T: class instance or namedtuple with temperature attributes
    :param L: class instance or namedtuple with gradient scale length attributes
    :param chi_r: float or array of floats for chi_radial
    :return: marfe density limit as float or array, depending on the shape of the inputs.
    """
    C2 = calc_C2(n)
    E_ion = calc_E_ion()
    Lz_t = Lz.t
    dLzdT = Lz.ddT

    #Lz_t, dLzdT = calc_Lz_gtedge(T.i.kev*1E3)
    #Lz_t = Lz_t * 1E-13  # now in Jm^3/s
    #dLzdT = dLzdT * 1E-13  # now in Jm^3/s

    nu = 5 / 2
    fz = calc_fz(n)
    f0 = calc_f0(n)
    f0c = calc_f0c(n)

    if False:
        print('MARFE Index Parameters')
        print(' chi_r = ', chi_r)
        print(' nu = ', nu)
        print(' L.T^-1 = ', 1/L.T)
        print(' L.n^-1 = ', 1/L.n)
        print(' C2 = ', C2)
        print(' fz = ', fz)
        print(' f0 = ', f0)
        print(' f0c = ', f0c)
        print(' E_ion = ', E_ion)
        print(' sv.ion = ', sv.ion)
        print(' sv.ion_ddT = ', sv.ion_ddT)
        print(' sv.cx = ', sv.cx)
        print(' sv.cx_ddT = ', sv.cx_ddT)
        print(' sv.el = ', sv.el)
        print(' sv.el_ddT = ', sv.el_ddT)
        print(' Lz_t = ', Lz_t)
        print(' dLzdT = ', dLzdT)
        print(' T.i.J = ',T.i.J)
        print(' n.e = ', n.e)
        print()

    term1 = chi_r * (nu * L.T**-2 + (C2 - 1.0) * L.T**-1 * L.n**-1)
    term2 = fz * ((nu + 1 - C2) * Lz_t / (T.i.J+T.e.J)/2 - dLzdT)
    term3 = f0 * (E_ion * sv.ion / T.i.J * (nu - (T.i.J+T.e.J)/2 / sv.ion * sv.ion_ddT))
    term4 = f0c * (3.0 / 2.0 * (sv.cx + sv.el) * (nu - 1.0 - (T.i.J+T.e.J)/2 * (sv.cx_ddT + sv.el_ddT) / (sv.cx + sv.el)))

    if False:
        print()
        print('term1 = ', term1)
        print('term2 = ', term2)
        print('term3 = ', term3)
        print('term4 = ', term4)
        print('term 2+3+4 = ', term2 + term3 + term4)

    n_marfe = term1 / \
              (
                      term2 +
                      term3 +
                      term4
              )
    return n_marfe, n.e / n_marfe


class Marfe:

    def __init__(self, inputs=None, core=None):
        """

        :param inputs:
        :param core:
        """
        if core is not None:
            if inputs is not None:
                print('ignoring \'inputs\' and attempting to use core instance')
            # use 'core' instance to create inputs for calc_n_marfe

            nn_mult = 1.0

            n_xpt = namedtuple('n', 'n e i C W Be')(
                namedtuple('nn', 's t tot')(
                    core.n.n.s[core.xpt_loc][0] * nn_mult,
                    core.n.n.t[core.xpt_loc][0] * nn_mult,
                    core.n.n.tot[core.xpt_loc][0] * nn_mult
                ),
                core.n.e[core.xpt_loc][0],
                core.n.i[core.xpt_loc][0],
                core.n.C[core.xpt_loc][0],
                core.n.C[core.xpt_loc][0],
                core.n.C[core.xpt_loc][0],
            )

            n_obmp = namedtuple('n', 'n e i C')(
                namedtuple('nn', 's t tot')(
                    core.n.n.s[core.obmp_loc][0] * nn_mult,
                    core.n.n.t[core.obmp_loc][0] * nn_mult,
                    core.n.n.tot[core.obmp_loc][0] * nn_mult,
                ),
                core.n.e[core.obmp_loc][0],
                core.n.i[core.obmp_loc][0],
                core.n.C[core.obmp_loc][0]
            )

            n_av = namedtuple('n', 'n e i C')(
                namedtuple('nn', 's t')(
                    (n_xpt.n.s + n_obmp.n.s) / 2,
                    (n_xpt.n.t + n_obmp.n.t) / 2
                ),
                (n_xpt.e + n_obmp.e) / 2,
                (n_xpt.i + n_obmp.i) / 2,
                (n_xpt.C + n_obmp.C) / 2
            )

            T_xpt = namedtuple('T', 'e i n')(
                namedtuple('Te', 'kev ev J')(
                    core.T.e.kev[core.xpt_loc][0],
                    core.T.e.ev[core.xpt_loc][0],
                    core.T.e.J[core.xpt_loc][0]
                ),
                namedtuple('Ti', 'kev ev J')(
                    core.T.i.kev[core.xpt_loc][0],
                    core.T.i.ev[core.xpt_loc][0],
                    core.T.i.J[core.xpt_loc][0]
                ),
                namedtuple('Tn', 's t')(
                    core.T.n.s.kev[core.xpt_loc][0],
                    core.T.n.t.kev[core.xpt_loc][0]
                )
            )

            T_obmp = namedtuple('T', 'e i n')(
                namedtuple('Te', 'kev ev J')(
                    core.T.e.kev[core.obmp_loc][0],
                    core.T.e.ev[core.obmp_loc][0],
                    core.T.e.J[core.obmp_loc][0]
                ),
                namedtuple('Ti', 'kev ev J')(
                    core.T.i.kev[core.obmp_loc][0],
                    core.T.i.ev[core.obmp_loc][0],
                    core.T.i.J[core.obmp_loc][0]
                ),
                namedtuple('Tn', 's t')(
                    core.T.n.s.kev[core.obmp_loc][0],
                    core.T.n.t.kev[core.obmp_loc][0]
                )
            )

            sv_xpt = namedtuple('sv', 'ion ion_ddT el el_ddT cx cx_ddT')(
                core.sv.ion.st[core.xpt_loc][0],
                core.sv.ion.d_dT.st[core.xpt_loc][0],
                core.sv.el.st[core.xpt_loc][0],
                core.sv.el.d_dT.st[core.xpt_loc][0],
                core.sv.cx.st[core.xpt_loc][0],
                core.sv.cx.d_dT.st[core.xpt_loc][0]
            )

            L_xpt = namedtuple('L', 'n T')(
                core.L.n.i[core.xpt_loc][0],
                core.L.T.i[core.xpt_loc][0]
            )

            # L_xpt = namedtuple('L', 'n T')(
            #     1.0/25.0,
            #     1.0/38.0
            # )

            # L_xpt = namedtuple('L', 'n T')(
            #     1.0/15.6,
            #     1.0/14.2
            # )

            L_obmp = namedtuple('L', 'n T')(
                core.L.n.i[core.obmp_loc][0],
                core.L.T.i[core.obmp_loc][0]
            )

            L_ibmp = namedtuple('L', 'n T')(
                core.L.n.i[core.ibmp_loc][0],
                core.L.T.i[core.ibmp_loc][0]
            )

            L_top = namedtuple('L', 'n T')(
                core.L.n.i[core.top_loc][0],
                core.L.T.i[core.top_loc][0]
            )

            L_av = namedtuple('L', 'n T')(
                (core.L.n.i[core.top_loc][0] + core.L.n.i[core.xpt_loc][0] + core.L.n.i[core.ibmp_loc][0] + core.L.n.i[core.obmp_loc][0])/4.0,
                (core.L.T.i[core.top_loc][0] + core.L.T.i[core.xpt_loc][0] + core.L.T.i[core.ibmp_loc][0] + core.L.T.i[core.obmp_loc][0])/4.0
            )

            xpt_weight = 1
            L_av2 = namedtuple('L', 'n T')(
                (xpt_weight * core.L.n.i[core.xpt_loc][0] + (1-xpt_weight) * core.L.n.i[core.obmp_loc][0])/2.0,
                (xpt_weight * core.L.T.i[core.xpt_loc][0] + (1-xpt_weight) * core.L.T.i[core.obmp_loc][0])/2.0
            )

            if False:
                print('L_xpt.n^-1 = ',1/L_xpt.n)
                print('L_xpt.T^-1 = ',1/L_xpt.T)
                print('L_obmp.n^-1 = ',1/L_obmp.n)
                print('L_obmp.T^-1 = ',1/L_obmp.T)
                print('L_top.n^-1 = ',1/L_top.n)
                print('L_top.T^-1 = ',1/L_top.T)
                print('L_ibmp.n^-1 = ',1/L_ibmp.n)
                print('L_ibmp.T^-1 = ',1/L_ibmp.T)
                # sys.exit()

            Lz_xpt = namedtuple('Lz', 't ddT')(
                core.Lz.t[core.xpt_loc][0],
                core.Lz.ddT.t[core.xpt_loc][0]
            )

            chi_r_xpt = core.chi.bohm[core.xpt_loc][0]
            chi_r_top = core.chi.bohm[core.top_loc][0]
            chi_r_ibmp = core.chi.bohm[core.ibmp_loc][0]
            chi_r_obmp = core.chi.bohm[core.obmp_loc][0]
            chi_r_av = (core.chi.bohm[core.xpt_loc][0] + core.chi.bohm[core.top_loc][0] + core.chi.bohm[core.ibmp_loc][0] + core.chi.bohm[core.obmp_loc][0])/4.0

            #self.n_marfe, self.MI = calc_n_marfe(n_xpt, sv_xpt, T_xpt, L_xpt, Lz_xpt, chi_r_xpt)
            self.n_marfe, self.MI = calc_n_marfe(n_xpt, sv_xpt, T_xpt, L_av2, Lz_xpt, chi_r_xpt)

        elif inputs is not None:
            # use 'inputs' to create inputs for calc_n_marfe
            n = inputs.n
            sv = inputs.sv
            T = inputs.T
            L = inputs.L
            Lz = inputs.Lz
            chi_r = inputs.chi_r
            self.n_marfe, self.MI = calc_n_marfe(n, sv, T, L, Lz, chi_r)
        else:
            raise Exception("No inputs provided")

        print()
        print('######################################')
        print()
        print('ni = ', n_xpt.i)
        print('ne = ', n_xpt.e)
        print('nn_xpt = ', n_xpt.n.tot)
        print('nn_obmp = ', n_obmp.n.tot)
        print('n_marfe = ',self.n_marfe)
        print()
        print('Ti(ev) = ', T_xpt.i.ev)
        print('Te(ev) = ', T_xpt.e.ev)
        print()
        print('L_n^-1 (xpt) = ',1/L_xpt.n)
        print('L_T^-1 (xpt) = ',1/L_xpt.T)
        print()
        print('L_n^-1 (obmp) = ',1/L_obmp.n)
        print('L_T^-1 (obmp) = ',1/L_obmp.T)
        print()
        print('L_n^-1 (av) = ',1/L_av.n)
        print('L_T^-1 (av) = ',1/L_av.T)
        print()
        print('nC / ne = ', n_xpt.W / n_xpt.e)
        print('nn / ne (xpt) = ', (n_xpt.n.s + n_xpt.n.t) / n_xpt.e)
        print('nn / ne (obmp) = ', (n_obmp.n.s + n_obmp.n.t) / n_obmp.e)
        print('MI = ', self.MI)
        print()
        print('######################################')
        print()
