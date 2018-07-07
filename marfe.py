#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:21:23 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from imp_rad import ImpRad
from collections import namedtuple
import sys
from exp_core_brnd import *


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

    print
    print 'Tn = ',Tn
    print 'nf = ',nf
    print 'T.e.kev = ',T.e.kev
    print 'Lz_C = ',Lz_C
    print 'dLzdT_C = ',dLzdT_C
    print
    return Lz_C, dLzdT_C


def calc_n_marfe(n, sv, T, L, Lz, chi_r):
    """
    Calculates the density limit for MARFE onset

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
    nu = 5 / 2
    fz = calc_fz(n)
    f0 = calc_f0(n)
    f0c = calc_f0c(n)

    print ' chi_r = ', chi_r
    print ' nu = ', nu
    print ' L.T = ', L.T
    print ' L.n = ', L.n
    print ' C2 = ', C2
    print ' fz = ', fz
    print ' f0 = ', f0
    print ' f0c = ', f0c
    print ' E_ion = ', E_ion
    print ' sv.ion = ', sv.ion
    print ' sv.ion_ddT = ', sv.ion_ddT
    print ' sv.cx = ', sv.cx
    print ' sv.cx_ddT = ', sv.cx_ddT
    print ' sv.el = ', sv.el
    print ' sv.el_ddT = ', sv.el_ddT
    print ' Lz_t = ', Lz_t
    print ' dLzdT = ', dLzdT
    print ' T.i.J = ',T.i.J
    print ' n.e = ', n.e
    print



    term1 = chi_r * (nu * L.T**-2 + (C2 - 1.0) * L.T**-1 * L.n**-1)
    term2 = fz * ((nu + 1 - C2) * Lz_t / T.i.J - dLzdT)
    term3 = f0 * (E_ion * sv.ion / T.i.J * (nu - T.i.J / sv.ion * sv.ion_ddT))
    term4 = f0c * (3.0 / 2.0 * (sv.cx + sv.el) * (nu - 1.0 - T.i.J * (sv.cx_ddT + sv.el_ddT) / (sv.cx + sv.el)))
    print
    print 'term1 = ', term1
    print 'term2 = ', term2
    print 'term3 = ', term3
    print 'term4 = ', term4
    print 'term 2+3+4 = ', term2 + term3 + term4

    n_marfe = term1 / \
              (
                      term2 +
                      term3 +
                      term4
              )

    print 'n_marfe = ', n_marfe
    print 'MI = ', n.e / n_marfe
    return n_marfe


class Marfe:
    """
    """
    def __init__(self, inputs=None, core=None):
        if core is not None:
            if inputs is not None:
                print 'ignoring \'inputs\' and attempting to use core instance'
            # use 'core' instance to create inputs for calc_n_marfe
            nn_dict = {}
            nn_dict['s'] = core.n.n.s
            nn_dict['t'] = core.n.n.t
            nn = namedtuple('nn', nn_dict.keys())(*nn_dict.values())

            n_dict = {}
            n_dict['n'] = nn
            n_dict['e'] = core.n.e
            n_dict['i'] = core.n.i
            n_dict['C'] = core.n.C
            n = namedtuple('n', n_dict.keys())(*n_dict.values())

            Te_dict = {}
            Te_dict['kev'] = core.T.e.kev
            Te_dict['J'] = core.T.e.J
            Te = namedtuple('Te', Te_dict.keys())(*Te_dict.values())

            Ti_dict = {}
            Ti_dict['kev'] = core.T.i.kev
            Ti_dict['J'] = core.T.i.J
            Ti = namedtuple('Ti', Ti_dict.keys())(*Ti_dict.values())

            Tn_dict = {}
            Tn_dict['s'] = core.T.n.s  # in kev
            Tn_dict['t'] = core.T.n.t  # in kev
            Tn = namedtuple('Tn', Tn_dict.keys())(*Tn_dict.values())

            T_dict = {}
            T_dict['e'] = Te
            T_dict['i'] = Ti
            T_dict['n'] = Tn
            T = namedtuple('T', T_dict.keys())(*T_dict.values())

            sv_dict = {}
            sv_dict['ion'] = core.sv.ion.st
            sv_dict['ion_ddT'] = core.sv.ion.d_dT.st
            sv_dict['el'] = core.sv.el.st
            sv_dict['el_ddT'] = core.sv.el.d_dT.st
            sv_dict['cx'] = core.sv.cx.st
            sv_dict['cx_ddT'] = core.sv.cx.d_dT.st
            sv = namedtuple('sv', sv_dict.keys())(*sv_dict.values())

            L_dict = {}
            L_dict['n'] = core.L.n.i
            L_dict['T'] = core.L.T.i
            L = namedtuple('L', L_dict.keys())(*L_dict.values())

            Lz_dict = {}
            Lz_dict['t'] = core.Lz.t
            Lz_dict['ddT'] = core.Lz.ddT.t
            Lz = namedtuple('Lz', Lz_dict.keys())(*Lz_dict.values())

            chi_r = core.chi_r

        elif inputs is not None:
            # use 'inputs' to create inputs for calc_n_marfe
            n = inputs.n
            sv = inputs.sv
            T = inputs.T
            L = inputs.L
            Lz = inputs.Lz
            chi_r = inputs.chi_r
        else:
            print 'you haven\'t provided any inputs'
            print 'stopping'
            sys.exit()

        self.n_marfe = calc_n_marfe(n, sv, T, L, Lz, chi_r)




if __name__ == '__main__':

    # Values for 92980.3600
    nn_dict = {}
    nn_dict['s'] = 0.0
    nn_dict['t'] = 0.29 * 1E20 * 1.2E-3
    nn = namedtuple('nn', nn_dict.keys())(*nn_dict.values())

    n_dict = {}
    n_dict['n'] = nn
    n_dict['e'] = 0.29 * 1E20
    n_dict['i'] = n_dict['e']
    n_dict['C'] = n_dict['e'] * 0.015
    n = namedtuple('n', n_dict.keys())(*n_dict.values())

    Te_dict = {}
    Te_dict['kev'] = np.array([0.03,0.05])
    Te_dict['ev'] = Te_dict['kev'] * 1E3
    Te_dict['J'] = Te_dict['kev'] * 1E3 * 1.6021E-19
    Te = namedtuple('Te', Te_dict.keys())(*Te_dict.values())

    Ti_dict = {}
    Ti_dict['kev'] = 0.08
    Ti_dict['ev'] = Ti_dict['kev'] * 1E3
    Ti_dict['J'] = Ti_dict['kev'] * 1E3 * 1.6021E-19
    Ti = namedtuple('Ti', Ti_dict.keys())(*Ti_dict.values())

    Tn_dict = {}
    Tn_dict['s'] = 0.002  # in kev
    Tn_dict['t'] = Ti_dict['kev']  # in kev
    Tn = namedtuple('Tn', Tn_dict.keys())(*Tn_dict.values())

    T_dict = {}
    T_dict['e'] = Te
    T_dict['i'] = Ti
    T_dict['n'] = Tn
    T = namedtuple('T', T_dict.keys())(*T_dict.values())

    sv_dict = {}
    sv_dict['ion'] = calc_svion_st(T)[0]
    sv_dict['ion_ddT'] = calc_svion_st(T)[1]
    sv_dict['el'] = calc_svel_st(T)[0]
    sv_dict['el_ddT'] = calc_svel_st(T)[1]
    sv_dict['cx'] = calc_svcx_st(T)[0]
    sv_dict['cx_ddT'] = calc_svcx_st(T)[1]
    sv = namedtuple('sv', sv_dict.keys())(*sv_dict.values())

    L_dict = {}
    L_dict['n'] = 1 / 17.9
    L_dict['T'] = 1 / 16.4
    L = namedtuple('L', L_dict.keys())(*L_dict.values())

    Lz_dict = {}
    Lz_dict['t'] = calc_Lz(n, T)[0]
    Lz_dict['ddT'] = calc_Lz(n, T)[1]
    Lz = namedtuple('Lz', Lz_dict.keys())(*Lz_dict.values())

    chi_r = 2.0

    calc_n_marfe(n, sv, T, L, Lz, chi_r)

























    # n_e = 0.46 * 1E20  # complete
    # n_i = n_e  # complete
    # n_C = 0.013 * n_e  # complete
    # nn_s = 0  # complete
    # nn_t = n_i * 1E-3  # complete
    #
    # T_units_nt = namedtuple('T_units', '')
    # T_nt = namedtuple('T', '')
    #
    # sv_ion = calc_svel_st(T)
    # sv_ion_ddT =
    # sv_el =
    # sv_el_ddT =
    # sv_cx =
    # sv_cx_ddT =
    #
    # T_e_kev =
    # T_e_ev = T_e_kev * 1E3
    # T_e_J = T_e_kev * 1E3 * 1.6021E-19
    #
    # T_i_kev =
    # T_i_ev = T_i_kev * 1E3
    # T_i_J = T_i_kev * 1E3 * 1.6021E-19
    #
    # T_n_kev =
    # T_n_ev = T_n_kev * 1E3
    # T_n_J = T_n_kev * 1E3 * 1.6021E-19
    #
    # L_n =
    # L_T =
    #
    # Lz_t =
    # Lz_ddT =
    #
    # chi_r = 2.0
    #
    # nn_nt = namedtuple('nn', 's t')
    # n_n = nn_nt(nn_s, nn_t)
    #
    # n_nt = namedtuple('n', 'n e i C')
    # n = n_nt(n_n, n_e, n_i, n_C)
    #
    # sv_nt = namedtuple('sv', 'ion ion_ddT el el_ddT cx cx_ddT')
    # sv = sv_nt(sv_ion, sv_ion_ddT, sv_el, sv_el_ddT, sv_cx, sv_cx_ddT)
    #
    # T_units_nt = namedtuple('T_units', 'kev ev J')
    # T_i = T_units_nt(T_i_kev, T_i_ev, T_i_J)
    # T_e = T_units_nt(T_e_kev, T_e_ev, T_e_J)
    # T_n = T_units_nt(T_n_kev, T_n_ev, T_n_J)
    #
    # T_nt = namedtuple('T', 'e i n')
    # T = T_nt(T_i, T_e, T_n)
    #
    # L_nt = namedtuple('L', 'n T')
    # L = L_nt(L_n, L_T)
    #
    # Lz_nt = namedtuple('Lz', 't ddT')
    # Lz = Lz_nt(Lz_t, Lz_ddT)
    #
    # inputs_nt = namedtuple('inputs', 'n sv T L Lz chi_r')
    # inputs = inputs_nt(n, sv, T, L, Lz, chi_r)


    # chi_r = 2.0
    # sv_ion = core.sv_ion
    # sv_cx = core.sv_cx
    # sv_el = core.sv_el
    # L_T = core.L_Ti_J
    # L_n = core.L_ni
    # T = core.Ti_J
    # fz = core.nC / core.ni  # ne?
    # f0 = core.n_n_total / core.ni
    # f0c = core.n_n_slow / core.ni
    # Lz = core.Lz_thermal
    # dLzdT = core.dLzdT_thermal
    # dsv_ion_dT = core.dsv_ion_dT
    # dsv_cxel_dT = core.dsv_cx_dT + core.dsv_el_dT



    # n_marfe_edge = np.where(core.psi_norm > 0.5, n_marfe, np.nan)
    # n_marfe_met = np.where((core.ni > n_marfe) & (n_marfe > 0) & (core.psi_norm > 0.5), core.ni - n_marfe, np.nan)
    #
    # marfe_fig1 = plt.figure(figsize=(13, 6))
    # ax1 = marfe_fig1.add_subplot(1, 3, 1)
    # ax1.axis('equal')
    # ax1.set_title(r'$n_{MARFE}$')
    # cs1 = ax1.contourf(core.R, core.Z, n_marfe_edge, 500)
    # ax1.plot(core.R[-1, :], core.Z[-1, :], lw=1, color='red')
    # ax1.plot(inp.wall_exp[:, 0], inp.wall_exp[:, 1], lw=1, color='black')
    # marfe_fig1.colorbar(cs1, ax=ax1)
    #
    # ax2 = marfe_fig1.add_subplot(1, 3, 2)
    # ax2.axis('equal')
    # ax2.set_title(r'$n_i$')
    # cs2 = ax2.contourf(core.R, core.Z, core.ni, 500)
    # ax2.plot(inp.wall_exp[:, 0], inp.wall_exp[:, 1], lw=1, color='black')
    # marfe_fig1.colorbar(cs2, ax=ax2)
    #
    # ax3 = marfe_fig1.add_subplot(1, 3, 3)
    # ax3.axis('equal')
    # ax3.set_title('$n_i-n_{MARFE}$ where greater')
    # cs3 = ax3.contourf(core.R, core.Z, n_marfe_met, 500)
    # ax3.plot(core.R[-1, :], core.Z[-1, :], lw=1, color='red')
    # ax3.plot(inp.wall_exp[:, 0], inp.wall_exp[:, 1], lw=1, color='black')
    # marfe_fig1.colorbar(cs3, ax=ax3)
    # plt.tight_layout()
