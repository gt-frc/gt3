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
from gtedge_LZ import calc_Lz_gtedge

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

    print ' chi_r = ', chi_r
    print ' nu = ', nu
    print ' L.T^-1 = ', 1/L.T
    print ' L.n^-1 = ', 1/L.n
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

    # term1 = chi_r * (nu * L.T**-2 + (C2 - 1.0) * L.T**-1 * L.n**-1)
    # term2 = fz * ((nu + 1 - C2) * Lz_t / T.i.J - dLzdT)
    # term3 = f0 * (E_ion * sv.ion / T.i.J * (nu - T.i.J / sv.ion * sv.ion_ddT))
    # term4 = f0c * (3.0 / 2.0 * (sv.cx + sv.el) * (nu - 1.0 - T.i.J * (sv.cx_ddT + sv.el_ddT) / (sv.cx + sv.el)))

    term1 = chi_r * (nu * L.T**-2 + (C2 - 1.0) * L.T**-1 * L.n**-1)
    term2 = fz * ((nu + 1 - C2) * Lz_t / (T.i.J+T.e.J)/2 - dLzdT)
    term3 = f0 * (E_ion * sv.ion / T.i.J * (nu - (T.i.J+T.e.J)/2 / sv.ion * sv.ion_ddT))
    term4 = f0c * (3.0 / 2.0 * (sv.cx + sv.el) * (nu - 1.0 - (T.i.J+T.e.J)/2 * (sv.cx_ddT + sv.el_ddT) / (sv.cx + sv.el)))
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
    print
    print 'MI = ', n.e / n_marfe
    return n_marfe, n.e / n_marfe


class Marfe:
    """
    """
    def __init__(self, inputs=None, core=None):
        if core is not None:
            if inputs is not None:
                print 'ignoring \'inputs\' and attempting to use core instance'
            # use 'core' instance to create inputs for calc_n_marfe

            n_xpt = namedtuple('n', 'n e i C')(
                namedtuple('nn', 's t')(
                    core.n.n.s,
                    core.n.n.t
                ),
                core.n.e,
                core.n.i,
                core.n.C
            )

            T_xpt = namedtuple('T', 'e i n')(
                namedtuple('Te', 'kev ev J')(
                    core.T.e.kev,
                    core.T.e.ev,
                    core.T.e.J
                ),
                namedtuple('Ti', 'kev ev J')(
                    core.T.i.kev,
                    core.T.i.ev,
                    core.T.i.J
                ),
                namedtuple('Tn', 's t')(
                    core.T.n.s,
                    core.T.n.t
                )
            )

            sv_xpt = namedtuple('sv', 'ion ion_ddT el el_ddT cx cx_ddT')(
                core.sv.ion.st,
                core.sv.ion.d_dT.st,
                core.sv.el.st,
                core.sv.el.d_dT.st,
                core.sv.cx.st,
                core.sv.cx.d_dT.st
            )

            L_xpt = namedtuple('L', 'n T')(
                core.L.n.i,
                core.L.T.i
            )

            Lz_xpt = namedtuple('Lz', 't ddT')(
                core.Lz.t,
                core.Lz.ddT.t
            )

            chi_r_xpt = core.chi_r

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

        self.n_marfe, self.MI = calc_n_marfe(n, sv, T, L, Lz, chi_r)
