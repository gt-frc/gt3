#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from RadialTransport.Functions.CalcQie import calc_qie
import numpy as np

def calc_Qe_int_method(r, n, T, en_src_nbi_e_tot, cool_rate):  # Piper Changes: Same as Qi changes.

    Qe = np.zeros(r.shape)
    qie = calc_qie(n, T, ion_species='D')

    Qe[0] = (en_src_nbi_e_tot[0] - qie[0] - cool_rate[0]) * (r[1] - r[0])
    Qe[1] = Qe[0] + (en_src_nbi_e_tot[1] - qie[1] - cool_rate[1]) * (r[1] - r[0])

    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the continuity equation, but different source term.
    for n in range(2, len(r)):
        Qe[n] = (r[n - 1] / r[n]) * Qe[n - 1] + (en_src_nbi_e_tot[n] - qie[n] - cool_rate[n]) * (r[n] - r[n - 1])

    return Qe