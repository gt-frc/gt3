#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from RadialTransport.Functions.CalcQie import calc_qie
import numpy as np
from math import exp

def calc_Qi_int_method(r, n, T, en_src_nbi_i_kept, cool_rate, iol_adjusted=False, E_orb=None):  # formerly qheat

    Qi = np.zeros(r.shape)

    qie = calc_qie(n, T, ion_species='D')

    # Boundary condition at the magnetic axis.
    # Only has the second term, since it's the center value. Also uses delta_r of the next point.
    # If not adjusted for IOL, en_src_nbi_kept = en_src_nbi_tot, so no need for an IOL check.
    Qi[0] = (en_src_nbi_i_kept[0] - cool_rate[0] - qie[0]) * (r[1] - r[0])
    Qi[1] = Qi[0] + (en_src_nbi_i_kept[1] - cool_rate[1] - qie[1]) * (r[1] - r[0])

    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the particle continuity equation, but different source term.
    for n in range(2, len(r)):
        if iol_adjusted:
            # Imported the exp() function for the thermal IOL attenuation.
            Qi[n] = (r[n - 1] / r[n]) * Qi[n - 1] * exp(-(E_orb[n] - E_orb[n - 1])) + (
                        en_src_nbi_i_kept[n] - cool_rate[n] - qie[n]) * (r[n] - r[n - 1])
        else:
            Qi[n] = (r[n - 1] / r[n]) * Qi[n - 1] + (en_src_nbi_i_kept[n] - cool_rate[n] - qie[n]) * (r[n] - r[n - 1])

    return Qi, qie
