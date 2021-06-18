#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import UnivariateSpline

def calc_return_cur(r, part_src_nbi_lost, gamma, ch_d, iol_adjusted=False, F_orb=None):
    # Piper Changes: This calculates the return current from fast and thermal IOL.
    # This is a new quantity that i want to calculate. Big part of my thesis.

    if iol_adjusted:
        diff_F_orb = UnivariateSpline(r, F_orb, k=1, s=0).derivative()(
            r)  # Need the derivative of F_orb with respect to r.
    else:
        diff_F_orb = np.zeros(r.shape)

    Jr_iol = np.zeros(r.shape)
    Jr_iol[0] = (part_src_nbi_lost[0] + (gamma[0] * diff_F_orb[0])) * (r[1] - r[0]) * ch_d
    Jr_iol[1] = Jr_iol[0] + (part_src_nbi_lost[1] + (gamma[1] * diff_F_orb[1])) * (r[1] - r[0]) * ch_d

    # Integral cylindrical form of the continuity of (current?) equation.
    for n in range(2, len(r)):
        Jr_iol[n] = (r[n - 1] / r[n]) * Jr_iol[n - 1] + (part_src_nbi_lost[n] + (gamma[n] * diff_F_orb[n])) * (
                    r[n] - r[n - 1]) * ch_d  # No factor of 2, because only 1 direction.

    return Jr_iol
