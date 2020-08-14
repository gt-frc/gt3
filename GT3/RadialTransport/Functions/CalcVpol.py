#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.constants import constants

e = constants.elementary_charge
z_c = 6  # atomic number of carbon

def calc_vpol(Er, vphi_j, Lp, T, n, z_d, B_t, B_p, vphi_k, vpol_k, z_k):
    vpol = (1.0/B_t) * (1.0/(e*z_d) * -Lp.i * T.i.J + vphi_j * B_p - Er)
    vpol_assum = vpol_k - (T.C.J / (e * B_t)) * (Lp.i - (Lp.C / z_c)) + (B_p / B_t) * (vphi_j - vphi_k)
    vpol_alt = vpol_k - 1.0 / (e * B_t) * (T.i.J * Lp.i - (T.C.J * Lp.C / z_c)) + (B_p / B_t) * (vphi_j - vphi_k)
    return vpol, vpol_assum, vpol_alt
