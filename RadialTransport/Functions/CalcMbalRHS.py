#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.constants import constants

E_phi = 0.04  # toroidal electrostatic potential

e = constants.elementary_charge

def calc_mbal_rhs(mom_src_ext, z, n, B_p, gamma):
    """ The function formerly known as y11 and y22
    """

    mbal_rhs = mom_src_ext + (z * e) * (n * E_phi + B_p * gamma)

    return mbal_rhs