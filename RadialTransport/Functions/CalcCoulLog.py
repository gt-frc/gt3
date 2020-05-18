#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import constants
from math import pi

eps_0 = constants.epsilon_0
e = constants.elementary_charge

def calc_coul_log(z1, z2, T_J, n2):
    """Calculate the coulomb logarithm.

    Reference Equation 1.36 in Dr. Stacey's Fusion Plasma Physics book

    :param z1:
    :param z2:
    :param T_J:
    :param n2:
    :return:
    """

    coul_log = np.log(12 * pi * np.sqrt((eps_0 * T_J)**3 / (n2 * (z2*e)**4 * (z1*e)**2)))
    return coul_log