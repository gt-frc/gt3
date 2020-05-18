#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from math import pi

from scipy.constants import constants

eps_0 = constants.epsilon_0
e = constants.elementary_charge


def calc_coul_log_j_k(z_j, z_k, T_j, n_k):
    # Coulomb logarithm calculation in GTEDGE PARAM subroutine.
    coul_log = np.log(12 * pi * (T_j**1.5) * ((eps_0/e)**1.5) / (np.sqrt(n_k) * (z_k**2.0) * z_j))
    return coul_log
