#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from math import sqrt, pi
import numpy as np

def calc_intrin_rot(M_orb, T_J, m):

    intrin_rot = 2 / sqrt(pi) * M_orb * np.sqrt(2 * T_J / m)
    return intrin_rot
