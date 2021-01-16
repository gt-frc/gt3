#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from math import pi

def calc_fsa(x, R, Z):


    R1 = R[:, :-1]
    R2 = np.roll(R[:, :-1], -1, axis=1)
    Z1 = Z[:, :-1]
    Z2 = np.roll(Z[:, :-1], -1, axis=1)
    x1 = x[:, :-1]
    x2 = np.roll(x[:, :-1], -1, axis=1)

    dl = np.sqrt((R2 - R1) ** 2 + (Z2 - Z1) ** 2)

    R_av = (R1 + R2)/2

    dA = dl * (2 * pi * R_av)

    x_av = (x1 + x2)/2

    fsa = np.sum(x_av * dA, axis=1) / np.sum(dA, axis=1)
    fsa[0] = x[0, 0]
    return fsa