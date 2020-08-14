#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

def calc_pwr_frac(ebeam):
    # calculate pwr fracs. These fits are only useful for DIII-D

    pwr_frac = np.zeros(3)
    pwr_frac[0] = (68 + 0.11 * ebeam) / 100
    pwr_frac[1] = (-159 + 6.53 * ebeam - 0.082 * (ebeam ** 2) + 0.00034 * (ebeam ** 3)) / 100
    pwr_frac[2] = (191 - 6.64 * ebeam + 0.082 * (ebeam ** 2) - 0.00034 * (ebeam ** 3)) / 100

    return pwr_frac
