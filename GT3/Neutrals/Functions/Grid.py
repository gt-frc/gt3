#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata

def grid(x, y, z, resX=100, resY=100):
    """Convert 3 column data to matplotlib grid"""
    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    X, Y = np.meshgrid(xi, yi)
    Z = griddata(np.column_stack((x, y)), z, (X, Y), method='linear')

    #interp = Rbf(x, y, z, function='linear')
    #Z = interp(X, Y)
    return X, Y, Z