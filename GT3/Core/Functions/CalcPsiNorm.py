#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata

def calc_psi_norm(R, Z, psi, xpt, axis_mag):
    # normalize psi
    # psi_interp = Rbf(R, Z, psi)
    # psi_min = psi_interp(axis_mag[0], axis_mag[1])
    #
    # psi_shifted = psi - psi_min  # set center to zero
    # psi_shifted_interp = Rbf(R, Z, psi_shifted)
    # psi_shifted_xpt = psi_shifted_interp(xpt[0], xpt[1])
    #
    # psi_norm = psi_shifted / psi_shifted_xpt

    psi_min = griddata(np.column_stack((R.flatten(), Z.flatten())),
                       psi.flatten(),
                       [axis_mag[0], axis_mag[1]],
                       method='cubic')

    psi_shifted = psi - psi_min  # set center to zero
    psi_shifted_xpt_l, psi_shifted_xpt_u = None, None

    if xpt[0] is not None:

        psi_shifted_xpt_l = griddata(np.column_stack((R.flatten(), Z.flatten())),
                               psi_shifted.flatten(),
                               [xpt[0][0], xpt[0][1]],
                               method='cubic')
    if xpt[1] is not None:
        psi_shifted_xpt_u = griddata(np.column_stack((R.flatten(), Z.flatten())),
                               psi_shifted.flatten(),
                               [xpt[1][0], xpt[1][1]],
                               method='cubic')
    psi_shifted_xpt = [psi_shifted_xpt_l, psi_shifted_xpt_u]
    if xpt[1] is None:
        psi_norm = psi_shifted / np.average(psi_shifted_xpt_l)
    elif xpt[0] is None:
        psi_norm = psi_shifted / np.average(psi_shifted_xpt_u)
    else:
        psi_norm = psi_shifted / np.average(psi_shifted_xpt)

    return psi_norm
