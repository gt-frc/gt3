#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata, interp1d
from shapely.geometry import Point, LineString

def calc_psi(rho, pts, psi_data):
    # move along line between m_axis and obmp and define psi values corresponding to rho values
    rho_line = LineString([Point(pts.axis.mag), Point(pts.obmp)])
    init_coords = np.zeros((0, 2))
    for i, rhoval in enumerate(rho[:,0]):
        pt_coords = np.asarray(rho_line.interpolate(rhoval, normalized=True).coords)[0]
        init_coords = np.vstack((init_coords, pt_coords))

    psi_vals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                             psi_data.psi.flatten(),
                             init_coords,
                             method='cubic')
    psi_norm_vals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                             psi_data.psi_norm.flatten(),
                             init_coords,
                             method='cubic')

    psi = interp1d(rho[:,0], psi_vals)(rho)
    psi_norm = interp1d(rho[:, 0], psi_norm_vals)(rho)

    # If rho goes all the way to 1 (the usual case) then force the last psi_norm value to be 1.
    # Otherwise things break later.
    if rho[-1,0] == 1.0:
        psi_norm[-1] = 1.0

    return psi, psi_norm

