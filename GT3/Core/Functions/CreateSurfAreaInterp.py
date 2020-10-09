#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from GT3.Core.Functions.DrawCoreLine import draw_core_line
from scipy.interpolate import UnivariateSpline
from shapely.geometry import LineString
from math import pi, sqrt


def create_surf_area_interp(rho2psinorm, psinorm2rho, psi_data, sep_pts, R0_a, a, sep_val):
    rho_vals = np.linspace(0, sep_val, 50, endpoint=True)
    r_vals = rho_vals * a
    psi_vals = rho2psinorm(rho_vals)

    surf_area = np.zeros(rho_vals.shape)
    # skip the first value. First value is zero.
    for i, psi_val in enumerate(psi_vals[1:]):
        if psi_val < 1:
            try:
                # try to get the actual length of the fs trace
                fs_line, fs_pts = draw_core_line(psi_data.R, psi_data.Z, psi_data.psi_norm, psi_val, sep_pts)
                surf_area[i + 1] = fs_line.length * 2 * pi * fs_pts.axis[0]
            except:
                # if too close to the seperatrix, it may return None. In this case, estimate perimeter using Ramanujan
                # TODO: Improve this estimation. Also, see comments in draw_core_line
                rho_val = psinorm2rho(psi_val)
                a_el = rho_val * a  # a of the ellipse
                b_el = rho_val * a * 1.5  # b of the ellipse, assumed kappa of 1.5. This needs to be improved

                surf_area[i + 1] = pi * (3*(a_el+b_el) - sqrt((3*a_el + b_el)*(a_el+3*b_el)))

        else:
            surf_area[i+1] = LineString(sep_pts).length * 2*pi*R0_a

    # For some reason, psi_vals will sometimes not be monotonically increasing, especially near the magnetic axis.
    # This prevents UnivariateSpline from working. To prevent this, we're just going to delete any points that are
    # lower than the previous point, prior to doing the fit. As long as the psi data isn't too wonky, this should be
    # fine.
    psinormvals_mi = []  # psivals that are monotonically increasing
    rvals_mi = []  # rvals corresponding to monotonically increasing psi
    rhovals_mi = []  # rhovals corresponding to monotonically increasing psi
    surfarea_mi = []  # surf_area corresponding to monotonically increasing psi
    for i, psi_val in enumerate(psi_vals):
        if i == 0:
            rvals_mi.append(0)
            rhovals_mi.append(0)
            psinormvals_mi.append(0)
            surfarea_mi.append(0)
        elif psi_val > psinormvals_mi[-1]:
            rvals_mi.append(r_vals[i])
            rhovals_mi.append(rho_vals[i])
            psinormvals_mi.append(psi_vals[i])
            surfarea_mi.append(surf_area[i])

    rvals_mi = np.asarray(rvals_mi)
    rhovals_mi = np.asarray(rhovals_mi)
    psinormvals_mi = np.asarray(psinormvals_mi)
    surfarea_mi = np.asarray(surfarea_mi)

    r2sa = UnivariateSpline(rvals_mi, surfarea_mi, k=3, s=0)
    rho2sa = UnivariateSpline(rhovals_mi, surfarea_mi, k=3, s=0)
    psinorm2sa = UnivariateSpline(psinormvals_mi, surfarea_mi, k=3, s=0)
    # r2sa = UnivariateSpline(r_vals, surf_area, k=3, s=0)
    # rho2sa = UnivariateSpline(rho_vals, surf_area, k=3, s=0)
    # psinorm2sa = UnivariateSpline(psi_vals, surf_area, k=3, s=0)

    return r2sa, rho2sa, psinorm2sa
