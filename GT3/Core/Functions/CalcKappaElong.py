#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from collections import namedtuple
from GT3.Core.Functions.DrawCoreLine import draw_core_line


def calc_kappa_elong(psi_data, sep_pts):
    # get kappa and elongation from the psi_norm=0.95 flux surface
    fs_line, fs_pts = draw_core_line(psi_data.R, psi_data.Z, psi_data.psi_norm, 0.95, sep_pts)
    elong_a = (fs_pts.top[1] - fs_pts.bot[1]) / (fs_pts.out[0] - fs_pts.inn[0])
    tri_a = (fs_pts.axis[0] - fs_pts.top[0]) / (fs_pts.out[0] - fs_pts.axis[0])

    # get kappa and elongation near the magnetic axis
    fs_line, fs_pts = draw_core_line(psi_data.R, psi_data.Z, psi_data.psi_norm, 0.05, sep_pts)

    # a small number isn't guaranteed to find a contour. If one is not found,
    # start near zero increase psi_norm until we find one
    if fs_line is None:
        for psi_val in np.linspace(0.05, 1, 100):
            fs_line, fs_pts = draw_core_line(psi_data.R, psi_data.Z, psi_data.psi_norm, psi_val, sep_pts)
            if fs_line is None:
                pass
            else:
                break

    elong_0 = (fs_pts.top[1] - fs_pts.bot[1]) / (fs_pts.out[0] - fs_pts.inn[0])
    tri_0 = 0

    kappa = namedtuple('kappa', 'axis sep')(elong_0, elong_a)
    tri = namedtuple('tri', 'axis sep')(tri_0, tri_a)
    return kappa, tri
