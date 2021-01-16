#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from contours.quad import QuadContourGenerator
from GT3.Core.Functions.IsClose import isclose
from GT3.Core.Functions.DrawCounterLine import draw_contour_line
from collections import namedtuple
from shapely.geometry import LineString
import numpy as np

def draw_core_line(R, Z, psinorm, psi_val, sep_pts):
    # create contour generator
    c = QuadContourGenerator.from_rectilinear(R[0], Z[:, 0], psinorm)

    # draw contours with psi_val
    contours = c.contour(psi_val)

    if len(contours) == 0 and isclose(psi_val, 0, abs_tol=1E-9):
        # This is probably the magnetic axis
        fs_line = None
        fs_inn_pt = None
        fs_out_pt = None
        fs_top_pt = None
        fs_bot_pt = None
        fs_axis = None
    elif len(contours) == 0 and not isclose(psi_val, 0, abs_tol=1E-9):
        # This either means that either:
        #    A) psi is a value not present in the psi data or
        #    B) psi is very close to the magnetic axis, but is too small for the contours
        #       package to give us a contour. This can happen if you you have an very fine
        #       radial mesh in the vicinity of the magnetic axis.
        # Passing back None for now, but should probably raise an exception in the future  # TODO
        fs_line = None
        fs_inn_pt = None
        fs_out_pt = None
        fs_top_pt = None
        fs_bot_pt = None
        fs_axis = None
    else:
        if len(contours) == 1:
            # then we're definitely dealing with a surface inside the seperatrix
            x, y = draw_contour_line(R, Z, psinorm, psi_val, 0)
        else:
            # we need to find which of the surfaces is inside the seperatrix
            for j, line in enumerate(contours):
                x, y = draw_contour_line(R, Z, psinorm, psi_val, j)

                if (np.amax(x) < np.amax(sep_pts[:, 0]) and
                    np.amin(x) > np.amin(sep_pts[:, 0]) and
                    np.amax(y) < np.amax(sep_pts[:, 1]) and
                    np.amin(y) > np.amin(sep_pts[:, 1])):
                    # then it's an internal flux surface
                    break

        pts = np.column_stack((x, y))
        fs_line = LineString(pts)
        fs_out_pt = pts[np.argmax(pts, axis=0)[0]]
        fs_inn_pt = pts[np.argmin(pts, axis=0)[0]]
        fs_top_pt = pts[np.argmax(pts, axis=0)[1]]
        fs_bot_pt = pts[np.argmin(pts, axis=0)[1]]
        fs_axis = [(fs_out_pt[0]+fs_inn_pt[0])/2, (fs_out_pt[1]+fs_inn_pt[1])/2]

    fs_pts = namedtuple('fs_pts', 'inn out top bot axis')(
        fs_inn_pt,
        fs_out_pt,
        fs_top_pt,
        fs_bot_pt,
        fs_axis
    )
    return fs_line, fs_pts