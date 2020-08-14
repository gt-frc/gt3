#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from contours.quad import QuadContourGenerator
from shapely.geometry import LineString
import numpy as np


def calc_core_lines_ntrl(core):
    c = QuadContourGenerator.from_rectilinear(core.psi_data.R[0], core.psi_data.Z[:, 0], core.psi_data.psi_norm)

    rhovals = np.linspace(0.7, 1, 5, endpoint=False)
    psivals = core.rho2psinorm(rhovals)

    core_lines_ntrl = []
    for i, psival in enumerate(psivals):
        contours = c.contour(psival)
        # determine how many surfaces have that psi value
        num_lines = len(contours)

        if num_lines == 1:
            # then we're definitely dealing with a surface inside the seperatrix
            core_lines_ntrl.append(LineString(contours[0]))
        else:
            # we need to find which of the surfaces is inside the seperatrix
            for contour in contours:
                x, y = contour[:, 0], contour[:, 1]
                if (
                        np.amax(x) < np.amax(np.asarray(core.lines.sep.coords)[:, 0]) and
                        np.amin(x) > np.amin(np.asarray(core.lines.sep.coords)[:, 0]) and
                        np.amax(y) < np.amax(np.asarray(core.lines.sep.coords)[:, 1]) and
                        np.amin(y) > np.amin(np.asarray(core.lines.sep.coords)[:, 1])
                ):
                    # then it's an internal flux surface
                    core_lines_ntrl.append(LineString(contour))
                    break

    return core_lines_ntrl
