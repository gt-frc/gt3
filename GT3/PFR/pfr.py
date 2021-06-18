#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:04:43 2018

@author: max
"""

import numpy as np
from shapely.geometry import LineString
from shapely.ops import polygonize, linemerge
import sys
from contours.quad import QuadContourGenerator
from GT3.Core.Functions.Cut import cut
from deprecation import deprecated

@deprecated(deprecated_in="0.0.3", removed_in="0.0.4", details="To be replaced with a more robust version at a future date.")
class Pfr:
    def __init__(self, inp, core):

        # R = inp.psirz_exp[:, 0].reshape(-1, 65)
        # Z = inp.psirz_exp[:, 1].reshape(-1, 65)

        self.pfr_lines(inp, core)
        self.pfr_nT(inp, core)

    def pfr_lines(self, inp, core):

        c = QuadContourGenerator.from_rectilinear(core.psi_data.R[0], core.psi_data.Z[:, 0], core.psi_data.psi_norm)
        contours = c.contour(0.999)

        if len(contours) == 1:
            # then we're definitely dealing with a surface inside the seperatrix
            raise ValueError("Did not find PFR flux surface. Stopping.")
        else:
            # we need to find the surface that is contained within the private flux region
            for j, contour in enumerate(contours):
                # if the contour is entirely below the x-point and contour extends to both the left and the
                # right of the x-point horizontally, then the contour is almost certainly the desired PFR contour
                if np.amax(contour[:, 1]) < core.pts.xpt[1] and \
                        np.amin(contour[:, 0]) < core.pts.xpt[0] < np.amax(contour[:, 0]):
                    # then it's a probably a pfr flux surface, might need to add additional checks later

                    # make sure the line goes from inboard to outboard
                    if contour[-1, 0] < contour[0, 0]:
                        contour = np.flipud(contour)

                    # find cut points
                    cut_pt_ib = LineString(contour).intersection(inp.wall_line)[0]
                    cut_pt_ob = LineString(contour).intersection(inp.wall_line)[1]

                    dist1 = LineString(contour).project(cut_pt_ib, normalized=True)
                    cutline_temp = cut(LineString(contour), dist1)[1]

                    # reverse line point order so we can reliably find the second intersection point
                    cutline_temp_rev = LineString(np.flipud(np.asarray(cutline_temp.coords)))

                    dist2 = cutline_temp_rev.project(cut_pt_ob, normalized=True)
                    cutline_final_rev = cut(cutline_temp_rev, dist2)[1]

                    # reverse again for final pfr flux line
                    pfr_flux_line = LineString(np.flipud(np.asarray(cutline_final_rev.coords)))

                    # add pfr_line intersection points on inboard side
                    # for some reason, union freaks out when I try to do inboard and outboard
                    # at the same time.
                    union = inp.wall_line.union(cut(LineString(contour), 0.5)[0])
                    result = [geom for geom in polygonize(union)][0]
                    inp.wall_line = LineString(result.exterior.coords)

                    # add pfr line intersection points on outboard side
                    union = inp.wall_line.union(cut(LineString(contour), 0.5)[1])
                    result = [geom for geom in polygonize(union)][0]
                    inp.wall_line = LineString(result.exterior.coords)

                    # cut out pfr section of wall line
                    wall_pts = np.asarray(inp.wall_line.xy).T

                    wall_start_pos = np.where((wall_pts == cut_pt_ob).all(axis=1))[0][0]
                    wall_line_rolled = LineString(np.roll(wall_pts, -wall_start_pos, axis=0))
                    wall_line_cut_pfr = cut(wall_line_rolled,
                                            wall_line_rolled.project(cut_pt_ib, normalized=True))[0]

                    # create LineString with pfr line and section of wall line
                    self.pfr_line = linemerge((pfr_flux_line, wall_line_cut_pfr))
                    break

    def pfr_nT(self, inp, core):

        pfr_pts = np.asarray(self.pfr_line.xy).T
        pfr_ni = np.zeros(len(pfr_pts)) + inp.pfr_ni_val
        pfr_ne = np.zeros(len(pfr_pts)) + inp.pfr_ne_val
        pfr_Ti = np.zeros(len(pfr_pts)) + inp.pfr_Ti_val
        pfr_Te = np.zeros(len(pfr_pts)) + inp.pfr_Te_val

        # TODO: May need to convert these temperatures to keV
        pts_ni_pfr = np.column_stack((pfr_pts, pfr_ni))
        pts_ne_pfr = np.column_stack((pfr_pts, pfr_ne))
        pts_Ti_pfr = np.column_stack((pfr_pts, pfr_Ti))
        pts_Te_pfr = np.column_stack((pfr_pts, pfr_Te))

        brnd_pfr_dict = {}
        brnd_pfr_dict['ni'] = pts_ni_pfr
        brnd_pfr_dict['ne'] = pts_ne_pfr
        brnd_pfr_dict['Ti'] = pts_Ti_pfr
        brnd_pfr_dict['Te'] = pts_Te_pfr

        # grid_x, grid_y = np.mgrid[1:2.5:500j, -1.5:1.5:500j]
        # ni_for_plot = griddata(self.ni_pts[:, :-1], self.ni_pts[:, -1], (grid_x, grid_y))
        # Ti_for_plot = griddata(self.Ti_kev_pts[:, :-1], self.Ti_kev_pts[:, -1], (grid_x, grid_y))
        # plt.contourf(grid_x, grid_y, np.log10(Ti_for_plot), 500)
        # plt.colorbar()
        # sys.exit()
        pass
