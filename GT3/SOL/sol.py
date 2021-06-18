#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:14:01 2018

@author: max
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, UnivariateSpline
from scipy.constants import elementary_charge
from shapely.geometry import Point, LineString
from shapely.ops import polygonize
from collections import namedtuple
from math import ceil
from GT3.Core.Functions.Cut import cut
from contours.quad import QuadContourGenerator
import sys


class Sol:
    def __init__(self, inp, core):
        # # NOTE: this assumes structured psi data on a regular grid. This way of getting R and Z is entirely inadeqaute for
        # # unstructured psi data
        # R = inp.psirz_exp[:, 0].reshape(-1, np.unique(inp.psirz_exp[:, 0]).size)
        # Z = inp.psirz_exp[:, 1].reshape(-1, np.unique(inp.psirz_exp[:, 1]).size)

        self.calc_sol_lines(inp, core)
        self.calc_sol_nT(inp, core)
        pass

    def calc_sol_lines(self, inp, core):
        c = QuadContourGenerator.from_rectilinear(core.psi_data.R[0], core.psi_data.Z[:, 0], core.psi_data.psi_norm)
        self.sol_lines = []
        self.sol_lines_cut = []

        #
        # first we need to make sure that the inp.sollines_psi_max specified in the input file doesn't
        # go outside of the first wall. For now, we will check for this as follows:
        #   1. draw the contours for psi_norm = inp.sollines_psi_max
        #   2. if there is only 1, make sure it only intersects the first wall no more than twice
        #   3. if there is only 1 and it intersects the first wall more than twice, raise an error
        #       that tells the user to reduce their sollines_psi_max value in the input file
        #   4. if there is more than 1 contour, then we need to determine which one is the correct one.
        #       for now, we will check that it's largest y value is higher than the magnetic axis,
        #       it's largest x value is to the right of the magnetic axis, it's lowest y value is lower
        #       than the magnetic axis, and it's lowest x value is to the left of the magnetic axis.
        #   5. Once the correct contour has been identified, do steps 2 and 3.
        #

        sollines_psi_max_contours = c.contour(inp.sollines_psi_max)
        num_lines = len(sollines_psi_max_contours)

        print('num_lines = ', num_lines)
        if num_lines == 1:
            # then this is probably the correct line. Check to see how many times it intersects
            # with the first wall
            num_wall_ints = len(LineString(sollines_psi_max_contours[0]).intersection(inp.wall_line))
            print('num_wall_ints = ', num_wall_ints)
            if num_wall_ints > 2:
                print('It looks like your sollines_psi_max value might be intersecting the wall.' \
                      'Try reducing it. Stopping.')
                raise Exception("It looks like your sollines_psi_max value might be intersecting the wall. Try reducing it.")
        else:
            for i, line in enumerate(sollines_psi_max_contours):
                max_x = np.amax(line[:, 0])
                min_x = np.amin(line[:, 0])
                max_y = np.amax(line[:, 1])
                min_y = np.amin(line[:, 1])
                if min_x < core.pts.axis.mag[0] < max_x and min_y < core.pts.axis.mag[1] < max_y:
                    # then this is probably the correct line. Check to see how many times it intersects
                    # with the first wall
                    num_wall_ints = len(LineString(line).intersection(inp.wall_line))
                    if num_wall_ints > 2:
                        print('It looks like your sollines_psi_max value might be intersecting the wall.' \
                              'Try reducing it. Stopping.')
                        raise Exception("It looks like your sollines_psi_max value might be intersecting the wall. Try reducing it.")
                    else:
                        break

        psi_pts = np.linspace(1, inp.sollines_psi_max, inp.num_sollines + 1, endpoint=True)[1:]

        for i, v in enumerate(psi_pts):
            num_lines = len(c.contour(v))
            if num_lines == 1:
                # this is probably the correct line
                self.sol_lines.append(LineString(c.contour(v)[0]))
            else:
                # we need to determine which contour to use
                for line in c.contour(v):
                    max_x = np.amax(line[:, 0])
                    min_x = np.amin(line[:, 0])
                    max_y = np.amax(line[:, 1])
                    min_y = np.amin(line[:, 1])
                    if min_x < core.pts.axis.mag[0] < max_x and min_y < core.pts.axis.mag[1] < max_y:
                        # then this is probably the correct line.
                        self.sol_lines.append(LineString(line))

        for line in self.sol_lines:
            # find intersection points with the wall
            int_pts = line.intersection(inp.wall_line)
            # cut line at intersection points
            cut_line = cut(line, line.project(int_pts[0], normalized=True))[1]
            cut_line = cut(cut_line, cut_line.project(int_pts[1], normalized=True))[0]
            self.sol_lines_cut.append(cut_line)

        # add wall intersection points from divertor legs and sol lines to wall_line.
        # This is necessary to prevent thousands of tiny triangles from forming if the
        # end of the flux line isn't exactly on top of the wall line.

        # add inboard seperatrix strike point

        # plt.plot(np.asarray(inp.wall_line)[:,0], np.asarray(inp.wall_line)[:,1])
        # plt.plot(np.asarray(core.lines.div.ib_long)[:,0], np.asarray(core.lines.div.ib_long)[:, 1])

        union = inp.wall_line.union(core.lines.div.ib_long)
        result = [geom for geom in polygonize(union)][0]
        inp.wall_line = LineString(result.exterior.coords)

        # add outboard seperatrix strike point
        union = inp.wall_line.union(core.lines.div.ob_long)
        result = [geom for geom in polygonize(union)][0]
        inp.wall_line = LineString(result.exterior.coords)

        # add sol line intersection points on inboard side
        # for some reason, union freaks out when I try to do inboard and outboard
        # at the same time.
        for num, line in enumerate(self.sol_lines):
            union = inp.wall_line.union(cut(line, 0.5)[0])
            result = [geom for geom in polygonize(union)][0]
            inp.wall_line = LineString(result.exterior.coords)

        # add sol line intersection points on outboard side
        for num, line in enumerate(self.sol_lines):
            union = inp.wall_line.union(cut(line, 0.5)[1])
            result = [geom for geom in polygonize(union)][0]
            inp.wall_line = LineString(result.exterior.coords)

    def calc_sol_nT(self, inp, core):
        # calculate spatial gradients for density and temperature along the seperatrix from dni/dr = dni/dpsi * dpsi/dr
        # specify the flux surface to get densities, temperatures, and their gradients
        sep_flx_surf = .98

        # calculate dni/dpsi and dTi/dpsi at the seperatrix
        # TODO: include both deuterium and tritium here
        ni_psi_fit = UnivariateSpline(core.psi.rho2psi(core.rho[:, 0]), core.n.i[:, 0], k=3, s=2.0)
        ne_psi_fit = UnivariateSpline(core.psi.rho2psi(core.rho[:, 0]), core.n.e[:, 0], k=3, s=2.0)
        Ti_psi_fit = UnivariateSpline(core.psi.rho2psi(core.rho[:, 0]), core.T.i.kev[:, 0], k=3, s=2.0)
        Te_psi_fit = UnivariateSpline(core.psi.rho2psi(core.rho[:, 0]), core.T.e.kev[:, 0], k=3, s=2.0)
        dni_dpsi_sep = ni_psi_fit.derivative()(sep_flx_surf)
        dTi_dpsi_sep = Ti_psi_fit.derivative()(sep_flx_surf)

        # calculate dpsidr everywhere (technically, we're calculating |dpsi/dr|. We don't care about the direction.
        dpsidR = np.abs(np.gradient(core.psi_data.psi_norm, core.psi_data.R[0, :], axis=1))
        dpsidZ = np.abs(np.gradient(core.psi_data.psi_norm, core.psi_data.Z[:, 0], axis=0))

        dpsidr = dpsidR + dpsidZ

        dpsidr_sep = griddata(np.column_stack((core.psi_data.R.flatten(), core.psi_data.Z.flatten())),
                              dpsidr.flatten(),
                              np.asarray(core.lines.sep.coords),
                              method='linear')

        # calculate dni/dr and dTi/dr at the seperatrix
        dnidr_sep = dni_dpsi_sep * dpsidr_sep
        dnedr_sep = dnidr_sep
        dTidr_sep = dTi_dpsi_sep * dpsidr_sep
        dTedr_sep = dTidr_sep

        num_sep_pts = len(np.asarray(core.lines.sep.coords))

        ni_sep = np.full(num_sep_pts, ni_psi_fit(sep_flx_surf))
        ne_sep = np.full(num_sep_pts, ne_psi_fit(sep_flx_surf))
        Ti_sep = np.full(num_sep_pts, Ti_psi_fit(sep_flx_surf)) * 1E3 * 1.6021E-19  # in Joules
        Te_sep = np.full(num_sep_pts, Te_psi_fit(sep_flx_surf)) * 1E3 * 1.6021E-19  # in Joules

        # calculate BT along the seperatrix
        BT_sep = inp.BT0 * core.pts.axis.mag[0] / np.asarray(core.lines.sep.coords)[:, 0]

        # remove a certain percentage of the seperatrix in the vicinity of the x-point
        # this fraction will be removed from both the inboard and the outboard sides of the x-point

        frac_to_remove = 0.1
        pts_to_remove = int(ceil(frac_to_remove * num_sep_pts))

        dnidr_sep_cut = dnidr_sep[pts_to_remove:-pts_to_remove]
        dnedr_sep_cut = dnedr_sep[pts_to_remove:-pts_to_remove]
        dTidr_sep_cut = dTidr_sep[pts_to_remove:-pts_to_remove]
        dTedr_sep_cut = dTedr_sep[pts_to_remove:-pts_to_remove]

        ni_sep_cut = ni_sep[pts_to_remove:-pts_to_remove]
        ne_sep_cut = ne_sep[pts_to_remove:-pts_to_remove]
        Ti_sep_cut = Ti_sep[pts_to_remove:-pts_to_remove]
        Te_sep_cut = Te_sep[pts_to_remove:-pts_to_remove]

        BT_sep_cut = BT_sep[pts_to_remove:-pts_to_remove]

        # define densities and temperatures along divertor legs
        # TODO: Specify these things in the input file
        ni_ib_wall = ni_sep_cut[0] * 10
        ni_ob_wall = ni_sep_cut[-1] * 10
        ni_ib = np.linspace(ni_ib_wall, ni_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        ni_ob = np.linspace(ni_sep_cut[-1], ni_ob_wall, inp.xi_ob_pts, endpoint=True)

        ne_ib_wall = ne_sep_cut[0] * 10
        ne_ob_wall = ne_sep_cut[-1] * 10
        ne_ib = np.linspace(ne_ib_wall, ne_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        ne_ob = np.linspace(ne_sep_cut[-1], ne_ob_wall, inp.xi_ob_pts, endpoint=True)

        Ti_ib_wall = Ti_sep_cut[0] / 4
        Ti_ob_wall = Ti_sep_cut[-1] / 4
        Ti_ib = np.linspace(Ti_ib_wall, Ti_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        Ti_ob = np.linspace(Ti_sep_cut[-1], Ti_ob_wall, inp.xi_ob_pts, endpoint=True)

        Te_ib_wall = Te_sep_cut[0] / 4
        Te_ob_wall = Te_sep_cut[-1] / 4
        Te_ib = np.linspace(Te_ib_wall, Te_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        Te_ob = np.linspace(Te_sep_cut[-1], Te_ob_wall, inp.xi_ob_pts, endpoint=True)

        print()
        print('#####################################')
        print(' divertor values')
        print(' ni_ib_wall = ', ni_ib_wall)
        print(' ni_ob_wall = ', ni_ob_wall)
        print(' ne_ib_wall = ', ne_ib_wall)
        print(' ne_ob_wall = ', ne_ob_wall)
        print(' Ti_ib_wall(ev) = ', Ti_ib_wall / 1.6021E-19)
        print(' Ti_ob_wall(ev) = ', Ti_ob_wall / 1.6021E-19)
        print(' Te_ib_wall(ev) = ', Te_ib_wall / 1.6021E-19)
        print(' Te_ob_wall(ev) = ', Te_ob_wall / 1.6021E-19)
        print('#####################################')
        print()

        # define density and temperature gradients along the inboard and outboard divertor legs
        dnidr_ib_wall = dnidr_sep_cut[0]
        dnidr_ob_wall = dnidr_sep_cut[-1]
        dnidr_ib = np.linspace(dnidr_ib_wall, dnidr_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        dnidr_ob = np.linspace(dnidr_sep_cut[-1], dnidr_ob_wall, inp.xi_ob_pts, endpoint=True)

        dnedr_ib_wall = dnedr_sep_cut[0]
        dnedr_ob_wall = dnedr_sep_cut[-1]
        dnedr_ib = np.linspace(dnedr_ib_wall, dnedr_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        dnedr_ob = np.linspace(dnedr_sep_cut[-1], dnedr_ob_wall, inp.xi_ob_pts, endpoint=True)

        dTidr_ib_wall = dTidr_sep_cut[0]
        dTidr_ob_wall = dTidr_sep_cut[-1]
        dTidr_ib = np.linspace(dTidr_ib_wall, dTidr_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        dTidr_ob = np.linspace(dTidr_sep_cut[-1], dTidr_ob_wall, inp.xi_ob_pts, endpoint=True)

        dTedr_ib_wall = dTedr_sep_cut[0]
        dTedr_ob_wall = dTedr_sep_cut[-1]
        dTedr_ib = np.linspace(dTedr_ib_wall, dTedr_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        dTedr_ob = np.linspace(dTedr_sep_cut[-1], dTedr_ob_wall, inp.xi_ob_pts, endpoint=True)

        BT_ib_wall = BT_sep_cut[0]
        BT_ob_wall = BT_sep_cut[-1]
        BT_ib = np.linspace(BT_ib_wall, BT_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        BT_ob = np.linspace(BT_sep_cut[-1], BT_ob_wall, inp.xi_ob_pts, endpoint=True)

        # combine inboard, seperatrix, and outboard points. These now comprise all the values in the xi direction.
        ni_xi = np.concatenate((ni_ib, ni_sep_cut, ni_ob))
        ne_xi = np.concatenate((ne_ib, ne_sep_cut, ne_ob))
        Ti_xi = np.concatenate((Ti_ib, Ti_sep_cut, Ti_ob))
        Te_xi = np.concatenate((Te_ib, Te_sep_cut, Te_ob))
        dnidr_xi = np.concatenate((dnidr_ib, dnidr_sep_cut, dnidr_ob))
        dnedr_xi = np.concatenate((dnedr_ib, dnedr_sep_cut, dnedr_ob))
        dTidr_xi = np.concatenate((dTidr_ib, dTidr_sep_cut, dTidr_ob))
        dTedr_xi = np.concatenate((dTedr_ib, dTedr_sep_cut, dTedr_ob))
        BT_xi = np.concatenate((BT_ib, BT_sep_cut, BT_ob))

        ib_leg_length = core.lines.div.ib.length
        ob_leg_length = core.lines.div.ob.length
        sep_length = core.lines.sep_closed.length
        ib_frac = ib_leg_length / (ib_leg_length + sep_length + ob_leg_length)
        sep_frac = sep_length / (ib_leg_length + sep_length + ob_leg_length)
        ob_frac = ob_leg_length / (ib_leg_length + sep_length + ob_leg_length)

        # specify the points along xi
        xi_ib_div = np.linspace(0,
                                frac_to_remove,
                                inp.xi_ib_pts,
                                endpoint=False)

        xi_sep = np.linspace(frac_to_remove,
                             1.0 - frac_to_remove,
                             len(dnidr_sep_cut),
                             endpoint=False)

        xi_ob_div = np.linspace(1.0 - frac_to_remove,
                                1.0,
                                inp.xi_ob_pts,
                                endpoint=True)

        xi_pts = np.concatenate((xi_ib_div, xi_sep, xi_ob_div))

        # model perpendicular particle and heat transport using Bohm Diffusion
        D_perp = Ti_xi / (16.0 * elementary_charge * BT_xi)
        Chi_perp = 5.0 * Ti_xi / (32.0 * elementary_charge * BT_xi)

        Gamma_perp = -D_perp * dnidr_xi
        Q_perp = -ni_xi * Chi_perp * dTidr_xi - \
                 3.0 * Ti_xi * D_perp * dnidr_xi

        delta_sol_n = D_perp * ni_xi / Gamma_perp
        delta_sol_T = Chi_perp / (Q_perp / (ni_xi * Ti_xi) - 3.0 * D_perp / delta_sol_n)
        delta_sol_E = 2 / 7 * delta_sol_T

        # now calculate densities and temperatures radially outward from the seperatrix for a distance
        # long enough that the wall is enclosed, so we can get densities and temperatures along the wall

        # first draw wall line through 2d strip model to get n, T along the line
        # we do this first so we can find the farthest point and make sure that we
        # make our SOL strip wide enough to go all the way to even the farthest point
        # on the wall
        wall_pts = np.asarray(inp.wall_line.xy).T
        # ib_int_pt = np.asarray(core.lines.div.ib.intersection(inp.wall_line).xy).T
        # ob_int_pt = core.lines.div.ob.intersection(inp.wall_line)
        ib_int_pt = np.asarray(core.pts.strike.ib)
        ob_int_pt = np.asarray(core.pts.strike.ob)

        wall_start_pos = np.where((wall_pts == ib_int_pt).all(axis=1))[0][0]
        wall_line_rolled = LineString(np.roll(wall_pts, -wall_start_pos, axis=0))

        wall_line_cut = cut(wall_line_rolled,
                            wall_line_rolled.project(Point(ob_int_pt), normalized=True))[0]

        # add points to wall line for the purpose of getting n, T along the wall. These points
        # won't be added to the main wall line or included in the triangulation.
        # for i, v in enumerate(np.linspace(0, 1, 300)):
        #    #interpolate along wall_line_cut to find point to add
        #    pt = wall_line_cut.interpolate(v, normalized=True)
        #    #add point to wall_line_cut
        #    union = wall_line_cut.union(pt)
        #    result = [geom for geom in polygonize(union)][0]
        #    wall_line_cut = LineString(result.exterior.coords)

        wall_nT_pts = np.asarray(wall_line_cut)
        num_wall_pts = len(wall_nT_pts)
        wall_pos_norm = np.zeros(num_wall_pts)
        wall_dist = np.zeros(num_wall_pts)

        for i, pt in enumerate(wall_nT_pts):
            wall_pt = Point(pt)
            sep_pt_pos = core.lines.ib2ob.project(Point(wall_pt), normalized=True)
            sep_pt = core.lines.ib2ob.interpolate(sep_pt_pos, normalized=True)
            wall_pos_norm[i] = wall_line_cut.project(wall_pt, normalized=True)
            wall_dist[i] = wall_pt.distance(sep_pt)

        r_max = np.amax(wall_dist) * 1.1  # the 1.1 multiplication just ensures that we go a bit farther than necessary
        twoptdiv_r_pts = 20

        r_pts = np.linspace(0, r_max, twoptdiv_r_pts)
        xi, r = np.meshgrid(xi_pts, r_pts)
        sol_ni = ni_xi * np.exp(-r / delta_sol_n)

        sol_ne = ne_xi * np.exp(-r / delta_sol_n)
        sol_Ti = Ti_xi * np.exp(-r / delta_sol_T)
        sol_Te = Te_xi * np.exp(-r / delta_sol_T)

        # set some minimum values
        min_n = 1E15
        min_T = 2.0 * 1.6021E-19
        sol_ni = np.where(sol_ni < min_n, min_n, sol_ni)
        sol_ne = np.where(sol_ne < min_n, min_n, sol_ne)
        sol_Ti = np.where(sol_Ti < min_T, min_T, sol_Ti)
        sol_Te = np.where(sol_Te < min_T, min_T, sol_Te)

        # draw sol lines through 2d strip model to get n, T along the lines
        sol_line_dist = np.zeros((len(xi_pts), len(self.sol_lines_cut)))
        sol_nT_pts = np.zeros((len(xi_pts), 2, len(self.sol_lines_cut)))
        for i, sol_line in enumerate(self.sol_lines_cut):
            for j, xi_val in enumerate(xi_pts):
                sol_pt = sol_line.interpolate(xi_val, normalized=True)
                sol_nT_pts[j, :, i] = np.asarray(sol_pt.xy).T
                sep_pt_pos = core.lines.ib2ob.project(sol_pt, normalized=True)
                sep_pt = core.lines.ib2ob.interpolate(sep_pt_pos, normalized=True)
                sol_line_dist[j, i] = sol_pt.distance(sep_pt)

        sol_line_ni = np.zeros((len(xi_pts), len(self.sol_lines_cut)))
        sol_line_ne = np.zeros((len(xi_pts), len(self.sol_lines_cut)))
        sol_line_Ti = np.zeros((len(xi_pts), len(self.sol_lines_cut)))
        sol_line_Te = np.zeros((len(xi_pts), len(self.sol_lines_cut)))

        pts_ni_sol = np.zeros((0, 3))
        pts_ne_sol = np.zeros((0, 3))
        pts_Ti_sol = np.zeros((0, 3))
        pts_Te_sol = np.zeros((0, 3))
        for i, sol_line in enumerate(self.sol_lines_cut):
            sol_line_ni[:, i] = griddata(np.column_stack((xi.flatten(), r.flatten())),
                                         sol_ni.flatten(),
                                         np.column_stack((np.linspace(0, 1, len(xi_pts)), sol_line_dist[:, i])),
                                         method='linear')
            sol_line_ne[:, i] = griddata(np.column_stack((xi.flatten(), r.flatten())),
                                         sol_ne.flatten(),
                                         np.column_stack((np.linspace(0, 1, len(xi_pts)), sol_line_dist[:, i])),
                                         method='linear')
            sol_line_Ti[:, i] = griddata(np.column_stack((xi.flatten(), r.flatten())),
                                         sol_Ti.flatten(),
                                         np.column_stack((np.linspace(0, 1, len(xi_pts)), sol_line_dist[:, i])),
                                         method='linear')
            sol_line_Te[:, i] = griddata(np.column_stack((xi.flatten(), r.flatten())),
                                         sol_Te.flatten(),
                                         np.column_stack((np.linspace(0, 1, len(xi_pts)), sol_line_dist[:, i])),
                                         method='linear')

            # append to SOL n, T arrays
            pts_ni_sol = np.vstack((pts_ni_sol, np.column_stack((sol_nT_pts[:, :, i], sol_line_ni[:, i]))))
            pts_ne_sol = np.vstack((pts_ne_sol, np.column_stack((sol_nT_pts[:, :, i], sol_line_ne[:, i]))))
            pts_Ti_sol = np.vstack((pts_Ti_sol, np.column_stack(
                (sol_nT_pts[:, :, i], sol_line_Ti[:, i] / 1.0E3 / 1.6021E-19))))  # converting back to kev
            pts_Te_sol = np.vstack((pts_Te_sol, np.column_stack(
                (sol_nT_pts[:, :, i], sol_line_Te[:, i] / 1.0E3 / 1.6021E-19))))  # converting back to kev

        self.sol_nT = namedtuple('sol_nT', 'ni ne Ti Te')(pts_ni_sol, pts_ne_sol, pts_Ti_sol, pts_Te_sol)

        wall_ni = griddata(np.column_stack((xi.flatten(), r.flatten())),
                           sol_ni.flatten(),
                           np.column_stack((wall_pos_norm, wall_dist)),
                           method='linear')
        wall_ne = griddata(np.column_stack((xi.flatten(), r.flatten())),
                           sol_ne.flatten(),
                           np.column_stack((wall_pos_norm, wall_dist)),
                           method='linear')
        wall_Ti = griddata(np.column_stack((xi.flatten(), r.flatten())),
                           sol_Ti.flatten(),
                           np.column_stack((wall_pos_norm, wall_dist)),
                           method='linear')
        wall_Te = griddata(np.column_stack((xi.flatten(), r.flatten())),
                           sol_Te.flatten(),
                           np.column_stack((wall_pos_norm, wall_dist)),
                           method='linear')

        # append to master arrays
        pts_ni_wall = np.column_stack((wall_nT_pts, wall_ni))
        pts_ne_wall = np.column_stack((wall_nT_pts, wall_ne))
        pts_Ti_wall = np.column_stack((wall_nT_pts, wall_Ti / 1.0E3 / 1.6021E-19))  # in kev
        pts_Te_wall = np.column_stack((wall_nT_pts, wall_Te / 1.0E3 / 1.6021E-19))  # in kev

        self.wall_nT = namedtuple('wall_nT', 'ni ne Ti Te')(pts_ni_wall, pts_ne_wall, pts_Ti_wall, pts_Te_wall)

        # # uncomment this if you're debugging
        # plt.contourf(xi, r, sol_Ti, 500)
        # plt.colorbar()
        # for i, v in enumerate(self.sol_lines_cut):
        #     plt.plot(xi_pts, sol_line_dist[:, i])
        # plt.plot(np.linspace(0, 1, num_wall_pts), wall_dist, color='black')
        # plt.show()
        # sys.exit()
