#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:14:01 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib._cntr as cntr
from scipy.interpolate import griddata, UnivariateSpline
from scipy.constants import elementary_charge
from shapely.geometry import Point, LineString
from shapely.ops import polygonize
from collections import namedtuple
from math import ceil
import sys

def draw_line(R, Z, array, val, pathnum):
    res = cntr.Cntr(R, Z, array).trace(val)[pathnum]
    x = res[:, 0]
    y = res[:, 1]
    return x, y

def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= 1.0:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p), normalized=True)
        if pd == distance:
            return [
                LineString(coords[:i+1]), 
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance, normalized=True)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]), 
                LineString([(cp.x, cp.y)] + coords[i:])]

class exp_sol_brnd():
    def __init__(self, inp, core):
        
        R = inp.psirz_exp[:, 0].reshape(-1, 65)
        Z = inp.psirz_exp[:, 1].reshape(-1, 65)
        
        self.calc_sol_lines(inp, R, Z, core)
        self.calc_sol_nT(inp, R, Z, core)
        pass
    
    def calc_sol_lines(self, inp, R, Z, core):
        # find value of psi at outside of what we're going to call the SOL
        self.sol_lines = []
        self.sol_lines_cut = []

        psi_pts = np.linspace(1, inp.sollines_psi_max, inp.num_sollines+1, endpoint=True)[1:]
        for i, v in enumerate(psi_pts):
            num_lines = int(len(cntr.Cntr(R, Z, core.psi_norm_raw).trace(v))/2)
            if num_lines==1:
                # then we're definitely dealing with a surface inside the seperatrix
                x, y = draw_line(R, Z, core.psi_norm_raw, v, 0)
                self.sol_lines.append(LineString(np.column_stack((x, y))))
            else:
                # TODO:
                pass
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
        union = inp.wall_line.union(core.ib_div_line)
        result = [geom for geom in polygonize(union)][0]
        inp.wall_line = LineString(result.exterior.coords)

        # add outboard seperatrix strike point
        union = inp.wall_line.union(core.ob_div_line)
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
            
    def calc_sol_nT(self, inp, R, Z, core):
        # calculate spatial gradients for density and temperature along the seperatrix from dni/dr = dni/dpsi * dpsi/dr
        # specify the flux surface to get densities, temperatures, and their gradients
        sep_flx_surf = 0.98

        # calculate dni/dpsi and dTi/dpsi at the seperatrix
        ni_psi_fit = UnivariateSpline(core.rho2psi(inp.ni_data[:, 0]), inp.ni_data[:, 1], k=3, s=2.0)
        ne_psi_fit = UnivariateSpline(core.rho2psi(inp.ne_data[:, 0]), inp.ne_data[:, 1], k=3, s=2.0)
        Ti_psi_fit = UnivariateSpline(core.rho2psi(inp.Ti_data[:, 0]), inp.Ti_data[:, 1], k=3, s=2.0)
        Te_psi_fit = UnivariateSpline(core.rho2psi(inp.Te_data[:, 0]), inp.Te_data[:, 1], k=3, s=2.0)
        dni_dpsi_sep = ni_psi_fit.derivative()(sep_flx_surf)
        dTi_dpsi_sep = Ti_psi_fit.derivative()(sep_flx_surf)

        # calculate dpsidr everywhere (technically, we're calculating |dpsi/dr|. We don't care about the direction.
        dpsidR = np.abs(np.gradient(core.psi_norm_raw, R[0, :], axis=1))
        dpsidZ = np.abs(np.gradient(core.psi_norm_raw, Z[:, 0], axis=0))

        dpsidr = dpsidR + dpsidZ

        dpsidr_sep = griddata(np.column_stack((R.flatten(), Z.flatten())),
                              dpsidr.flatten(),
                              (core.main_sep_pts[:,0], core.main_sep_pts[:,1]),
                              method='linear')

        # calculate dni/dr and dTi/dr at the seperatrix
        dnidr_sep = dni_dpsi_sep * dpsidr_sep
        dnedr_sep = dnidr_sep
        dTidr_sep = dTi_dpsi_sep * dpsidr_sep
        dTedr_sep = dTidr_sep

        num_sep_pts = len(core.main_sep_pts)

        ni_sep = np.full(num_sep_pts, ni_psi_fit(sep_flx_surf))
        ne_sep = np.full(num_sep_pts, ne_psi_fit(sep_flx_surf))
        Ti_sep = np.full(num_sep_pts, Ti_psi_fit(sep_flx_surf)) * 1E3 * 1.6021E-19  # in Joules
        Te_sep = np.full(num_sep_pts, Te_psi_fit(sep_flx_surf)) * 1E3 * 1.6021E-19  # in Joules

        #calculate BT along the seperatrix
        BT_sep = inp.BT0 * core.m_axis[0] / core.main_sep_pts[:,0]

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
        ni_ib_wall = ni_sep_cut[0] * 1.5
        ni_ob_wall = ni_sep_cut[-1] * 1.5
        ni_ib = np.linspace(ni_ib_wall, ni_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        ni_ob = np.linspace(ni_sep_cut[-1], ni_ob_wall, inp.xi_ob_pts, endpoint=True)
        
        ne_ib_wall = ne_sep_cut[0] * 1.5
        ne_ob_wall = ne_sep_cut[-1] * 1.5
        ne_ib = np.linspace(ne_ib_wall, ne_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        ne_ob = np.linspace(ne_sep_cut[-1], ne_ob_wall, inp.xi_ob_pts, endpoint=True)
        
        Ti_ib_wall = Ti_sep_cut[0] * 1.5
        Ti_ob_wall = Ti_sep_cut[-1] * 1.5
        Ti_ib = np.linspace(Ti_ib_wall, Ti_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        Ti_ob = np.linspace(Ti_sep_cut[-1], Ti_ob_wall, inp.xi_ob_pts, endpoint=True)

        Te_ib_wall = Te_sep_cut[0] * 1.5
        Te_ob_wall = Te_sep_cut[-1] * 1.5
        Te_ib = np.linspace(Te_ib_wall, Te_sep_cut[0], inp.xi_ib_pts, endpoint=False)
        Te_ob = np.linspace(Te_sep_cut[-1], Te_ob_wall, inp.xi_ob_pts, endpoint=True)

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
        
        ib_leg_length = core.ib_div_line_cut.length
        ob_leg_length = core.ob_div_line_cut.length
        sep_length = core.main_sep_line_closed.length
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
        delta_sol_T = Chi_perp / (Q_perp/(ni_xi*Ti_xi) - 3.0*D_perp/delta_sol_n)
        delta_sol_E = 2/7*delta_sol_T

        # now calculate densities and temperatures radially outward from the seperatrix for a distance
        # long enough that the wall is enclosed, so we can get densities and temperatures along the wall
        # r_max~0.5 is enough for DIII-D. May not be enough for ITER.
        # TODO: calculate the maximum wall distance before setting this parameter

        r_max = 0.5
        twoptdiv_r_pts = 20
        
        r_pts = np.linspace(0, r_max, twoptdiv_r_pts)
        xi, r = np.meshgrid(xi_pts, r_pts)
        sol_ni = ni_xi * np.exp(-r/delta_sol_n)
        sol_ne = ne_xi * np.exp(-r/delta_sol_n)
        sol_Ti = Ti_xi * np.exp(-r/delta_sol_T)
        sol_Te = Te_xi * np.exp(-r/delta_sol_T)

        # set some minimum values

        min_n = 1E15
        min_T = 2.0*1.6021E-19
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
                sep_pt_pos = core.entire_sep_line.project(sol_pt, normalized=True)
                sep_pt = core.entire_sep_line.interpolate(sep_pt_pos, normalized=True)
                sol_line_dist[j, i] = sol_pt.distance(sep_pt)
        
        sol_line_ni = np.zeros((len(xi_pts), len(self.sol_lines_cut)))
        sol_line_ne = np.zeros((len(xi_pts), len(self.sol_lines_cut)))
        sol_line_Ti = np.zeros((len(xi_pts), len(self.sol_lines_cut)))
        sol_line_Te = np.zeros((len(xi_pts), len(self.sol_lines_cut)))

        pts_ni_sol = np.zeros((0,3))
        pts_ne_sol = np.zeros((0,3))
        pts_Ti_sol = np.zeros((0,3))
        pts_Te_sol = np.zeros((0,3))
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
            pts_Ti_sol = np.vstack((pts_Ti_sol, np.column_stack((sol_nT_pts[:, :, i], sol_line_Ti[:, i]/1.0E3/1.6021E-19))))  # converting back to kev
            pts_Te_sol = np.vstack((pts_Te_sol, np.column_stack((sol_nT_pts[:, :, i], sol_line_Te[:, i]/1.0E3/1.6021E-19))))  # converting back to kev

        sol_nT_dict = {}
        sol_nT_dict['ni'] = pts_ni_sol
        sol_nT_dict['ne'] = pts_ne_sol
        sol_nT_dict['Ti'] = pts_Ti_sol
        sol_nT_dict['Te'] = pts_Te_sol
        self.sol_nT = namedtuple('sol_nT', sol_nT_dict.keys())(*sol_nT_dict.values())

        # draw wall line through 2d strip model to get n, T along the line
        wall_pts = np.asarray(inp.wall_line.xy).T
        ib_int_pt = np.asarray(core.ib_div_line.intersection(inp.wall_line).xy).T
        ob_int_pt = core.ob_div_line.intersection(inp.wall_line)
        wall_start_pos = np.where((wall_pts == ib_int_pt).all(axis=1))[0][0]
        wall_line_rolled = LineString(np.roll(wall_pts, -wall_start_pos, axis=0))
        wall_line_cut = cut(wall_line_rolled, 
                            wall_line_rolled.project(ob_int_pt, normalized=True))[0]
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
            sep_pt_pos = core.entire_sep_line.project(Point(wall_pt), normalized=True)
            sep_pt = core.entire_sep_line.interpolate(sep_pt_pos, normalized=True)
            wall_pos_norm[i] = wall_line_cut.project(wall_pt, normalized=True)
            wall_dist[i] = wall_pt.distance(sep_pt)
        
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

        wall_nT_dict = {}
        wall_nT_dict['ni'] = pts_ni_wall
        wall_nT_dict['ne'] = pts_ne_wall
        wall_nT_dict['Ti'] = pts_Ti_wall
        wall_nT_dict['Te'] = pts_Te_wall
        self.wall_nT = namedtuple('wall_nT', wall_nT_dict.keys())(*wall_nT_dict.values())

        # plt.contourf(xi, r, sol_Ti, 500)
        # plt.colorbar()
        # for i, v in enumerate(self.sol_lines_cut):
        #     plt.plot(xi_pts, sol_line_dist[:, i])
        # plt.plot(np.linspace(0, 1, num_wall_pts), wall_dist, color='black')
        # plt.show()
        # sys.exit()
