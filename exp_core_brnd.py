#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 13:22:31 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib._cntr as cntr
from scipy.interpolate import griddata, UnivariateSpline, interp2d
from scipy.constants import elementary_charge
from shapely.geometry import Point, LineString
from shapely.ops import polygonize, linemerge
import sys
from math import atan2, pi, ceil, sin, cos


def draw_contour_line(R, Z, array, val, pathnum):
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


def draw_core_line(R, Z, psi, psi_val, sep_pts):
    num_lines = int(len(cntr.Cntr(R, Z, psi).trace(psi_val))/2)
    if num_lines == 1:
        # then we're definitely dealing with a surface inside the seperatrix
        x, y = draw_contour_line(R, Z, psi, psi_val, 0)
    else:
        # we need to find which of the surfaces is inside the seperatrix
        for j, line in enumerate(cntr.Cntr(R, Z, psi).trace(psi_val)[:num_lines]):
            x, y = draw_contour_line(R, Z, psi, psi_val, j)
            if (np.amax(x) < np.amax(sep_pts[:, 0]) and
                np.amin(x) > np.amin(sep_pts[:, 0]) and
                np.amax(y) < np.amax(sep_pts[:, 1]) and
                np.amin(y) > np.amin(sep_pts[:, 1])):
                # then it's an internal flux surface
                break
    pts = np.column_stack((x, y))
    line = LineString(pts)
    out_pt = pts[np.argmax(pts, axis=0)[0]]
    in_pt = pts[np.argmin(pts, axis=0)[0]]
    top_pt = pts[np.argmax(pts, axis=0)[1]]
    bot_pt = pts[np.argmin(pts, axis=0)[1]]
    fs_axis = [(out_pt[0]+in_pt[0])/2, (out_pt[1]+in_pt[1])/2]
    return line, fs_axis


class exp_core_brnd():
    def __init__(self, inp, ntrl_switch):
        
        R_psi = inp.psirz_exp[:, 0].reshape(-1, 65)
        Z_psi = inp.psirz_exp[:, 1].reshape(-1, 65)
        psi = inp.psirz_exp[:, 2].reshape(-1, 65)
        B_pol_raw = np.sqrt((np.gradient(psi, axis=1)/R_psi)**2 + (-np.gradient(psi, axis=0)/R_psi)**2)
        self.sep_lines(inp, R_psi, Z_psi, psi)
        # self.core_lines_main(inp, R_psi, Z_psi, psi)
        self.core_main(inp, R_psi, Z_psi, psi, B_pol_raw)
        if ntrl_switch == 2:
            self.core_lines_ntrl(inp, R_psi, Z_psi, psi)
            self.core_nT_ntrl(inp, R_psi, Z_psi, psi)
        self.xsec(inp)
    
    def update_ntrl_data(self, n_n_slow, n_n_thermal, izn_rate_slow, izn_rate_thermal):
        self.n_n_slow = n_n_slow
        self.n_n_thermal = n_n_thermal
        self.n_n_total = n_n_slow + n_n_thermal
        self.izn_rate_slow = izn_rate_slow
        self.izn_rate_thermal = izn_rate_thermal
        self.izn_rate_total = izn_rate_slow + izn_rate_thermal

    def update_Lz_data(self, Lz_slow, dLzdT_slow, Lz_thermal, dLzdT_thermal):
        self.Lz_slow = Lz_slow
        self.dLzdT_slow = dLzdT_slow
        self.Lz_thermal = Lz_thermal
        self.dLzdT_thermal = dLzdT_thermal

    def update_imprad_data(self, n_n_slow, n_n_thermal, izn_rate_slow, izn_rate_thermal):
        self.n_n_slow = n_n_slow
        self.n_n_thermal = n_n_thermal
        self.n_n_total = n_n_slow + n_n_thermal
        self.izn_rate_slow = izn_rate_slow
        self.izn_rate_thermal = izn_rate_thermal
        self.izn_rate_total = izn_rate_slow + izn_rate_thermal
    
    def sep_lines(self, inp, R, Z, psi):
        # find x-point location  
        dpsidR = np.gradient(psi, R[0, :], axis=1)
        dpsidZ = np.gradient(psi, Z[:, 0], axis=0)
        d2psidR2 = np.gradient(dpsidR, R[0, :], axis=1)
        d2psidZ2 = np.gradient(dpsidZ, Z[:, 0], axis=0)
        
        # find line(s) where dpsidR=0
        self.dpsidR_0 = cntr.Cntr(R, Z, dpsidR).trace(0.0)
        # find line(s) where dpsidZ=0
        self.dpsidZ_0 = cntr.Cntr(R, Z, dpsidZ).trace(0.0)
    
        for i, path1 in enumerate(self.dpsidR_0):
            for j, path2 in enumerate(self.dpsidZ_0):
                try:
                    # find intersection points between curves for dpsidR=0 and dpsidZ=0
                    ints = LineString(path1).intersection(LineString(path2))
                    # if there is only one intersection ('Point'), then we're probably not
                    # dealing with irrelevant noise in psi
                    if ints.type == 'Point':
                        # check if local maximum or minimum
                        d2psidR2_pt = griddata(np.column_stack((R.flatten(), Z.flatten())), 
                                               d2psidR2.flatten(), 
                                               [ints.x, ints.y], 
                                               method='cubic')
                        d2psidZ2_pt = griddata(np.column_stack((R.flatten(), Z.flatten())), 
                                               d2psidZ2.flatten(), 
                                               [ints.x, ints.y], 
                                               method='cubic')
                        
                        if d2psidR2_pt > 0 and d2psidZ2_pt > 0:
                            # we've found the magnetic axis
                            self.m_axis = np.array([ints.x, ints.y])
                        elif d2psidR2_pt < 0 and d2psidZ2_pt < 0:
                            # we've found a magnet. Do nothing.
                            pass
                        elif ints.y<0:
                            # we've probably found our x-point, although this isn't super robust
                            # and obviously only applies to a single-diverted, lower-null configuration
                            # TODO: make this more robust, I could easily see this failing on some shots
                            self.xpt = np.array([ints.x, ints.y])
                            
                        # uncomment this line when debugging
                        # print list(ints.coords), d2psidR2(ints.x, ints.y), d2psidZ2(ints.x, ints.y)
                except:
                    pass
        
        # normalize psi
        psi_shift = psi + abs(np.amin(psi))  # set center to zero
        psi_shift_xpt = griddata(np.column_stack((R.flatten(), Z.flatten())), 
                                 psi_shift.flatten(), 
                                 self.xpt, 
                                 method='cubic')
        # psi_shift_xpt = interp2d(R, Z, psi_shift, kind='linear')(xpt[0], xpt[1])  # get new value at sep
        self.psi_norm_raw = psi_shift / psi_shift_xpt
        
        # create lines for seperatrix and divertor legs of seperatrix
        num_lines = int(len(cntr.Cntr(R, Z, self.psi_norm_raw).trace(1.0))/2)
        if num_lines == 1:
            # in this case, the contour points that matplotlib returned constitute
            # a single line from inboard divertor to outboard divertor. We need to
            # add in the x-point in at the appropriate locations and split into a
            # main and a lower seperatrix line, each of which will include the x-point.
            x_psi , y_psi = draw_contour_line(R, Z, self.psi_norm_raw, 1.0, 0)
            
            loc1 = np.argmax(y_psi > self.xpt[1])
            loc2 = len(y_psi) - np.argmin(y_psi[::-1] < self.xpt[1])
    
            x_psi = np.insert(x_psi, (loc1, loc2), self.xpt[0])
            y_psi = np.insert(y_psi, (loc1, loc2), self.xpt[1])
            
            psi_1_pts = np.column_stack((x_psi, y_psi))
            self.main_sep_pts = psi_1_pts[loc1:loc2+1, :]
            self.main_sep_line = LineString(self.main_sep_pts[:-1])
            self.main_sep_line_closed = LineString(self.main_sep_pts)

            # get the inboard and outboard divertor legs seperately. This is so that
            # everything that includes the x-point can start with the x-point, which
            # elliminates the risk of tiny triangles in the vicinity of the x-point
            self.inboard_div_sep = np.flipud(psi_1_pts[:loc1+1])
            self.outboard_div_sep = psi_1_pts[loc2+1:]

            # cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.inboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ib_div_line = line
            self.ib_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]
            # self.ib_div_line_cut = line
            # TODO: add point to wall line
            
            # cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.outboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ob_div_line = line
            self.ob_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]
            
            ib_div_pts = np.flipud(np.asarray(self.ib_div_line_cut.xy).T)
            sep_pts = np.asarray(self.main_sep_line.xy).T
            ob_div_pts = np.asarray(self.ob_div_line_cut.xy).T
            
            entire_sep_pts = np.vstack((ib_div_pts, sep_pts[1:, :], ob_div_pts))
            self.entire_sep_line = LineString(entire_sep_pts)

            # these are useful later
            self.obmp_pt = self.main_sep_pts[np.argmax(self.main_sep_pts, axis=0)[0]]
            self.ibmp_pt = self.main_sep_pts[np.argmin(self.main_sep_pts, axis=0)[0]]
            self.top_pt = self.main_sep_pts[np.argmax(self.main_sep_pts, axis=0)[1]]
            self.bot_pt = self.main_sep_pts[np.argmin(self.main_sep_pts, axis=0)[1]]
            self.geo_axis = [(self.obmp_pt[0] + self.ibmp_pt[0])/2, (self.obmp_pt[1] + self.ibmp_pt[1])/2]
            self.R0_a = self.geo_axis[0]
            # TODO: Is this how a is actually defined?
            # a is used by nbeams. I'm not sure it's used anywhere else.
            self.a = self.obmp_pt[0] - self.R0_a
            # TODO: add point to wall line
                
        elif num_lines==2:
            # in this case, we have a lower seperatrix trace (line 0), and a main
            # seperatrix trace (line 1).
            
            # first do lower seperatrix line
            x_psi, y_psi = draw_contour_line(R, Z, self.psi_norm_raw, 1.0, 0)
            loc = np.argmax(x_psi > self.xpt[0])
            
            x_psi = np.insert(x_psi, loc, self.xpt[0])
            y_psi = np.insert(y_psi, loc, self.xpt[1])
            psi_1_pts = np.column_stack((x_psi, y_psi))
            
            self.inboard_div_sep = np.flipud(psi_1_pts[:loc+1])
            self.outboard_div_sep = psi_1_pts[loc+1:]
            
            #cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.inboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ib_div_line = line
            self.ib_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]

            #cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.outboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ob_div_line = line
            self.ob_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]
            #TODO: add point to wall line

            #now to main seperatrix line
            x_psi, y_psi = draw_contour_line(R, Z, self.psi_norm_raw, 1.0, 1)
            self.main_sep_pts = np.insert(np.column_stack((x_psi, y_psi)), 0, self.xpt, axis=0)
            self.main_sep_line = LineString(self.main_sep_pts[:-1])
            self.main_sep_line_closed = LineString(self.main_sep_pts)
            
            entire_sep_pts = np.vstack((ib_div_pts, sep_pts[1:, :], ob_div_pts))
            self.entire_sep_line = LineString(entire_sep_pts)
            #now clean up the lines by removing any points that are extremely close
            #to the x-point 
            #TODO: 

    def core_lines_main(self, inp, R, Z, psi):                    
        #define lines to be used for the main computational grid                    
        self.core_main_lines = []
        psi_pts_main = np.concatenate((np.linspace(0, 0.8, 20, endpoint=False), np.linspace(0.8, 1.0, 20, endpoint=False)))
        for i, v in enumerate(psi_pts_main):
            num_lines = int(len(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v))/2)
            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x, y = draw_contour_line(R, Z, self.psi_norm_raw, v, 0)
                self.core_main_lines.append(LineString(np.column_stack((x[:-1], y[:-1]))))
            else:
                #we need to find which of the surfaces is inside the seperatrix
                for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)[:num_lines]):
                #for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)):
                    x, y = draw_contour_line(R, Z, self.psi_norm_raw, v, j)
                    if (np.amax(x) < np.amax(self.main_sep_pts[:, 0]) and \
                        np.amin(x) > np.amin(self.main_sep_pts[:, 0]) and \
                        np.amax(y) < np.amax(self.main_sep_pts[:, 1]) and \
                        np.amin(y) > np.amin(self.main_sep_pts[:, 1])):
                        #then it's an internal flux surface
                        self.core_main_lines.append(LineString(np.column_stack((x[:-1], y[:-1]))))
                        break

    def core_lines_ntrl(self, inp, R, Z, psi):
        #define lines to be used in the neutrals calculation
        self.core_ntrl_lines = []
        psi_pts_ntrl = np.linspace(inp.edge_rho_ntrl, 1, inp.rhopts_edge_ntrl, endpoint=False)
        for i, v in enumerate(psi_pts_ntrl):
            num_lines = int(len(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v))/2)
            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x, y = draw_contour_line(R, Z, self.psi_norm_raw, v, 0)
                self.core_lines.append(LineString(np.column_stack((x[:-1], y[:-1]))))
            else:
                #we need to find which of the surfaces is inside the seperatrix
                for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)[:num_lines]):
                #for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)):
                    x, y = draw_contour_line(R, Z, self.psi_norm_raw, v, j)
                    if (np.amax(x) < np.amax(self.main_sep_pts[:, 0]) and \
                        np.amin(x) > np.amin(self.main_sep_pts[:, 0]) and \
                        np.amax(y) < np.amax(self.main_sep_pts[:, 1]) and \
                        np.amin(y) > np.amin(self.main_sep_pts[:, 1])):
                        #then it's an internal flux surface
                        self.core_ntrl_lines.append(LineString(np.column_stack((x[:-1], y[:-1]))))
                        break

    def core_main(self, inp, R_psi, Z_psi, psi, B_pol_raw):
        #define rho points
        try:
            rho1d = np.concatenate((np.linspace(0, inp.edge_rho, inp.rhopts_core, endpoint=False), 
                                    np.linspace(inp.edge_rho, 1, inp.rhopts_edge, endpoint=True)))
        except:
            try:
                rho1d = np.linspace(0, 1, inp.rhopts)
            except:
                print 'rho parameters not defined. Using 100 evenly spaced rho values'
                rho1d = np.linspace(0, 1, 100)
                
        #define theta points
        def atan3(y, x):
            result = atan2(y, x)
            if result<0:
                result = result + 2*pi
            return result
        
        #these theta markers correspond to the obmp, top, ibmp, bot, and obmp+2pi, respectively
        theta_marker = np.zeros(5)
        theta_marker[0] = atan2((self.obmp_pt[1]-self.geo_axis[1]), (self.obmp_pt[0]-self.geo_axis[0]))
        theta_marker[1] = atan3((self.top_pt[1]-self.geo_axis[1]), (self.top_pt[0]-self.geo_axis[0]))
        theta_marker[2] = atan3((self.ibmp_pt[1]-self.geo_axis[1]), (self.ibmp_pt[0]-self.geo_axis[0]))
        theta_marker[3] = atan3((self.xpt[1]-self.geo_axis[1]), (self.xpt[0]-self.geo_axis[0]))
        theta_marker[4] = theta_marker[0] + 2*pi
        
        try:
            min_delta_theta = 2*pi/inp.thetapts_approx
        except:
            print 'thetapts_approx not defined. Setting to 30'
            min_delta_theta = 2*pi/30
        
        theta1d = np.zeros(0)
        for i in range(4):
            quad_pts = ceil((theta_marker[i+1] - theta_marker[i])/min_delta_theta)
            quad_theta = np.linspace(theta_marker[i], theta_marker[i+1], quad_pts)
            theta1d = np.concatenate((theta1d, quad_theta))
        
        self.thetapts = len(theta1d)
        self.theta, self.rho = np.meshgrid(theta1d, rho1d)
        self.r = self.rho * self.a #TODO: revisit what 'r' means in 2D.      
        #fill in parameters that are only functions of rho
        try:
            self.ni = UnivariateSpline(inp.ni_data[:, 0], inp.ni_data[:, 1], k=3, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=3, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.Ti_kev = UnivariateSpline(inp.Ti_data[:, 0], inp.Ti_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.Te_kev = UnivariateSpline(inp.Te_data[:, 0], inp.Te_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass        
        
        try:
            E_r_fit = UnivariateSpline(inp.er_data[:, 0], inp.er_data[:, 1], k=5, s=2.0)
            self.E_r = E_r_fit(self.rho)
            self.E_pot = np.zeros(self.rho.shape)
            for i, rhoval in enumerate(rho1d):
                self.E_pot[i] = E_r_fit.integral(rhoval, 1.0)
        except AttributeError:
            pass
        
        try:
            self.fracz = UnivariateSpline(inp.fracz_data[:, 0], inp.fracz_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.fracz = np.zeros(self.rho.shape) + 0.025            
        self.nC = self.ne * self.fracz        
        self.z_0 = self.nC*6.0**2 / self.ni
        self.z_eff = (self.ni*1.0**2 + self.nC*6.0**2) / self.ne
        try:
            self.fz1 = UnivariateSpline(inp.fz1_data[:, 0], inp.fz1_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.q = UnivariateSpline(inp.q_data[:, 0], inp.q_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.vpolC = UnivariateSpline(inp.vpolC_data[:, 0], inp.vpolC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.vpolD = UnivariateSpline(inp.vpolD_data[:, 0], inp.vpolD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.vtorC = UnivariateSpline(inp.vtorC_data[:, 0], inp.vtorC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.vtorD = UnivariateSpline(inp.vtorD_data[:, 0], inp.vtorD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass
        
        try:
            self.zbar2 = UnivariateSpline(inp.zbar2_data[:, 0], inp.zbar2_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            pass

        #get parameters that depend on both rho and theta
        self.R = np.zeros(self.rho.shape)
        self.Z = np.zeros(self.rho.shape)
        self.psi = np.zeros(self.rho.shape)
        self.psi_norm = np.zeros(self.rho.shape)
        
        #move along line between m_axis and obmp and define psi values corresponding to rho values
        rho_line = LineString([Point(self.m_axis), Point(self.obmp_pt)])
        init_coords = np.zeros((0, 2))
        for i, rhoval in enumerate(rho1d):
            pt_coords = np.asarray(rho_line.interpolate(rhoval, normalized=True).coords)[0]
            init_coords = np.vstack((init_coords, pt_coords))

        psi_vals = griddata(np.column_stack((R_psi.flatten(), Z_psi.flatten())), 
                            psi.flatten(), 
                            init_coords, 
                            method='linear')
        
        psi_norm_vals = griddata(np.column_stack((R_psi.flatten(), Z_psi.flatten())), 
                            self.psi_norm_raw.flatten(), 
                            init_coords, 
                            method='linear')

        for i, (psi_val, psi_norm_val) in enumerate(zip(psi_vals, psi_norm_vals)): 
            self.psi[i] = psi_val
            self.psi_norm[i] = psi_norm_val
                 
            fs_line, fs_axis = draw_core_line(R_psi, Z_psi, self.psi_norm_raw, psi_norm_val, self.main_sep_pts)

            for j, thetaval in enumerate(theta1d):
                if psi_norm_val<1.0:
                    thetaline = LineString([Point(fs_axis), 
                                            Point([3.0*cos(thetaval)+fs_axis[0], 
                                                   3.0*sin(thetaval)+fs_axis[1]])])
                    int_pt = fs_line.intersection(thetaline)
                else:
                    if thetaval == theta_marker[3]:
                        int_pt = Point(self.xpt)

                    else:
                        thetaline = LineString([Point(self.geo_axis), 
                                                Point([3.0*cos(thetaval)+self.geo_axis[0], 
                                                       3.0*sin(thetaval)+self.geo_axis[1]])])
                        int_pt = self.main_sep_line_closed.intersection(thetaline)

                self.R[i, j] = int_pt.x
                self.Z[i, j] = int_pt.y               

        self.R[0] = self.m_axis[0]        
        self.Z[0] = self.m_axis[1]     
        self.psi_norm[0] = 0
        self.B_p = griddata(np.column_stack((R_psi.flatten(), Z_psi.flatten())), 
                                    B_pol_raw.flatten(), 
                                    (self.R, self.Z), 
                                    method='linear') 

        self.B_t = inp.BT0 * self.m_axis[0] / self.R
        self.B_tot = np.sqrt(self.B_p**2 + self.B_t**2)
        self.f_phi = self.B_t/self.B_tot
        #plt.axis('equal')
        #for i, (Rvals, Zvals) in enumerate(zip(self.R, self.Z)):
        #    plt.plot(Rvals, Zvals, lw=0.5)

        #calculate gradients and gradient scale lengths
        #get quantities on a fairly fine R, Z grid for the purpose of taking gradients, etc.
        R_temp, Z_temp = np.meshgrid(np.linspace(0.98*self.ibmp_pt[0], 1.02*self.obmp_pt[0], 500), 
                                    np.linspace(1.02*self.top_pt[1], 1.02*self.bot_pt[1], 500))

        ni_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())), 
                           self.ni.flatten(), 
                           (R_temp, Z_temp), 
                           method='cubic', 
                           fill_value = self.ni[-1, 0])
        ne_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())), 
                           self.ne.flatten(), 
                           (R_temp, Z_temp), 
                           method='cubic', 
                           fill_value = self.ne[-1, 0])
        Ti_kev_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())), 
                           self.Ti_kev.flatten(), 
                           (R_temp, Z_temp), 
                           method='cubic', 
                           fill_value = self.Ti_kev[-1, 0])
        Te_kev_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())), 
                           self.Te_kev.flatten(), 
                           (R_temp, Z_temp), 
                           method='cubic', 
                           fill_value = self.Te_kev[-1, 0])
        
        dnidr_temp = -1.0*(np.abs(np.gradient(ni_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(ni_grid, R_temp[0, :], axis=0)))
        dnedr_temp = -1.0*(np.abs(np.gradient(ne_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(ne_grid, R_temp[0, :], axis=0)))
        dTidr_temp = -1.0*(np.abs(np.gradient(Ti_kev_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(Ti_kev_grid, R_temp[0, :], axis=0)))
        dTedr_temp = -1.0*(np.abs(np.gradient(Te_kev_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(Te_kev_grid, R_temp[0, :], axis=0)))
        
        self.dni_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())), 
                              dnidr_temp.flatten(), 
                              (self.R, self.Z), 
                              method='cubic')
        self.dne_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())), 
                              dnedr_temp.flatten(), 
                              (self.R, self.Z), 
                              method='cubic')
        self.dTi_kev_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())), 
                              dTidr_temp.flatten(), 
                              (self.R, self.Z), 
                              method='cubic')
        self.dTe_kev_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())), 
                              dTedr_temp.flatten(), 
                              (self.R, self.Z), 
                              method='cubic')

        self.Ti_J = self.Ti_kev * 1.6021E-16
        self.Ti_ev = self.Ti_kev * 1E3
        self.Te_J = self.Te_kev * 1.6021E-16
        self.Te_ev = self.Te_kev * 1E3
        self.dTi_J_dr = self.dTi_kev_dr * 1.6021E-16
        self.dTi_ev_dr = self.dTi_kev_dr * 1E3
        self.dTe_J_dr = self.dTe_kev_dr * 1.6021E-16
        self.dTe_ev_dr = self.dTe_kev_dr * 1E3
        self.L_ni = -self.dni_dr / self.ni
        self.L_ne = -self.dne_dr / self.ne
        self.L_Ti_J = -self.dTi_J_dr / self.Ti_J
        self.L_Te_J = -self.dTe_J_dr / self.Te_J
        #plt.axis('equal')
        #plt.contourf(self.R, self.Z, self.dnidr, 500)
        #plt.colorbar()
        #sys.exit()
            
        #create neutrals-related variables. These will remain zero unless set by exp_neutpy_prep or read_ntrl_data modules
        self.n_n_slow = np.zeros(self.rho.shape)
        self.n_n_thermal = np.zeros(self.rho.shape)
        self.n_n_total = np.zeros(self.rho.shape)
        
        self.izn_rate_slow = np.zeros(self.rho.shape)
        self.izn_rate_thermal = np.zeros(self.rho.shape)
        self.izn_rate_total = np.zeros(self.rho.shape)

        self.T_n_slow = np.full(self.rho.shape, 0.002)
        self.T_n_thermal = self.Ti_kev

        # create Lz-related variables. These will remain zero unless set by the ImpRad module
        self.Lz_slow = np.zeros(self.rho.shape)
        self.dLzdT_slow = np.zeros(self.rho.shape)
        self.Lz_thermal = np.zeros(self.rho.shape)
        self.dLzdT_thermal = np.zeros(self.rho.shape)

    def core_nT_ntrl(self, inp, R, Z, psi):
        #CREATE ARRAYS OF POINTS, DENSITIES AND TEMPERATURES FOR THE NEUTRALS CALCULATION
        
        #Master arrays that will contain all the points we'll use to get n, T
        #throughout the plasma chamber via 2-D interpolation
        self.ni_pts = np.zeros((0, 3), dtype='float')
        self.ne_pts = np.zeros((0, 3), dtype='float')
        self.Ti_kev_pts = np.zeros((0, 3), dtype='float')
        self.Te_kev_pts = np.zeros((0, 3), dtype='float')
        
        ##########################################
        #Calculate n, T throughout the core plasma using radial profile input files, uniform on flux surface
        ni = UnivariateSpline(inp.ni_data[:, 0], inp.ni_data[:, 1], k=5, s=2.0)
        ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=5, s=2.0)
        Ti_kev = UnivariateSpline(inp.Ti_data[:, 0], inp.Ti_data[:, 1], k=5, s=2.0)
        Te_kev = UnivariateSpline(inp.Te_data[:, 0], inp.Te_data[:, 1], k=5, s=2.0)
        
        #get approximate rho values associated with the psi values we're using
        #draw line between magnetic axis and the seperatrix at the outboard midplane
        rho_line = LineString([Point(self.m_axis), Point(self.obmp_pt)])
        rho_pts = np.concatenate((np.linspace(0, 0.95, 20, endpoint=False), 
                                 np.linspace(0.95, 1, 50, endpoint=False)), axis=0)
        
        thetapts = np.linspace(0, 1, 100, endpoint=False)
        for i, rho in enumerate(rho_pts): 
            #get n, T information at the point by interpolating the rho-based input file data
            ni_val = ni(rho)
            ne_val = ne(rho)
            Ti_kev_val = Ti_kev(rho)
            Te_kev_val = Te_kev(rho)
            #get R, Z coordinates of each point along the rho_line
            pt_coords = np.asarray(rho_line.interpolate(rho, normalized=True).coords)[0]

            #get psi value at that point
            psi_val = griddata(np.column_stack((R.flatten(), Z.flatten())), 
                                 self.psi_norm_raw.flatten(), 
                                 pt_coords, 
                                 method='linear')
            #map this n, T data to every point on the corresponding flux surface
            num_lines = int(len(cntr.Cntr(R, Z, self.psi_norm_raw).trace(psi_val))/2)

            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x, y = draw_contour_line(R, Z, self.psi_norm_raw, psi_val, 0)
                surf = LineString(np.column_stack((x, y)))
            else:
                #we need to find which of the surfaces is inside the seperatrix
                for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(psi_val)[:num_lines]):
                    #for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)):
                    x, y = draw_contour_line(R, Z, self.psi_norm_raw, psi_val, j)
                    if (np.amax(x) < np.amax(self.main_sep_pts[:, 0]) and \
                        np.amin(x) > np.amin(self.main_sep_pts[:, 0]) and \
                        np.amax(y) < np.amax(self.main_sep_pts[:, 1]) and \
                        np.amin(y) > np.amin(self.main_sep_pts[:, 1])):
                        #then it's an internal flux surface
                        surf = LineString(np.column_stack((x, y)))
                        break
            
            for j, theta_norm in enumerate(thetapts):
                pt = np.asarray(surf.interpolate(theta_norm, normalized=True).coords).T
                self.ni_pts = np.vstack((self.ni_pts, np.append(pt, ni_val)))
                self.ne_pts = np.vstack((self.ne_pts, np.append(pt, ne_val)))
                self.Ti_kev_pts = np.vstack((self.Ti_kev_pts, np.append(pt, Ti_kev_val)))
                self.Te_kev_pts = np.vstack((self.Te_kev_pts, np.append(pt, Te_kev_val)))

        #Do seperatrix separately so we don't accidentally assign the input n, T data to the divertor legs
        self.ni_sep_val = ni(1.0)
        self.ne_sep_val = ne(1.0)
        self.Ti_kev_sep_val = Ti_kev(1.0)
        self.Te_kev_sep_val = Te_kev(1.0)
        self.Ti_J_sep_val = self.Ti_kev_sep_val * 1.0E3 * 1.6021E-19
        self.Te_J_sep_val = self.Te_kev_sep_val * 1.0E3 * 1.6021E-19
        for j, theta_norm in enumerate(thetapts): 
            pt = np.asarray(self.main_sep_line.interpolate(theta_norm, normalized=False).coords, dtype='float').T
            self.ni_pts = np.vstack((self.ni_pts, np.append(pt, self.ni_sep_val)))
            self.ne_pts = np.vstack((self.ne_pts, np.append(pt, self.ne_sep_val)))
            self.Ti_kev_pts = np.vstack((self.Ti_kev_pts, np.append(pt, self.Ti_kev_sep_val)))
            self.Te_kev_pts = np.vstack((self.Te_kev_pts, np.append(pt, self.Te_kev_sep_val)))

    def xsec(self, inp):
        #Fusion Reactivity calculation  
        def calc_sigv_fus(mode='dd'):   
            def sigv(Ti, mode): #function takes T in kev
                if mode=='dt':
                    B_G = 34.3827    
                    m_rc2 = 1124656
                    
                    C1 = 1.17302E-9
                    C2 = 1.51361E-2
                    C3 = 7.51886E-2
                    C4 = 4.60643E-3
                    C5 = 1.35000E-2
                    C6 = -1.06750E-4
                    C7 = 1.36600E-5
                
                    theta = Ti/(1.0-(Ti*(C2+Ti*(C4+Ti*C6)))/(1.0+Ti*(C3+Ti*(C5+Ti*C7))))
                    xi = (B_G**2.0/(4.0*theta))**(1.0/3.0)
                    sigv = C1 * theta * np.sqrt(xi/(m_rc2 * Ti**3.0)) * np.exp(-3.0*xi)
                    sigv = sigv/1.0E6 #convert from cm^3/s to m^3/s
                    
                elif mode=='dd':
                    
                    B_G = 31.3970 
                    m_rc2 = 937814
                    
                    #first for the D(d, p)T reaction
                    C1_1 = 5.65718E-12
                    C2_1 = 3.41267E-3
                    C3_1 = 1.99167E-3
                    C4_1 = 0.0
                    C5_1 = 1.05060E-5
                    C6_1 = 0.0
                    C7_1 = 0.0
                
                    theta_1 = Ti/(1.0-(Ti*(C2_1+Ti*(C4_1+Ti*C6_1)))/(1.0+Ti*(C3_1+Ti*(C5_1+Ti*C7_1))))
                    xi_1 = (B_G**2.0/(4.0*theta_1))**(1.0/3.0)
                    sigv_1 = C1_1 * theta_1 * np.sqrt(xi_1/(m_rc2 * Ti**3.0)) * np.exp(-3.0*xi_1)
                    
                    #then for the D(d, n)He3 reaction
                    
                    C1_2 = 5.43360E-12
                    C2_2 = 5.85778E-3
                    C3_2 = 7.68222E-3
                    C4_2 = 0.0
                    C5_2 = -2.96400E-6
                    C6_2 = 0.0
                    C7_2 = 0.0
                
                    theta_2 = Ti/(1.0-(Ti*(C2_2+Ti*(C4_2+Ti*C6_2)))/(1.0+Ti*(C3_2+Ti*(C5_2+Ti*C7_2))))
                    xi_2 = (B_G**2.0/(4.0*theta_2))**(1.0/3.0)
                    sigv_2 = C1_2 * theta_2 * np.sqrt(xi_2/(m_rc2 * Ti**3.0)) * np.exp(-3.0*xi_2)                
                    
                    sigv = (0.5*sigv_1 + 0.5*sigv_2) / 1.0E6  # convert from cm^3/s to m^3/s
                return sigv
            
            # create logspace over the relevant temperature range
            # (bosch hale technically only valid over 0.2 - 100 kev)
            Ti_range = np.logspace(-1, 2, 1000)  # values in kev
            sigv_fus_range = sigv(Ti_range, mode='dd')  # in m^3/s
            sigv_fus_interp = UnivariateSpline(Ti_range*1.0E3*1.6021E-19, sigv_fus_range, s=0)  # converted to Joules
            self.sv_fus = sigv_fus_interp(self.Ti_J)
            self.dsv_fus_dT = sigv_fus_interp.derivative()(self.Ti_J)
            self.dsv_fus_dT_eq9 = sigv_fus_interp.derivative()(5.0E2*1.6021E-19)
        
        def calc_sigv_ion():
            # TODO: configure so it can use any of the cross section libraries
            # currently using the Stacey-Thomas cross sections
            T_exps_fit = np.array([-1, 0, 1, 2, 3, 4, 5])
            sigv_exps_fit = np.array([-2.8523E+01, -1.7745E+01, -1.3620E+01, 
                                            -1.3097E+01, -1.3301E+01, -1.3301E+01, -1.3301E+01])
            interp1 = UnivariateSpline(T_exps_fit, sigv_exps_fit, s=0)
            
            T_exps_range = np.linspace(-1, 5, 1000)
            sigv_vals_range = 10.0**interp1(T_exps_range) # in m^3/s
            
            T_vals_range = np.logspace(-1, 5, 1000)*1.6021E-19 # in joules
            interp2 = UnivariateSpline(T_vals_range, sigv_vals_range, s=0)
            
            self.sv_ion = interp2(self.Ti_J)
            self.dsv_ion_dT = interp2.derivative()(self.Ti_J)

        def calc_svel():
        
            tint = np.array([-1, 0, 1, 2, 3])
            tnnt = np.array([0, 1, 2])

            elast = np.array([[-1.3569E+01, -1.3337E+01, -1.3036E+01, -1.3569E+01, -1.3337E+01], 
                              [-1.3036E+01, -1.3337E+01, -1.3167E+01, -1.3046E+01, -1.3036E+01], 
                              [-1.3046E+01, -1.2796E+01, -1.3036E+01, -1.3046E+01, -1.2796E+01]])
                
            interp1 = interp2d(tint, tnnt, elast)
            
            Ti_exps = np.linspace(-1, 3, 100)
            Tn_exps = np.linspace( 0, 2, 100)
            svel_vals = 10.0**(interp1(Ti_exps, Tn_exps)) # in m^3/s
            
            Ti_vals = np.logspace(-1, 3, 100)*1.6021E-19 # in joules
            Tn_vals = np.logspace( 0, 2, 100)*1.6021E-19 # in joules

            dsvel_dTi_vals = np.gradient(svel_vals, Ti_vals, axis=0)

            Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)
            
            Ti_mod = np.where(self.Ti_ev>1E3, 1E3 * 1.6021E-19, self.Ti_ev*1.6021E-19)
            Tn_mod = np.zeros(Ti_mod.shape) + 2.0*1.6021E-19

            self.sv_el = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       svel_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)
            self.dsv_el_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       dsvel_dTi_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)
         
        def calc_svcxi_st():
            
            tint = np.array([-1, 0, 1, 2, 3])
            tnnt = np.array([0, 1, 2])
            
            cx = np.array([[-1.4097E+01, -1.3921E+01, -1.3553E+01, -1.4097E+01, -1.3921E+01], 
                           [-1.3553E+01, -1.3921E+01, -1.3824E+01, -1.3538E+01, -1.3553E+01], 
                           [-1.3538E+01, -1.3432E+01, -1.3553E+01, -1.3538E+01, -1.3432E+01]])
                
            interp1 = interp2d(tint, tnnt, cx)
            
            Ti_exps = np.linspace(-1, 3, 100)
            Tn_exps = np.linspace( 0, 2, 100)
            svcx_vals = 10.0**(interp1(Ti_exps, Tn_exps)) # in m^3/s
            
            Ti_vals = np.logspace(-1, 3, 100)*1.6021E-19 # in joules
            Tn_vals = np.logspace( 0, 2, 100)*1.6021E-19 # in joules

            dsvcx_dTi_vals = np.gradient(svcx_vals, Ti_vals, axis=0)

            Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)
            
            Ti_mod = np.where(self.Ti_ev>1E3, 1E3 * 1.6021E-19, self.Ti_ev*1.6021E-19)
            Tn_mod = np.zeros(Ti_mod.shape) + 2.0*1.6021E-19

            self.sv_cx = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       svcx_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)
            
            self.dsv_cx_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       dsvcx_dTi_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)

        def calc_svrec_st(): 
            # TODO: check this calculation. -MH
            znint = np.array([16, 18, 20, 21, 22])
            Tint = np.array([-1, 0, 1, 2, 3])
        
            rec = np.array([[-1.7523E+01, -1.6745E+01, -1.5155E+01, -1.4222E+01, -1.3301E+01], 
                            [-1.8409E+01, -1.8398E+01, -1.8398E+01, -1.7886E+01, -1.7000E+01], 
                            [-1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01], 
                            [-2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01], 
                            [-2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01]])
            
            interp1 = interp2d(znint, Tint, rec)
            
            zni_exps = np.linspace(16, 22, 100)
            Ti_exps = np.linspace(-1, 3, 100)
            svrec_vals = 10.0**(interp1(zni_exps, Ti_exps)) #in m^3/s
            
            zni_vals = np.logspace(16, 22, 100)
            Ti_vals = np.logspace(-1, 3, 100)*1.6021E-19 #in joules
            
            dsvrec_dTi_vals = np.gradient(svrec_vals, Ti_vals, axis=0)
                
            zni_vals2d, Ti_vals2d = np.meshgrid(zni_vals, Ti_vals)
            
            zni_mod = np.where(self.ni>1E22, 1E22, self.ni)
            zni_mod = np.where(self.ni<1E16, 1E16, zni_mod)
            Ti_mod = np.where(self.Ti_ev>1E3, 1E3 * 1.6021E-19, self.Ti_ev*1.6021E-19)
            Ti_mod = np.where(self.Ti_ev<1E-1, 1E-1 * 1.6021E-19, Ti_mod)
                
            self.sv_rec = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())), 
                                       svrec_vals.flatten(), 
                                       (zni_mod, Ti_mod), 
                                       method='linear', rescale=False)
            
            self.dsv_rec_dT = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())), 
                                       dsvrec_dTi_vals.flatten(), 
                                       (zni_mod, Ti_mod), 
                                       method='linear', rescale=False)

        calc_sigv_fus()
        calc_sigv_ion()
        calc_svel()
        calc_svcxi_st()
        #calc_svrec_st()
