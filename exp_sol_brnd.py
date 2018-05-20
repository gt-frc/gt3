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
import sys

def draw_line(R,Z,array,val,pathnum):
    res = cntr.Cntr(R,Z,array).trace(val)[pathnum]
    x = res[:,0]
    y = res[:,1]
    return x,y

def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= 1.0:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p),normalized=True)
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance,normalized=True)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]

class exp_sol_brnd():
    def __init__(self,inp,core):
        
        R         = inp.psirz_exp[:,0].reshape(-1,65)
        Z         = inp.psirz_exp[:,1].reshape(-1,65)
        
        self.sol_lines(inp,R,Z,core)
        self.sol_nT(inp,R,Z,core)
        pass
    
    def sol_lines(self,inp,R,Z,core):
        #find value of psi at outside of what we're going to call the SOL
        self.sol_lines = []
        self.sol_lines_cut = []

        sol_width_obmp = 0.02
        psi_pts = np.linspace(1,inp.sollines_psi_max,inp.num_sollines+1,endpoint=True)[1:]
        for i,v in enumerate(psi_pts):
            num_lines = int(len(cntr.Cntr(R,Z,core.psi_norm).trace(v))/2)
            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x,y = draw_line(R,Z,core.psi_norm,v,0)
                self.sol_lines.append(LineString(np.column_stack((x,y))))
            else:
                #TODO: 
                pass
        for line in self.sol_lines:
            #find intersection points with the wall
            int_pts = line.intersection(inp.wall_line)
            #cut line at intersection points
            cut_line = cut(line,line.project(int_pts[0],normalized=True))[1]
            cut_line = cut(cut_line,cut_line.project(int_pts[1],normalized=True))[0]
            self.sol_lines_cut.append(cut_line)
            
        #add wall intersection points from divertor legs and sol lines to wall_line.
        #This is necessary to prevent thousands of tiny triangles from forming if the 
        #end of the flux line isn't exactly on top of the wall line.

        #add inboard seperatrix strike point
        union = inp.wall_line.union(core.ib_div_line)
        result = [geom for geom in polygonize(union)][0]
        inp.wall_line = LineString(result.exterior.coords)

        #add outboard seperatrix strike point
        union = inp.wall_line.union(core.ob_div_line)
        result = [geom for geom in polygonize(union)][0]
        inp.wall_line = LineString(result.exterior.coords)  
        
        #add sol line intersection points on inboard side
        #for some reason, union freaks out when I try to do inboard and outboard
        #at the same time.
        for num,line in enumerate(self.sol_lines):
            union = inp.wall_line.union(cut(line,0.5)[0])    
            result = [geom for geom in polygonize(union)][0]
            inp.wall_line = LineString(result.exterior.coords)

        #add sol line intersection points on outboard side            
        for num,line in enumerate(self.sol_lines):
            union = inp.wall_line.union(cut(line,0.5)[1])    
            result = [geom for geom in polygonize(union)][0]
            inp.wall_line = LineString(result.exterior.coords)
            
    def sol_nT(self,inp,R,Z,core):
        #draw core lines in the psi_norm range in which we want to evaluate d__/dr for calculating bohm
        #diffusion into the SOL. These lines won't be used for anything else. 
        
        ni     = UnivariateSpline(inp.ni_data[:,0],inp.ni_data[:,1],k=5,s=2.0)
        Ti_kev = UnivariateSpline(inp.Ti_data[:,0],inp.Ti_data[:,1],k=5,s=2.0)
        
        gradient_pts = np.zeros((0,2))
        gradient_pts_ni = np.zeros(0)
        gradient_pts_Ti = np.zeros(0)
        
        psi_pts = np.linspace(0.95,1.0,20,endpoint=False)
        for i,psi_val in enumerate(psi_pts):
            num_lines = int(len(cntr.Cntr(R,Z,core.psi_norm).trace(psi_val))/2)
            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x,y = draw_line(R,Z,core.psi_norm,psi_val,0)
            else:
                #we need to find which of the surfaces is inside the seperatrix
                for j,line in enumerate(cntr.Cntr(R,Z,core.psi_norm).trace(psi_val)[:num_lines]):
                #for j,line in enumerate(cntr.Cntr(R,Z,self.psi_norm).trace(v)):
                    x,y = draw_line(R,Z,core.psi_norm,psi_val,j)
                    if (np.amax(x) < np.amax(core.main_sep_pts[:,0]) and \
                        np.amin(x) > np.amin(core.main_sep_pts[:,0]) and \
                        np.amax(y) < np.amax(core.main_sep_pts[:,1]) and \
                        np.amin(y) > np.amin(core.main_sep_pts[:,1])):
                        #then it's an internal flux surface
                        break
            #TODO: double check if radial profile data is given in rho or psi_norm. Might need to adjust this.
            gradient_pts_ni = np.concatenate((gradient_pts_ni, np.zeros(len(x)) + ni(psi_val)))
            gradient_pts_Ti = np.concatenate((gradient_pts_Ti, np.zeros(len(x)) + Ti_kev(psi_val)))
            gradient_pts = np.vstack((gradient_pts,np.column_stack((x,y))))

        gradient_pts = np.vstack((gradient_pts,core.main_sep_pts))
        gradient_pts_ni = np.concatenate((gradient_pts_ni, np.zeros(len(core.main_sep_pts)) + core.ni_sep_val))
        gradient_pts_Ti = np.concatenate((gradient_pts_Ti, np.zeros(len(core.main_sep_pts)) + core.Ti_kev_sep_val))

        #get quantities on a fairly fine R,Z grid for the purpose of taking gradients, etc.
        R_temp,Z_temp = np.meshgrid(np.linspace(0.98*core.ibmp_pt[0],1.02*core.obmp_pt[0],500),
                                    np.linspace(1.02*core.top_pt[1],1.02*core.bot_pt[1],500))

        ni_grid = griddata(gradient_pts,
                           gradient_pts_ni,
                           (R_temp,Z_temp),
                           method='linear',
                           fill_value = core.ni_sep_val)
        
        Ti_grid = griddata(gradient_pts,
                           gradient_pts_Ti*1.0E3*1.6021E-19,
                           (R_temp,Z_temp),
                           method='linear',
                           fill_value = core.Ti_J_sep_val)
        
        dnidr = -1.0*(np.abs(np.gradient(ni_grid,Z_temp[:,0],axis=1)) + np.abs(np.gradient(ni_grid,R_temp[0,:],axis=0)))
        dTidr = -1.0*(np.abs(np.gradient(Ti_grid,Z_temp[:,0],axis=1)) + np.abs(np.gradient(Ti_grid,R_temp[0,:],axis=0)))

        #Calculate n,T in SOL using Bohm diffusion, core data from radial profile input files, and input
        #divertor target densities and temperatures (replace with 2-pt divertor model later)

        #draw core line just inside the seperatrix (seperatrix would be too noisy)
        psi_val = 0.98
        num_lines = int(len(cntr.Cntr(R,Z,core.psi_norm).trace(psi_val))/2)
        if num_lines==1:
            #then we're definitely dealing with a surface inside the seperatrix
            x,y = draw_line(R,Z,core.psi_norm,psi_val,0)
        else:
            #we need to find which of the surfaces is inside the seperatrix
            for j,line in enumerate(cntr.Cntr(R,Z,core.psi_norm).trace(psi_val)[:num_lines]):
            #for j,line in enumerate(cntr.Cntr(R,Z,psi_norm).trace(v)):
                x,y = draw_line(R,Z,core.psi_norm,psi_val,j)
                if (np.amax(x) < np.amax(core.main_sep_pts[:,0]) and \
                    np.amin(x) > np.amin(core.main_sep_pts[:,0]) and \
                    np.amax(y) < np.amax(core.main_sep_pts[:,1]) and \
                    np.amin(y) > np.amin(core.main_sep_pts[:,1])):
                    #then it's an internal flux surface
                    break
                
        #Get densities, temperatures, and other quantities along the flux surface we just drew
        #note densities and temperatures on seperatrix are assumed to be constant for all theta and are
        #obtained above, i.e. self.ni_sep_val, etc.
        dnidr_sep_raw = griddata(np.column_stack((R_temp.flatten(),Z_temp.flatten())),
                             dnidr.flatten(),
                             np.column_stack((x,y)),
                             method='cubic'
                             )
    
        dTidr_sep_raw = griddata(np.column_stack((R_temp.flatten(),Z_temp.flatten())),
                             dTidr.flatten(),
                             np.column_stack((x,y)),
                             method='cubic'
                             )

        BT_sep_raw     = griddata(np.column_stack((R_temp.flatten(),Z_temp.flatten())),
                             (core.m_axis[0]*inp.BT0/R_temp).flatten(),
                             np.column_stack((x,y)),
                             method='cubic'
                             )
        
        dnedr_sep_raw = dnidr_sep_raw
        
        dTedr_sep_raw = dTidr_sep_raw
        
        #norm factor used to divide by the order of magnitude to facilitate easier smoothing
        ni_norm_factor = 1.0#*10**(int(np.log10(np.average(dnidr_sep)))-1)
        ne_norm_factor = 1.0#*10**(int(np.log10(np.average(dnedr_sep)))-1)
        Ti_norm_factor = 1.0#*10**(int(np.log10(np.average(dTidr_sep)))-1)
        Te_norm_factor = 1.0#*10**(int(np.log10(np.average(dTedr_sep)))-1)

        ni_sep = np.zeros(inp.xi_sep_pts) + core.ni_sep_val
        ne_sep = np.zeros(inp.xi_sep_pts) + core.ne_sep_val
        Ti_sep = np.zeros(inp.xi_sep_pts) + core.Ti_J_sep_val
        Te_sep = np.zeros(inp.xi_sep_pts) + core.Te_J_sep_val
        
        dnidr_sep = UnivariateSpline(np.linspace(0,1,len(dnidr_sep_raw)),
                                            dnidr_sep_raw/ni_norm_factor,
                                            k=5,
                                            s=0.0)(np.linspace(inp.ib_trim_off,1.0-inp.ob_trim_off,inp.xi_sep_pts))*ni_norm_factor

        dnedr_sep = UnivariateSpline(np.linspace(0,1,len(dnedr_sep_raw)),
                                            dnedr_sep_raw/ne_norm_factor,
                                            k=5,
                                            s=0.0)(np.linspace(inp.ib_trim_off,1.0-inp.ob_trim_off,inp.xi_sep_pts))*ne_norm_factor

        dTidr_sep = UnivariateSpline(np.linspace(0,1,len(dTidr_sep_raw)),
                                            dTidr_sep_raw/Ti_norm_factor,
                                            k=5,
                                            s=0.0)(np.linspace(inp.ib_trim_off,1.0-inp.ob_trim_off,inp.xi_sep_pts))*Ti_norm_factor
                                            
        dTedr_sep = UnivariateSpline(np.linspace(0,1,len(dTedr_sep_raw)),
                                            dTedr_sep_raw/Te_norm_factor,
                                            k=5,
                                            s=0.0)(np.linspace(inp.ib_trim_off,1.0-inp.ob_trim_off,inp.xi_sep_pts))*Te_norm_factor
                                            
        BT_sep    = UnivariateSpline(np.linspace(0,1,len(dnidr_sep_raw)),
                                            BT_sep_raw,
                                            k=5,
                                            s=0.0)(np.linspace(inp.ib_trim_off,1.0-inp.ob_trim_off,inp.xi_sep_pts))
        
        ni_ib_wall = core.ni_sep_val * 1.0
        ni_ob_wall = core.ni_sep_val * 1.0
        ni_ib = np.linspace(ni_ib_wall,core.ni_sep_val,inp.xi_ib_pts,endpoint=False)
        ni_ob = np.linspace(core.ni_sep_val,ni_ob_wall,inp.xi_ob_pts,endpoint=True)     
        
        ne_ib_wall = core.ne_sep_val * 1.0
        ne_ob_wall = core.ne_sep_val * 1.0
        ne_ib = np.linspace(ne_ib_wall,core.ne_sep_val,inp.xi_ib_pts,endpoint=False)
        ne_ob = np.linspace(core.ne_sep_val,ne_ob_wall,inp.xi_ob_pts,endpoint=True) 
        
        Ti_ib_wall = core.Ti_J_sep_val * 1.0
        Ti_ob_wall = core.Ti_J_sep_val * 1.0
        Ti_ib = np.linspace(Ti_ib_wall,core.Ti_J_sep_val,inp.xi_ib_pts,endpoint=False)
        Ti_ob = np.linspace(core.Ti_J_sep_val,Ti_ob_wall,inp.xi_ob_pts,endpoint=True)    

        Te_ib_wall = core.Te_J_sep_val * 1.0
        Te_ob_wall = core.Te_J_sep_val * 1.0
        Te_ib = np.linspace(Te_ib_wall,core.Te_J_sep_val,inp.xi_ib_pts,endpoint=False)
        Te_ob = np.linspace(core.Te_J_sep_val,Te_ob_wall,inp.xi_ob_pts,endpoint=True)    

        dnidr_ib_wall = dnidr_sep[0]
        dnidr_ob_wall = dnidr_sep[-1]
        dnidr_ib = np.linspace(dnidr_ib_wall,dnidr_sep[0],inp.xi_ib_pts,endpoint=False)
        dnidr_ob = np.linspace(dnidr_sep[-1],dnidr_ob_wall,inp.xi_ob_pts,endpoint=True)
        
        dnedr_ib_wall = dnedr_sep[0]
        dnedr_ob_wall = dnedr_sep[-1]
        dnedr_ib = np.linspace(dnedr_ib_wall,dnedr_sep[0],inp.xi_ib_pts,endpoint=False)
        dnedr_ob = np.linspace(dnedr_sep[-1],dnedr_ob_wall,inp.xi_ob_pts,endpoint=True)
        
        dTidr_ib_wall = dTidr_sep[0]
        dTidr_ob_wall = dTidr_sep[-1]
        dTidr_ib = np.linspace(dTidr_ib_wall,dTidr_sep[0],inp.xi_ib_pts,endpoint=False)
        dTidr_ob = np.linspace(dTidr_sep[-1],dTidr_ob_wall,inp.xi_ob_pts,endpoint=True)
        
        dTedr_ib_wall = dTedr_sep[0]
        dTedr_ob_wall = dTedr_sep[-1]
        dTedr_ib = np.linspace(dTedr_ib_wall,dTedr_sep[0],inp.xi_ib_pts,endpoint=False)
        dTedr_ob = np.linspace(dTedr_sep[-1],dTedr_ob_wall,inp.xi_ob_pts,endpoint=True)
        
        BT_ib_wall = BT_sep[0]
        BT_ob_wall = BT_sep[-1]
        BT_ib = np.linspace(BT_ib_wall,BT_sep[0],inp.xi_ib_pts,endpoint=False)
        BT_ob = np.linspace(BT_sep[-1],BT_ob_wall,inp.xi_ob_pts,endpoint=True) 
        
        ni_xi    = np.concatenate((ni_ib,ni_sep,ni_ob))
        ne_xi    = np.concatenate((ne_ib,ne_sep,ne_ob))
        Ti_xi    = np.concatenate((Ti_ib,Ti_sep,Ti_ob))
        Te_xi    = np.concatenate((Te_ib,Te_sep,Te_ob))
        dnidr_xi = np.concatenate((dnidr_ib,dnidr_sep,dnidr_ob))
        dnedr_xi = np.concatenate((dnedr_ib,dnedr_sep,dnedr_ob))
        dTidr_xi = np.concatenate((dTidr_ib,dTidr_sep,dTidr_ob))
        dTedr_xi = np.concatenate((dTedr_ib,dTedr_sep,dTedr_ob))
        BT_xi    = np.concatenate((BT_ib,BT_sep,BT_ob))
        
        ib_leg_length = core.ib_div_line_cut.length
        ob_leg_length = core.ob_div_line_cut.length
        sep_length    = core.main_sep_line_closed.length
        ib_frac  = ib_leg_length / (ib_leg_length + sep_length + ob_leg_length)
        sep_frac = sep_length    / (ib_leg_length + sep_length + ob_leg_length)
        ob_frac  = ob_leg_length / (ib_leg_length + sep_length + ob_leg_length)
  

        xi_ib_div = np.linspace(0,
                                ib_frac+sep_frac*inp.ib_trim_off,
                                inp.xi_ib_pts,
                                endpoint=False)
        
        xi_sep    = np.linspace(ib_frac+sep_frac*inp.ib_trim_off,
                                ib_frac+sep_frac*inp.ib_trim_off + sep_frac-(inp.ib_trim_off + inp.ob_trim_off),
                                inp.xi_sep_pts,
                                endpoint=False)

        xi_ob_div    = np.linspace(ib_frac+sep_frac*inp.ib_trim_off + sep_frac-(inp.ib_trim_off + inp.ob_trim_off),
                                1,
                                inp.xi_ob_pts,
                                endpoint=True)
        
        xi_pts = np.concatenate((xi_ib_div,xi_sep,xi_ob_div))
                                        
        #model perpendicular particle and heat transport using Bohm Diffusion
        D_perp   =       Ti_xi / (16.0 * elementary_charge * BT_xi)
        Chi_perp = 5.0 * Ti_xi / (32.0 * elementary_charge * BT_xi)
        
        Gamma_perp = -D_perp * dnidr_xi
        Q_perp     = -ni_xi * Chi_perp * dTidr_xi - \
                     3.0 * Ti_xi * D_perp * dnidr_xi
        
        delta_sol_n = D_perp * ni_xi / Gamma_perp
        delta_sol_T = Chi_perp / \
                    (Q_perp/(ni_xi*Ti_xi) \
                     - 3.0*D_perp/delta_sol_n)
        delta_sol_E = 2/7*delta_sol_T
        #plt.plot(delta_sol_n)

        
        #plt.axis('equal')
        #plt.contourf(R_temp,
        #             Z_temp,
        #             dnidr,
        #             500)
        #plt.colorbar()
        #sys.exit()
        #plt.plot(xi_pts,ni_xi)
        #plt.plot(dTidr_sep_raw)
        #plt.plot(xi_pts,dTidr_sep_smooth)
        #plt.plot(xi_pts,np.nan_to_num(-Ti_J_sep_val/dTidr_sep_smooth))
        #plt.plot(xi_pts,BT_sep)
        #plt.plot(xi_pts,BT_sep_smooth)
        #plt.plot(xi_pts,D_perp,label='D_perp')
        #plt.plot(xi_pts,Chi_perp,label='Chi_perp')
        #plt.plot(xi_pts,Gamma_perp,label='Gamma_perp')
        #plt.plot(xi_pts,Q_perp,label='Q_perp')
        #plt.plot(xi_pts,-ni_sep_val * Chi_perp * dTidr_sep_smooth,label='Q term 1')
        #plt.plot(xi_pts,3.0 * Ti_J_sep_val * D_perp * dnidr_sep_smooth,label='Q term 2')
        #plt.plot(xi_pts,Q_perp,label='Q_perp')
        #plt.plot(xi_pts,3.0*D_perp*ni_sep_val*Ti_J_sep_val/delta_sol_n,label='term2')
        #plt.plot(xi_pts,delta_sol_n,label='delta_sol_n')
        #plt.plot(xi_pts,delta_sol_T,label='delta_sol_T')
        #plt.plot(xi_pts,delta_sol_E,label='delta_sol_E')
        #plt.legend()
        #plt.plot()
        #sys.exit()
        #delta_n_xi_ib  = np.array([0])
        #delta_n_xi_ib  = np.array([0])
        
        #pts = np.concatenate((np.array([0.0]),
        #                      np.linspace(ib_frac,ib_frac+sep_frac,xi_sep_pts),
        #                      np.array([1.0])))
        #vals = np.concatenate((delta_n_xi_ib,delta_sol_n_trim,delta_n_xi_ib))

        #delta_n_xi_sep = griddata(pts,
        #         vals,
        #         xi_pts,
        #         method='linear',
        #         fill_value='np.nan')


        r_max  = 0.45
        twoptdiv_r_pts = 20
        
        r_pts  = np.linspace(0,r_max,twoptdiv_r_pts)
        xi,r = np.meshgrid(xi_pts,r_pts)
        sol_ni =   ni_xi * np.exp(-r/delta_sol_n)
        sol_ne =   ne_xi * np.exp(-r/delta_sol_n)
        sol_Ti =   Ti_xi * np.exp(-r/delta_sol_T)
        sol_Te =   Te_xi * np.exp(-r/delta_sol_T)

        
        #draw sol lines through 2d strip model to get n,T along the lines
        sol_line_dist = np.zeros((len(xi_pts),len(self.sol_lines_cut)))
        sol_nT_pts = np.zeros((len(xi_pts),2,len(self.sol_lines_cut)))
        for i,sol_line in enumerate(self.sol_lines_cut):
            for j, xi_val in enumerate(xi_pts):
                sol_pt = sol_line.interpolate(xi_val,normalized=True)
                sol_nT_pts[j,:,i] = np.asarray(sol_pt.xy).T
                sep_pt_pos = core.entire_sep_line.project(sol_pt,normalized=True)
                sep_pt = core.entire_sep_line.interpolate(sep_pt_pos,normalized=True)
                sol_line_dist[j,i] = sol_pt.distance(sep_pt)
        
        sol_line_ni = np.zeros((len(xi_pts),len(self.sol_lines_cut)))
        sol_line_ne = np.zeros((len(xi_pts),len(self.sol_lines_cut)))
        sol_line_Ti = np.zeros((len(xi_pts),len(self.sol_lines_cut)))
        sol_line_Te = np.zeros((len(xi_pts),len(self.sol_lines_cut)))
        for i,sol_line in enumerate(self.sol_lines_cut):
            sol_line_ni[:,i] = griddata(np.column_stack((xi.flatten(),r.flatten())),
                                        sol_ni.flatten(),
                                        np.column_stack((np.linspace(0,1,len(xi_pts)),sol_line_dist[:,i])),
                                        method='linear')
            sol_line_ne[:,i] = griddata(np.column_stack((xi.flatten(),r.flatten())),
                                        sol_ne.flatten(),
                                        np.column_stack((np.linspace(0,1,len(xi_pts)),sol_line_dist[:,i])),
                                        method='linear')
            sol_line_Ti[:,i] = griddata(np.column_stack((xi.flatten(),r.flatten())),
                                        sol_Ti.flatten(),
                                        np.column_stack((np.linspace(0,1,len(xi_pts)),sol_line_dist[:,i])),
                                        method='linear')
            sol_line_Te[:,i] = griddata(np.column_stack((xi.flatten(),r.flatten())),
                                        sol_Te.flatten(),
                                        np.column_stack((np.linspace(0,1,len(xi_pts)),sol_line_dist[:,i])),
                                        method='linear')
        
        
        #append to master arrays
        #for i,line in enumerate(self.sol_lines_cut):
            pts_ni_sol = np.column_stack((sol_nT_pts[:,:,i],sol_line_ni[:,i]))
            pts_ne_sol = np.column_stack((sol_nT_pts[:,:,i],sol_line_ne[:,i]))
            pts_Ti_sol = np.column_stack((sol_nT_pts[:,:,i],sol_line_Ti[:,i]/1.0E3/1.6021E-19)) #converting back to kev
            pts_Te_sol = np.column_stack((sol_nT_pts[:,:,i],sol_line_Te[:,i]/1.0E3/1.6021E-19)) #converting back to kev
            
            core.ni_pts     = np.vstack((core.ni_pts,pts_ni_sol))
            core.ne_pts     = np.vstack((core.ne_pts,pts_ne_sol))
            core.Ti_kev_pts = np.vstack((core.Ti_kev_pts,pts_Ti_sol))
            core.Te_kev_pts = np.vstack((core.Te_kev_pts,pts_Te_sol))

        #draw wall line through 2d strip model to get n,T along the line
        wall_pts = np.asarray(inp.wall_line.xy).T
        ib_int_pt = np.asarray(core.ib_div_line.intersection(inp.wall_line).xy).T
        ob_int_pt = core.ob_div_line.intersection(inp.wall_line)
        wall_start_pos = np.where((wall_pts==ib_int_pt).all(axis=1))[0][0]
        wall_line_rolled = LineString(np.roll(wall_pts,-wall_start_pos,axis=0))
        wall_line_cut = cut(wall_line_rolled, 
                            wall_line_rolled.project(ob_int_pt,normalized=True))[0]
        #add points to wall line for the purpose of getting n,T along the wall. These points
        #won't be added to the main wall line or included in the triangulation.
        #for i,v in enumerate(np.linspace(0,1,300)):
        #    #interpolate along wall_line_cut to find point to add
        #    pt = wall_line_cut.interpolate(v,normalized=True)
        #    #add point to wall_line_cut
        #    union = wall_line_cut.union(pt)
        #    result = [geom for geom in polygonize(union)][0]
        #    wall_line_cut = LineString(result.exterior.coords)
        
        wall_nT_pts = np.asarray(wall_line_cut)
        num_wall_pts = len(wall_nT_pts)
        wall_pos_norm = np.zeros(num_wall_pts)
        wall_dist = np.zeros(num_wall_pts)

        for i,pt in enumerate(wall_nT_pts):
            wall_pt = Point(pt)
            sep_pt_pos = core.entire_sep_line.project(Point(wall_pt),normalized=True)
            sep_pt = core.entire_sep_line.interpolate(sep_pt_pos,normalized=True)
            wall_pos_norm[i] = wall_line_cut.project(wall_pt,normalized=True)
            wall_dist[i] = wall_pt.distance(sep_pt)
        
        wall_ni = griddata(np.column_stack((xi.flatten(),r.flatten())),
                           sol_ni.flatten(),
                           np.column_stack((wall_pos_norm,wall_dist)),
                           method='linear')
        wall_ne = griddata(np.column_stack((xi.flatten(),r.flatten())),
                           sol_ne.flatten(),
                           np.column_stack((wall_pos_norm,wall_dist)),
                           method='linear')
        wall_Ti = griddata(np.column_stack((xi.flatten(),r.flatten())),
                           sol_Ti.flatten(),
                           np.column_stack((wall_pos_norm,wall_dist)),
                           method='linear')
        wall_Te = griddata(np.column_stack((xi.flatten(),r.flatten())),
                           sol_Te.flatten(),
                           np.column_stack((wall_pos_norm,wall_dist)),
                           method='linear')
        
        #set minimum wall densities and temperatures
        #TODO: this needs to be more robust

        
        wall_ni[wall_ni < inp.wall_ni_min] = inp.wall_ni_min
        wall_ne[wall_ne < inp.wall_ne_min] = inp.wall_ne_min
        wall_Ti[wall_Ti < inp.wall_Ti_min] = inp.wall_Ti_min
        wall_Te[wall_Te < inp.wall_Te_min] = inp.wall_Te_min
        
        #append to master arrays
        pts_ni_wall = np.column_stack((wall_nT_pts,wall_ni))
        pts_ne_wall = np.column_stack((wall_nT_pts,wall_ne))
        pts_Ti_wall = np.column_stack((wall_nT_pts,wall_Ti)) #in kev
        pts_Te_wall = np.column_stack((wall_nT_pts,wall_Te)) #in kev
        
        core.ni_pts     = np.vstack((core.ni_pts,pts_ni_wall))
        core.ne_pts     = np.vstack((core.ne_pts,pts_ne_wall))
        core.Ti_kev_pts = np.vstack((core.Ti_kev_pts,pts_Ti_wall))
        core.Te_kev_pts = np.vstack((core.Te_kev_pts,pts_Te_wall))
        
        #plt.contourf(xi,r,np.log10(sol_ni),500)
        #plt.colorbar()
        #for i,v in enumerate(self.sol_lines_cut):
        #    plt.plot(xi_pts,sol_line_dist[:,i])
        #plt.plot(np.linspace(0,1,num_wall_pts),wall_dist,color='black')
        #sys.exit()
        