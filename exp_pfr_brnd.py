#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:04:43 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib._cntr as cntr
from scipy.interpolate import griddata, UnivariateSpline
from scipy.constants import elementary_charge
from shapely.geometry import Point, LineString
from shapely.ops import polygonize,linemerge
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

class exp_pfr_brnd():
    def __init__(self,inp,core):
        
        R = inp.psirz_exp[:,0].reshape(-1,65)
        Z = inp.psirz_exp[:,1].reshape(-1,65)
        
        self.pfr_lines(inp,R,Z,core)
        self.pfr_nT(inp,core)
        
    def pfr_lines(self,inp,R,Z,core):
        num_lines = int(len(cntr.Cntr(R,Z,core.psi_norm_raw).trace(0.999))/2)
        if num_lines==1:
            #then we're definitely dealing with a surface inside the seperatrix
            print 'Did not find PFR flux surface. Stopping.'
            sys.exit()
        else:
            #we need to find the surface that is contained within the private flux region
            for j,line in enumerate(cntr.Cntr(R,Z,core.psi_norm_raw).trace(0.99)[:num_lines]):
            #for j,line in enumerate(cntr.Cntr(R,Z,self.psi_norm_raw).trace(v)):
                x,y = draw_line(R,Z,core.psi_norm_raw,0.99,j)
                if (np.amax(y) < np.amin(core.main_sep_pts[:,1])):
                    #then it's a pfr flux surface, might need to add additional checks later
                    pfr_line_raw = LineString(np.column_stack((x,y)))
                    #find cut points
                    cut_pt1 = pfr_line_raw.intersection(inp.wall_line)[0]
                    dist1   = pfr_line_raw.project(cut_pt1,normalized=True)
                    cutline_temp  = cut(pfr_line_raw,dist1)[1]
                    
                    #reverse line point order so we can reliably find the second intersection point
                    cutline_temp_rev = LineString(np.flipud(np.asarray(cutline_temp.xy).T))
                    
                    cut_pt2 = cutline_temp_rev.intersection(inp.wall_line)
                    dist2   = cutline_temp_rev.project(cut_pt2,normalized=True)
                    cutline_final_rev  = cut(cutline_temp_rev,dist2)[1] 

                    #reverse again for final pfr flux line
                    pfr_flux_line = LineString(np.flipud(np.asarray(cutline_final_rev.xy).T))
                    
                    #add pfr_line intersection points on inboard side
                    #for some reason, union freaks out when I try to do inboard and outboard
                    #at the same time.
                    union = inp.wall_line.union(cut(pfr_line_raw,0.5)[0])    
                    result = [geom for geom in polygonize(union)][0]
                    inp.wall_line = LineString(result.exterior.coords)
            
                    #add pfr line intersection points on outboard side   
                    union = inp.wall_line.union(cut(pfr_line_raw,0.5)[1])    
                    result = [geom for geom in polygonize(union)][0]
                    inp.wall_line = LineString(result.exterior.coords)

                    #cut out pfr section of wall line
                    wall_pts = np.asarray(inp.wall_line.xy).T

                    #ib_int_pt = np.asarray(self.ib_div_line.intersection(inp.wall_line).xy).T
                    #ob_int_pt = self.ob_div_line.intersection(inp.wall_line)
                    wall_start_pos = np.where((wall_pts==cut_pt2).all(axis=1))[0][0]
                    wall_line_rolled = LineString(np.roll(wall_pts,-wall_start_pos,axis=0))
                    wall_line_cut_pfr = cut(wall_line_rolled, 
                                             wall_line_rolled.project(cut_pt1,normalized=True))[0]
                    
                    #create LineString with pfr line and section of wall line
                    self.pfr_line = linemerge((pfr_flux_line,wall_line_cut_pfr))
                    break     

    def pfr_nT(self,inp,core):
        
        pfr_pts = np.asarray(self.pfr_line.xy).T
        pfr_ni = np.zeros(len(pfr_pts)) + inp.pfr_ni_val
        pfr_ne = np.zeros(len(pfr_pts)) + inp.pfr_ne_val
        pfr_Ti = np.zeros(len(pfr_pts)) + inp.pfr_Ti_val
        pfr_Te = np.zeros(len(pfr_pts)) + inp.pfr_Te_val
        
        pts_ni_pfr = np.column_stack((pfr_pts,pfr_ni))
        pts_ne_pfr = np.column_stack((pfr_pts,pfr_ne))
        pts_Ti_pfr = np.column_stack((pfr_pts,pfr_Ti))
        pts_Te_pfr = np.column_stack((pfr_pts,pfr_Te))
        
        core.ni_pts     = np.vstack((core.ni_pts,pts_ni_pfr))
        core.ne_pts     = np.vstack((core.ne_pts,pts_ne_pfr))
        core.Ti_kev_pts = np.vstack((core.Ti_kev_pts,pts_Ti_pfr))
        core.Te_kev_pts = np.vstack((core.Te_kev_pts,pts_Te_pfr))
        
        #grid_x, grid_y = np.mgrid[1:2.5:500j, -1.5:1.5:500j]
        #ni_for_plot = griddata(self.ni_pts[:,:-1],self.ni_pts[:,-1],(grid_x,grid_y))
        #Ti_for_plot = griddata(self.Ti_kev_pts[:,:-1],self.Ti_kev_pts[:,-1],(grid_x,grid_y))
        #plt.contourf(grid_x,grid_y,np.log10(Ti_for_plot),500)
        #plt.colorbar()
        #sys.exit()
        pass