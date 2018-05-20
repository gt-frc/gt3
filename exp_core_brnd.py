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

class exp_core_brnd():
    def __init__(self,inp):
        
        R         = inp.psirz_exp[:,0].reshape(-1,65)
        Z         = inp.psirz_exp[:,1].reshape(-1,65)
        psi       = inp.psirz_exp[:,2].reshape(-1,65)        
        
        self.sep_lines(inp,R,Z,psi)
        self.core_lines(inp,R,Z,psi)
        self.core_nT(inp,R,Z,psi)
    
    def sep_lines(self,inp,R,Z,psi):
        #find x-point location  
        dpsidR = np.gradient(psi,R[0,:],axis=1)
        dpsidZ = np.gradient(psi,Z[:,0],axis=0)
        d2psidR2 = np.gradient(dpsidR,R[0,:],axis=1)
        d2psidZ2 = np.gradient(dpsidZ,Z[:,0],axis=0)
        
        #find line(s) where dpsidR=0
        self.dpsidR_0 = cntr.Cntr(R,Z,dpsidR).trace(0.0)
        #find line(s) where dpsidZ=0
        self.dpsidZ_0 = cntr.Cntr(R,Z,dpsidZ).trace(0.0)
    
        for i,path1 in enumerate(self.dpsidR_0):
            for j,path2 in enumerate(self.dpsidZ_0):
                try:
                    #find intersection points between curves for dpsidR=0 and dpsidZ=0
                    ints = LineString(path1).intersection(LineString(path2))
                    #if there is only one intersection ('Point'), then we're probably not
                    #dealing with irrelevant noise in psi
                    if ints.type=='Point':
                        #check if local maximum or minimum
                        d2psidR2_pt = griddata(np.column_stack((R.flatten(),Z.flatten())),
                                 d2psidR2.flatten(),
                                 [ints.x,ints.y],
                                 method='cubic')
                        d2psidZ2_pt = griddata(np.column_stack((R.flatten(),Z.flatten())),
                                 d2psidZ2.flatten(),
                                 [ints.x,ints.y],
                                 method='cubic')
                        
                        if d2psidR2_pt>0 and d2psidZ2_pt>0:
                            #we've found the magnetic axis
                            self.m_axis = np.array([ints.x,ints.y])
                        elif d2psidR2_pt<0 and d2psidZ2_pt<0:
                            #we've found a magnet. Do nothing.
                            pass
                        elif ints.y<0:
                            #we've probably found our x-point, although this isn't super robust
                            #and obviously only applies to a single-diverted, lower-null configuration
                            #TODO: make this more robust, I could easily see this failing on some shots
                            self.xpt = np.array([ints.x,ints.y])
                            
                        #uncomment this line when debugging
                        #print list(ints.coords),d2psidR2(ints.x,ints.y),d2psidZ2(ints.x,ints.y)
                except:
                    pass
        
        #normalize psi
        psi_shift = psi + abs(np.amin(psi)) #set center to zero
        psi_shift_xpt = griddata(np.column_stack((R.flatten(),Z.flatten())),
                                 psi_shift.flatten(),
                                 self.xpt,
                                 method='cubic')
        #psi_shift_xpt = interp2d(R,Z,psi_shift,kind='linear')(xpt[0],xpt[1]) #get new value at sep
        self.psi_norm = psi_shift / psi_shift_xpt
        
        #create lines for seperatrix and divertor legs of seperatrix
        num_lines = int(len(cntr.Cntr(R,Z,self.psi_norm).trace(1.0))/2)
        if num_lines==1:
            #in this case, the contour points that matplotlib returned constitute
            #a single line from inboard divertor to outboard divertor. We need to
            #add in the x-point in at the appropriate locations and split into a
            #main and a lower seperatrix line, each of which will include the x-point.
            x_psi,y_psi = draw_line(R,Z,self.psi_norm,1.0,0)
            
            
            loc1 = np.argmax(y_psi>self.xpt[1])
            loc2 = len(y_psi) - np.argmin(y_psi[::-1]<self.xpt[1])
    
            x_psi = np.insert(x_psi, (loc1,loc2), self.xpt[0])
            y_psi = np.insert(y_psi, (loc1,loc2), self.xpt[1])
            
            psi_1_pts = np.column_stack((x_psi,y_psi))
            self.main_sep_pts = psi_1_pts[loc1:loc2+1,:]
            self.main_sep_line = LineString(self.main_sep_pts[:-1])
            self.main_sep_line_closed = LineString(self.main_sep_pts)

            #get the inboard and outboard divertor legs seperately. This is so that
            #everything that includes the x-point can start with the x-point, which
            #elliminates the risk of tiny triangles in the vicinity of the x-point
            self.inboard_div_sep = np.flipud(psi_1_pts[:loc1+1])
            self.outboard_div_sep = psi_1_pts[loc2+1:]

            #cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.inboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ib_div_line = line
            self.ib_div_line_cut = cut(line,line.project(int_pt,normalized=True))[0]
            #self.ib_div_line_cut = line
            #TODO: add point to wall line
            
            
            #cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.outboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ob_div_line = line
            self.ob_div_line_cut = cut(line,line.project(int_pt,normalized=True))[0]
            
            ib_div_pts = np.flipud(np.asarray(self.ib_div_line_cut.xy).T)
            sep_pts    = np.asarray(self.main_sep_line.xy).T
            ob_div_pts = np.asarray(self.ob_div_line_cut.xy).T
            
            entire_sep_pts = np.vstack((ib_div_pts,sep_pts[1:,:],ob_div_pts))
            self.entire_sep_line = LineString(entire_sep_pts)

            #TODO: add point to wall line
                
        elif num_lines==2:
            #in this case, we have a lower seperatrix trace (line 0), and a main
            #seperatrix trace (line 1).
            
            #first do lower seperatrix line
            x_psi,y_psi = draw_line(R,Z,self.psi_norm,1.0,0)
            loc = np.argmax(x_psi>self.xpt[0])
            
            x_psi = np.insert(x_psi, loc, self.xpt[0])
            y_psi = np.insert(y_psi, loc, self.xpt[1])
            psi_1_pts = np.column_stack((x_psi,y_psi))
            
            self.inboard_div_sep = np.flipud(psi_1_pts[:loc+1])
            self.outboard_div_sep = psi_1_pts[loc+1:]
            
            #cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.inboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ib_div_line = line
            self.ib_div_line_cut = cut(line,line.project(int_pt,normalized=True))[0]

            #cut inboard line at the wall and add intersection point to wall_line
            line = LineString(self.outboard_div_sep)
            int_pt = line.intersection(inp.wall_line)
            self.ob_div_line = line
            self.ob_div_line_cut = cut(line,line.project(int_pt,normalized=True))[0]
            #TODO: add point to wall line
    
            #now to main seperatrix line
            x_psi,y_psi = draw_line(R,Z,self.psi_norm,1.0,1)
            self.main_sep_pts = np.insert(np.column_stack((x_psi,y_psi)),0,self.xpt,axis=0)
            self.main_sep_line = LineString(self.main_sep_pts[:-1])
            self.main_sep_line_closed = LineString(self.main_sep_pts)
            
            entire_sep_pts = np.vstack((ib_div_pts,sep_pts[1:,:],ob_div_pts))
            self.entire_sep_line = LineString(entire_sep_pts)
            #now clean up the lines by removing any points that are extremely close
            #to the x-point 
            #TODO: 

    def core_lines(self,inp,R,Z,psi):
        self.core_lines = []
        #psi_pts = np.concatenate((np.linspace(0,0.8,5,endpoint=False),np.linspace(0.8,1.0,4,endpoint=False)))
        psi_pts = np.linspace(inp.corelines_begin,1,inp.num_corelines,endpoint=False)
        for i,v in enumerate(psi_pts):
            num_lines = int(len(cntr.Cntr(R,Z,self.psi_norm).trace(v))/2)
            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x,y = draw_line(R,Z,self.psi_norm,v,0)
                self.core_lines.append(LineString(np.column_stack((x[:-1],y[:-1]))))
            else:
                #we need to find which of the surfaces is inside the seperatrix
                for j,line in enumerate(cntr.Cntr(R,Z,self.psi_norm).trace(v)[:num_lines]):
                #for j,line in enumerate(cntr.Cntr(R,Z,self.psi_norm).trace(v)):
                    x,y = draw_line(R,Z,self.psi_norm,v,j)
                    if (np.amax(x) < np.amax(self.main_sep_pts[:,0]) and \
                        np.amin(x) > np.amin(self.main_sep_pts[:,0]) and \
                        np.amax(y) < np.amax(self.main_sep_pts[:,1]) and \
                        np.amin(y) > np.amin(self.main_sep_pts[:,1])):
                        #then it's an internal flux surface
                        self.core_lines.append(LineString(np.column_stack((x[:-1],y[:-1]))))
                        break

    def core_nT(self,inp,R,Z,psi):
        
        #Master arrays that will contain all the points we'll use to get n,T
        #throughout the plasma chamber via 2-D interpolation
        self.ni_pts = np.zeros((0,3),dtype='float')
        self.ne_pts = np.zeros((0,3),dtype='float')
        self.Ti_kev_pts = np.zeros((0,3),dtype='float')
        self.Te_kev_pts = np.zeros((0,3),dtype='float')
        
        ##########################################
        #Calculate n,T throughout the core plasma using radial profile input files, uniform on flux surface
        
        ni     = UnivariateSpline(inp.ni_data[:,0],inp.ni_data[:,1],k=5,s=2.0)
        ne     = UnivariateSpline(inp.ne_data[:,0],inp.ne_data[:,1],k=5,s=2.0)
        Ti_kev = UnivariateSpline(inp.Ti_data[:,0],inp.Ti_data[:,1],k=5,s=2.0)
        Te_kev = UnivariateSpline(inp.Te_data[:,0],inp.Te_data[:,1],k=5,s=2.0)
        
        #get approximate rho values associated with the psi values we're using
        #draw line between magnetic axis and the seperatrix at the outboard midplane
        self.obmp_pt = self.main_sep_pts[np.argmax(self.main_sep_pts,axis=0)[0]]
        self.ibmp_pt = self.main_sep_pts[np.argmin(self.main_sep_pts,axis=0)[0]]
        self.top_pt  = self.main_sep_pts[np.argmax(self.main_sep_pts,axis=0)[1]]
        self.bot_pt  = self.main_sep_pts[np.argmin(self.main_sep_pts,axis=0)[1]]

        rho_line = LineString([Point(self.m_axis),Point(self.obmp_pt)])
        #for several points on the rho line specified above:
        
        #To get smooth gradients for use in the SOL calculation, you need around
        # 50-100 radial points in the far edge and around 100 or so theta points
        # TODO: There is almost certainly a faster way to get these gradients.
        rho_pts = np.concatenate((np.linspace(0, 0.95, 20, endpoint=False), 
                                 np.linspace(0.95, 1, 50, endpoint=False)),axis=0)
        
        thetapts = np.linspace(0,1,100,endpoint=False)
        for i,rho in enumerate(rho_pts): 
            #get n,T information at the point by interpolating the rho-based input file data
            ni_val = ni(rho)
            ne_val = ne(rho)
            Ti_kev_val = Ti_kev(rho)
            Te_kev_val = Te_kev(rho)
            #get R,Z coordinates of each point along the rho_line
            pt_coords = np.asarray(rho_line.interpolate(rho,normalized=True).coords)[0]

            #get psi value at that point
            psi_val = griddata(np.column_stack((R.flatten(),Z.flatten())),
                                 self.psi_norm.flatten(),
                                 pt_coords,
                                 method='linear')
            #map this n,T data to every point on the corresponding flux surface
            num_lines = int(len(cntr.Cntr(R,Z,self.psi_norm).trace(psi_val))/2)

            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x,y = draw_line(R,Z,self.psi_norm,psi_val,0)
                surf = LineString(np.column_stack((x,y)))
            else:
                #we need to find which of the surfaces is inside the seperatrix
                for j,line in enumerate(cntr.Cntr(R,Z,self.psi_norm).trace(psi_val)[:num_lines]):
                    #for j,line in enumerate(cntr.Cntr(R,Z,self.psi_norm).trace(v)):
                    x,y = draw_line(R,Z,self.psi_norm,psi_val,j)
                    if (np.amax(x) < np.amax(self.main_sep_pts[:,0]) and \
                        np.amin(x) > np.amin(self.main_sep_pts[:,0]) and \
                        np.amax(y) < np.amax(self.main_sep_pts[:,1]) and \
                        np.amin(y) > np.amin(self.main_sep_pts[:,1])):
                        #then it's an internal flux surface
                        surf = LineString(np.column_stack((x,y)))
                        break
            
            for j,theta_norm in enumerate(thetapts):
                pt = np.asarray(surf.interpolate(theta_norm,normalized=True).coords).T
                self.ni_pts = np.vstack((self.ni_pts,np.append(pt,ni_val)))
                self.ne_pts = np.vstack((self.ne_pts,np.append(pt,ne_val)))
                self.Ti_kev_pts = np.vstack((self.Ti_kev_pts,np.append(pt,Ti_kev_val)))
                self.Te_kev_pts = np.vstack((self.Te_kev_pts,np.append(pt,Te_kev_val)))

        #Do seperatrix separately so we don't accidentally assign the input n,T data to the divertor legs
        self.ni_sep_val = ni(1.0)
        self.ne_sep_val = ne(1.0)
        self.Ti_kev_sep_val = Ti_kev(1.0)
        self.Te_kev_sep_val = Te_kev(1.0)
        self.Ti_J_sep_val = self.Ti_kev_sep_val * 1.0E3 * 1.6021E-19
        self.Te_J_sep_val = self.Te_kev_sep_val * 1.0E3 * 1.6021E-19
        for j,theta_norm in enumerate(thetapts): 
            pt = np.asarray(self.main_sep_line.interpolate(theta_norm,normalized=False).coords,dtype='float').T
            self.ni_pts = np.vstack((self.ni_pts,np.append(pt,self.ni_sep_val)))
            self.ne_pts = np.vstack((self.ne_pts,np.append(pt,self.ne_sep_val)))
            self.Ti_kev_pts = np.vstack((self.Ti_kev_pts,np.append(pt,self.Ti_kev_sep_val)))
            self.Te_kev_pts = np.vstack((self.Te_kev_pts,np.append(pt,self.Te_kev_sep_val)))





