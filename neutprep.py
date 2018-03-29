# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:54:44 2017

@author: max
"""
import numpy as np
from helpers import getangle, getpoint, cut, isinline, getangle3ptsdeg
from scipy import interpolate
from math import pi, cos, sin, sqrt, degrees
from shapely.geometry import LineString, Polygon, Point, MultiPoint
from fipy import * #Gmsh2D, CellVariable, DiffusionTerm, ImplicitSourceTerm, Viewer
import matplotlib.pyplot as plt
import os
from subprocess import call
import re
import sys
from neutpy import neutpy,neutpyplot

class neutprep():
    """Prepare input information for neutpy and run neutpy
    
    Dependencies:
        Triangle (http://www.cs.cmu.edu/%7Equake/triangle.html)
    
    Methods:
        solhalomesh
        plasmamesh
        solnT
        pfrnT
        prep_gtneut_input
        prep_neutpy_input
        
    Attributes:
        sol_pts
        sol_lim_pts
        sol_flx_line
        id_int_pt
        od_int_pt
        sep_line
        sol_flx_line
        sol_segs
        divsep_pts
        divsep_segs
        plasma_pts
        plasma_segs
        plasma_lim_pts
        plasma_line
        plasma_param
        
    """
    def __init__(self,brnd,param):
        self.skipsolhalo = 1
        self.solhalomesh(brnd,param)
        self.plasmamesh(brnd,param)
        self.solnT(brnd,param)
        self.pfrnT(brnd,param)
        self.prep_neutpy_input(brnd,param)
        self.ntrl2grid(param,brnd)
        
    def solhalomesh(self,brnd,p):
        '''make the mesh to give to gtneut'''
        EF_ID = 0.05 #Expansion factor inner divertor
        EF_OD = 0.05 #Expansion factor outer divertor
        
        #####################################################################
        ## DEFINE SOL FLUX LINES
        #####################################################################
        num_flx_lines = 1
        num_flx_pts = 200
    
        flx_line_count = 0
        self.sol_pts = np.zeros((num_flx_lines*num_flx_pts,2))
        self.sol_lim_pts = np.zeros((0,2))
        self.sol_flx_line = np.zeros((num_flx_pts,0))
        
        #for each flux line:
        for count,wscale in enumerate(np.linspace(0.1,0.99,num_flx_lines)):
            ## ADD IN POINTS BASED ON CALCULATED WIDTHS FOR SOME POINTS IN INNER DIV REGION IN REVERSE ORDER
            
            #draw a line from the x-point down into the inner divertor and find the intersection point
            #with the first wall
            id_line_long = LineString([p.xpt,getpoint(p.xpt,p.xtheta1,2.0,'s')])
            id_int_pt_s = id_line_long.intersection(p.lim_line)
            self.id_int_pt = [id_int_pt_s.coords.xy[0][0],id_int_pt_s.coords.xy[1][0]]
            
            #create a bunch (i.e. 1000) of points along the line and get the outward normal
            #angle for all of them (will be the same for all of them).
            id_seppts = np.zeros((1000,2))
            for i,val in enumerate(np.linspace(0,1,1000)):
                id_seppts[i] = np.asarray(id_line_long.interpolate(val, normalized=True).coords.xy).T[0]

            onmls = np.ones(1000)*p.xtheta1-pi/2
            
            thetavals = np.linspace(0, 1, 1000)
            conv_dist_id = 0.5
            condition = np.where(thetavals<=(1-conv_dist_id),1.0,wscale/conv_dist_id)
            condition = np.where(thetavals>=(conv_dist_id))
            #condition2 = np.where(thetavals>=(wscale)**1.5)
            seppts1 = np.flipud(id_seppts[condition])
            onmls1 = np.flipud(onmls[condition])
            #id_solpts1 = getpoint(seppts1,onmls1,condition*EF_ID/cos(3*pi/2-p.xtheta1),'n')
            id_solpts1 = getpoint(seppts1,onmls1,wscale*EF_ID,'n')

            ## ADD IN POINTS BASED ON CALCULATED WIDTHS FOR SOME POINTS THETA IN [0,5pi/4] IN REVERSE ORDER
            seppts = np.column_stack((brnd.R[-1],brnd.Z[-1]))
            self.sep_line = LineString(seppts)
            onmls = getangle(np.roll(seppts,-1,axis=0),np.roll(seppts,+1,axis=0))-pi/2 #outward normal angle for each seperatrix point
            onmls = np.where(onmls<=0,onmls + 2*pi,onmls)
            
            thetavals = np.linspace(0, 2*pi, p.thetapts)
            condition = np.where(thetavals<=(5+(1-wscale**1.5))*pi/4)
            seppts1 = np.flipud(seppts[condition])
            onmls1 = np.flipud(onmls[condition])
            solpts1 = getpoint(seppts1,onmls1,(wscale**1)*0.04,'n')
            
            ## ADD IN POINTS BASED ON CALCULATED WIDTHS FOR THETA IN [7pi/4,2pi]
            condition = np.where(thetavals>=(6+(wscale**1.5))*pi/4-pi/200)
            seppts2 = np.flipud(seppts[condition])
            onmls2 = np.flipud(onmls[condition])
            solpts2 = getpoint(seppts2,onmls2,(wscale**1)*0.04,'n')
            
            ## ADD IN POINTS BASED ON CALCULATED WIDTHS FOR SOME POINTS IN OUTER DIV REGION
            od_line_long = LineString([p.xpt,getpoint(p.xpt,p.xtheta4,2.0,'s')])
            od_int_pt_s = od_line_long.intersection(p.lim_line)

            self.od_int_pt = [od_int_pt_s.coords.xy[0][0],od_int_pt_s.coords.xy[1][0]] #used later       
            
            #dist = od_line_long.project(od_line_long.intersection(p.lim_line))
            #od_line = cut(od_line_long,dist)[0]
            od_seppts = np.zeros((1000,2))
            for i,val in enumerate(np.linspace(0,1,1000)):
                pt = od_line_long.interpolate(val, normalized=True)
                od_seppts[i,0] = pt.coords.xy[0][0]
                od_seppts[i,1] = pt.coords.xy[1][0]
    
            #onmls = np.linspace(p.xtheta4 - 3*pi/2,0,1000)
            onmls = np.ones(1000)*p.xtheta4 - 3*pi/2
            
            thetavals = np.linspace(0, 1, 1000)
            conv_dist_id = 0.3
            condition = np.where(thetavals>=wscale**0.1)
            condition = np.where(thetavals>=(conv_dist_id))
            seppts1 = od_seppts[condition]
            onmls1 = onmls[condition]
            #od_solpts1 = getpoint(seppts1,onmls1,0.99*wscale*EF_OD/sin(2*pi-p.xtheta4),'n')
            od_solpts1 = getpoint(seppts1,onmls1,wscale*EF_OD,'n')
            
            ## COMBINE ALL THE POINTS AND INTERPOLATE
            solpts = np.vstack((id_solpts1,solpts1,solpts2,od_solpts1))
            tck,u     = interpolate.splprep( [solpts[:,0],solpts[:,1]], s = 0) #try setting k=5 at some point
            sol_x,sol_y = interpolate.splev( np.linspace( 0, 1, num_flx_pts), tck,der = 0)
            solpts = np.column_stack((sol_x,sol_y))
            self.sol_flx_line = np.column_stack((self.sol_flx_line,solpts))
            
            
            ## MAKE A LINE AND TRIM IT TO ONLY BE WITHIN PLASMA CHAMBER
            sol_line_long = LineString(solpts)
            
            x1, y1 = sol_line_long.xy
            x2, y2 = p.lim_line.xy
            
            #plt.plot(x1,y1)
            #plt.plot(x2,y2)
            
            int_pt1 = sol_line_long.intersection(p.lim_line)[0]
            int_pt2 = sol_line_long.intersection(p.lim_line)[1]
                
            dist1 = sol_line_long.project(int_pt1, normalized=True)
            sol_line = cut(sol_line_long,dist1)[1]
            
            dist2 = sol_line.project(int_pt2, normalized=True)
            sol_line = cut(sol_line,dist2)[0]
            
            #GET COORDINATES FOR PLOTTING FLUX LINES
            sol_x, sol_y = np.asarray(sol_line.coords.xy[0]), np.asarray(sol_line.coords.xy[1])    
            ## INTERPOLATE ALONG FLUX LINES AND GET POINT COORDINATES
            for i,val in enumerate(np.linspace(0,1,num_flx_pts)):
                pt = sol_line.interpolate(val, normalized=True)
                self.sol_pts[i+flx_line_count*num_flx_pts,0] = pt.coords.xy[0][0]
                self.sol_pts[i+flx_line_count*num_flx_pts,1] = pt.coords.xy[1][0]
                if val==0 or val==1:
                    self.sol_lim_pts = np.vstack((self.sol_lim_pts,self.sol_pts[i+flx_line_count*num_flx_pts]))
    
            flx_line_count +=1
          
        #####################################################################       
        ## CREATE SEGMENTS FOR TRIANGULATION       
        #####################################################################       
        line1segs = np.column_stack((np.linspace(1,num_flx_pts-1,num_flx_pts-1),np.linspace(2,num_flx_pts,num_flx_pts-1)))
        self.sol_segs = line1segs
        for i in np.linspace(1,num_flx_lines-1,num_flx_lines-1):
            self.sol_segs = np.vstack((self.sol_segs,line1segs + i*num_flx_pts))
            
        #####################################################################       
        ## CREATE INNER AND OUTER DIVERTOR LEGS OF THE SEPERATRIX     
        #####################################################################            
        self.divsep_pts = np.vstack((self.id_int_pt,p.xpt,self.od_int_pt))
        self.divsep_segs = np.array([[1,2],[2,3]])
        
        if self.skipsolhalo == 1:
            self.sol_pts = np.zeros((0,2))
            self.sol_segs = np.zeros((0,2))
            self.divsep_pts = np.zeros((0,2))
            self.divsep_segs = np.zeros((0,2))

    def plasmamesh(self,brnd,p):
        '''make the mesh to give to gtneut'''
        radmesh_start = 0.8
        # GET POINTS FOR TRIANGULATION
        R_tri = np.delete(brnd.R,-1,axis=1)[np.delete(brnd.r,-1,axis=1)>=radmesh_start*p.a].flatten()
        Z_tri = np.delete(brnd.Z,-1,axis=1)[np.delete(brnd.r,-1,axis=1)>=radmesh_start*p.a].flatten()

        self.plasma_pts = np.column_stack((R_tri,Z_tri))
    
        # GET SEGMENTS FOR TRIANGULATION
        line1seg = np.column_stack((
                                    np.linspace(1,brnd.R.shape[1]-1,brnd.R.shape[1]-1),
                                    np.roll(np.linspace(1,brnd.R.shape[1]-1,brnd.R.shape[1]-1),-1,axis=0)
                                    ))

        rmeshnum = (brnd.r[:,0] >= 0.8*p.a).sum()
        
        self.plasma_segs = line1seg
        for i,val in enumerate(np.linspace(1,rmeshnum-1,rmeshnum-1)):
            self.plasma_segs = np.vstack((self.plasma_segs,line1seg+val*(brnd.r.shape[1]-1)))
        
        # ADD POINTS FOR INBOARD AND OUTBOARD DIVERTOR LEGS
        self.plasma_pts = np.vstack((self.plasma_pts,self.id_int_pt,self.od_int_pt))


        # ADD SEGMENTS FOR INBOARD AND OUTBOARD DIVERTOR LEGS
        #find xpt in plasma_pts
        #xpt_loc = 0
        #for i,val in enumerate(plasma_pts):
        #    if val[0]==1.5 and val[1]==-1.2:
        #        xpt_loc = i+1
        #        break
    

        xpt_loc = len(self.plasma_pts) - 2 - (p.thetapts-1)/4
        id_int_pt_loc = self.plasma_pts.shape[0]-1
        od_int_pt_loc = self.plasma_pts.shape[0]
        
        self.plasma_segs = np.vstack((self.plasma_segs,[id_int_pt_loc,xpt_loc],[od_int_pt_loc,xpt_loc]))
        
        # GET POINTS THAT NEED TO BE ADDED TO LIM_LINE
        self.plasma_lim_pts = np.vstack((self.id_int_pt,self.od_int_pt))
        
        #CREATE SHAPELY LINE THAT CONSISTS OF POINTS ON THE INNER SIDE OF THE MESH.
        #USED LATER WHEN MAKING GTNEUT INPUT FILE (GTNEUT NEEDS TO KNOW WHICH CELLS
        #BORDER ON THE PLASMA REGION, WHICH IS CHECKED BY SEEING IF A POINT ON THE
        #FACES OF A CELL LIE WITHIN OR VERY NEAR TO THE LINE WE'RE CREATING RIGHT 
        #NOW)
        R_tri_in = R_tri[:(brnd.R.shape[1]-1)]
        Z_tri_in = Z_tri[:(brnd.Z.shape[1]-1)]
        plasma_in_pts = np.vstack((np.column_stack((R_tri_in,Z_tri_in)),np.column_stack((R_tri_in,Z_tri_in))[0]))
        self.plasma_line = LineString(plasma_in_pts)
        self.plasma_param = np.column_stack((brnd.R.flatten(),brnd.Z.flatten(),brnd.ni.flatten(),brnd.ne.flatten(),brnd.Ti_kev.flatten(),brnd.Te_kev.flatten()))

    def solnT(self,brnd,p):         
        #####################################################################       
        ## POPULATE DENSITY AND TEMPERATURE VALUES FOR THE SOL AND HALO REGIONS  
        ##################################################################### 
    
        ## THESE WILL GET COMBINED WITH AN UNSTRUCTURED LIST OF DENSITIES AND TEMPERATURES
        ## AT EVERY OTHER VERTEX IN THE PROBLEM AND USED IN A 2D INTERPOLATION
        ## TO GET THE VALUES AT THE MIDDLES OF THE TRIANGLES GENERATED FOR GTNEUT
        
        ####################################################################
        # ADD INT POINTS TO LIM_LINE (DELETE THIS LATER OR PUT IT IN A BETTER PLACE)
        ####################################################################
        
        all_pts_limx, all_pts_limy = p.lim_line.coords.xy
        
        all_pts_lim = np.column_stack((all_pts_limx,all_pts_limy))
        for testpt in [self.id_int_pt,self.od_int_pt]:
            for indx,point in enumerate(all_pts_lim):
                if indx <= len(all_pts_lim)-2:
                    pt1 = all_pts_lim[indx]
                    pt2 = all_pts_lim[indx+1]
                    dist_p1_p2 = sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)
                    dist_p1_testpt = sqrt((pt1[0]-testpt[0])**2 + (pt1[1]-testpt[1])**2)
                    dist_testpt_p2 = sqrt((testpt[0]-pt2[0])**2 + (testpt[1]-pt2[1])**2)
                    if abs(dist_p1_p2 - (dist_p1_testpt+dist_testpt_p2)) < 1E-15:
                        all_pts_lim2 = np.insert(all_pts_lim,indx+1,testpt,axis=0)
                        all_pts_lim = all_pts_lim2
                        break    
        lim_line = LineString(all_pts_lim)
        limx, limy = lim_line.coords.xy
        
        ####################################################################
        # FIND POSITION OF XPOINT IN SEP_LINE COORDINATES AND ROTATE TO START GOING
        # AROUND THE SEPERATRIX COUNTER-CLOCKWISE
        ####################################################################
    
        sepx, sepy = self.sep_line.coords.xy
        sepx = np.where(sepx[:-1]<1E-10,0,sepx[:-1])
        sepy = np.where(sepy[:-1]<1E-10,0,sepy[:-1])

        xpt_pos = np.where(np.isclose(sepy, p.xpt[1],atol=5e-03)) and np.where(np.isclose(sepx, p.xpt[0],atol=5e-03))[0][0]

        sepx = np.roll(sepx,-xpt_pos)
        sepy = np.roll(sepy,-xpt_pos)
        
        #make new sep_line with x-point repeated for use in the bdry condition
        sep_line2 = LineString(np.vstack((np.column_stack((sepx,sepy)),p.xpt)))
        #start making mesh
        mesh_info = ""
        pt_count = 1
        for i, (xval, yval) in enumerate(zip(sepx, sepy)):
            mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(xval) + ', ' + str(yval) + ', 0, 1/10};')
            pt_count+=1
            
        ####################################################################
        # ADD IN LIM_LINE POINTS BETWEEN THE SEP-LIM INTERSECTION POINTS
        # ON THE INBOARD AND OUTBOARD SIDES
        ####################################################################

        #rotate lim_line to start at ib_int_pt
        start_pos = 0
        limx, limy = lim_line.coords.xy
        limx = np.delete(limx,-1)
        limy = np.delete(limy,-1)
        
        #for i, (xval, yval) in enumerate(zip(limx, limy)):
        #    if xval == self.id_int_pt[0] and yval == self.id_int_pt[1]:
        #        start_pos = i
        #        break
        start_pos = np.where(np.isclose(limx, self.id_int_pt[0],atol=5e-04)) and np.where(np.isclose(limy, self.id_int_pt[1],atol=5e-04))[0][0]
        limx = np.roll(limx,-start_pos)
        limy = np.roll(limy,-start_pos) 
        
        #make new lim_line for use in the bdry condition
        lim_line_new = LineString(np.column_stack((limx,limy)))
        dist = lim_line_new.project(Point(self.od_int_pt[0],self.od_int_pt[1]),normalized=True)
        lim_line_cut = cut(lim_line_new,dist)[0]
        limx, limy = lim_line_cut.coords.xy
            
        for i, (xval, yval) in enumerate(zip(limx, limy)):
            mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(xval) + ', ' + str(yval) + ', 0, 1/10};')
            pt_count+=1
            
        ####################################################################
        # CREATE LINES FOR MESH_INFO
        ####################################################################
        line_count = 1
        for i in np.linspace(1,len(sepx)-1,len(sepx)-1):
            mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(i)) + ', ' + str(int(i+1)) + '};')
            line_count+=1
            
        #return to the xpoint
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(len(sepx))) + ', ' + str(int(1)) + '};')
        line_count+=1
        
        #line from xpt to first wall on inboard side
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(1)) + ', ' + str(int(len(sepx)+1)) + '};')
        line_count+=1
        
        #lines around the first wall to the outboard side
        for i in np.linspace(int(len(sepx)+1),int(len(sepx)+1) + len(limx)-2,len(limx)-1):
            mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(i)) + ', ' + str(int(i+1)) + '};')
            line_count+=1
        
        #return to the xpoint on the outboard side
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(line_count-1)) + ', ' + str(int(1)) + '};')
     
       
        ####################################################################
        # CREATE LINELOOP AND PLANE SURFACE
        ####################################################################
        lineloop = "Line Loop(1) = {1"
        for i in np.linspace(2,line_count,line_count-1):
            lineloop = lineloop + ", " + str(int(i))
        lineloop = lineloop + "};"
        mesh_info = mesh_info + lineloop
        mesh_info = mesh_info + "Plane Surface(1) = {1};\n"
        # MESH_INFO COMPLETE!
        ####################################################################
        ## CREATE MESH OBJECT
        ####################################################################
        m = Gmsh2D(mesh_info)
        
        ####################################################################
        ## CREATE LINES FOR BOUNDARY CONDITIONS
        ####################################################################
        #already made lines for sep and first wall
        #now make them for inner and outer divertor legs
        idiv_line = LineString([(p.xpt[0],p.xpt[1]),(self.id_int_pt[0], self.id_int_pt[1])])
        odiv_line = LineString([(p.xpt[0],p.xpt[1]),(self.od_int_pt[0], self.od_int_pt[1])])
        
        var = CellVariable(mesh=m)
        
        ####################################################################
        ## LOOP OVER QUANTITIES TO CALCULATE (ni, ne, Ti, Te)
        ####################################################################
        self.sol_param = np.column_stack((m.cellCenters.value[0],m.cellCenters.value[1]))
        for j in [0,1,2,3]:
            
            mask_array = np.zeros(m.numberOfCells)
            bcval = np.zeros(m.numberOfCells)
            for i,val in enumerate(mask_array):
                face1num = m.cellFaceIDs.data[0,i]
                face2num = m.cellFaceIDs.data[1,i]
                face3num = m.cellFaceIDs.data[2,i]
                
                face1x = m.faceCenters.value[0,face1num]
                face2x = m.faceCenters.value[0,face2num]
                face3x = m.faceCenters.value[0,face3num]
                
                face1y = m.faceCenters.value[1,face1num]
                face2y = m.faceCenters.value[1,face2num]
                face3y = m.faceCenters.value[1,face3num]
                
                p1 = Point(face1x, face1y)
                p2 = Point(face2x, face2y)
                p3 = Point(face3x, face3y)
                
                ptinsol = sep_line2.distance(p1) < 1e-8 or sep_line2.distance(p2) < 1e-8 or sep_line2.distance(p3) < 1e-8
                ptinidiv = idiv_line.distance(p1) < 1e-8 or idiv_line.distance(p2) < 1e-8 or idiv_line.distance(p3) < 1e-8
                ptinodiv = odiv_line.distance(p1) < 1e-8 or odiv_line.distance(p2) < 1e-8 or odiv_line.distance(p3) < 1e-8
                ptinwall = lim_line_cut.distance(p1) < 1e-8 or lim_line_cut.distance(p2) < 1e-8 or lim_line_cut.distance(p3) < 1e-8
                
                if ptinsol or ptinidiv or ptinodiv or ptinwall:
                    mask_array[i] = 1
                
                if j==0:
                    solbc = brnd.ni[-1,0]*1E-19
                elif j==1:
                    solbc = brnd.ne[-1,0]*1E-19
                elif j==2:
                    solbc = brnd.Te_kev[-1,0]
                else:
                    solbc = brnd.Te_kev[-1,0]
                
                if ptinsol:
                    bcval[i] = solbc
                if ptinidiv:
                    bcval[i] = solbc
                if ptinodiv:
                    bcval[i] = solbc
                if ptinwall:
                    bcval[i] = solbc/100.0
                    
            mask = (CellVariable(mesh=m,value=mask_array)==1)
            
            largeValue = 1e+10
            value = CellVariable(mesh=m,value=bcval)
            eqn = DiffusionTerm(coeff=1.0) - ImplicitSourceTerm(largeValue * mask) + largeValue * mask * value - 200.0*(var - 0.06)*var 
            print ("solving for quantity: ",j)
            eqn.solve(var)
            
            if j==0 or j == 1:
                var_final = var.value*1E19
            else:
                var_final = var.value
    
            self.sol_param = np.column_stack((self.sol_param,var_final))
    
        #ADD DENSITIES AND TEMPERATURES ALONG THE BOUNDARY TO sol_param
        #FOR INCLUSION IN THE 2D INTERPOLATION
        nibc = np.asarray(limx)*0 + brnd.ni[-1,0]
        nebc = np.asarray(limx)*0 + brnd.ne[-1,0]
        Tibc = np.asarray(limx)*0 + brnd.Ti_kev[-1,0]
        Tebc = np.asarray(limx)*0 + brnd.Te_kev[-1,0]
        
        solbc_array = np.column_stack((limx,limy,nibc,nebc,Tibc,Tebc))
        self.sol_param = np.vstack((self.sol_param,solbc_array))
        
            
        viewer = Viewer(vars=var, datamin=0., datamax=np.amax(var_final))
        viewer.plot()
        plt.show()
    def pfrnT(self,brnd,p):         
        #start making mesh
        pt_count = 1
        mesh_info = ""
        mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(p.xpt[0]) + ', ' + str(p.xpt[1]) + ', 0, 1/20};')
        pt_count+=1
        mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(self.id_int_pt[0]) + ', ' + str(self.id_int_pt[1]) + ', 0, 1/30};')
        pt_count+=1
        mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(self.od_int_pt[0]) + ', ' + str(self.od_int_pt[1]) + ', 0, 1/30};')
        pt_count+=1
            
        ####################################################################
        # CREATE LINES FOR MESH_INFO
        ####################################################################
        line_count = 1
        
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(1)) + ', ' + str(int(2)) + '};')
        line_count+=1
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(2)) + ', ' + str(int(3)) + '};')
        line_count+=1
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(3)) + ', ' + str(int(1)) + '};')
        line_count+=1
     
       
        ####################################################################
        # CREATE LINELOOP AND PLANE SURFACE
        ####################################################################
        lineloop = "Line Loop(1) = {1"
        for i in np.linspace(2,line_count,line_count-1):
            lineloop = lineloop + ", " + str(int(i))
        lineloop = lineloop + "};"
        mesh_info = mesh_info + lineloop
        mesh_info = mesh_info + "Plane Surface(1) = {1};\n"
        # MESH_INFO COMPLETE!
        
        ####################################################################
        ## CREATE MESH OBJECT
        ####################################################################
        m = Gmsh2D(mesh_info)
        
        ####################################################################
        ## CREATE LINES FOR BOUNDARY CONDITIONS
        ####################################################################
        #already made lines for sep and first wall
        #now make them for inner and outer divertor legs
        idiv_line = LineString([(p.xpt[0],p.xpt[1]),(self.id_int_pt[0], self.id_int_pt[1])])
        odiv_line = LineString([(p.xpt[0],p.xpt[1]),(self.od_int_pt[0], self.od_int_pt[1])])
        bottom_line = LineString([(self.id_int_pt[0], self.id_int_pt[1]),(self.od_int_pt[0], self.od_int_pt[1])])
        
        var = CellVariable(mesh=m)
        
        ####################################################################
        ## LOOP OVER QUANTITIES TO CALCULATE (ni, ne, Ti, Te)
        ####################################################################
        self.pfr_param = np.column_stack((m.cellCenters.value[0],m.cellCenters.value[1]))
        for j in [0,1,2,3]:
            
            mask_array = np.zeros(m.numberOfCells)
            bcval = np.zeros(m.numberOfCells)
            for i,val in enumerate(mask_array):
                face1num = m.cellFaceIDs.data[0,i]
                face2num = m.cellFaceIDs.data[1,i]
                face3num = m.cellFaceIDs.data[2,i]
                
                face1x = m.faceCenters.value[0,face1num]
                face2x = m.faceCenters.value[0,face2num]
                face3x = m.faceCenters.value[0,face3num]
                
                face1y = m.faceCenters.value[1,face1num]
                face2y = m.faceCenters.value[1,face2num]
                face3y = m.faceCenters.value[1,face3num]
                
                p1 = Point(face1x, face1y)
                p2 = Point(face2x, face2y)
                p3 = Point(face3x, face3y)
                
                #ptinsol = sep_line2.distance(p1) < 1e-8 or sep_line2.distance(p2) < 1e-8 or sep_line2.distance(p3) < 1e-8
                ptinidiv = idiv_line.distance(p1) < 1e-8 or idiv_line.distance(p2) < 1e-8 or idiv_line.distance(p3) < 1e-8
                ptinodiv = odiv_line.distance(p1) < 1e-8 or odiv_line.distance(p2) < 1e-8 or odiv_line.distance(p3) < 1e-8
                ptinwall = bottom_line.distance(p1) < 1e-8 or bottom_line.distance(p2) < 1e-8 or bottom_line.distance(p3) < 1e-8
                
                if ptinidiv or ptinodiv or ptinwall:
                    mask_array[i] = 1
                
                if j==0:
                    solbc = brnd.ni[-1,0]*1E-19
                elif j==1:
                    solbc = brnd.ne[-1,0]*1E-19
                elif j==2:
                    solbc = brnd.Ti_kev[-1,0]
                else:
                    solbc = brnd.Te_kev[-1,0]
                
                if ptinidiv:
                    bcval[i] = solbc
                if ptinodiv:
                    bcval[i] = solbc
                if ptinwall:
                    bcval[i] = solbc/100.0
                    
            mask = (CellVariable(mesh=m,value=mask_array)==1)
            
            largeValue = 1e+10
            value = CellVariable(mesh=m,value=bcval)
            eqn = DiffusionTerm() - ImplicitSourceTerm(largeValue * mask) + largeValue * mask * value 
            print ("solving for quantity: ",j)
            eqn.solve(var)
            
            if j==0 or j == 1:
                var_final = var.value*1E19
            else:
                var_final = var.value
    
            self.pfr_param = np.column_stack((self.pfr_param,var_final))
    
        #ADD DENSITIES AND TEMPERATURES ALONG THE BOUNDARY TO sol_param
        #FOR INCLUSION IN THE 2D INTERPOLATION
        #nibc = np.asarray(limx)*0 + solnibc
        #nebc = np.asarray(limx)*0 + solnebc
        #Tibc = np.asarray(limx)*0 + solTibc
        #Tebc = np.asarray(limx)*0 + solTebc
        
        #solbc_array = np.column_stack((limx,limy,nibc,nebc,Tibc,Tebc))
        #sol_param = np.vstack((sol_param,solbc_array))
    
            
        viewer = Viewer(vars=var, datamin=0., datamax=np.amax(var_final))
        viewer.plot()
        plt.show()
    def prep_gtneut_input(self,brnd,p):
        ####################################################################
        ####################################################################
        ## CREATE VERSION OF FIRST WALL THAT INCLUDES FLUX LINE INT. PTS.
        ####################################################################
        ####################################################################
        #this code is also in sol_halo and I think one other place.
        # Should probably write a helper function at some point
        
        all_pts_limx, all_pts_limy = p.lim_line.coords.xy
        all_pts_lim = np.column_stack((all_pts_limx,all_pts_limy))
        
        if self.skipsolhalo==1:
            pts2add = np.zeros((0,2))
        else:
            pts2add = np.vstack((self.plasma_lim_pts,self.sol_lim_pts))
            
            for testpt in pts2add:
                for indx,point in enumerate(all_pts_lim):
                    if indx <= len(all_pts_lim)-2:
                        pt1 = all_pts_lim[indx]
                        pt2 = all_pts_lim[indx+1]
                        dist_p1_p2 = sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)
                        dist_p1_testpt = sqrt((pt1[0]-testpt[0])**2 + (pt1[1]-testpt[1])**2)
                        dist_testpt_p2 = sqrt((testpt[0]-pt2[0])**2 + (testpt[1]-pt2[1])**2)
                        if abs(dist_p1_p2 - (dist_p1_testpt+dist_testpt_p2)) < 1E-15:
                            all_pts_lim2 = np.insert(all_pts_lim,indx+1,testpt,axis=0)
                            all_pts_lim = all_pts_lim2
                            break
        
        lim_line = LineString(all_pts_lim)
        limx, limy = lim_line.coords.xy
        self.lim_pts = np.column_stack((limx, limy))
        self.lim_segs = np.column_stack((
                                    np.linspace(1,self.lim_pts.shape[0],self.lim_pts.shape[0]),
                                    np.roll(np.linspace(1,self.lim_pts.shape[0],self.lim_pts.shape[0]),-1,axis=0)
                                    ))
        
        ####################################################################
        ####################################################################
        ## MAKE ARRAY OF DENSITIES AND TEMPERATURES FOR INTERPOLATION
        ####################################################################
        ####################################################################
        tri_param = np.vstack((self.plasma_param,self.sol_param,self.pfr_param))
    
        ####################################################################
        ####################################################################
        ## COMBINE all_pts, segments, holes, AND attributes FOR TRIANGLE
        ####################################################################
        ####################################################################
    
        #vertices = np.vstack((plasma_pts,sol_pts,lim_pts))
        vertices = np.vstack((self.plasma_pts,self.sol_pts,self.lim_pts,self.divsep_pts))
        vert_number = np.linspace(1,vertices.shape[0],vertices.shape[0])
        vertices = np.column_stack((vert_number,vertices))    
        
        #segments = np.vstack((plasma_segs,sol_segs+plasma_pts.shape[0], lim_segs+plasma_pts.shape[0]+sol_pts.shape[0]))
        segments = np.vstack((
                              self.plasma_segs, 
                              self.sol_segs    +self.plasma_pts.shape[0],
                              self.lim_segs    +self.plasma_pts.shape[0] + self.sol_pts.shape[0],
                              self.divsep_segs +self.plasma_pts.shape[0] + self.sol_pts.shape[0] + self.lim_pts.shape[0]
                              ))
        seg_number = np.linspace(1,segments.shape[0],segments.shape[0])
        seg_attributes = seg_number * 0
        segments = np.column_stack((seg_number,segments,seg_attributes))
        
    
        ####################################################################
        ####################################################################
        ## OUTPUT .poly FILE AND RUN TRIANGLE PROGRAM
        ####################################################################
        ####################################################################
        open('mil_mesh.poly', 'w').close()
        outfile = open('mil_mesh.poly','ab')
        filepath = os.path.realpath(outfile.name)
        np.savetxt(outfile,np.array([vertices.shape[0],2,0,0])[None],fmt='%i %i %i %i')
        np.savetxt(outfile,vertices,fmt='%i %f %f')
        np.savetxt(outfile,np.array([segments.shape[0],0])[None],fmt='%i %i')
        np.savetxt(outfile,segments.astype(int),fmt='%i %i %i %i')
        np.savetxt(outfile,np.array([1])[None],fmt='%i')
        np.savetxt(outfile,np.array([1,p.R0_a,p.Z0])[None],fmt='%i %f %f')
        np.savetxt(outfile,np.array([0])[None],fmt='%i')
        outfile.close()
        
        #http://www.cs.cmu.edu/%7Equake/triangle.switch.html
        call(["triangle", "-pqn",filepath])
        
        ####################################################################
        ####################################################################
        ## READ TRIANGLE OUTPUT
        ####################################################################
        ####################################################################

        ## DECLARE FILE PATHS
        nodepath = os.getcwd() + '/mil_mesh.1.node'
        elepath = os.getcwd() + '/mil_mesh.1.ele'
        neighpath = os.getcwd() + '/mil_mesh.1.neigh'

        ## GET NODE DATA
        with open(nodepath, 'r') as node:
            #dummy = next(mil_mesh)
            nodecount = re.findall(r'\d+', next(node))
            nNodes = int(nodecount[0])
            nodenum = np.zeros(nNodes)
            nodesx = np.zeros(nNodes)
            nodesy = np.zeros(nNodes)
            
            for i in range (0,nNodes):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(node)) 
                nodenum[i] = int(data1[0])
                nodesx[i] = data1[1]
                nodesy[i] = data1[2]
                
        ## GET TRIANGLE DATA
        with open(elepath, 'r') as tri_file:
            tricount = re.findall(r'\d+', next(tri_file))
            nTri = int(tricount[0])
            triangles = np.zeros((nTri,3))
            tri_regions = np.zeros(nTri)
            for i in range (0,nTri):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(tri_file))
                triangles[i,0] = data1[1]
                triangles[i,1] = data1[2]
                triangles[i,2] = data1[3]
                #tri_regions[i] = data1[4]
        triangles = triangles.astype('int')
        tri_regions = tri_regions.astype('int')

        ## GET NEIGHBOR DATA
        with open(neighpath, 'r') as neigh_file:
            neighcount = re.findall(r'\d+', next(neigh_file))
            nNeigh = int(neighcount[0])
            neighbors = np.zeros((nNeigh,3))
            for i in range (0,nNeigh):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(neigh_file))
                neighbors[i,0] = data1[1]
                neighbors[i,1] = data1[2]
                neighbors[i,2] = data1[3]
        neighbors = neighbors.astype('int')
        
        ## REARRANGE TRIANGLES TO CONFORM TO GTNEUT CONVENTION      
        triangles = np.fliplr(triangles) #triangle vertices are given counterclockwise, but we want clockwise
        neighbors = np.fliplr(neighbors) #neighbor 1 is opposite vertex 1, so also counterclockwise 
    
        y=np.zeros(3)
        for i,tri in enumerate(triangles):
            # Find lowest value of y component of vertices
            y[0] = nodesy[tri[0]-1]
            y[1] = nodesy[tri[1]-1]
            y[2] = nodesy[tri[2]-1]
            miny = np.amin(y)
            miny_count = np.sum(y == miny)
            if miny_count == 1:
                #identify position of minimum and roll array accordingly
                miny_index = np.where(y==miny)[0][0]
            else:
                #identify which points are the two minima and determine
                #which of them is farthest to the left (or right if I change it)
                miny_index = np.where(y==miny)[0][1] #change this 1 to a zero to choose the rightmost of the two bottom vertices
            triangles[i] = np.roll(triangles[i],-1*miny_index)
            neighbors[i] = np.roll(neighbors[i],-1*miny_index-2) # the -2 is because the side 1 is opposite vertex 1. We want side 1 to start at vertex 1

        ## GET VALUES TO ORIENT THE FIRST CELL WHEN PLOTTING
        point1_x = nodesx[triangles[0,0]-1]
        point1_y = nodesx[triangles[0,0]-1]
        point2_x = nodesx[triangles[0,1]-1]
        point2_y = nodesx[triangles[0,1]-1]
        point3_x = nodesx[triangles[0,2]-1]
        point3_y = nodesx[triangles[0,2]-1]
        
        center_x = (point1_x + point2_x + point3_x) / 3
        center_y = (point1_y + point2_y + point3_y) / 3

        ## CALCULATE ANGLE BY WHICH TO ROTATE THE FIRST CELL WHEN PLOTTING
        theta0 = degrees(getangle([point1_x,point1_y],[point3_x,point3_y]))
        
        
        ## CALCULATE MID POINTS OF TRIANGLES, AS WELL AS MIDPOINTS FOR EACH FACE     
        ptsx = np.zeros((nTri,3))
        ptsy = np.zeros((nTri,3))
        #for index,tri in ndenumerate(triangles):
        for i in range(0,nTri):
            ptsx[i,0] = nodesx[triangles[i,0]-1]
            ptsy[i,0] = nodesy[triangles[i,0]-1]
            ptsx[i,1] = nodesx[triangles[i,1]-1]
            ptsy[i,1] = nodesy[triangles[i,1]-1]
            ptsx[i,2] = nodesx[triangles[i,2]-1]
            ptsy[i,2] = nodesy[triangles[i,2]-1]
        
        mid_x = np.mean(ptsx,axis=1)
        mid_y = np.mean(ptsy,axis=1)
        midpts = np.column_stack((mid_x,mid_y))
        
        side1_midx = (ptsx[:,0] + ptsx[:,1])/2
        side2_midx = (ptsx[:,1] + ptsx[:,2])/2
        side3_midx = (ptsx[:,2] + ptsx[:,0])/2
        
        side1_midy = (ptsy[:,0] + ptsy[:,1])/2
        side2_midy = (ptsy[:,1] + ptsy[:,2])/2
        side3_midy = (ptsy[:,2] + ptsy[:,0])/2
        
        side1_midpt = np.column_stack((side1_midx,side1_midy))
        side2_midpt = np.column_stack((side2_midx,side2_midy))
        side3_midpt = np.column_stack((side3_midx,side3_midy))
        
        #COMBINE POINTS FOR THE PLASMA, SOL, AND DIVERTOR REGIONS
        #first fill in plasma cells
        plasmacells = np.zeros((1,2))
        pcellnum = nTri + 1
        pcellcount = 0
        for index,nei in enumerate(neighbors):
            #for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index],self.plasma_line)
            side2inline = isinline(side2_midpt[index],self.plasma_line)
            side3inline = isinline(side3_midpt[index],self.plasma_line)
            
            if side1inline or side2inline or side3inline:
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                if nb == 1: #cell has one plasma border
                    
                    #create plasma cell
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index + 1
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1
                elif nb == 2: #cell has two plasma borders (this will probably never happen. It would require a local concavity in the inner-most meshed flux surface)
                    #create plasma cell #1
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index + 1
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1   
                    
                    #create plasma cell #2
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index + 1
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1
        plasmacells = np.delete(plasmacells,-1,0)
        plasmacells = plasmacells.astype('int')
        
        #now fill in wall cells
        wallcells = np.zeros((1,2))
        wcellnum = pcellnum #was already advanced in the plasmacell loop. Don't add 1.
        wcellcount = 0
        for index, nei in enumerate(neighbors):
            #for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index],p.lim_line)
            side2inline = isinline(side2_midpt[index],p.lim_line)
            side3inline = isinline(side3_midpt[index],p.lim_line)
            
            if side1inline or side2inline or side3inline:
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                if nb == 1: #cell has one wall border
                    #create wall cell
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index + 1
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1
                elif nb == 2: #cell has two wall borders (This can easily happen because the wall has many concave points.)
                    #create wall cell #1
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index + 1
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1   
                    
                    #create wall cell #2
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index + 1
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1
        wallcells = np.delete(wallcells,-1,0)
        wallcells = wallcells.astype('int')
        
        ## POPULATE CELL DENSITIES AND TEMPERATURES
        #create array of all points in plasma, sol, id, and od
        #tri_param = np.vstack((plasma_param,sol_param,id_param,od_param))
        
        ni_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,2],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)
        ne_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,3],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)
        Ti_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,4],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)
        Te_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,5],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)
                
        ni_tri = np.where(ni_tri<1.0E16,1.0E16,ni_tri)
        ne_tri = np.where(ni_tri<1.0E16,1.0E16,ne_tri)
        Ti_tri = np.where(ni_tri<0.002,0.002,Ti_tri)
        Te_tri = np.where(ni_tri<0.002,0.002,Te_tri)
        ## CALCULATE DENSITIES AND TEMPERATURES TO ASSIGN TO PLASMA CORE
        ni_core = 1.5E19 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,2],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)
        ne_core = 1.5E19 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,3],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)
        Ti_core = 0.849 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,4],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)
        Te_core = 0.849 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,5],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)    
        
        ## CALCULATE LENGTHS OF SIDES
        lsides = np.zeros((nTri,3))
        for i in range (0,nTri):
            lsides[i,0] = sqrt((ptsx[i,0]-ptsx[i,1])**2 + (ptsy[i,0]-ptsy[i,1])**2)
            lsides[i,1] = sqrt((ptsx[i,1]-ptsx[i,2])**2 + (ptsy[i,1]-ptsy[i,2])**2)
            lsides[i,2] = sqrt((ptsx[i,2]-ptsx[i,0])**2 + (ptsy[i,2]-ptsy[i,0])**2)
        
        ## CALCULATE CELL ANGLES
        angles = np.zeros((nTri,3))
        for i in range (0,nTri):
            p1 = np.array([ptsx[i,0],ptsy[i,0]])
            p2 = np.array([ptsx[i,1],ptsy[i,1]])
            p3 = np.array([ptsx[i,2],ptsy[i,2]])
            angles[i,0] = getangle3ptsdeg(p1,p2,p3)
            angles[i,1] = getangle3ptsdeg(p2,p3,p1)
            angles[i,2] = getangle3ptsdeg(p3,p1,p2)
        
        ## WRITE GTNEUT INPUT FILE!    
        f = open(os.getcwd() + '/toneut_generated','w')
        f.write(' $inp')
        f.write('\n' + ' nCells = ' + str(nTri) + ' nPlasmReg = ' + str(pcellcount) + ' nWallSegm = ' + str(wcellcount))
        for i in range(0,nTri):
            f.write('\n'+'iType( ' + str(i+1) + ' ) = 0 nsides( ' + str(i+1) + ' ) = 3')
            f.write('\n'+'adjCell(1, '+str(i+1)+') = '+str(neighbors[i,0])+' lside(1, '+str(i+1)+ ') = '+str(lsides[i,0])+' angle(1, '+str(i+1)+') = '+str(angles[i,0]))
            f.write('\n'+'adjCell(2, '+str(i+1)+') = '+str(neighbors[i,1])+' lside(2, '+str(i+1)+ ') = '+str(lsides[i,1])+' angle(2, '+str(i+1)+') = '+str(angles[i,1]))
            f.write('\n'+'adjCell(3, '+str(i+1)+') = '+str(neighbors[i,2])+' lside(3, '+str(i+1)+ ') = '+str(lsides[i,2])+' angle(3, '+str(i+1)+') = '+str(angles[i,2]))
        for i,pcell in enumerate(plasmacells):
            f.write('\n'+'iType(' + str(pcell[0]) + ') = 1 nsides(' + str(pcell[0]) + ') = 1 adjCell(1, ' + str(pcell[0]) + ') = ' + str(pcell[1]))
        for i,wcell in enumerate(wallcells):
            f.write('\n'+'iType(' + str(wcell[0]) + ') = 2 nsides(' + str(wcell[0]) + ') = 1 adjCell(1, ' + str(wcell[0]) + ') = ' + str(wcell[1]))
        for i in range(0,nTri):
            f.write('\n'+' elecTemp(' + str(i+1) + ') = ' + str(Te_tri[i]) + ' elecDens(' + str(i+1) + ') = ' + str(ne_tri[i]))
        for i in range(nTri,nTri+pcellcount):
            f.write('\n'+' elecTemp(' + str(i+1) + ') = ' + str(Te_core) + ' elecDens(' + str(i+1) + ') = ' + str(ne_core)) 
        for i in range(0,nTri):
            f.write('\n'+' ionTemp(' + str(i+1) + ') = ' + str(Ti_tri[i]) + ' ionDens(' + str(i+1) + ') = ' + str(ni_tri[i])) 
        for i in range(nTri,nTri+pcellcount):
            f.write('\n'+' ionTemp(' + str(i+1) + ') = ' + str(Ti_core) + ' ionDens(' + str(i+1) + ') = ' + str(ni_core)) 
        for i in range(0,wcellcount):
            f.write('\n'+' g_ex(' + str(i+1) + ') = 1.000000e+19') 
        ## general end of file stuff
        f.write('\n'+'zion =  1    aion =  2    aneut =  2    eneut = 0.002')
        f.write('\n'+'iquad=2')
        f.write('\n'+'icosn=0')
        f.write('\n'+'nph=21')
        f.write('\n'+'idbug =   1')
        f.write('\n'+'scalFact=1.')
        f.write('\n'+'prntOrdr =  -1')
        f.write('\n'+'i_e0=3')
        f.write('\n'+'iatdat=1')
        f.write('\n'+'leh0=1')
        f.write('\n'+'ifrstcol = 1')
        f.write('\n'+'ifjsv=1')
        f.write('\n'+'irefl=1')
        f.write('\n'+'Rwall= '+str(wcellcount)+'*0')
        f.write('\n'+'awall= '+str(wcellcount)+'*12.0')
        f.write('\n'+'zwall= '+str(wcellcount)+'*6.0')
        f.write('\n'+'twall= '+str(wcellcount)+'*2e-3')
        f.write('\n'+'fwabsorb=0.58 13*0.0 15*0.58 13*0.0')
        f.write('\n'+'idp=1')
        f.write('\n'+'isparsitr=5')
        f.write('\n'+'nd0=0')
        f.write('\n'+'neitr = 1')
        f.write('\n'+'$end' + '\n')
        f.close()
        
        f = open(os.getcwd() + '/gtneut_cells','w')
        f.write('\n' + ' nCells = ' + str(nTri) + ' nPlasmReg = ' + str(pcellcount) + ' nWallSegm = ' + str(wcellcount))
        for i in range(0,nTri):
            cell_p = (lsides[i,0] + lsides[i,1] + lsides[i,2])/2
            cell_area = sqrt(cell_p*(cell_p-lsides[i,0])*(cell_p-lsides[i,1])*(cell_p-lsides[i,2]))
            f.write('\n'+'cell: ' + str(i+1) + ' center: ' + str(midpts[i,0]) + ' ' + str(midpts[i,1]) + ' area: ' + str(cell_area) + ' angles: ' + str(angles[i,0]) +' '+ str(angles[i,1]) +' '+ str(angles[i,2]))
        f.close()
        
    def prep_neutpy_input(self,brnd,p):
        ####################################################################
        ## CREATE VERSION OF FIRST WALL THAT INCLUDES FLUX LINE INT. PTS.
        ####################################################################
        ####################################################################
        #this code is also in sol_halo and I think one other place.
        # Should probably write a helper function at some point
        
        all_pts_limx, all_pts_limy = p.lim_line.coords.xy
        all_pts_lim = np.column_stack((all_pts_limx,all_pts_limy))
        pts2add = np.vstack((self.plasma_lim_pts,self.sol_lim_pts))
        for testpt in pts2add:
            for indx,point in enumerate(all_pts_lim):
                if indx <= len(all_pts_lim)-2:
                    pt1 = all_pts_lim[indx]
                    pt2 = all_pts_lim[indx+1]
                    dist_p1_p2 = sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)
                    dist_p1_testpt = sqrt((pt1[0]-testpt[0])**2 + (pt1[1]-testpt[1])**2)
                    dist_testpt_p2 = sqrt((testpt[0]-pt2[0])**2 + (testpt[1]-pt2[1])**2)
                    if abs(dist_p1_p2 - (dist_p1_testpt+dist_testpt_p2)) < 1E-15:
                        all_pts_lim2 = np.insert(all_pts_lim,indx+1,testpt,axis=0)
                        all_pts_lim = all_pts_lim2
                        break    
        lim_line = LineString(all_pts_lim)
        limx, limy = lim_line.coords.xy
        self.lim_pts = np.column_stack((limx, limy))
        self.lim_segs = np.column_stack((
                                    np.linspace(1,self.lim_pts.shape[0],self.lim_pts.shape[0]),
                                    np.roll(np.linspace(1,self.lim_pts.shape[0],self.lim_pts.shape[0]),-1,axis=0)
                                    ))
        
        ####################################################################
        ####################################################################
        ## MAKE ARRAY OF DENSITIES AND TEMPERATURES FOR INTERPOLATION
        ####################################################################
        ####################################################################
        tri_param = np.vstack((self.plasma_param,self.sol_param,self.pfr_param))
    
        ####################################################################
        ####################################################################
        ## COMBINE all_pts, segments, holes, AND attributes FOR TRIANGLE
        ####################################################################
        ####################################################################
    
        #vertices = np.vstack((plasma_pts,sol_pts,lim_pts))
        vertices = np.vstack((self.plasma_pts,self.sol_pts,self.lim_pts,self.divsep_pts))
        vert_number = np.linspace(1,vertices.shape[0],vertices.shape[0])
        vertices = np.column_stack((vert_number,vertices))    
        
        #segments = np.vstack((plasma_segs,sol_segs+plasma_pts.shape[0], lim_segs+plasma_pts.shape[0]+sol_pts.shape[0]))
        segments = np.vstack((
                              self.plasma_segs, 
                              self.sol_segs    +self.plasma_pts.shape[0],
                              self.lim_segs    +self.plasma_pts.shape[0] + self.sol_pts.shape[0],
                              self.divsep_segs +self.plasma_pts.shape[0] + self.sol_pts.shape[0] + self.lim_pts.shape[0]
                              ))
        seg_number = np.linspace(1,segments.shape[0],segments.shape[0])
        seg_attributes = seg_number * 0
        segments = np.column_stack((seg_number,segments,seg_attributes))
        
    
        ####################################################################
        ####################################################################
        ## OUTPUT .poly FILE AND RUN TRIANGLE PROGRAM
        ####################################################################
        ####################################################################
        open('mil_mesh.poly', 'w').close()
        outfile = open('mil_mesh.poly','ab')
        filepath = os.path.realpath(outfile.name)
        np.savetxt(outfile,np.array([vertices.shape[0],2,0,0])[None],fmt='%i %i %i %i')
        np.savetxt(outfile,vertices,fmt='%i %f %f')
        np.savetxt(outfile,np.array([segments.shape[0],0])[None],fmt='%i %i')
        np.savetxt(outfile,segments.astype(int),fmt='%i %i %i %i')
        np.savetxt(outfile,np.array([1])[None],fmt='%i')
        np.savetxt(outfile,np.array([1,p.R0_a,p.Z0])[None],fmt='%i %f %f')
        np.savetxt(outfile,np.array([0])[None],fmt='%i')
        outfile.close()
        call(["triangle", "-pqa0.03n",filepath])
        ####################################################################
        ####################################################################
        ## READ TRIANGLE OUTPUT
        ####################################################################
        ####################################################################

        ## DECLARE FILE PATHS
        nodepath = os.getcwd() + '/mil_mesh.1.node'
        elepath = os.getcwd() + '/mil_mesh.1.ele'
        neighpath = os.getcwd() + '/mil_mesh.1.neigh'

        ## GET NODE DATA
        with open(nodepath, 'r') as node:
            #dummy = next(mil_mesh)
            nodecount = re.findall(r'\d+', next(node))
            nNodes = int(nodecount[0])
            nodenum = np.zeros(nNodes)
            nodesx = np.zeros(nNodes)
            nodesy = np.zeros(nNodes)
            
            for i in range (0,nNodes):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(node)) 
                nodenum[i] = int(data1[0])
                nodesx[i] = data1[1]
                nodesy[i] = data1[2]
                
        ## GET TRIANGLE DATA
        with open(elepath, 'r') as tri_file:
            tricount = re.findall(r'\d+', next(tri_file))
            nTri = int(tricount[0])
            triangles = np.zeros((nTri,3))
            tri_regions = np.zeros(nTri)
            for i in range (0,nTri):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(tri_file))
                triangles[i,0] = data1[1]
                triangles[i,1] = data1[2]
                triangles[i,2] = data1[3]
                #tri_regions[i] = data1[4]
        triangles = triangles.astype('int')
        tri_regions = tri_regions.astype('int')

        ## GET NEIGHBOR DATA
        with open(neighpath, 'r') as neigh_file:
            neighcount = re.findall(r'\d+', next(neigh_file))
            nNeigh = int(neighcount[0])
            neighbors = np.zeros((nNeigh,3))
            for i in range (0,nNeigh):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(neigh_file))
                neighbors[i,0] = data1[1]
                neighbors[i,1] = data1[2]
                neighbors[i,2] = data1[3]
        neighbors = neighbors.astype('int')

        ## REARRANGE TRIANGLES TO CONFORM TO GTNEUT CONVENTION      
        triangles = np.fliplr(triangles) #triangle vertices are given counterclockwise, but we want clockwise
        neighbors = np.fliplr(neighbors) #neighbor 1 is opposite vertex 1, so also counterclockwise 
    
        y=np.zeros(3)
        for i,tri in enumerate(triangles):
            # Find lowest value of y component of vertices
            y[0] = nodesy[tri[0]-1]
            y[1] = nodesy[tri[1]-1]
            y[2] = nodesy[tri[2]-1]
            miny = np.amin(y)
            miny_count = np.sum(y == miny)
            if miny_count == 1:
                #identify position of minimum and roll array accordingly
                miny_index = np.where(y==miny)[0][0]
            else:
                #identify which points are the two minima and determine
                #which of them is farthest to the left (or right if I change it)
                miny_index = np.where(y==miny)[0][1] #change this 1 to a zero to choose the rightmost of the two bottom vertices
            triangles[i] = np.roll(triangles[i],-1*miny_index)
            neighbors[i] = np.roll(neighbors[i],-1*miny_index-2) # the -2 is because the side 1 is opposite vertex 1. We want side 1 to start at vertex 1

        ## GET VALUES TO ORIENT THE FIRST CELL WHEN PLOTTING
        point1_x = nodesx[triangles[0,0]-1]
        point1_y = nodesy[triangles[0,0]-1]
        point2_x = nodesx[triangles[0,1]-1]
        point2_y = nodesy[triangles[0,1]-1]
        point3_x = nodesx[triangles[0,2]-1]
        point3_y = nodesy[triangles[0,2]-1]
        
        cell1_ctr_x = (point1_x + point2_x + point3_x) / 3
        cell1_ctr_y = (point1_y + point2_y + point3_y) / 3

        ## CALCULATE ANGLE BY WHICH TO ROTATE THE FIRST CELL WHEN PLOTTING
        cell1_theta0 = degrees(getangle([point3_x,point3_y],[point1_x,point1_y]))

        ## CALCULATE MID POINTS OF TRIANGLES, AS WELL AS MIDPOINTS FOR EACH FACE     
        ptsx = np.zeros((nTri,3))
        ptsy = np.zeros((nTri,3))
        #for index,tri in ndenumerate(triangles):
        for i in range(0,nTri):
            ptsx[i,0] = nodesx[triangles[i,0]-1]
            ptsy[i,0] = nodesy[triangles[i,0]-1]
            ptsx[i,1] = nodesx[triangles[i,1]-1]
            ptsy[i,1] = nodesy[triangles[i,1]-1]
            ptsx[i,2] = nodesx[triangles[i,2]-1]
            ptsy[i,2] = nodesy[triangles[i,2]-1]
        
        mid_x = np.mean(ptsx,axis=1)
        mid_y = np.mean(ptsy,axis=1)
        self.midpts = np.column_stack((mid_x,mid_y))
        
        side1_midx = (ptsx[:,0] + ptsx[:,1])/2
        side2_midx = (ptsx[:,1] + ptsx[:,2])/2
        side3_midx = (ptsx[:,2] + ptsx[:,0])/2
        
        side1_midy = (ptsy[:,0] + ptsy[:,1])/2
        side2_midy = (ptsy[:,1] + ptsy[:,2])/2
        side3_midy = (ptsy[:,2] + ptsy[:,0])/2
        
        side1_midpt = np.column_stack((side1_midx,side1_midy))
        side2_midpt = np.column_stack((side2_midx,side2_midy))
        side3_midpt = np.column_stack((side3_midx,side3_midy))
        
        #COMBINE POINTS FOR THE PLASMA, SOL, AND DIVERTOR REGIONS
        #first fill in plasma cells
        plasmacells = np.zeros((1,2))
        pcellnum = nTri + 1
        pcellcount = 0
        for index,nei in enumerate(neighbors):
            #for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index],self.plasma_line)
            side2inline = isinline(side2_midpt[index],self.plasma_line)
            side3inline = isinline(side3_midpt[index],self.plasma_line)
            
            if side1inline or side2inline or side3inline:
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                if nb == 1: #cell has one plasma border
                    
                    #create plasma cell
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index + 1
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1
                elif nb == 2: #cell has two plasma borders (this will probably never happen. It would require a local concavity in the inner-most meshed flux surface)
                    #create plasma cell #1
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index + 1
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1   
                    
                    #create plasma cell #2
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index + 1
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1
        plasmacells = np.delete(plasmacells,-1,0)
        plasmacells = plasmacells.astype('int')
        
        #now fill in wall cells
        wallcells = np.zeros((1,2))
        wcellnum = pcellnum #was already advanced in the plasmacell loop. Don't add 1.
        wcellcount = 0
        for index, nei in enumerate(neighbors):
            #for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index],p.lim_line)
            side2inline = isinline(side2_midpt[index],p.lim_line)
            side3inline = isinline(side3_midpt[index],p.lim_line)
            
            if side1inline or side2inline or side3inline:
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                if nb == 1: #cell has one wall border
                    #create wall cell
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index + 1
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1
                elif nb == 2: #cell has two wall borders (This can easily happen because the wall has many concave points.)
                    #create wall cell #1
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index + 1
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1   
                    
                    #create wall cell #2
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index + 1
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1
        wallcells = np.delete(wallcells,-1,0)
        wallcells = wallcells.astype('int')
        
        ## POPULATE CELL DENSITIES AND TEMPERATURES
        #create array of all points in plasma, sol, id, and od
        #tri_param = np.vstack((plasma_param,sol_param,id_param,od_param))
        
        ni_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,2],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)
        ne_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,3],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)
        Ti_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,4],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)
        Te_tri = interpolate.griddata(np.column_stack((tri_param[:,0],tri_param[:,1])), tri_param[:,5],(mid_x, mid_y),method='linear',fill_value=0,rescale=True)

        ni_tri[ni_tri<1.0E16] = 1.0E16
        ne_tri[ne_tri<1.0E16] = 1.0E16
        Ti_tri[Ti_tri<0.002] = 0.002
        Te_tri[Te_tri<0.002] = 0.002
                
        ## CALCULATE DENSITIES AND TEMPERATURES TO ASSIGN TO PLASMA CORE
        ni_core = 1.5E19 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,2],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)
        ne_core = 1.5E19 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,3],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)
        Ti_core = 0.849 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,4],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)
        Te_core = 0.849 #interpolate.griddata(np.column_stack((param[:,0],param[:,1])), param[:,5],(inner_obmp[0], inner_obmp[1]),method='linear',fill_value=0,rescale=True)    
        
        ## CALCULATE LENGTHS OF SIDES
        lsides = np.zeros((nTri,3))
        for i in range (0,nTri):
            lsides[i,0] = sqrt((ptsx[i,0]-ptsx[i,1])**2 + (ptsy[i,0]-ptsy[i,1])**2)
            lsides[i,1] = sqrt((ptsx[i,1]-ptsx[i,2])**2 + (ptsy[i,1]-ptsy[i,2])**2)
            lsides[i,2] = sqrt((ptsx[i,2]-ptsx[i,0])**2 + (ptsy[i,2]-ptsy[i,0])**2)
        
        ## CALCULATE CELL ANGLES
        angles = np.zeros((nTri,3))
        for i in range (0,nTri):
            p1 = np.array([ptsx[i,0],ptsy[i,0]])
            p2 = np.array([ptsx[i,1],ptsy[i,1]])
            p3 = np.array([ptsx[i,2],ptsy[i,2]])
            angles[i,0] = getangle3ptsdeg(p1,p2,p3)
            angles[i,1] = getangle3ptsdeg(p2,p3,p1)
            angles[i,2] = getangle3ptsdeg(p3,p1,p2)
        
        ## WRITE NEUTPY INPUT FILE!    
        f = open(os.getcwd() + '/neutpy_in_generated','w')
        f.write('nCells = ' + str(nTri) + ' nPlasmReg = ' + str(pcellcount) + ' nWallSegm = ' + str(wcellcount))
        for i in range(0,nTri):
            f.write('\n'+'iType(' + str(i) + ') = 0 nSides(' + str(i) + ') = 3 ' + 'adjCell('+str(i)+') = '+', '.join(map(str, neighbors[i,:]-1)))
        f.write('\n')
        f.write('\n#lsides and angles for normal cells')
        for i in range(0,nTri):
            f.write('\n'+'lsides(' + str(i) + ') = '+', '.join(map(str, lsides[i,:]))+' angles(' + str(i) + ') = '+', '.join(map(str, angles[i,:])))
        f.write('\n')
        f.write('\n#densities and temperatures for normal cells')
        for i in range(0,nTri):
            f.write('\n'+'elecTemp('+str(i)+') = '+str(Te_tri[i]) +' elecDens(' + str(i) + ') = '+str(ne_tri[i])+' ionTemp('+str(i)+') = '+str(Ti_tri[i]) +' ionDens(' + str(i) + ') = '+str(ni_tri[i]))
        f.write('\n')
        f.write('\n#wall cells')
        for i,wcell in enumerate(wallcells):
            f.write('\n'+'iType('+str(wcell[0]-1)+') = 2 nSides('+str(wcell[0]-1)+') = 1 adjCell('+str(wcell[0]-1)+') = '+str(wcell[1]-1)+' zwall('+str(wcell[0]-1)+') = 6 awall('+str(wcell[0]-1)+') = 12 twall('+str(wcell[0]-1)+') = 0.002 f_abs('+str(wcell[0]-1)+') = 0.0 s_ext('+str(wcell[0]-1)+') = 1e+19') 
        f.write('\n')
        f.write('\n#plasma core and vacuum cells')
        for i,pcell in enumerate(plasmacells):
            f.write('\n'+'iType(' + str(pcell[0]-1) + ') = 1 nSides(' + str(pcell[0]-1) + ') = 1 adjCell(1, ' + str(pcell[0]-1) + ') = ' + str(pcell[1]-1) + ' twall(' + str(pcell[0]-1) + ') = 5000  alb_s(' + str(pcell[0]-1) + ') = 0  alb_t(' + str(pcell[0]-1) + ') = 0  s_ext(' + str(pcell[0]-1) + ') = 0 ')
        f.write('\n')
        f.write('\n#general parameters')
        f.write('\nzion = 1 ')
        f.write('\naion = 2 ')
        f.write('\naneut = 2 ')
        f.write('\ntslow = 0.002 ')
        f.write('\n')
        f.write('\n#cross section and reflection model parameters')
        f.write('\nxsec_ioni = janev')
        f.write('\nxsec_ione = janev')
        f.write('\nxsec_cx = janev')
        f.write('\nxsec_el = janev')
        f.write('\nxsec_eln = stacey_thomas')
        f.write('\nxsec_rec = stacey_thomas')
        f.write('\nrefmod_e = stacey')
        f.write('\nrefmod_n = stacey')
        f.write('\n')
        f.write('\n#transmission coefficient parameters')
        f.write('\nint_method = midpoint')
        f.write('\nphi_int_pts = 10')
        f.write('\nxi_int_pts = 10')
        f.write('\n')
        f.write('\n#make a bickley-naylor interpolated lookup file. (y or n)')
        f.write('\nmake_bn_int = n')
        f.write('\n')
        f.write('\n#extra (optional) arguments for plotting')
        f.write('\ncell1_ctr_x  = ' + str(cell1_ctr_x))
        f.write('\ncell1_ctr_y  = ' + str(cell1_ctr_y))
        f.write('\ncell1_theta0 = ' + str(cell1_theta0))
        f.write('\n')
        f.close()

        
        #create dictionary to pass to neutpy
        toneutpy={}
        toneutpy["nCells"]       = nTri
        toneutpy["nPlasmReg"]    = pcellcount
        toneutpy["nWallSegm"]    = wcellcount
        toneutpy["aneut"]        = 2
        toneutpy["zion"]         = 1
        toneutpy["aion"]         = 2
        toneutpy["tslow"]        = 0.002
        toneutpy["int_method"]   = 'midpoint'
        toneutpy["phi_int_pts"]  = 10
        toneutpy["xi_int_pts"]   = 10
        toneutpy["xsec_ioni"]    = 'janev'
        toneutpy["xsec_ione"]    = 'janev'
        toneutpy["xsec_cx"]      = 'janev'
        toneutpy["xsec_rec"]     = 'stacey_thomas'
        toneutpy["xsec_el"]      = 'janev'
        toneutpy["xsec_eln"]     = 'stacey_thomas'
        toneutpy["refmod_e"]     = 'stacey'
        toneutpy["refmod_n"]     = 'stacey'
        
        toneutpy["iType"]        = np.asarray([0]*nTri + [1]*pcellcount + [2]*wcellcount)
        toneutpy["nSides"]       = np.asarray([3]*nTri + [1]*(pcellcount + wcellcount))
        toneutpy["zwall"]        = np.asarray([0]*(nTri+pcellcount) + [6]*wcellcount)
        toneutpy["awall"]        = np.asarray([0]*(nTri+pcellcount) + [12]*wcellcount)
        toneutpy["elecTemp"]     = Te_tri[:nTri]
        toneutpy["ionTemp"]      = Ti_tri[:nTri]
        toneutpy["elecDens"]     = ne_tri[:nTri]
        toneutpy["ionDens"]      = ni_tri[:nTri]
        toneutpy["twall"]        = np.asarray([0]*nTri + [5000]*pcellcount + [0.002]*wcellcount)
        toneutpy["f_abs"]        = np.asarray([0]*(nTri+pcellcount) + [0]*wcellcount)
        toneutpy["alb_s"]        = np.asarray([0]*nTri + [0]*pcellcount + [0]*wcellcount)
        toneutpy["alb_t"]        = np.asarray([0]*nTri + [0]*pcellcount + [0]*wcellcount)
        toneutpy["s_ext"]        = np.asarray([0]*nTri + [0]*pcellcount + [1.0E19]*wcellcount)
        
        toneutpy["adjCell"]      = neighbors-1
        toneutpy["lsides"]       = lsides
        toneutpy["angles"]       = angles
        toneutpy["cell1_ctr_x"]  = cell1_ctr_x
        toneutpy["cell1_ctr_y"]  = cell1_ctr_y
        toneutpy["cell1_theta0"] = cell1_theta0
        
        self.neutpy_inst = neutpy(inarrs=toneutpy)
        plot = neutpyplot(self.neutpy_inst)
        
    def ntrl2grid(self,inp,brnd):
        """
        
        """
        self.nn = interpolate.griddata(self.midpts, 
                                       self.neutpy_inst.cell_nn,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape)
        self.izn_rate = interpolate.griddata(self.midpts, 
                                       self.neutpy_inst.cell_izn_rate,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape)        
        #DOING A FLUX SURFACE AVERAGE OF THE NEUTRALS DATA.
        #FTR, I DON'T ACTUALLY THINK THIS IS A GOOD WAY TO HANDLE THIS. - MH
        nn_1D_col = np.sum(self.nn * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        self.nn_1D = np.repeat(nn_1D_col.reshape(-1,1),inp.thetapts,axis=1)

        izn_rate_1D_col = np.sum(self.izn_rate * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        self.izn_rate_1D = np.repeat(izn_rate_1D_col.reshape(-1,1),inp.thetapts,axis=1)
        
        nn_plot = plt.figure(figsize=(6,4))
        ax1 = nn_plot.add_subplot(1,1,1)
        ax1.set_title('FSA Neutral Densities')
        ax1.set_xlabel(r'$\rho$')
        ax1.set_ylabel(r'$m^{-3}$')
        ax1.plot(brnd.rho[:,0],self.nn_1D[:,0],label='neutral densities')
        ax1.legend()
        
        izn_rate_plot = plt.figure(figsize=(6,4))
        ax1 = izn_rate_plot.add_subplot(1,1,1)
        ax1.set_title('FSA Ionization Rate')
        ax1.set_xlabel(r'$\rho$')
        ax1.set_ylabel(r'???')
        ax1.plot(brnd.rho[:,0],self.izn_rate_1D[:,0],label='izn rate')
        ax1.legend()
