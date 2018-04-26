# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:54:44 2017

@author: max
"""
from __future__ import division
import numpy as np
from helpers import getangle, getpoint, cut, isinline, getangle3ptsdeg,makeline
from scipy import interpolate
from math import pi, cos, sin, sqrt, degrees
from shapely.ops import polygonize
from shapely.geometry import LineString, Polygon, Point, MultiPoint, LinearRing
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
        in_strike
        out_strike
        sep_line
        sol_flx_line
        sol_segs
        divsep_pts
        divsep_segs
        plasma_pts
        plasma_segs
        plasma_lim_pts
        core_line
        plasma_param
        
    """
    def __init__(self,inp,brnd):
        sys.dont_write_bytecode = True
        
        #check if specified neutpy_outfile exists. If so, read in and skip everything else.
        outfile_found=0
        try:
            for root, subdirs, files in os.walk(os.getcwd()):
                for filename in files:
                    if filename == inp.neut_outfile:
                        outfile_found = 1
                        os.path.join(root,filename)
                        self.neutfile_loc = os.path.join(root,filename)
            if outfile_found==0:
                self.neutfile_loc = os.getcwd() + '/' + inp.neut_outfile
                print
                print '#########################################################'
                print 'Neutrals file',inp.neut_outfile,'not found.'
                print 'Proceeding with neutrals calculation. Results'
                print 'will be saved as',inp.neut_outfile,'in the main directory'
                print '#########################################################'
                print
        except AttributeError:
            #neut_outfile wasn't specified in input file.
            #Assign default value of neut_outfile.dat and print notice 
            #so user can decide whether to modify input file
            self.neutfile_loc = os.getcwd() + '/neut_outfile.dat'
            print
            print '#########################################################'
            print 'Neutrals output file not specified in main input file.'
            print 'GT3 will generate one and save as neut_outfile.dat in'
            print 'the main directory. To save time in the future, include'
            print 'the following line somewhere in the input file. Note that'
            print 'you can modify the filename; it just has to be specified.'
            print 
            print 'neut_outfile = neut_outfile.dat'
            print
            print '#########################################################'
            print
        
        if outfile_found==1:
            self.read_outfile(inp,brnd)
            self.ntrl2grid(inp,brnd)
        else:
            self.plasmamesh(inp,brnd)
            self.solmesh(inp,brnd)
            self.solnT(inp,brnd)
            self.pfrnT(inp,brnd)
            self.prep_neutpy_input(inp,brnd)
            self.ntrl2grid(inp,brnd)

    def read_outfile(self,inp,brnd):
        neutdata = np.loadtxt(self.neutfile_loc,skiprows=1)
        self.midpts = np.column_stack((neutdata[:,0],neutdata[:,1]))
        self.nn_s_raw       = neutdata[:,2]
        self.nn_t_raw       = neutdata[:,3]
        self.nn_raw         = neutdata[:,4]
        self.iznrate_s_raw  = neutdata[:,5]
        self.iznrate_t_raw  = neutdata[:,6]
        self.iznrate_raw    = neutdata[:,7]

    def plasmamesh(self,inp,brnd):
        #IF PARAMETERS ARE SPECIFIED IN INPUT, CREATE NEW R, Z, ni, etc. ARRAYS
        #FOR NEUTRALS CALCULATION
        
        r_interp       = interpolate.interp2d(brnd.rho,brnd.theta,brnd.r,kind='linear')
        rho_interp     = interpolate.interp2d(brnd.rho,brnd.theta,brnd.rho,kind='linear')
        theta_interp   = interpolate.interp2d(brnd.rho,brnd.theta,brnd.theta,kind='linear')
        R_interp       = interpolate.interp2d(brnd.rho,brnd.theta,brnd.R,kind='linear')
        Z_interp       = interpolate.interp2d(brnd.rho,brnd.theta,brnd.Z,kind='linear')
        ni_interp      = interpolate.interp2d(brnd.rho,brnd.theta,brnd.ni,kind='linear')
        ne_interp      = interpolate.interp2d(brnd.rho,brnd.theta,brnd.ne,kind='linear')
        Ti_kev_interp  = interpolate.interp2d(brnd.rho,brnd.theta,brnd.Ti_kev,kind='linear')
        Te_kev_interp  = interpolate.interp2d(brnd.rho,brnd.theta,brnd.Te_kev,kind='linear')
        
        r      = r_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        rho    = rho_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        theta  = theta_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        R      = R_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        Z      = Z_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        ni     = ni_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        ne     = ne_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        Ti_kev = Ti_kev_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T
        Te_kev = Te_kev_interp(np.linspace(inp.ntrl_rho_start,1,inp.ntrl_rpts),
                                    np.linspace(0,2*pi,inp.ntrl_thetapts)).T

        # GET POINTS FOR TRIANGULATION
        self.plasma_pts = np.column_stack((R[:,:-1].flatten(),Z[:,:-1].flatten()))
        
        self.plasma_segs = np.zeros((0,2))
        #for i in range( (ntrl_rho[:,0] >= radmesh_start).sum() ):
        for i,v in enumerate(rho[:,0]):
            new_segs = np.column_stack((
                                        np.arange(inp.ntrl_thetapts-1),
                                        np.roll(np.arange(inp.ntrl_thetapts-1),-1)
                                        )) + (inp.ntrl_thetapts-1) * i
            self.plasma_segs = np.vstack((self.plasma_segs,new_segs))
        
        #calculate location of the x-point in plasma_pts, in case it's necessary
        self.xpt_loc = len(self.plasma_pts) - (inp.ntrl_thetapts-1)/4

        #shapely LinearRing used later to identify which sides of 
        #which cells border on core plasma region
        self.core_line = LineString(self.plasma_pts[:inp.ntrl_thetapts-1])
        self.core_ring = LinearRing(self.plasma_pts[:inp.ntrl_thetapts-1])
        
        #shapely LinearRing of seperatrix. Used later.
        self.sep_line = LineString(self.plasma_pts[len(self.plasma_pts)-inp.ntrl_thetapts+1:])
        self.sep_ring = LinearRing(self.plasma_pts[len(self.plasma_pts)-inp.ntrl_thetapts+1:])

        #parameters get combined with other parameters later for one 
        #big interpolation over the plasma chamber
        self.plasma_param = np.column_stack((
                                            R.flatten(),
                                            Z.flatten(),
                                            ni.flatten(),
                                            ne.flatten(),
                                            Ti_kev.flatten(),
                                            Te_kev.flatten()
                                            ))
        
        #also create the lim_pts array. This might get moved somewhere else
        #TODO: figure out where to put this.
        
    def solmesh(self,inp,brnd):
        #for now, all this does is calculate the strike point locations
        self.in_strike  = np.asarray(makeline(inp.xpt,10.0,inp.xtheta1)[1].intersection(inp.lim_line).xy)
        self.out_strike = np.asarray(makeline(inp.xpt,10.0,inp.xtheta4)[1].intersection(inp.lim_line).xy)
        
        #add inner strike point
        union = inp.lim_line.union(makeline(inp.xpt,10.0,inp.xtheta1)[1])
        result = [geom for geom in polygonize(union)][0]
        
        #add outer strike point
        union = result.union(makeline(inp.xpt,10.0,inp.xtheta4)[1])
        result = [geom for geom in polygonize(union)][0]
        
        self.wall_x, self.wall_y = result.exterior.coords.xy
        
        #delete the repeated point. It causes problems for FiPy
        self.wall_x = np.delete(self.wall_x,-1)
        self.wall_y = np.delete(self.wall_y,-1)

        self.wall_pts = np.column_stack((self.wall_x,self.wall_y))
        self.wall_line = LineString(self.wall_pts)
        self.wall_ring = LinearRing(self.wall_pts)
        #create wall segments for triangulation later
        self.wall_segs = np.column_stack((
                                        np.arange(len( self.wall_pts )),
                                        np.roll(np.arange(len( self.wall_pts )),-1)
                                        ))

    def solnT(self,inp,brnd):         
        #####################################################################       
        ## POPULATE DENSITY AND TEMPERATURE VALUES FOR THE SOL AND HALO REGIONS  
        ##################################################################### 
    
        ## THESE WILL GET COMBINED WITH AN UNSTRUCTURED LIST OF DENSITIES AND TEMPERATURES
        ## AT EVERY OTHER VERTEX IN THE PROBLEM AND USED IN A 2D INTERPOLATION
        ## TO GET THE VALUES AT THE MIDDLES OF THE TRIANGLES GENERATED FOR GTNEUT
        
        
        # FIND POSITION OF XPOINT IN SEP_LINE COORDINATES AND ROTATE TO START GOING
        # AROUND THE SEPERATRIX COUNTER-CLOCKWISE
    
        sepx, sepy = self.sep_line.coords.xy
        xpt_pos = int(len(sepx) - (inp.ntrl_thetapts-1)/4)
        self.sep_pts = np.roll(np.column_stack((sepx,sepy)),-xpt_pos,axis=0)

        #make new sep_line with x-point repeated for use in the bdry condition
        sep_line2 = LineString(np.vstack((self.sep_pts,self.sep_pts[0])))
        
        #rotate wall points to start at inner strike point
        in_str_loc = np.where((self.wall_pts[:,0] == self.in_strike[0]) & 
                             (self.wall_pts[:,1] == self.in_strike[1]))[0][0]
        wall_pts_rot = np.roll(self.wall_pts,-in_str_loc,axis=0)

        #make new wall_line for use in the bdry condition
        wall_line = LineString(wall_pts_rot)
        
        #cut the new wall line at the outboard strike point
        dist = wall_line.project(Point(self.out_strike[0],self.out_strike[1]),normalized=True)
        wall_line_cut = cut(wall_line,dist)[0]
        
        #get the final wall points for creating the FiPy input string
        wallx, wally = wall_line_cut.coords.xy
        wall_pts = np.column_stack((wallx,wally))
        
        #start making mesh. Start with points
        mesh_info = ""
        pt_count = 1
        for i, (xval, yval) in enumerate(zip(self.sep_pts[:,0],self.sep_pts[:,1])):
            mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(xval) + ', ' + str(yval) + ', 0, 1/10};')
            pt_count+=1
            
        for i, (xval, yval) in enumerate(zip(wall_pts[:,0],wall_pts[:,1])):
            mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(xval) + ', ' + str(yval) + ', 0, 1/10};')
            pt_count+=1

        #define lines for the mesh. Start with the seperatrix
        line_count = 1
        for i in np.linspace(1,len(sepx)-1,len(sepx)-1):
            mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(i)) + ', ' + str(int(i+1)) + '};')
            line_count+=1

        #return to the xpoint
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(len(sepx))) + ', ' + str(int(1)) + '};')
        line_count+=1
        
        #line from xpt to inboard side strike point
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(1)) + ', ' + str(int(len(sepx)+1)) + '};')
        line_count+=1
        
        #lines around the first wall to the outboard strike point
        for i in np.linspace(int(len(sepx)+1),int(len(sepx)+1) + len(wallx)-2,len(wallx)-1):
            mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(i)) + ', ' + str(int(i+1)) + '};')
            line_count+=1
        
        #return to the xpoint on the outboard side
        mesh_info = mesh_info + ('Line(' + str(line_count) + ') = {' + str(int(line_count-1)) + ', ' + str(int(1)) + '};')
       
        #create lineloop and plane surface
        lineloop = "Line Loop(1) = {1"
        for i in np.linspace(2,line_count,line_count-1):
            lineloop = lineloop + ", " + str(int(i))
        lineloop = lineloop + "};"
        mesh_info = mesh_info + lineloop
        mesh_info = mesh_info + "Plane Surface(1) = {1};\n"

        # Mesh info complete. Create mesh object.
        m = Gmsh2D(mesh_info)
        
        ####################################################################
        ## CREATE LINES FOR BOUNDARY CONDITIONS
        ####################################################################
        #already made lines for sep and first wall
        #now make them for inner and outer divertor legs
        idiv_line = LineString([self.sep_pts[0],self.in_strike])
        odiv_line = LineString([self.sep_pts[0],self.out_strike])
        
        var = CellVariable(mesh=m)

        ####################################################################
        ## LOOP OVER QUANTITIES TO CALCULATE (ni, ne, Ti, Te)
        ####################################################################
        self.sol_param = np.column_stack((m.cellCenters.value[0],m.cellCenters.value[1]))
        for j in [0,1,2,3]:
            mask_array  = np.zeros(m.numberOfCells)
            bcval       = np.zeros(m.numberOfCells)
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
                ptinwall = wall_line_cut.distance(p1) < 1e-8 or wall_line_cut.distance(p2) < 1e-8 or wall_line_cut.distance(p3) < 1e-8
                
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
        nibc = np.asarray(wallx)*0 + brnd.ni[-1,0]
        nebc = np.asarray(wallx)*0 + brnd.ne[-1,0]
        Tibc = np.asarray(wallx)*0 + brnd.Ti_kev[-1,0]
        Tebc = np.asarray(wallx)*0 + brnd.Te_kev[-1,0]
        
        solbc_array = np.column_stack((wallx,wally,nibc,nebc,Tibc,Tebc))
        self.sol_param = np.vstack((self.sol_param,solbc_array))
        
            
        viewer = Viewer(vars=var, datamin=0., datamax=np.amax(var_final))
        viewer.plot()
        plt.show()
        
    def pfrnT(self,inp,brnd):         
        #start making mesh
        pt_count = 1
        mesh_info = ""
        mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(self.sep_pts[0,0]) + ', ' + str(self.sep_pts[0,1]) + ', 0, 1/20};')
        pt_count+=1
        mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(self.in_strike[0][0]) + ', ' + str(self.in_strike[1][0]) + ', 0, 1/30};')
        pt_count+=1
        mesh_info = mesh_info + ('Point(' + str(pt_count) + ') = {' + str(self.out_strike[0][0]) + ', ' + str(self.out_strike[1][0]) + ', 0, 1/30};')
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
        lineloop = "Line Loop(1) = {1, 2, 3};"
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
        idiv_line = LineString([(inp.xpt[0],inp.xpt[1]),(self.in_strike[0], self.in_strike[1])])
        odiv_line = LineString([(inp.xpt[0],inp.xpt[1]),(self.out_strike[0], self.out_strike[1])])
        bottom_line = LineString([(self.in_strike[0], self.in_strike[1]),(self.out_strike[0], self.out_strike[1])])
        
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

        
    def prep_neutpy_input(self,inp,brnd):

        # MAKE ARRAY OF DENSITIES AND TEMPERATURES FOR INTERPOLATION
        tri_param = np.vstack((self.plasma_param,self.sol_param,self.pfr_param))
    
        # COMBINE all_pts, segments, holes, AND attributes FOR TRIANGLE
        vertices = np.vstack((self.plasma_pts,self.wall_pts))
        vert_number = np.arange(len(vertices))
        vertices = np.column_stack((vert_number,vertices))    
        
        segments = np.vstack((
                              self.plasma_segs, 
                              self.wall_segs + self.plasma_pts.shape[0],
                              ))

        seg_number = np.linspace(1,segments.shape[0],segments.shape[0])
        seg_attributes = seg_number * 0
        segments = np.column_stack((seg_number,segments,seg_attributes))
    
        ## OUTPUT .poly FILE AND RUN TRIANGLE PROGRAM
        open('mil_mesh.poly', 'w').close()
        outfile = open('mil_mesh.poly','ab')
        filepath = os.path.realpath(outfile.name)
        np.savetxt(outfile,np.array([vertices.shape[0],2,0,0])[None],fmt='%i %i %i %i')
        np.savetxt(outfile,vertices,fmt='%i %f %f')
        np.savetxt(outfile,np.array([segments.shape[0],0])[None],fmt='%i %i')
        np.savetxt(outfile,segments.astype(int),fmt='%i %i %i %i')
        np.savetxt(outfile,np.array([1])[None],fmt='%i')
        np.savetxt(outfile,np.array([1,inp.R0_a,inp.Z0])[None],fmt='%i %f %f')
        np.savetxt(outfile,np.array([0])[None],fmt='%i')
        outfile.close()
        
        #construct options to pass to triangle, as specified in input file
        #refer to https://www.cs.cmu.edu/~quake/triangle.html
        tri_options = '-p'
        try:
            tri_options = tri_options + 'q' + inp.tri_min_angle
        except:
            pass
        
        try:
            tri_options = tri_options + 'a' + inp.tri_min_area
        except:
            pass
        
        tri_options = tri_options + 'nz'
        
        #call triangle
        try:
            call([inp.triangle_loc+'triangle', tri_options,filepath])
        except AttributeError:
            try:
                call(['triangle', tri_options,filepath])
            except:
                print 'triangle could not be found. Stopping.'
                sys.exit

        ## READ TRIANGLE OUTPUT

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
            y[0] = nodesy[tri[0]]
            y[1] = nodesy[tri[1]]
            y[2] = nodesy[tri[2]]
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
        point1_x = nodesx[triangles[0,0]]
        point1_y = nodesy[triangles[0,0]]
        point2_x = nodesx[triangles[0,1]]
        point2_y = nodesy[triangles[0,1]]
        point3_x = nodesx[triangles[0,2]]
        point3_y = nodesy[triangles[0,2]]
        
        cell1_ctr_x = (point1_x + point2_x + point3_x) / 3
        cell1_ctr_y = (point1_y + point2_y + point3_y) / 3

        ## CALCULATE ANGLE BY WHICH TO ROTATE THE FIRST CELL WHEN PLOTTING
        cell1_theta0 = degrees(getangle([point3_x,point3_y],[point1_x,point1_y]))    

        ## GET VALUES TO ORIENT THE FIRST CELL WHEN PLOTTING
        point1_x = nodesx[triangles[0,0]]
        point1_y = nodesy[triangles[0,0]]
        point2_x = nodesx[triangles[0,1]]
        point2_y = nodesy[triangles[0,1]]
        point3_x = nodesx[triangles[0,2]]
        point3_y = nodesy[triangles[0,2]]
        
        cell1_ctr_x = (point1_x + point2_x + point3_x) / 3
        cell1_ctr_y = (point1_y + point2_y + point3_y) / 3

        ## CALCULATE ANGLE BY WHICH TO ROTATE THE FIRST CELL WHEN PLOTTING
        cell1_theta0 = degrees(getangle([point3_x,point3_y],[point1_x,point1_y]))

        ## CALCULATE MID POINTS OF TRIANGLES, AS WELL AS MIDPOINTS FOR EACH FACE     
        ptsx = np.zeros((nTri,3))
        ptsy = np.zeros((nTri,3))
        #for index,tri in ndenumerate(triangles):
        for i in range(0,nTri):
            ptsx[i,0] = nodesx[triangles[i,0]]
            ptsy[i,0] = nodesy[triangles[i,0]]
            ptsx[i,1] = nodesx[triangles[i,1]]
            ptsy[i,1] = nodesy[triangles[i,1]]
            ptsx[i,2] = nodesx[triangles[i,2]]
            ptsy[i,2] = nodesy[triangles[i,2]]
        
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
        pcellnum = nTri
        pcellcount = 0
        for index,nei in enumerate(neighbors):
            #for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index],self.core_ring)
            side2inline = isinline(side2_midpt[index],self.core_ring)
            side3inline = isinline(side3_midpt[index],self.core_ring)
            
            if side1inline or side2inline or side3inline:
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                if nb == 1: #cell has one plasma border
                    
                    #create plasma cell
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1
                elif nb == 2: #cell has two plasma borders (this will probably never happen. It would require a local concavity in the inner-most meshed flux surface)
                    #create plasma cell #1
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index
                    plasmacells = np.vstack((plasmacells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1   
                    
                    #create plasma cell #2
                    plasmacells[pcellcount,0] = pcellnum
                    plasmacells[pcellcount,1] = index
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
            side1inline = isinline(side1_midpt[index],self.wall_ring)
            side2inline = isinline(side2_midpt[index],self.wall_ring)
            side3inline = isinline(side3_midpt[index],self.wall_ring)
            
            if side1inline or side2inline or side3inline:
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                if nb == 1: #cell has one wall border
                    #create wall cell
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1
                elif nb == 2: #cell has two wall borders (This can easily happen because the wall has many concave points.)
                    #create wall cell #1
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index
                    wallcells = np.vstack((wallcells,[0,0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1   
                    
                    #create wall cell #2
                    wallcells[wcellcount,0] = wcellnum
                    wallcells[wcellcount,1] = index
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
            f.write('\n'+'iType(' + str(i) + ') = 0 nSides(' + str(i) + ') = 3 ' + 'adjCell('+str(i)+') = '+', '.join(map(str, neighbors[i,:])))
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
            f.write('\n'+'iType('+str(wcell[0])+') = 2 nSides('+str(wcell[0])+') = 1 adjCell('+str(wcell[0])+') = '+str(wcell[1])+' zwall('+str(wcell[0])+') = 6 awall('+str(wcell[0])+') = 12 twall('+str(wcell[0])+') = 0.002 f_abs('+str(wcell[0])+') = 0.0 s_ext('+str(wcell[0])+') = 1e+19') 
        f.write('\n')
        f.write('\n#plasma core and vacuum cells')
        for i,pcell in enumerate(plasmacells):
            f.write('\n'+'iType(' + str(pcell[0]) + ') = 1 nSides(' + str(pcell[0]) + ') = 1 adjCell(1, ' + str(pcell[0]) + ') = ' + str(pcell[1]) + ' twall(' + str(pcell[0]) + ') = 5000  alb_s(' + str(pcell[0]) + ') = 0  alb_t(' + str(pcell[0]) + ') = 0  s_ext(' + str(pcell[0]) + ') = 0 ')
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
        
        toneutpy["adjCell"]      = neighbors
        toneutpy["lsides"]       = lsides
        toneutpy["angles"]       = angles
        toneutpy["cell1_ctr_x"]  = cell1_ctr_x
        toneutpy["cell1_ctr_y"]  = cell1_ctr_y
        toneutpy["cell1_theta0"] = cell1_theta0
        
        self.neutpy_inst = neutpy(inarrs=toneutpy)
        plot = neutpyplot(self.neutpy_inst)
        self.nn_s_raw = self.neutpy_inst.cell_nn_s
        self.nn_t_raw = self.neutpy_inst.cell_nn_t
        self.nn_raw = self.nn_s_raw + self.nn_t_raw
        
        self.iznrate_s_raw = self.neutpy_inst.cell_izn_rate_s
        self.iznrate_t_raw = self.neutpy_inst.cell_izn_rate_t
        self.iznrate_raw = self.iznrate_s_raw + self.iznrate_t_raw
        
        #create output file
        #the file contains R,Z coordinates and then the values of several calculated parameters
        #at each of those points.
        
        f = open(self.neutfile_loc,'w')
        f.write(('{:^18s}'*8).format('R','Z','n_n_slow','n_n_thermal','n_n_total','izn_rate_slow','izn_rate_thermal','izn_rate_total'))
        for i,pt in enumerate(self.midpts):
            f.write(('\n'+'{:>18.5f}'*2+'{:>18.5E}'*6).format(
                                        self.midpts[i,0],
                                        self.midpts[i,1],
                                        self.nn_s_raw[i],
                                        self.nn_t_raw[i],
                                        self.nn_raw[i],
                                        self.iznrate_s_raw[i],
                                        self.iznrate_t_raw[i],
                                        self.iznrate_raw[i]))
        f.close()
        


    def ntrl2grid(self,inp,brnd):
        """
        
        """
        self.nn_s = interpolate.griddata(self.midpts, 
                                       self.nn_s_raw,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape)
        self.nn_t = interpolate.griddata(self.midpts, 
                                       self.nn_t_raw,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape)
        self.nn = interpolate.griddata(self.midpts, 
                                       self.nn_raw,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape)
        self.izn_rate_s = interpolate.griddata(self.midpts, 
                                       self.iznrate_s_raw,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape)     
        self.izn_rate_t = interpolate.griddata(self.midpts, 
                                       self.iznrate_t_raw,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape) 
        self.izn_rate = interpolate.griddata(self.midpts, 
                                       self.iznrate_raw,
                                       np.column_stack((brnd.R.flatten(),brnd.Z.flatten()))
                                       ).reshape(brnd.rho.shape) 

        #DOING A FLUX SURFACE AVERAGE OF THE NEUTRALS DATA.
        #FTR, I DON'T ACTUALLY THINK THIS IS A GOOD WAY TO HANDLE THIS. - MH
        nn_s_1D_col  = np.sum(self.nn_s * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        nn_t_1D_col  = np.sum(self.nn_t * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        nn_1D_col    = np.sum(self.nn   * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        self.nn_s_1D = np.repeat(nn_s_1D_col.reshape(-1,1),inp.thetapts,axis=1)
        self.nn_t_1D = np.repeat(nn_t_1D_col.reshape(-1,1),inp.thetapts,axis=1)
        self.nn_1D   = np.repeat(nn_1D_col.reshape(-1,1),inp.thetapts,axis=1)

        izn_rate_s_1D_col  = np.sum(self.izn_rate_s * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        izn_rate_t_1D_col  = np.sum(self.izn_rate_t * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        izn_rate_1D_col    = np.sum(self.izn_rate   * (brnd.L_seg*brnd.R),axis=1) / (brnd.L_r[:,0]*brnd.R0[:,0])
        self.izn_rate_s_1D = np.repeat(izn_rate_s_1D_col.reshape(-1,1),inp.thetapts,axis=1)
        self.izn_rate_t_1D = np.repeat(izn_rate_t_1D_col.reshape(-1,1),inp.thetapts,axis=1)
        self.izn_rate_1D   = np.repeat(izn_rate_1D_col.reshape(-1,1),inp.thetapts,axis=1)
        
        #nn_plot = plt.figure(figsize=(6,4))
        #ax1 = nn_plot.add_subplot(1,1,1)
        #ax1.set_title('FSA Neutral Densities')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'$m^{-3}$')
        #ax1.plot(brnd.rho[:,0],self.nn_1D[:,0],label='neutral densities')
        #ax1.legend()
        
        #izn_rate_plot = plt.figure(figsize=(6,4))
        #ax1 = izn_rate_plot.add_subplot(1,1,1)
        #ax1.set_title('FSA Ionization Rate')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.izn_rate_1D[:,0],label='izn rate')
        #ax1.legend()
