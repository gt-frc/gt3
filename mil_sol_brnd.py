#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:03:04 2018

@author: max
"""

class mil_sol_brnd():
    def __init__(self,inp,brnd):
        self.solmesh(inp,brnd)
        self.solnT(inp,brnd)
    
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