#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:05:53 2018

@author: max
"""

class mil_pfr_brnd():
    def __init__(self):
        pass
    
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
