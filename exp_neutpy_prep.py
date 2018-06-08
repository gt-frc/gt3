#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 14:27:42 2018

@author: max
"""
from __future__ import division
import numpy as np
from shapely.geometry import LinearRing, Point, LineString
from neutpy import neutpy, neutpyplot
from math import degrees, sqrt, acos, pi
from scipy.interpolate import griddata
from collections import namedtuple
import sys
import os
import re
import time
from subprocess import call
import matplotlib.pyplot as plt

def getangle(p1, p2):
    if isinstance(p1, Point) and isinstance(p2, Point):
        p1 = [p1.coords.xy[0][0], p1.coords.xy[1][0]]
        p2 = [p2.coords.xy[0][0], p2.coords.xy[1][0]]
    p1 = np.asarray(p1)
    p1 = np.reshape(p1, (-1, 2))
    p2 = np.asarray(p2)
    p2 = np.reshape(p2, (-1, 2))
    theta = np.arctan2(p1[:, 1]-p2[:, 1], p1[:, 0]-p2[:, 0])
    theta_mod = np.where(theta<0, theta+pi, theta) #makes it so the angle is always measured counterclockwise from the horizontal
    return theta

def getangle3ptsdeg(p1, p2, p3):
    a = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    b = sqrt((p2[0]-p3[0])**2+(p2[1]-p3[1])**2)
    c = sqrt((p1[0]-p3[0])**2+(p1[1]-p3[1])**2) 
    #print ("a, b, c = ", a, b, c)
    #print ("line2 = ", (c**2 - a**2 - b**2)/(-2*a*b))
    theta = degrees(acos((c**2 - a**2 - b**2)/(-2*a*b))) #returns degree in radians
    return theta

def isinline(pt, line):
    pt_s = Point(pt)
    dist = line.distance(pt_s)
    if dist < 1E-6:
        return True
    else:
        return False

class exp_neutpy_prep():
    def __init__(self, inp, core, sol, pfr):
        self.triangle_prep(inp, core, sol, pfr)
        self.read_triangle(core)
        self.run_neutpy(inp, core)
        
    def triangle_prep(self, inp, core, sol, pfr):
        sol_pol_pts = inp.core_thetapts_ntrl + inp.ib_thetapts_ntrl + inp.ob_thetapts_ntrl
        
        #GET POINTS FOR TRIANGULATION
        #main seperatrix
        sep_pts = np.zeros((inp.core_thetapts_ntrl, 2))
        for i, v in enumerate(np.linspace(0, 1, inp.core_thetapts_ntrl, endpoint=False)):
            sep_pts[i] = np.asarray(core.main_sep_line.interpolate(v, normalized=True).xy).T[0]
        
        #inboard divertor leg
        ib_div_pts = np.zeros((inp.ib_thetapts_ntrl, 2))
        for i, v in enumerate(np.linspace(0, 1, inp.ib_thetapts_ntrl, endpoint=True)): #skipping the x-point (point 0)
            ib_div_pts[i] = np.asarray(core.ib_div_line_cut.interpolate(v, normalized=True).xy).T[0]
        
        #outboard divertor leg
        ob_div_pts = np.zeros((inp.ob_thetapts_ntrl, 2))
        for i, v in enumerate(np.linspace(0, 1, inp.ob_thetapts_ntrl, endpoint=True)): #skipping the x-point (point 0)
            ob_div_pts[i] = np.asarray(core.ob_div_line_cut.interpolate(v, normalized=True).xy).T[0]
                
        #core
        core_pts = np.zeros((inp.core_thetapts_ntrl*len(core.core_ntrl_lines), 2))
        for num, line in enumerate(core.core_ntrl_lines):
            for i, v in enumerate(np.linspace(0, 1, inp.core_thetapts_ntrl, endpoint=False)):
                core_pts[num*inp.core_thetapts_ntrl + i] = np.asarray(line.interpolate(v, normalized=True).xy).T[0]
            
        self.core_ring = LinearRing(core_pts[:inp.core_thetapts_ntrl])
        
        
        #sol
        sol_pts = np.zeros((sol_pol_pts*len(sol.sol_lines_cut), 2))
        for num, line in enumerate(sol.sol_lines_cut):
            for i, v in enumerate(np.linspace(0, 1, sol_pol_pts, endpoint=True)):
                sol_pts[num*sol_pol_pts + i] = np.asarray(line.interpolate(v, normalized=True).xy).T[0]         

        #wall        
        wall_pts = np.asarray(inp.wall_line.coords)[:-1]
        self.wall_ring = LinearRing(wall_pts)

        all_pts = np.vstack((sep_pts, 
                             ib_div_pts, 
                             ob_div_pts, 
                             core_pts, 
                             sol_pts, 
                             wall_pts))
        
        #CREATE SEGMENTS FOR TRIANGULATION
        #WHEN DOING WALL, CHECK EACH POINT TO SEE IF IT HAS ALREADY BEEN
        #CREATED. IF SO, USE THE NUMBER OF THAT POINT AND DELETE THE WALL
        #VERSION OF IT IN THE ALL_PTS ARRAY.

        sep_segs    = np.column_stack((np.arange(inp.core_thetapts_ntrl), 
                                       np.roll(np.arange(inp.core_thetapts_ntrl), -1)))
        
        ib_div_segs = np.column_stack((np.arange(inp.ib_thetapts_ntrl), 
                                       np.roll(np.arange(inp.ib_thetapts_ntrl), -1)))[:-1]
        
        ob_div_segs = np.column_stack((np.arange(inp.ob_thetapts_ntrl), 
                                       np.roll(np.arange(inp.ob_thetapts_ntrl), -1)))[:-1]

        core_segs   = np.zeros((0, 2), dtype='int')
        for i, v in enumerate(core.core_ntrl_lines):
            new_segs = np.column_stack((np.arange(inp.core_thetapts_ntrl), 
                                        np.roll(np.arange(inp.core_thetapts_ntrl), -1))) \
                                        + inp.core_thetapts_ntrl * i
            core_segs = np.vstack((core_segs, new_segs))

        sol_segs    = np.zeros((0, 2), dtype='int')
        for i, v in enumerate(sol.sol_lines):
            new_segs = np.column_stack((np.arange(sol_pol_pts), 
                                        np.roll(np.arange(sol_pol_pts), -1)))[:-1] \
                                        + sol_pol_pts * i
            sol_segs = np.vstack((sol_segs, new_segs))
        
        wall_segs   = np.column_stack((np.arange(len(wall_pts)), 
                                       np.roll(np.arange(len(wall_pts)), -1)))
        
        all_segs = np.vstack((sep_segs, 
                              ib_div_segs + len(sep_segs), 
                              ob_div_segs + len(ib_div_segs) + len(sep_segs) + 1, 
                              core_segs + len(ob_div_segs) + len(ib_div_segs) + len(sep_segs) + 1 + 1, 
                              sol_segs + len(core_segs) + len(ob_div_segs) + len(ib_div_segs) + len(sep_segs)  + 1 + 1, 
                              wall_segs + len(sol_segs) + len(core_segs) + len(ob_div_segs) + len(ib_div_segs) + len(sep_segs) + 1 + 1 + inp.num_sollines
                              ))

        all_pts_unique = np.unique(all_pts, axis=0)

        #CLEANUP
        #NOTE: this process will result in a segments array that looks fairly chaotic, 
        #but will ensure that the triangulation goes smoothly.
        
        #loop over each point in all_segs
        #look up the point's coordinates in all_pts
        #find the location of those coordinates in all_pts_unique
        #put that location in the corresponding location in all_segs_unique
        
        all_segs_unique = np.zeros(all_segs.flatten().shape, dtype='int')
        for i, pt in enumerate(all_segs.flatten()):
            pt_coords = all_pts[pt]
            loc_unique = np.where((all_pts_unique==pt_coords).all(axis=1))[0][0]
            all_segs_unique[i] = loc_unique
        all_segs_unique = all_segs_unique.reshape(-1, 2)  

        ## OUTPUT .poly FILE AND RUN TRIANGLE PROGRAM
        open('./outputs/exp_mesh.poly', 'w').close()
        outfile = open('./outputs/exp_mesh.poly', 'ab')
        filepath = os.path.realpath(outfile.name)
        np.savetxt(outfile, 
                   np.array([all_pts_unique.shape[0], 2, 0, 0])[None], 
                   fmt='%i %i %i %i')
        np.savetxt(outfile, 
                   np.column_stack((np.arange(len(all_pts_unique)), 
                                    all_pts_unique)), 
                   fmt='%i %f %f')
        np.savetxt(outfile, 
                   np.array([all_segs_unique.shape[0], 0])[None], 
                   fmt='%i %i')
        np.savetxt(outfile, 
                   np.column_stack((np.arange(len(all_segs_unique)), 
                                    all_segs_unique, 
                                    np.zeros(len(all_segs_unique), dtype='int'))), 
                   fmt='%i %i %i %i')
        np.savetxt(outfile, 
                   np.array([1])[None], 
                   fmt='%i')
        np.savetxt(outfile, 
                   np.array([1, core.m_axis[0], core.m_axis[1]])[None], 
                   fmt='%i %f %f')
        np.savetxt(outfile, 
                   np.array([0])[None], 
                   fmt='%i')
        outfile.close()
        
        #construct options to pass to triangle, as specified in input file
        #refer to https://www.cs.cmu.edu/~quake/triangle.html
        
        tri_options = '-p'
        try:
            tri_options = tri_options + 'q' + str(inp.tri_min_angle)
        except:
            pass
        
        try:
            tri_options = tri_options + 'a' + str(inp.tri_min_area)
        except:
            pass
        
        tri_options = tri_options + 'nz'
        #call triangle
        try:
            call(['triangle', tri_options, filepath])
        except AttributeError:
            try:
                call(['triangle', tri_options, filepath])
            except:
                print 'triangle could not be found. Stopping.'
                sys.exit
      
    def read_triangle(self, core):
        ## READ TRIANGLE OUTPUT

        ## DECLARE FILE PATHS
        nodepath = os.getcwd() + '/outputs/exp_mesh.1.node'
        elepath = os.getcwd() + '/outputs/exp_mesh.1.ele'
        neighpath = os.getcwd() + '/outputs/exp_mesh.1.neigh'

        ## GET NODE DATA
        with open(nodepath, 'r') as node:
            #dummy = next(mil_mesh)
            nodecount = re.findall(r'\d+', next(node))
            nNodes = int(nodecount[0])
            nodenum = np.zeros(nNodes)
            nodesx = np.zeros(nNodes)
            nodesy = np.zeros(nNodes)
            
            for i in range (0, nNodes):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(node)) 
                nodenum[i] = int(data1[0])
                nodesx[i] = data1[1]
                nodesy[i] = data1[2]
                
        ## GET TRIANGLE DATA
        with open(elepath, 'r') as tri_file:
            tricount = re.findall(r'\d+', next(tri_file))
            nTri = int(tricount[0])
            print 'number of triangles = ', nTri
            triangles = np.zeros((nTri, 3))
            tri_regions = np.zeros(nTri)
            for i in range (0, nTri):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(tri_file))
                triangles[i, 0] = data1[1]
                triangles[i, 1] = data1[2]
                triangles[i, 2] = data1[3]
                #tri_regions[i] = data1[4]
        triangles = triangles.astype('int')
        tri_regions = tri_regions.astype('int')

        ## GET NEIGHBOR DATA
        with open(neighpath, 'r') as neigh_file:
            neighcount = re.findall(r'\d+', next(neigh_file))
            nNeigh = int(neighcount[0])
            neighbors = np.zeros((nNeigh, 3))
            for i in range (0, nNeigh):
                data1 = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', next(neigh_file))
                neighbors[i, 0] = data1[1]
                neighbors[i, 1] = data1[2]
                neighbors[i, 2] = data1[3]
        neighbors = neighbors.astype('int')

        ## REARRANGE TRIANGLES TO CONFORM TO GTNEUT CONVENTION      
        triangles = np.fliplr(triangles) #triangle vertices are given counterclockwise, but we want clockwise
        neighbors = np.fliplr(neighbors) #neighbor 1 is opposite vertex 1, so also counterclockwise 
        
        y=np.zeros(3)
        for i, tri in enumerate(triangles):
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
            triangles[i] = np.roll(triangles[i], -1*miny_index)
            neighbors[i] = np.roll(neighbors[i], -1*miny_index-2) # the -2 is because the side 1 is opposite vertex 1. We want side 1 to start at vertex 1
        
        ## GET VALUES TO ORIENT THE FIRST CELL WHEN PLOTTING
        point1_x = nodesx[triangles[0, 0]]
        point1_y = nodesy[triangles[0, 0]]
        point2_x = nodesx[triangles[0, 1]]
        point2_y = nodesy[triangles[0, 1]]
        point3_x = nodesx[triangles[0, 2]]
        point3_y = nodesy[triangles[0, 2]]
        
        cell1_ctr_x = (point1_x + point2_x + point3_x) / 3
        cell1_ctr_y = (point1_y + point2_y + point3_y) / 3

        ## CALCULATE ANGLE BY WHICH TO ROTATE THE FIRST CELL WHEN PLOTTING
        cell1_theta0 = degrees(getangle([point3_x, point3_y], [point1_x, point1_y]))    

        ## GET VALUES TO ORIENT THE FIRST CELL WHEN PLOTTING
        point1_x = nodesx[triangles[0, 0]]
        point1_y = nodesy[triangles[0, 0]]
        point2_x = nodesx[triangles[0, 1]]
        point2_y = nodesy[triangles[0, 1]]
        point3_x = nodesx[triangles[0, 2]]
        point3_y = nodesy[triangles[0, 2]]
        
        cell1_ctr_x = (point1_x + point2_x + point3_x) / 3
        cell1_ctr_y = (point1_y + point2_y + point3_y) / 3

        ## CALCULATE ANGLE BY WHICH TO ROTATE THE FIRST CELL WHEN PLOTTING
        cell1_theta0 = degrees(getangle([point3_x, point3_y], [point1_x, point1_y]))

        ## CALCULATE MID POINTS OF TRIANGLES, AS WELL AS MIDPOINTS FOR EACH FACE     
        ptsx = np.zeros((nTri, 3))
        ptsy = np.zeros((nTri, 3))
        #for index, tri in ndenumerate(triangles):
        for i in range(0, nTri):
            ptsx[i, 0] = nodesx[triangles[i, 0]]
            ptsy[i, 0] = nodesy[triangles[i, 0]]
            ptsx[i, 1] = nodesx[triangles[i, 1]]
            ptsy[i, 1] = nodesy[triangles[i, 1]]
            ptsx[i, 2] = nodesx[triangles[i, 2]]
            ptsy[i, 2] = nodesy[triangles[i, 2]]
        
        mid_x = np.mean(ptsx, axis=1)
        mid_y = np.mean(ptsy, axis=1)
        self.midpts = np.column_stack((mid_x, mid_y))
        
        side1_midx = (ptsx[:, 0] + ptsx[:, 1])/2
        side2_midx = (ptsx[:, 1] + ptsx[:, 2])/2
        side3_midx = (ptsx[:, 2] + ptsx[:, 0])/2
        
        side1_midy = (ptsy[:, 0] + ptsy[:, 1])/2
        side2_midy = (ptsy[:, 1] + ptsy[:, 2])/2
        side3_midy = (ptsy[:, 2] + ptsy[:, 0])/2
        
        side1_midpt = np.column_stack((side1_midx, side1_midy))
        side2_midpt = np.column_stack((side2_midx, side2_midy))
        side3_midpt = np.column_stack((side3_midx, side3_midy))
        
        #COMBINE POINTS FOR THE PLASMA, SOL, AND DIVERTOR REGIONS
        #first fill in plasma cells
        plasmacells = np.zeros((1, 2))
        pcellnum = nTri
        pcellcount = 0
        
        for index, nei in enumerate(neighbors):
            #for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index], self.core_ring)
            side2inline = isinline(side2_midpt[index], self.core_ring)
            side3inline = isinline(side3_midpt[index], self.core_ring)

            if side1inline or side2inline or side3inline:
                
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                #print 'nb = ', nb
                #print side1_midpt[index]
                #print side2_midpt[index]
                #print side3_midpt[index]
                if nb == 1: #cell has one plasma border
                    
                    #create plasma cell
                    plasmacells[pcellcount, 0] = pcellnum
                    plasmacells[pcellcount, 1] = index
                    plasmacells = np.vstack((plasmacells, [0, 0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1
                elif nb == 2: 
                    #cell has two plasma borders (this will probably never happen. It would require a local 
                    #concavity in the inner-most meshed flux surface)
                    #create plasma cell #1
                    plasmacells[pcellcount, 0] = pcellnum
                    plasmacells[pcellcount, 1] = index
                    plasmacells = np.vstack((plasmacells, [0, 0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1   
                    
                    #create plasma cell #2
                    plasmacells[pcellcount, 0] = pcellnum
                    plasmacells[pcellcount, 1] = index
                    plasmacells = np.vstack((plasmacells, [0, 0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    #get ready for next run
                    pcellnum +=1
                    pcellcount +=1
        plasmacells = np.delete(plasmacells, -1, 0)
        plasmacells = plasmacells.astype('int')

        #now fill in wall cells
        wallcells = np.zeros((1, 6))
        wcellnum = pcellnum #was already advanced in the plasmacell loop. Don't add 1.
        wcellcount = 0
        #print np.asarray(self.wall_ring.xy).T
        #plt.scatter(np.asarray(self.wall_ring.xy).T[:, 0], np.asarray(self.wall_ring.xy).T[:, 1], marker='o', s=0.5)
        
        for index, nei in enumerate(neighbors):
            #for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index], self.wall_ring)
            side2inline = isinline(side2_midpt[index], self.wall_ring)
            side3inline = isinline(side3_midpt[index], self.wall_ring)
            
            if side1inline or side2inline or side3inline:
                #print index, nei, side1inline, side2inline, side3inline
                nb = (nei == -1).sum() #count number of times -1 occurs in nei
                if nb == 1: #cell has one wall border
                    #identify the side that is the wall cell
                    sidenum = np.where(np.asarray([side1inline, side2inline, side3inline]))[0][0]
                    if sidenum==0:
                        pt = side1_midpt[index]
                    elif sidenum==1:
                        pt = side2_midpt[index]
                    elif sidenum==2:
                        pt = side3_midpt[index]
                    
                    #create wall cell
                    wallcells[wcellcount, 0] = wcellnum
                    wallcells[wcellcount, 1] = index
                    wallcells[wcellcount, 2] = griddata(core.ni_pts[:, :-1], core.ni_pts[:, -1], pt, method='nearest', rescale=True)
                    wallcells[wcellcount, 3] = griddata(core.ne_pts[:, :-1], core.ne_pts[:, -1], pt, method='nearest', rescale=True)
                    wallcells[wcellcount, 4] = griddata(core.Ti_kev_pts[:, :-1], core.Ti_kev_pts[:, -1], pt, method='nearest', rescale=True)
                    wallcells[wcellcount, 5] = griddata(core.Te_kev_pts[:, :-1], core.Te_kev_pts[:, -1], pt, method='nearest', rescale=True)
                    wallcells = np.vstack((wallcells, [0, 0, 0, 0, 0, 0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1
                elif nb == 2: #cell has two wall borders (This can easily happen because the wall has many concave points.)
                    #create wall cell #1
                    wallcells[wcellcount, 0] = wcellnum
                    wallcells[wcellcount, 1] = index
                    wallcells = np.vstack((wallcells, [0, 0, 0, 0, 0, 0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1   
                    
                    #create wall cell #2
                    wallcells[wcellcount, 0] = wcellnum
                    wallcells[wcellcount, 1] = index
                    wallcells = np.vstack((wallcells, [0, 0, 0, 0, 0, 0]))
                    #update neighbors
                    nei[np.argmax(nei==-1)] = wcellnum
                    #get ready for next run
                    wcellnum +=1
                    wcellcount +=1
        wallcells = np.delete(wallcells, -1, 0)
        wallcells = wallcells.astype('int')
        ## POPULATE CELL DENSITIES AND TEMPERATURES
        #create array of all points in plasma, sol, id, and od
        #tri_param = np.vstack((plasma_param, sol_param, id_param, od_param))
        
        ni_tri = griddata(core.ni_pts[:, :-1], 
                          core.ni_pts[:, -1], 
                          (mid_x, mid_y), 
                          method='linear', 
                          fill_value=0, 
                          rescale=True)
        ne_tri = griddata(core.ne_pts[:, :-1], 
                          core.ne_pts[:, -1], 
                          (mid_x, mid_y), 
                          method='linear', 
                          fill_value=0, 
                          rescale=True)
        Ti_tri = griddata(core.Ti_kev_pts[:, :-1], 
                          core.Ti_kev_pts[:, -1], 
                          (mid_x, mid_y), 
                          method='linear', 
                          fill_value=0, 
                          rescale=True)
        Te_tri = griddata(core.Te_kev_pts[:, :-1], 
                          core.Te_kev_pts[:, -1], 
                          (mid_x, mid_y), 
                          method='linear', 
                          fill_value=0, 
                          rescale=True)

        #ni_tri[ni_tri<1.0E16] = 1.0E16
        #ne_tri[ne_tri<1.0E16] = 1.0E16
        #Ti_tri[Ti_tri<0.002] = 0.002
        #Te_tri[Te_tri<0.002] = 0.002
        
        ## CALCULATE LENGTHS OF SIDES
        lsides = np.zeros((nTri, 3))
        for i in range (0, nTri):
            lsides[i, 0] = sqrt((ptsx[i, 0]-ptsx[i, 1])**2 + (ptsy[i, 0]-ptsy[i, 1])**2)
            lsides[i, 1] = sqrt((ptsx[i, 1]-ptsx[i, 2])**2 + (ptsy[i, 1]-ptsy[i, 2])**2)
            lsides[i, 2] = sqrt((ptsx[i, 2]-ptsx[i, 0])**2 + (ptsy[i, 2]-ptsy[i, 0])**2)
        
        ## CALCULATE CELL ANGLES
        angles = np.zeros((nTri, 3))
        for i in range (0, nTri):
            p1 = np.array([ptsx[i, 0], ptsy[i, 0]])
            p2 = np.array([ptsx[i, 1], ptsy[i, 1]])
            p3 = np.array([ptsx[i, 2], ptsy[i, 2]])
            angles[i, 0] = getangle3ptsdeg(p1, p2, p3)
            angles[i, 1] = getangle3ptsdeg(p2, p3, p1)
            angles[i, 2] = getangle3ptsdeg(p3, p1, p2)
        
        ## WRITE NEUTPY INPUT FILE!    
        f = open(os.getcwd() + '/neutpy_in_generated', 'w')
        f.write('nCells = ' + str(nTri) + ' nPlasmReg = ' + str(pcellcount) + ' nWallSegm = ' + str(wcellcount))
        for i in range(0, nTri):
            f.write('\n'+'iType(' + str(i) + ') = 0 nSides(' + str(i) + ') = 3 ' + 'adjCell('+str(i)+') = '+', '.join(map(str, neighbors[i, :])))
        f.write('\n')
        f.write('\n#lsides and angles for normal cells')
        for i in range(0, nTri):
            f.write('\n'+'lsides(' + str(i) + ') = '+', '.join(map(str, lsides[i, :]))+' angles(' + str(i) + ') = '+', '.join(map(str, angles[i, :])))
        f.write('\n')
        f.write('\n#densities and temperatures for normal cells')
        for i in range(0, nTri):
            f.write('\n'+'elecTemp('+str(i)+') = '+str(Te_tri[i]) +' elecDens(' + str(i) + ') = '+str(ne_tri[i])+' ionTemp('+str(i)+') = '+str(Ti_tri[i]) +' ionDens(' + str(i) + ') = '+str(ni_tri[i]))
        f.write('\n')
        f.write('\n#wall cells')
        for i, wcell in enumerate(wallcells):
            f.write('\n'+'iType('+str(wcell[0])+') = 2 nSides('+str(wcell[0])+') = 1 adjCell('+str(wcell[0])+') = '+str(wcell[1])+' zwall('+str(wcell[0])+') = 6 awall('+str(wcell[0])+') = 12 twall('+str(wcell[0])+') = '+str(wcell[4])+' f_abs('+str(wcell[0])+') = 0.0 s_ext('+str(wcell[0])+') = 1.0E19') 
        f.write('\n')
        f.write('\n#plasma core and vacuum cells')
        for i, pcell in enumerate(plasmacells):
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
        self.toneutpy={}
        self.toneutpy["nCells"]       = nTri
        self.toneutpy["nPlasmReg"]    = pcellcount
        self.toneutpy["nWallSegm"]    = wcellcount
        self.toneutpy["aneut"]        = 2
        self.toneutpy["zion"]         = 1
        self.toneutpy["aion"]         = 2
        self.toneutpy["tslow"]        = 0.002
        self.toneutpy["int_method"]   = 'midpoint'
        self.toneutpy["phi_int_pts"]  = 10
        self.toneutpy["xi_int_pts"]   = 10
        self.toneutpy["xsec_ioni"]    = 'janev'
        self.toneutpy["xsec_ione"]    = 'janev'
        self.toneutpy["xsec_cx"]      = 'janev'
        self.toneutpy["xsec_rec"]     = 'stacey_thomas'
        self.toneutpy["xsec_el"]      = 'janev'
        self.toneutpy["xsec_eln"]     = 'stacey_thomas'
        self.toneutpy["refmod_e"]     = 'stacey'
        self.toneutpy["refmod_n"]     = 'stacey'
        
        self.toneutpy["iType"]        = np.asarray([0]*nTri + [1]*pcellcount + [2]*wcellcount)
        self.toneutpy["nSides"]       = np.asarray([3]*nTri + [1]*(pcellcount + wcellcount))
        self.toneutpy["zwall"]        = np.asarray([0]*(nTri+pcellcount) + [6]*wcellcount)
        self.toneutpy["awall"]        = np.asarray([0]*(nTri+pcellcount) + [12]*wcellcount)
        self.toneutpy["elecTemp"]     = Te_tri[:nTri]
        self.toneutpy["ionTemp"]      = Ti_tri[:nTri]
        self.toneutpy["elecDens"]     = ne_tri[:nTri]
        self.toneutpy["ionDens"]      = ni_tri[:nTri]
        self.toneutpy["twall"]        = np.asarray([0]*nTri + [5000]*pcellcount + [0.002]*wcellcount)
        self.toneutpy["f_abs"]        = np.asarray([0]*(nTri+pcellcount) + [0]*wcellcount)
        self.toneutpy["alb_s"]        = np.asarray([0]*nTri + [0]*pcellcount + [0]*wcellcount)
        self.toneutpy["alb_t"]        = np.asarray([0]*nTri + [0]*pcellcount + [0]*wcellcount)
        self.toneutpy["s_ext"]        = np.asarray([0.0]*nTri + [0.0]*pcellcount + [0.0]*wcellcount)
        
        self.toneutpy["adjCell"]      = neighbors
        self.toneutpy["lsides"]       = lsides
        self.toneutpy["angles"]       = angles
        self.toneutpy["cell1_ctr_x"]  = cell1_ctr_x
        self.toneutpy["cell1_ctr_y"]  = cell1_ctr_y
        self.toneutpy["cell1_theta0"] = cell1_theta0
        
    def run_neutpy(self, inp, core):
        time0 = time.time()
        self.neutpy_inst = neutpy(inarrs=self.toneutpy)
        time1 = time.time()
        print 'neutpy time = ', time1-time0
        plot = neutpyplot(self.neutpy_inst)
        self.n_n_slow = self.neutpy_inst.cell_nn_s
        self.n_n_thermal = self.neutpy_inst.cell_nn_t
        self.n_n_total = self.nn_s_raw + self.nn_t_raw
        
        self.izn_rate_slow = self.neutpy_inst.cell_izn_rate_s
        self.izn_rate_thermal = self.neutpy_inst.cell_izn_rate_t
        self.izn_rate_total = self.iznrate_s_raw + self.iznrate_t_raw
        
        #create output file
        #the file contains R, Z coordinates and then the values of several calculated parameters
        #at each of those points.
        
        f = open(inp.neut_outfile, 'w')
        f.write(('{:^18s}'*8).format('R', 'Z', 'n_n_slow', 'n_n_thermal', 'n_n_total', 'izn_rate_slow', 'izn_rate_thermal', 'izn_rate_total'))
        for i, pt in enumerate(self.midpts):
            f.write(('\n'+'{:>18.5f}'*2+'{:>18.5E}'*6).format(
                                        self.midpts[i, 0], 
                                        self.midpts[i, 1], 
                                        self.n_n_slow[i], 
                                        self.n_n_thermal[i],
                                        self.izn_rate_slow[i], 
                                        self.izn_rate_thermal[i]))
        f.close()

        # update neutral data in core instance
        dict = {}
        dict['R'] = self.midpts[:, 0]
        dict['Z'] = self.midpts[:, 1]
        dict['n_n_slow'] = self.n_n_slow,
        dict['n_n_thermal'] = self.n_n_thermal
        dict['izn_rate_slow'] = self.izn_rate_slow
        dict['izn_rate_thermal'] = self.izn_rate_thermal

        data = namedtuple('data', dict.keys())(*dict.values())
        core.update_ntrl_data(data)
