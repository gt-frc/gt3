#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 14:27:42 2018

@author: max
"""
from __future__ import division
import numpy as np
from neutpy import neutpy
from neutpy_tools import NeutpyTools
from math import degrees, sqrt
from scipy.interpolate import griddata
from collections import namedtuple
import sys
import os
import re
from subprocess import Popen, PIPE
from SOL.sol import Sol
from PFR.pfr import Pfr
from Functions.CalcCoreLinesNtrl import calc_core_lines_ntrl
from Functions.CreateTriangleInfile import create_triangle_infile
from Functions.CreateTriangleOpts import create_triangle_opts
from Functions.CreateTriPts import create_tri_pts
from Functions.CreateTriSegs import create_tri_segs
from Functions.GetAngle import getangle
from Functions.GetAngle3ptsDeg import getangle3ptsdeg
from Functions.IsInLine import isinline


class Neutrals:
    def __init__(self, inp, core):
        # Try to read in specified neutrals data file. If it's not there, then prepare inputs for and run neutpy

        try:
            ntrl_data = np.loadtxt(inp.neutfile_loc, delimiter=',', skiprows=1)
            self.data = namedtuple('data', 'R Z n_n_slow n_n_thermal izn_rate_slow izn_rate_thermal')(
                ntrl_data[:, 1],
                ntrl_data[:, 2],
                ntrl_data[:, 3],
                ntrl_data[:, 4],
                ntrl_data[:, 6],
                ntrl_data[:, 7]
            )
        except:
            # instantiate sol and pfr
            sol = Sol(inp, core)
            pfr = Pfr(inp, core)

            # assemble lines for neutrals calculation
            lines = namedtuple('lines', 'core sep sol pfr ib_div ob_div wall')(
                calc_core_lines_ntrl(core),
                core.lines.sep,
                sol.sol_lines_cut,
                pfr.pfr_line,
                core.lines.div.ib,
                core.lines.div.ob,
                inp.wall_line
            )

            # assemble density and temperature points and values for neutrals calculation
            nT = namedtuple('nT', 'core sol wall')(
                namedtuple('core_nT', 'ni ne Ti Te')(
                    np.column_stack((core.R.flatten(), core.Z.flatten(), core.n.i.flatten())),
                    np.column_stack((core.R.flatten(), core.Z.flatten(), core.n.e.flatten())),
                    np.column_stack((core.R.flatten(), core.Z.flatten(), core.T.i.kev.flatten())),
                    np.column_stack((core.R.flatten(), core.Z.flatten(), core.T.e.kev.flatten()))
                ),
                sol.sol_nT,
                sol.wall_nT
            )

            # get points from those lines
            pts, core_ring = create_tri_pts(inp, lines)

            # get segments from those points
            pts_unique, segs_unique = create_tri_segs(inp, lines, pts)

            # create the triangle input file from the points and segments
            triangle_infile = create_triangle_infile(pts_unique, segs_unique, core)

            # specify triangle options based on input file specs
            triangle_opts = create_triangle_opts(inp)


            # run triangle
            print 'running triangle'

            # set the name of the executable based on the operating system
            if os.name == 'nt':
                triangle_name = 'triangle.exe'
            elif os.name == 'posix':
                triangle_name = 'triangle'
            else:
                print 'Not sure what os you\'re running. If mac, you might need to add some code \
                        to the neutpy_prep module to help it find and run triangle.'
                sys.exit()
            try:
                p = Popen([triangle_name, triangle_opts, triangle_infile], stdin=PIPE, stdout=PIPE).wait()
                #call(['triangle', triangle_opts, triangle_infile])
            except:
                try:
                    p = Popen([inp.triangle_loc, triangle_opts, triangle_infile], stdin=PIPE, stdout=PIPE).wait()
                    #call([inp.triangle_loc, triangle_opts, triangle_infile])
                except:
                    print 'Unable to find triangle executable. Stopping.'
                    sys.exit(status=1)

            midpts, toneutpy = self.create_neutpy_input(inp, core, lines, nT, core_ring)

            # run neutpy
            self.neutpy_inst = neutpy(inarrs=toneutpy)

            # output cell values (this is also done by neutpy tools in the next step, but stored somewhere else)
            create_cell_outfile(self.neutpy_inst, inp.neut_cell_file, midpts)

            # instantiating NeutpyTools and create neutrals data files. To make neutral densities and related
            # charts, you can use neutpy tools directly by passing it a namedtuple of the output files
            ntools = NeutpyTools(self.neutpy_inst)
            ntools.create_flux_outfile(inp.neut_face_file)
            ntools.create_cell_outfile(inp.neut_cell_file)

            self.data = namedtuple('data', 'R Z n_n_slow n_n_thermal izn_rate_slow izn_rate_thermal')(
                midpts[:, 0],
                midpts[:, 1],
                self.neutpy_inst.nn.s,
                self.neutpy_inst.nn.t,
                self.neutpy_inst.izn_rate.s,
                self.neutpy_inst.izn_rate.t
            )

        try:
            core.update_ntrl_data(self.data)
        except:
            print 'unable to update values in core instance.'
            pass


    @staticmethod
    def create_neutpy_input(inp, core, lines, nT, core_ring):

        # Assemble global density and temperature data
        ni_global = np.vstack((nT.core.ni,
                               nT.sol.ni,
                               nT.wall.ni))

        ne_global = np.vstack((nT.core.ne,
                               nT.sol.ne,
                               nT.wall.ne))

        Ti_global = np.vstack((nT.core.Ti,
                               nT.sol.Ti,
                               nT.wall.Ti))

        Te_global = np.vstack((nT.core.Te,
                               nT.sol.Te,
                               nT.wall.Te))

        # READ TRIANGLE OUTPUT
        # DECLARE FILE PATHS
        nodepath = os.getcwd() + '/outputs/exp_mesh.1.node'
        elepath = os.getcwd() + '/outputs/exp_mesh.1.ele'
        neighpath = os.getcwd() + '/outputs/exp_mesh.1.neigh'

        # GET NODE DATA
        with open(nodepath, 'r') as node:
            # dummy = next(mil_mesh)
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

        # GET TRIANGLE DATA
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
                # tri_regions[i] = data1[4]
        triangles = triangles.astype('int')
        tri_regions = tri_regions.astype('int')

        # GET NEIGHBOR DATA
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

        # REARRANGE TRIANGLES TO CONFORM TO GTNEUT CONVENTION
        triangles = np.fliplr(triangles)  # triangle vertices are given counterclockwise, but we want clockwise
        neighbors = np.fliplr(neighbors)  # neighbor 1 is opposite vertex 1, so also counterclockwise

        y=np.zeros(3)
        for i, tri in enumerate(triangles):
            # Find lowest value of y component of vertices
            y[0] = nodesy[tri[0]]
            y[1] = nodesy[tri[1]]
            y[2] = nodesy[tri[2]]
            miny = np.amin(y)
            miny_count = np.sum(y == miny)
            if miny_count == 1:
                # identify position of minimum and roll array accordingly
                miny_index = np.where(y==miny)[0][0]
            else:
                # identify which points are the two minima and determine
                # which of them is farthest to the left (or right if I change it)
                miny_index = np.where(y==miny)[0][1] #change this 1 to a zero to choose the rightmost of the two bottom vertices
            triangles[i] = np.roll(triangles[i], -1*miny_index)
            neighbors[i] = np.roll(neighbors[i], -1*miny_index-2) # the -2 is because the side 1 is opposite vertex 1. We want side 1 to start at vertex 1

        # GET VALUES TO ORIENT THE FIRST CELL WHEN PLOTTING
        point1_x = nodesx[triangles[0, 0]]
        point1_y = nodesy[triangles[0, 0]]
        point2_x = nodesx[triangles[0, 1]]
        point2_y = nodesy[triangles[0, 1]]
        point3_x = nodesx[triangles[0, 2]]
        point3_y = nodesy[triangles[0, 2]]

        cell1_ctr_x = (point1_x + point2_x + point3_x) / 3
        cell1_ctr_y = (point1_y + point2_y + point3_y) / 3

        # CALCULATE ANGLE BY WHICH TO ROTATE THE FIRST CELL WHEN PLOTTING
        cell1_theta0 = degrees(getangle([point3_x, point3_y], [point1_x, point1_y]))

        # CALCULATE MID POINTS OF TRIANGLES, AS WELL AS MIDPOINTS FOR EACH FACE
        ptsx = np.zeros((nTri, 3))
        ptsy = np.zeros((nTri, 3))
        for i in range(0, nTri):
            ptsx[i, 0] = nodesx[triangles[i, 0]]
            ptsy[i, 0] = nodesy[triangles[i, 0]]
            ptsx[i, 1] = nodesx[triangles[i, 1]]
            ptsy[i, 1] = nodesy[triangles[i, 1]]
            ptsx[i, 2] = nodesx[triangles[i, 2]]
            ptsy[i, 2] = nodesy[triangles[i, 2]]

        mid_x = np.mean(ptsx, axis=1)
        mid_y = np.mean(ptsy, axis=1)
        midpts = np.column_stack((mid_x, mid_y))

        # get side midpoints
        side1_midx = (ptsx[:, 0] + ptsx[:, 1])/2
        side2_midx = (ptsx[:, 1] + ptsx[:, 2])/2
        side3_midx = (ptsx[:, 2] + ptsx[:, 0])/2

        side1_midy = (ptsy[:, 0] + ptsy[:, 1])/2
        side2_midy = (ptsy[:, 1] + ptsy[:, 2])/2
        side3_midy = (ptsy[:, 2] + ptsy[:, 0])/2

        side1_midpt = np.column_stack((side1_midx, side1_midy))
        side2_midpt = np.column_stack((side2_midx, side2_midy))
        side3_midpt = np.column_stack((side3_midx, side3_midy))

        # COMBINE POINTS FOR THE PLASMA, SOL, AND DIVERTOR REGIONS
        # first fill in plasma cells
        plasmacells = np.zeros((1, 2))
        pcellnum = nTri
        pcellcount = 0

        for index, nei in enumerate(neighbors):
            # for each face of the cell, find the mid-point and check if it falls along the innermost flux surface being used
            side1inline = isinline(side1_midpt[index], core_ring)
            side2inline = isinline(side2_midpt[index], core_ring)
            side3inline = isinline(side3_midpt[index], core_ring)

            if side1inline or side2inline or side3inline:
                # count number of times -1 occurs in nei
                nb = (nei == -1).sum()

                if nb == 1:  # cell has one plasma border

                    # create plasma cell
                    plasmacells[pcellcount, 0] = pcellnum
                    plasmacells[pcellcount, 1] = index
                    plasmacells = np.vstack((plasmacells, [0, 0]))
                    # update neighbors
                    nei[np.argmax(nei == -1)] = pcellnum
                    # get ready for next run
                    pcellnum += 1
                    pcellcount += 1
                elif nb == 2:
                    # cell has two plasma borders (this will probably never happen. It would require a local
                    # concavity in the inner-most meshed flux surface)
                    # create plasma cell #1
                    plasmacells[pcellcount, 0] = pcellnum
                    plasmacells[pcellcount, 1] = index
                    plasmacells = np.vstack((plasmacells, [0, 0]))
                    # update neighbors
                    nei[np.argmax(nei == -1)] = pcellnum
                    # get ready for next run
                    pcellnum +=1
                    pcellcount +=1

                    # create plasma cell #2
                    plasmacells[pcellcount, 0] = pcellnum
                    plasmacells[pcellcount, 1] = index
                    plasmacells = np.vstack((plasmacells, [0, 0]))
                    # update neighbors
                    nei[np.argmax(nei==-1)] = pcellnum
                    # get ready for next run
                    pcellnum +=1
                    pcellcount +=1

        plasmacells = np.delete(plasmacells, -1, 0)
        plasmacells = plasmacells.astype('int')

        # now fill in wall cells
        wallcells = np.zeros((1, 6))
        wcellnum = pcellnum  # was already advanced in the plasmacell loop. Don't add 1.
        wcellcount = 0

        for index, nei in enumerate(neighbors):
            # for each face of the cell, find the mid-point and check if it falls in line
            side1inline = isinline(side1_midpt[index], inp.wall_line)
            side2inline = isinline(side2_midpt[index], inp.wall_line)
            side3inline = isinline(side3_midpt[index], inp.wall_line)

            if side1inline or side2inline or side3inline:
                # print index, nei, side1inline, side2inline, side3inline
                nb = (nei == -1).sum()  # count number of times -1 occurs in nei
                if nb == 1:  # cell has one wall border
                    # identify the side that is the wall cell
                    sidenum = np.where(np.asarray([side1inline, side2inline, side3inline]))[0][0]
                    if sidenum == 0:
                        pt = side1_midpt[index]
                    elif sidenum == 1:
                        pt = side2_midpt[index]
                    elif sidenum == 2:
                        pt = side3_midpt[index]

                    # create wall cell
                    wallcells[wcellcount, 0] = wcellnum
                    wallcells[wcellcount, 1] = index
                    wallcells[wcellcount, 2] = griddata(ni_global[:, :2], ni_global[:, 2], pt, method='nearest', rescale=True)
                    wallcells[wcellcount, 3] = griddata(ne_global[:, :2], ne_global[:, 2], pt, method='nearest', rescale=True)
                    wallcells[wcellcount, 4] = griddata(Ti_global[:, :2], Ti_global[:, 2], pt, method='nearest', rescale=True)
                    wallcells[wcellcount, 5] = griddata(Te_global[:, :2], Te_global[:, 2], pt, method='nearest', rescale=True)
                    wallcells = np.vstack((wallcells, [0, 0, 0, 0, 0, 0]))
                    # update neighbors
                    nei[np.argmax(nei == -1)] = wcellnum
                    # get ready for next run
                    wcellnum +=1
                    wcellcount +=1
                elif nb == 2:  # cell has two wall borders (This can easily happen because the wall has many concave points.)
                    # create wall cell #1
                    wallcells[wcellcount, 0] = wcellnum
                    wallcells[wcellcount, 1] = index
                    wallcells = np.vstack((wallcells, [0, 0, 0, 0, 0, 0]))
                    # update neighbors
                    nei[np.argmax(nei == -1)] = wcellnum
                    # get ready for next run
                    wcellnum += 1
                    wcellcount += 1

                    # create wall cell #2
                    wallcells[wcellcount, 0] = wcellnum
                    wallcells[wcellcount, 1] = index
                    wallcells = np.vstack((wallcells, [0, 0, 0, 0, 0, 0]))
                    # update neighbors
                    nei[np.argmax(nei == -1)] = wcellnum
                    # get ready for next run
                    wcellnum += 1
                    wcellcount += 1
        wallcells = np.delete(wallcells, -1, 0)
        wallcells = wallcells.astype('int')

        # POPULATE CELL DENSITIES AND TEMPERATURES
        # create array of all points in plasma, sol, id, and od
        # tri_param = np.vstack((plasma_param, sol_param, id_param, od_param))

        ni_tri = griddata(ni_global[:, :2],
                          ni_global[:, 2],
                          (mid_x, mid_y),
                          method='nearest',
                          fill_value=0,
                          rescale=True)
        ne_tri = griddata(ne_global[:, :2],
                          ne_global[:, 2],
                          (mid_x, mid_y),
                          method='nearest',
                          fill_value=0,
                          rescale=True)
        Ti_tri = griddata(Ti_global[:, :2],
                          Ti_global[:, 2],
                          (mid_x, mid_y),
                          method='nearest',
                          fill_value=0,
                          rescale=True)
        Te_tri = griddata(Te_global[:, :2],
                          Te_global[:, 2],
                          (mid_x, mid_y),
                          method='nearest',
                          fill_value=0,
                          rescale=True)

        # ni_tri[ni_tri<1.0E16] = 1.0E16
        # ne_tri[ne_tri<1.0E16] = 1.0E16
        # Ti_tri[Ti_tri<0.002] = 0.002
        # Te_tri[Te_tri<0.002] = 0.002

        # CALCULATE LENGTHS OF SIDES
        lsides = np.zeros((nTri, 3))
        for i in range (0, nTri):
            lsides[i, 0] = sqrt((ptsx[i, 0]-ptsx[i, 1])**2 + (ptsy[i, 0]-ptsy[i, 1])**2)
            lsides[i, 1] = sqrt((ptsx[i, 1]-ptsx[i, 2])**2 + (ptsy[i, 1]-ptsy[i, 2])**2)
            lsides[i, 2] = sqrt((ptsx[i, 2]-ptsx[i, 0])**2 + (ptsy[i, 2]-ptsy[i, 0])**2)

        # CALCULATE CELL ANGLES
        angles = np.zeros((nTri, 3))
        for i in range (0, nTri):
            p1 = np.array([ptsx[i, 0], ptsy[i, 0]])
            p2 = np.array([ptsx[i, 1], ptsy[i, 1]])
            p3 = np.array([ptsx[i, 2], ptsy[i, 2]])
            angles[i, 0] = getangle3ptsdeg(p1, p2, p3)
            angles[i, 1] = getangle3ptsdeg(p2, p3, p1)
            angles[i, 2] = getangle3ptsdeg(p3, p1, p2)

        # create dictionary to pass to neutpy
        toneutpy={}
        toneutpy["nCells"] = nTri
        toneutpy["nPlasmReg"] = pcellcount
        toneutpy["nWallSegm"] = wcellcount
        toneutpy["aneut"] = 2
        toneutpy["zion"] = 1
        toneutpy["aion"] = 2
        toneutpy["tslow"] = 0.002
        toneutpy["int_method"] = 'quad'
        toneutpy["phi_int_pts"] = 10
        toneutpy["xi_int_pts"] = 10
        toneutpy["xsec_ioni"] = 'degas'
        toneutpy["xsec_ione"] = 'degas'
        toneutpy["xsec_cx"] = 'degas'
        toneutpy["xsec_rec"] = 'degas'
        toneutpy["xsec_el"] = 'stacey_thomas'
        toneutpy["xsec_eln"] = 'stacey_thomas'
        toneutpy["refmod_e"] = 'stacey'
        toneutpy["refmod_n"] = 'stacey'

        toneutpy["iType"] = np.asarray([0]*nTri + [1]*pcellcount + [2]*wcellcount)
        toneutpy["nSides"] = np.asarray([3]*nTri + [1]*(pcellcount + wcellcount))
        toneutpy["zwall"] = np.asarray([0]*(nTri+pcellcount) + [6]*wcellcount)
        toneutpy["awall"] = np.asarray([0]*(nTri+pcellcount) + [12]*wcellcount)
        toneutpy["elecTemp"] = Te_tri[:nTri]
        toneutpy["ionTemp"] = Ti_tri[:nTri]
        toneutpy["elecDens"] = ne_tri[:nTri]
        toneutpy["ionDens"] = ni_tri[:nTri]
        toneutpy["twall"] = np.asarray([0]*nTri + [5000]*pcellcount + [0.002]*wcellcount)
        toneutpy["f_abs"] = np.asarray([0]*(nTri+pcellcount) + [0]*wcellcount)
        toneutpy["alb_s"] = np.asarray([0]*nTri + [0]*pcellcount + [0]*wcellcount)
        toneutpy["alb_t"] = np.asarray([0]*nTri + [0]*pcellcount + [0]*wcellcount)
        toneutpy["s_ext"] = np.asarray([0.0]*nTri + [0.0]*pcellcount + [0.0]*wcellcount)

        toneutpy["adjCell"] = neighbors
        toneutpy["lsides"] = lsides
        toneutpy["angles"] = angles
        toneutpy["cell1_ctr_x"] = cell1_ctr_x
        toneutpy["cell1_ctr_y"] = cell1_ctr_y
        toneutpy["cell1_theta0"] = cell1_theta0

        # write neutpy input file
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
        # f.write('\nint_method = midpoint')
        f.write('\nint_method = quad')
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

        return midpts, toneutpy
