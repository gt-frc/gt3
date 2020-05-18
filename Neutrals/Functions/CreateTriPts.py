#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from shapely.geometry import LinearRing
from collections import namedtuple
import numpy as np

def create_tri_pts(inp, lines):
    sol_pol_pts = inp.core_thetapts_ntrl + inp.ib_thetapts_ntrl + inp.ob_thetapts_ntrl

    # GET POINTS FOR TRIANGULATION
    # main seperatrix
    sep_pts = np.zeros((inp.core_thetapts_ntrl, 2))
    for i, v in enumerate(np.linspace(0, 1, inp.core_thetapts_ntrl, endpoint=False)):
        sep_pts[i] = np.asarray(lines.sep.interpolate(v, normalized=True).xy).T[0]

    # inboard divertor leg
    ib_div_pts = np.zeros((inp.ib_thetapts_ntrl, 2))
    for i, v in enumerate(np.linspace(0, 1, inp.ib_thetapts_ntrl, endpoint=True)):  # skipping the x-point (point 0)
        ib_div_pts[i] = np.asarray(lines.ib_div.interpolate(v, normalized=True).xy).T[0]

    # outboard divertor leg
    ob_div_pts = np.zeros((inp.ob_thetapts_ntrl, 2))
    for i, v in enumerate(np.linspace(0, 1, inp.ob_thetapts_ntrl, endpoint=True)):  # skipping the x-point (point 0)
        ob_div_pts[i] = np.asarray(lines.ob_div.interpolate(v, normalized=True).xy).T[0]

    # core
    core_pts = np.zeros((inp.core_thetapts_ntrl*len(lines.core), 2))
    for num, line in enumerate(lines.core):
        for i, v in enumerate(np.linspace(0, 1, inp.core_thetapts_ntrl, endpoint=False)):
            core_pts[num*inp.core_thetapts_ntrl + i] = np.asarray(line.interpolate(v, normalized=True).xy).T[0]

    core_ring = LinearRing(core_pts[:inp.core_thetapts_ntrl])

    # sol
    sol_pts = np.zeros((sol_pol_pts*len(lines.sol), 2))
    for num, line in enumerate(lines.sol):
        for i, v in enumerate(np.linspace(0, 1, sol_pol_pts, endpoint=True)):
            sol_pts[num*sol_pol_pts + i] = np.asarray(line.interpolate(v, normalized=True).xy).T[0]

    # wall
    wall_pts = np.asarray(inp.wall_line.coords)[:-1]

    pts_dict = {}
    pts_dict['core'] = core_pts
    pts_dict['sep'] = sep_pts
    pts_dict['sol'] = sol_pts
    #pts_dict['pfr'] = pfr_pts
    pts_dict['ib_div'] = ib_div_pts
    pts_dict['ob_div'] = ob_div_pts
    pts_dict['wall'] = wall_pts
    pts = namedtuple('pts', pts_dict.keys())(*pts_dict.values())

    return pts, core_ring