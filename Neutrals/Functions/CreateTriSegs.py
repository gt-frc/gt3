#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

def create_tri_segs(inp, lines, pts):

    # CREATE SEGMENTS FOR TRIANGULATION
    # WHEN DOING WALL, CHECK EACH POINT TO SEE IF IT HAS ALREADY BEEN
    # CREATED. IF SO, USE THE NUMBER OF THAT POINT AND DELETE THE WALL
    # VERSION OF IT IN THE ALL_PTS ARRAY.
    sol_pol_pts = inp.core_thetapts_ntrl + inp.ib_thetapts_ntrl + inp.ob_thetapts_ntrl

    sep_segs = np.column_stack((np.arange(inp.core_thetapts_ntrl),
                                np.roll(np.arange(inp.core_thetapts_ntrl), -1)))

    ib_div_segs = np.column_stack((np.arange(inp.ib_thetapts_ntrl),
                                   np.roll(np.arange(inp.ib_thetapts_ntrl), -1)))[:-1]

    ob_div_segs = np.column_stack((np.arange(inp.ob_thetapts_ntrl),
                                   np.roll(np.arange(inp.ob_thetapts_ntrl), -1)))[:-1]

    core_segs = np.zeros((0, 2), dtype='int')
    for i, v in enumerate(lines.core):
        new_segs = np.column_stack((np.arange(inp.core_thetapts_ntrl),
                                    np.roll(np.arange(inp.core_thetapts_ntrl), -1))) \
                                    + inp.core_thetapts_ntrl * i
        core_segs = np.vstack((core_segs, new_segs))

    sol_segs = np.zeros((0, 2), dtype='int')
    for i, v in enumerate(lines.sol):
        new_segs = np.column_stack((np.arange(sol_pol_pts),
                                    np.roll(np.arange(sol_pol_pts), -1)))[:-1] \
                                    + sol_pol_pts * i
        sol_segs = np.vstack((sol_segs, new_segs))

    wall_segs = np.column_stack((np.arange(len(pts.wall)),
                                 np.roll(np.arange(len(pts.wall)), -1)))

    all_segs = np.vstack((sep_segs,
                          ib_div_segs + len(sep_segs),
                          ob_div_segs + len(ib_div_segs) + len(sep_segs) + 1,
                          core_segs + len(ob_div_segs) + len(ib_div_segs) + len(sep_segs) + 1 + 1,
                          sol_segs + len(core_segs) + len(ob_div_segs) + len(ib_div_segs) + len(sep_segs) + 1 + 1,
                          wall_segs + len(sol_segs) + len(core_segs) + len(ob_div_segs) + len(ib_div_segs) + len(sep_segs) + 1 + 1 + inp.num_sollines
                          ))

    # CLEANUP
    # NOTE: this process will result in a segments array that looks fairly chaotic,
    # but will ensure that the triangulation goes smoothly.

    all_pts = np.vstack((pts.sep,
                         pts.ib_div,
                         pts.ob_div,
                         pts.core,
                         pts.sol,
                         pts.wall))

    all_pts_unique = np.unique(all_pts, axis=0)

    # Steps:
    #   loop over each point in all_segs
    #   look up the point's coordinates in all_pts
    #   find the location of those coordinates in all_pts_unique
    #   put that location in the corresponding location in all_segs_unique

    all_segs_unique = np.zeros(all_segs.flatten().shape, dtype='int')
    for i, pt in enumerate(all_segs.flatten()):
        pt_coords = all_pts[pt]
        loc_unique = np.where((all_pts_unique == pt_coords).all(axis=1))[0][0]
        all_segs_unique[i] = loc_unique

    all_segs_unique = all_segs_unique.reshape(-1, 2)

    return all_pts_unique, all_segs_unique
