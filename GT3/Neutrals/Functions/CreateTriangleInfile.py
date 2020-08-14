#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import os


def create_triangle_infile(all_pts_unique, all_segs_unique, core):
    # OUTPUT .poly FILE AND RUN TRIANGLE PROGRAM
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
               np.array([1, core.pts.axis.mag[0], core.pts.axis.mag[1]])[None],
               fmt='%i %f %f')
    np.savetxt(outfile,
               np.array([0])[None],
               fmt='%i')
    outfile.close()

    return filepath
