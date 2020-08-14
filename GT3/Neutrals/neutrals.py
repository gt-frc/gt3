#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 14:27:42 2018

@author: max
"""
from __future__ import division
import numpy as np
from neutpy import neutrals
from collections import namedtuple



class Neutrals:
    def __init__(self, inp, core):
        # Try to read in specified neutrals data file. If it's not there, then prepare inputs for and run neutpy
        NeutralData = namedtuple('NeutralsData', 'R Z n_n_slow n_n_thermal izn_rate_slow izn_rate_thermal')

        try:
            ntrl_data = np.loadtxt(inp.neutfile_loc, delimiter=',', skiprows=1)
            self.data = NeutralData(ntrl_data[:, 1],
                ntrl_data[:, 2],
                ntrl_data[:, 3],
                ntrl_data[:, 4],
                ntrl_data[:, 6],
                ntrl_data[:, 7])
        except:
            # Run NeutPy
            npi = neutrals()
            npi.from_gt3(core, inp)
            self.data = NeutralData(npi.midpts[:, 0],
                                    npi.midpts[:, 1],
                                    npi.nn_s_raw,
                                    npi.nn_t_raw,
                                    npi.iznrate_s_raw,
                                    npi.iznrate_t_raw)
        try:
            core.update_ntrl_data(self.data)
        except:
            print 'unable to update values in core instance.'
            pass
