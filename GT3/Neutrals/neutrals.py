#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 14:27:42 2018

@author: max
"""
from __future__ import division
from neutpy import neutrals
from collections import namedtuple
import json
from GT3.Core.Processors import NumpyEncoder



class Neutrals:

    def __init__(self, inp, core):

        if abs(1.0 - core.sep_val) > .0001:
            print "The separatrix value has been overwritten. Cannot run Neutrals calculation"
            return
        # Try to read in specified neutrals data file. If it's not there, then prepare inputs for and run neutpy
        self.NeutralDataNT = namedtuple('NeutralsData', 'R Z n_n_slow n_n_thermal izn_rate_slow izn_rate_thermal')
        self.inp = inp
        self.core = core
        try:
            with open(inp.neutfile_loc, "r") as f:
                ntrl_data = json.load(f)
                self.data = self.NeutralDataNT(ntrl_data['R'],
                                        ntrl_data['Z'],
                                        ntrl_data['nn_s_raw'],
                                        ntrl_data['nn_t_raw'],
                                        ntrl_data['izn_rate_slow'],
                                        ntrl_data['izn_rate_thermal'])
        except:
            # Run NeutPy
            print "Neutrals data not found. Running NeutPy"
            self.npi = neutrals()
            self.npi.from_gt3(core, inp)
            self.data = self.NeutralDataNT(self.npi.midpts[:, 0],
                                    self.npi.midpts[:, 1],
                                    self.npi.nn_s_raw,
                                    self.npi.nn_t_raw,
                                    self.npi.iznrate_s_raw,
                                    self.npi.iznrate_t_raw)
            # Save data
            self._save_data()
        try:
            core.update_ntrl_data(self.data)
        except:
            print 'unable to update values in core instance.'
            pass

    def reRun(self):
        print ("Manually re-running NeutPy")
        self.npi = neutrals()
        self.npi.from_gt3(self.core, self.inp)
        self.data = self.NeutralDataNT(self.npi.midpts[:, 0],
                                self.npi.midpts[:, 1],
                                self.npi.nn_s_raw,
                                self.npi.nn_t_raw,
                                self.npi.iznrate_s_raw,
                                self.npi.iznrate_t_raw)
        self._save_data()

    def _save_data(self):
        try:
            out_dict = {}
            with open(self.inp.neutfile_loc, "w") as f:
                out_dict['R'] = self.data.R
                out_dict['Z'] = self.data.Z
                out_dict['nn_s_raw'] = self.data.n_n_slow
                out_dict['nn_t_raw'] = self.data.n_n_thermal
                out_dict['izn_rate_slow'] = self.data.izn_rate_slow
                out_dict['izn_rate_thermal'] = self.data.izn_rate_thermal
                json.dump(out_dict, f, indent=4, cls=NumpyEncoder)
        except Exception as e:
            print "Unable to save NeutPy data to file: %s" % str(e)


