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
import yaml
from GT3.Core.Processors import NumpyEncoder



class Neutrals:
    def __init__(self, inp, core):
        # Try to read in specified neutrals data file. If it's not there, then prepare inputs for and run neutpy
        NeutralData = namedtuple('NeutralsData', 'R Z n_n_slow n_n_thermal izn_rate_slow izn_rate_thermal')

        try:
            with open(inp.neutfile_loc, "r") as f:
                ntrl_data = yaml.safe_load(f)
                self.data = NeutralData(ntrl_data['R'],
                                        ntrl_data['Z'],
                                        ntrl_data['nn_s_raw'],
                                        ntrl_data['nn_t_raw'],
                                        ntrl_data['izn_rate_slow'],
                                        ntrl_data['izn_rate_thermal'])
        except:
            # Run NeutPy
            print "Neutrals data not found. Running NeutPy"
            npi = neutrals()
            npi.from_gt3(core, inp)
            self.data = NeutralData(npi.midpts[:, 0],
                                    npi.midpts[:, 1],
                                    npi.nn_s_raw,
                                    npi.nn_t_raw,
                                    npi.iznrate_s_raw,
                                    npi.iznrate_t_raw)
            # Save data
            try:
                out_dict = {}
                with open(inp.neutfile_loc, "w") as f:
                    out_dict['R'] = self.data.R
                    out_dict['Z'] = self.data.Z
                    out_dict['nn_s_raw'] = self.data.n_n_slow
                    out_dict['nn_t_raw'] = self.data.n_n_thermal
                    out_dict['izn_rate_slow'] = self.data.izn_rate_slow
                    out_dict['izn_rate_thermal'] = self.data.izn_rate_thermal
                    json.dump(out_dict, f, indent=4, cls=NumpyEncoder)
            except:
                print "Unable to save NeutPy data to file"
        try:
            core.update_ntrl_data(self.data)
        except:
            print 'unable to update values in core instance.'
            pass
