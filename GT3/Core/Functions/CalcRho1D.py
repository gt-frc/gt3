#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np


def calc_rho1d(edge_rho=None, rhopts_core=None, rhopts_edge=None, rhopts=None):
    # define rho points
    try:
        rho1d = np.concatenate((np.linspace(0, edge_rho, rhopts_core, endpoint=False),
                                np.linspace(edge_rho, 1, rhopts_edge, endpoint=True)))
    except:
        try:
            rho1d = np.linspace(0, 1, rhopts)
        except:
            print 'rho parameters not defined. Using 100 evenly spaced rho values'
            rho1d = np.linspace(0, 1, 100)
    return rho1d
