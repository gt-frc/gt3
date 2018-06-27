#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
"""
import sys
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from exp_neutpy_prep import Neutrals
from thermaliol import iol_calc
from read_infile import read_infile
from imp_rad import ImpRad
from exp_core_brnd import exp_core_brnd
from exp_sol_brnd import exp_sol_brnd
from exp_pfr_brnd import exp_pfr_brnd
from mil_core_brnd import mil_core_brnd
from mil_sol_brnd import mil_sol_brnd
from mil_pfr_brnd import mil_pfr_brnd
from beamdep import beamdep
from dens_lim import dens_lim
from marfe import Marfe
from rad_trans import rad_trans


class gt3:
    """GT3 calculates various tokamak-related quantities

    Methods:
        allthethings
        justneutrals
        plotstuff

    Attributes:

    External Dependencies:
        Triangle:       Used to create the triangular mesh for neutpy. Source
                        code and documentation can be found at
                        https://www.cs.cmu.edu/~quake/triangle.html. Can be
                        installed from Ubuntu repositories as well.
        Neutpy:
        nbeams:
        adpack:
    """
    def __init__(self, shotlabel=None, mode=None):
        sys.dont_write_bytecode = True 
        # Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel
        self.inp = read_infile(self.shotlabel)
        self.core = exp_core_brnd(self.inp) if self.inp.exp_inp else mil_core_brnd(self.inp)

        if mode == 'coreonly':
            pass
        elif mode == 'thermaliol':
            self.iol = iol_calc(self.inp, self.core)
        elif mode == 'fulliol':
            self.nbi = beamdep(self.inp, self.core)
            self.iol = iol_calc(self.inp, self.core)
        elif mode == 'imp':
            self.imp = ImpRad(self.inp, self.core)
        elif mode == 'ntrls':
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'ntrlsandiol':
            self.ntrl = Neutrals(self.inp, self.core)
            self.iol = iol_calc(self.inp, self.core)
        elif mode == 'nbi':
            self.nbi = beamdep(self.inp, self.core)
        elif mode == 'therm_instab':
            self.nbi = beamdep(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.dl = dens_lim(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = marfe(self.inp, self.core, self.imp)
        elif mode == 'allthethings':
            self.nbi = beamdep(self.inp, self.core)
            self.iol = iol_calc(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.dl = dens_lim(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            self.iol = iol_calc(self.inp, self.core)
            self.nbi = beamdep(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = rad_trans(self.inp,self.core,self.iol,self.ntrl,self.nbi)
