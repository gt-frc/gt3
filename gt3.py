#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
"""
import sys
from exp_neutpy_prep import Neutrals
from iol import IOL
from read_infile import ReadInfile
from imp_rad import ImpRad
from exp_core_brnd import ExpCoreBrnd
from mil_core_brnd import MilCoreBrnd
from beamdep import BeamDeposition
from dens_lim import dens_lim
from marfe import Marfe
from radial_transport import RadialTransport


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
        self.inp = ReadInfile(self.shotlabel)
        self.core = ExpCoreBrnd(self.inp) if self.inp.exp_inp else MilCoreBrnd(self.inp)

        if mode == 'coreonly':
            pass
        elif mode == 'thermaliol':
            self.iol = IOL(self.inp, self.core)
        elif mode == 'fulliol':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
        elif mode == 'imp':
            self.imp = ImpRad(self.inp, self.core)
        elif mode == 'ntrls':
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'ntrlsandiol':
            self.ntrl = Neutrals(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
        elif mode == 'nbi':
            self.nbi = BeamDeposition(self.inp, self.core)
        elif mode == 'therm_instab':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.dl = dens_lim(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = Marfe(self.inp, self.core, self.imp)
        elif mode == 'allthethings':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.dl = dens_lim(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = Marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            #self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = RadialTransport(self.inp, self.core, self.iol, self.nbi)
