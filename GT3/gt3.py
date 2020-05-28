#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
from Neutrals.Neutrals import Neutrals
from IOL.iol import IOL
from SOL.sol import Sol
from ReadInFIle.read_infile import ReadInfile
from ImpRadiation.ImpurityRadiation import ImpRad
from Core.Core import Core
from BeamDeposition.BeamDeposition import BeamDeposition
from DensityLimit.dens_lim import DensityLimit
from Marfe.marfe import Marfe
from RadialTransport.radial_transport import RadialTransport


class gt3:

    def __init__(self, shotlabel=None, mode=None, iolFlag=True, neutFlag=True, debugRT=False):
        sys.dont_write_bytecode = True
        # Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel
        self.inp = ReadInfile(self.shotlabel)
        self.core = Core(self.inp)
        self.iolFlag = iolFlag
        self.neutFlag = neutFlag
        self.debugRT = debugRT

        if mode == 'coreonly':
            pass

        if mode == 'coreandsol':
            self.sol = Sol(self.inp, self.core)
        elif mode == 'thermaliol':
            self.iol = IOL(self.inp, self.core)
        elif mode == 'fulliol':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
        elif mode == 'imp':
            self.imp = ImpRad(core=self.core)
        elif mode == 'ntrls':
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'ntrlsandiol':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'nbi':
            self.nbi = BeamDeposition(self.inp, self.core)
        elif mode == 'marfe_denlim':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(core=self.core)
        elif mode == 'marfe':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.mar = Marfe(core=self.core)
        elif mode == 'allthethings':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = Marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            self.sol = Sol(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = RadialTransport(self.inp, self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag,
                                          debugFlag=self.debugRT)
