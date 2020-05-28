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

    def __init__(self, preparedInput = None, shotlabel=None, mode=None, iolFlag=True, neutFlag=True, debugRT=False):
        sys.dont_write_bytecode = True
        # Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel
        if preparedInput:
            self.inp = preparedInput
        else:
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
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            self.sol = Sol(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag,
                                          debugFlag=self.debugRT)

    def run_SOL(self):
        self.sol = Sol(self.inp, self.core)
        return self

    def run_IOL(self):
        self.iol = IOL(self.inp, self.core)
        return self

    def run_NBI(self):
        self.nbi = BeamDeposition(self.inp, self.core)
        return self

    def run_impurities(self):
        self.imp = ImpRad(core=self.core)
        return self

    def run_neutrals(self):
        self.ntrl = Neutrals(self.inp, self.core)
        return self

    def run_density_limit(self):
        self.dl = DensityLimit(self.core, self.nbi)
        return self

    def run_marf(self):
        self.mar = Marfe(self.inp, self.core)
        return self

    def run_radial_transport(self):
        if not self.iol: self.run_IOL()
        if not self.nbi: self.run_NBI()
        self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag,
                                      debugFlag=self.debugRT)
        return self
