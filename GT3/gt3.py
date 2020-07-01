#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
from Neutrals.Neutrals import Neutrals
from IOL.iol import IOL
from SOL.sol import Sol
from ReadInFIle.read_infile import ReadInfile
from ImpRadiation.ImpurityRadiation import ImpRad
from Core.Core import Core
from BeamDeposition import BeamDeposition
from DensityLimit.dens_lim import DensityLimit
from Marfe.marfe import Marfe
from RadialTransport.radial_transport import RadialTransport


class gt3:

    def __init__(self, preparedInput = None, shotlabel=None, mode=None, iolFlag=True, neutFlag=True, verbose=False):
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
        self.verbose = verbose

        if mode == 'coreonly':
            pass
        if mode == 'coreandsol':
            self.sol = Sol(self.inp, self.core)
        elif mode == 'thermaliol':
            self.iol = IOL(self.inp, self.core)
        elif mode == 'fulliol':
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol)
        elif mode == 'imp':
            self.imp = ImpRad(core=self.core)
        elif mode == 'ntrls':
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'ntrlsandiol':
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol)
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
            self.nbi = BeamDeposition(self.inp, self.core, self.iol)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.mar = Marfe(core=self.core)
        elif mode == 'allthethings':
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(self.inp, self.core)
        elif mode == 'radialtrans':
            self.sol = Sol(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag,
                                          debugFlag=self.verbose)

    def run_SOL(self):
        self.sol = Sol(self.inp, self.core)
        return self

    def run_IOL(self):
        self.iol = IOL(self.inp, self.core)
        return self

    def run_NBI(self):
        self.iol = IOL(self.inp, self.core)
        self.nbi = BeamDeposition(self.inp, self.core, self.iol)
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
        try:
            self.iol
        except AttributeError:
            print ("IOL module not run. Running now...")
            self.run_IOL()
        try:

            self.nbi
        except AttributeError:
            print ("MBI module not run. Running now...")
            self.run_NBI()

        self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag,
                                      debugFlag=self.verbose)
        return self
