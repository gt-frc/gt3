#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from GT3.IOL import IOL
from GT3.SOL import Sol
from GT3.ReadInFIle import ReadInfile
from GT3.ImpRadiation import ImpRad
from GT3.Core import Core
from GT3.BeamDeposition import BeamDeposition
from GT3.DensityLimit import DensityLimit
from GT3.Marfe import Marfe
from RadialTransport.radial_transport import RadialTransport

try:
    from GT3.Neutrals import Neutrals
except ImportError:
    pass

class gt3:

    def __init__(self, inputFile=None, preparedInput = None, mode="coreonly", **kwargs):
        sys.dont_write_bytecode = True
        # Create shotlabel as an attribute of plasma class
        if "iolFlag" in kwargs:
            iolFlag = kwargs['iolFlag']
        else:
            iolFlag = True
        if "neutFlag" in kwargs:
            neutFlag = kwargs['neutFlag']
        else:
            neutFlag = True
        if "verbose" in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = False
        if inputFile:
            self.inputFile = inputFile
        if preparedInput:
            self.inp = preparedInput
        else:
            self.inp = ReadInfile(self.inputFile)
        self.core = Core(self.inp)
        self.iolFlag = iolFlag
        self.neutFlag = neutFlag
        self.verbose = verbose
        self.beamPowerFracOverride = None

        try:
            import neutpy
            self.neutpyLoaded = True
        except ImportError:
            self.neutpyLoaded = False

        if mode == 'coreonly':
            pass


        if mode == 'coreandsol':
            self.sol = Sol(self.inp, self.core)
        elif mode == 'thermaliol':
            self.iol = IOL(self.inp, self.core)
        elif mode == 'fulliol':
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
        elif mode == 'imp':
            self.imp = ImpRad(core=self.core)
        elif mode == 'ntrls':
            self._run_neutpy()
        elif mode == 'ntrlsandiol':
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
        elif mode == 'nbi':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core, pwrFracOverride=self.beamPowerFracOverride)
        elif mode == 'marfe_denlim':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core, pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(core=self.core)
        elif mode == 'marfe':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
            self.imp = ImpRad(core=self.core)
            self.mar = Marfe(core=self.core)
        elif mode == 'allthethings':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol,  pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol,  pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
            self.sol = Sol(self.inp, self.core)
            self._run_neutpy()
            self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag)

    def _run_neutpy(self):
        if self.neutpyLoaded:
            self.ntrl = Neutrals(self.inp, self.core)
        else:
            print "NeutPy is not loaded. Cannot run Neutrals calculation"

    def override_NBI_Pwrfrac(self, frac):
        if isinstance(frac, list):
            self.beamPowerFracOverride = frac
        else:
            print "Please provide the NBI power fraction override as a list"

    def run_SOL(self):
        self.sol = Sol(self.inp, self.core)
        return self

    def run_IOL(self):
        self.iol = IOL(self.inp, self.core)
        return self

    def run_NBI(self, reRun=False):
        try:
            self.iol
        except AttributeError:
            print ("IOL module not run. Running now...")
            self.run_IOL()
        if self.iolFlag:
            self.nbi = BeamDeposition(self.inp, self.core, self.iol, reRun=reRun,
                                      pwrFracOverride=self.beamPowerFracOverride)
        else:
            self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
        return self

    def run_impurities(self):
        self.imp = ImpRad(core=self.core)
        return self

    def run_neutrals(self):
        self._run_neutpy()
        return self

    def run_density_limit(self):
        if self.iolFlag:
            self.nbi = BeamDeposition(self.inp, self.core, self.iol,  pwrFracOverride=self.beamPowerFracOverride)
        else:
            self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
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
            print ("NBI module not run. Running now...")
            self.run_NBI()
        try:
            self.ntrl
        except AttributeError:
            print ("Neutrals module not run. Running now...")
            self.run_neutrals()

        try:
            self.imp
        except AttributeError:
            print ("Impurity radiation module not run. Running now...")
            self.imp = ImpRad(z=6, core=self.core)
        self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag)
        return self

    def disable_IOL(self):
        self.iolFlag = False
        print ("Re-running Radial Transport without IOL")
        try:
            self.rtrans
        except:
            self.run_radial_transport()
        return self
