#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import sys
from GT3.IOL import IOL
from GT3.SOL import Sol
from GT3.ReadInFIle import ReadInfile
from GT3.ImpRadiation import ImpRad
from GT3.Core import Core
from GT3.BeamDeposition import BeamDeposition
from GT3.DensityLimit import DensityLimit
from GT3.Marfe import Marfe
from GT3.RadialTransport import RadialTransport

try:
    from GT3.Neutrals import Neutrals
except ImportError:
    pass
class gt3:

    debug=False

    def __init__(self, inputFile=None, preparedInput = None, mode="coreonly", **kwargs):
        sys.dont_write_bytecode = True
        # Create shotlabel as an attribute of plasma class
        self.iolFlag = kwargs.get("iolFlag", True)
        self.neutFlag = kwargs.get('neutFlag', True)
        self.verbose = kwargs.get('verbose', False)
        self.debug = True if kwargs.get('debug') else False

        if inputFile:
            self.inputFile = inputFile
        if preparedInput:
            self.inp = preparedInput
        else:
            self.inp = ReadInfile(self.inputFile, **kwargs)
        self.core = Core(self.inp, debug=self.debug, **kwargs)
        self.beamPowerFracOverride = None
        self.ntrl_cpu_override = False

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
            self.iol = IOL(self.inp, self.core, **kwargs)
        elif mode == 'fulliol':
            self.iol = IOL(self.inp, self.core, **kwargs)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
        elif mode == 'imp':
            self.imp = ImpRad(core=self.core)
        elif mode == 'ntrls':
            self._run_neutpy()
        elif mode == 'ntrlsandiol':
            self.iol = IOL(self.inp, self.core, **kwargs)
            self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
        elif mode == 'nbi':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core, **kwargs)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core, pwrFracOverride=self.beamPowerFracOverride)
        elif mode == 'marfe_denlim':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core, **kwargs)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core, pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(core=self.core)
        elif mode == 'marfe':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core, **kwargs)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol, pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
            self.imp = ImpRad(core=self.core)
            self.mar = Marfe(core=self.core)
        elif mode == 'allthethings':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core, **kwargs)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol,  pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
            self._run_neutpy()
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            if self.iolFlag:
                self.iol = IOL(self.inp, self.core, **kwargs)
                self.nbi = BeamDeposition(self.inp, self.core, self.iol,  pwrFracOverride=self.beamPowerFracOverride)
            else:
                self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
            self.sol = Sol(self.inp, self.core)
            self._run_neutpy()
            self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag)

    def _run_neutpy(self, reRun=False, **kwargs):
        if self.neutpyLoaded:
            self.ntrl = Neutrals(self.inp, self.core, cpus=self.ntrl_cpu_override, **kwargs)
            if reRun:
                self.ntrl.reRun(cpus=self.ntrl_cpu_override,  **kwargs)
        else:
            print("NeutPy is not loaded. Cannot run Neutrals calculation")


    def override_NBI_Pwrfrac(self, frac):
        if isinstance(frac, list):
            self.beamPowerFracOverride = frac
        else:
            print("Please provide the NBI power fraction override as a list")

    def run_SOL(self):
        self.sol = Sol(self.inp, self.core)
        return self

    def run_IOL(self, *args, **kwargs):
        self.iol = IOL(self.inp, self.core, **kwargs)
        return self

    def run_NBI(self, reRun=False, **kwargs):
        try:
            self.iol
        except AttributeError:
            print ("IOL module not run. Running now...")
            self.run_IOL(**kwargs)
        if self.iolFlag:
            self.nbi = BeamDeposition(self.inp, self.core, self.iol, reRun=reRun,
                                      pwrFracOverride=self.beamPowerFracOverride)
        else:
            self.nbi = BeamDeposition(self.inp, self.core,  pwrFracOverride=self.beamPowerFracOverride)
        return self

    def run_impurities(self):
        self.imp = ImpRad(core=self.core)
        return self

    def run_neutrals(self, reRun=False, **kwargs):
        self._run_neutpy(reRun=reRun, **kwargs)
        return self

    def override_ntrl_cpus(self, num):
        self.ntrl_cpu_override = num
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

    def run_radial_transport(self, nbiReRun=False, ntrlReRun=False, debug=False, *args, **kwargs):
        if kwargs.get("neutpy_iterate"):
            ntrlReRun=True
        if self.iolFlag:
            try:
                self.iol
            except AttributeError:
                print("IOL module not run. Running now...")
                self.run_IOL(**kwargs)
        try:
            self.nbi
        except AttributeError:
            print("NBI module not run. Running now...")
            self.run_NBI(reRun=nbiReRun)

        try:
            self.imp
        except AttributeError:
            print ("Impurity radiation module not run. Running now...")
            self.imp = ImpRad(z=6, core=self.core, neutFlag=self.neutFlag)

        if self.neutFlag:
            try:
                self.ntrl
            except AttributeError:
                print ("Neutrals module not run. Running now...")
                self.run_neutrals(reRun=ntrlReRun, **kwargs)

        self.rtrans = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag, **kwargs)
        if self.rtrans.reRun:
            kwargs['neutpy_D'] = self.rtrans.D_i.val[-1]
            kwargs['neutpy_chi'] = self.rtrans.chi.i.chi0[-1]
            kwargs['neutpy_q'] = self.rtrans.Q.D.diff.val[-1]
            kwargs['neutpy_gamma'] = self.rtrans.gamma.D.diff.val[-1]
            print(kwargs)
            self.run_neutrals(reRun=True, **kwargs)
            self.rtrans_iter = RadialTransport(self.core, self.iol, self.nbi, self.iolFlag, self.neutFlag, **kwargs)
            return self
        else:
            return self


    def disable_IOL(self):
        self.iolFlag = False
        return self

    def disable_neutrals(self):
        self.neutFlag = False
        self.core.izn_rate.t.set_to_zeros()
        self.core.izn_rate.tot.set_to_zeros()
        self.core.izn_rate.s.set_to_zeros()
        self.core.n.n.t.set_to_zeros()
        self.core.n.n.s.set_to_zeros()
        self.core.n.n.tot.set_to_zeros()
        return self