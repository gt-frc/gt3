#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from Test.DebugCore import debugCore
from Test.DebugInp import debugInp
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy.constants import physical_constants
import os, sys
from Functions.Beam import Beam
import json
import yaml
from pathos.multiprocessing import ProcessPool as Pool
from pathos.multiprocessing import cpu_count

m_d = physical_constants['deuteron mass'][0]
z_d = 1


class BeamDeposition:
    beam_result = None  # type: list[Beam]
    beams = None  # type: list[Beam]

    BeamResult = namedtuple('BeamResult', 'rtang beamWidth gaussR beamE beamA beamP rho Hofrho_D1 Hofrho_D2 '
                                          'Hofrho_D3 angle_D1 angle_D2 angle_D3 shine')
    BeamConfiguration = namedtuple('BeamConfiguration', 'name rho r2vol dVdr dVdrho shaf_shift kappa_vals a R0'
                                                        ' Te Ti ne Z_eff rtang beamWidth gaussR beamE beamA '
                                                        'beamP iol coCurr verbose timing iolFlag '
                                                        'beamResult pwrFracOverride')

    BeamSources = namedtuple('BeamSources', 'Snbi Qnbi Mnbi')

    def __init__(self, inp, core, iol=None, reRun=False, pwrFracOverride=None):
        sys.dont_write_bytecode = True
        # If deposition profile is provided as an input, use that.
        # Note, this is designed to read in the dP/dr(r) quantity, not dP/dV(r).
        # When integrated over r, this should equal P_beam, or at least the amount
        # that is deposited in the plasma.

        self.beams_space = np.linspace(0., core.sep_val, 50)
        """The rho space used by the Beams module on [0., 1.] or [0., sep_val]"""

        self.pwrFracOverride = pwrFracOverride
        """Allows us to override the beam power fraction for sensitivity studies"""

        if not reRun:
            # nbi_vals = calc_nbi_vals(inp, core)
            # Look for NBI output file
            try:
                self.load_beams(inp, core, iol)
            except BaseException as e:
                print "Failed to load beams data. Running beams module..."
                self.call_beams(inp, core, iol)
        else:
            self.call_beams(inp, core, iol)

        try:
            # Combine into one beam.

            self.combined_beam_src_total = self.BeamSources(np.array([a.part_src.total for a in self.beam_result]).sum(axis=0),
                                          np.array([a.en_src.total for a in self.beam_result]).sum(axis=0),
                                          np.array([a.mom_src.total for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_lost = self.BeamSources(np.array([a.part_src.lost for a in self.beam_result]).sum(axis=0),
                                          np.array([a.en_src.lost for a in self.beam_result]).sum(axis=0),
                                          np.array([a.mom_src.lost for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_kept = self.BeamSources(np.array([a.part_src.kept for a in self.beam_result]).sum(axis=0),
                                          np.array([a.en_src.kept for a in self.beam_result]).sum(axis=0),
                                          np.array([a.mom_src.kept for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_dens_total = self.BeamSources(np.array([a.part_src_dens.total for a in self.beam_result]).sum(axis=0),
                                          np.array([a.en_src_dens.total for a in self.beam_result]).sum(axis=0),
                                          np.array([a.mom_src_dens.total for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_dens_lost = self.BeamSources(np.array([a.part_src_dens.lost for a in self.beam_result]).sum(axis=0),
                                          np.array([a.en_src_dens.lost for a in self.beam_result]).sum(axis=0),
                                          np.array([a.mom_src_dens.lost for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_dens_kept = self.BeamSources(np.array([a.part_src_dens.kept for a in self.beam_result]).sum(axis=0),
                                          np.array([a.en_src_dens.kept for a in self.beam_result]).sum(axis=0),
                                          np.array([a.mom_src_dens.kept for a in self.beam_result]).sum(axis=0))
        except Exception as e:
            raise Exception("NBI module failed to create main beam profiles with error: %s" % str(e))

    def load_beams(self, inp, core, iol):
        try:
            self.beam_result = []
            with open(os.path.join(os.getcwd(), inp.beams_out_json), "r") as f:
                l = yaml.safe_load(f)
                for a in l.keys():
                    beam = l[a]
                    self.beam_result.append(Beam(self.BeamConfiguration('name',
                                                                        self.beams_space,
                                                                        core.r2vol,
                                                                        core.dVdr,
                                                                        core.dVdrho,
                                                                        core.shaf_shift,
                                                                        core.kappa_vals,
                                                                        core.a,
                                                                        core.R0_a,
                                                                        UnivariateSpline(core.rho[:, 0],
                                                                                         core.T_fsa.e.kev)(
                                                                            self.beams_space),
                                                                        UnivariateSpline(core.rho[:, 0],
                                                                                         core.T_fsa.i.kev)(
                                                                            self.beams_space),
                                                                        UnivariateSpline(core.rho[:, 0], core.n_fsa.e)(
                                                                            self.beams_space),
                                                                        UnivariateSpline(core.rho[:, 0],
                                                                                         core.z_eff_fsa)(
                                                                            self.beams_space),
                                                                        beam[0],
                                                                        beam[1],
                                                                        beam[2],
                                                                        beam[3],
                                                                        beam[4],
                                                                        beam[5],
                                                                        iol,
                                                                        coCurr=True,
                                                                        verbose=False,
                                                                        timing=False,
                                                                        iolFlag=True,
                                                                        beamResult=beam,
                                                                        pwrFracOverride=self.pwrFracOverride)))

        except IOError as e:
            raise Exception("Beams output file IO error: %s" % str(e))
        except Exception as e:
            raise Exception("Failed to load beams output file: %s" % str(e))

    def call_beams(self, inp, core, iol):
        # Call signature: def calcHofRho(rho, r2vol, dVdr, rtang, shaf_shift, kappa_vals, beamHeight, beamWidth, gaussR, Te, TC, R0, ne, beamE, Zeff, a)

        try:
            try:
                f = open(os.path.join(os.getcwd(), inp.beams_json), "r")
            except IOError as error:
                raise IOError("Could not open Beams JSON. Error: %s" % str(error))
            newdata = yaml.safe_load(f)
            f.close()
            if type(newdata) is not dict:
                raise TypeError("Beams JSON data not loaded as dictionary")
            else:
                beam_config_list = []
                for n in newdata.keys():
                    d = newdata[n]
                    if float(d['beamP']) <= 0. or float(d['beamE']) < .01:
                        continue
                    try:
                        beam_config = self.BeamConfiguration(n,
                                                             self.beams_space,
                                                             core.r2vol,
                                                             core.dVdr,
                                                             core.dVdrho,
                                                             core.shaf_shift,
                                                             core.kappa_vals,
                                                             core.a,
                                                             core.R0_g,
                                                             UnivariateSpline(core.rho[:, 0], core.T_fsa.e.kev)(
                                                                 self.beams_space),
                                                             UnivariateSpline(core.rho[:, 0], core.T_fsa.i.kev)(
                                                                 self.beams_space),
                                                             UnivariateSpline(core.rho[:, 0], core.n_fsa.e)(
                                                                 self.beams_space),
                                                             UnivariateSpline(core.rho[:, 0], core.z_eff_fsa)(
                                                                 self.beams_space),
                                                             d['rtang'],
                                                             d['beamWidth'],
                                                             d['gaussR'],
                                                             d['beamE'],
                                                             d['beamA'],
                                                             d['beamP'],
                                                             iol,
                                                             coCurr=True,
                                                             verbose=False,
                                                             timing=False,
                                                             iolFlag=True,
                                                             beamResult=None,
                                                             pwrFracOverride=self.pwrFracOverride,
                                                             )
                        beam_config_list.append(beam_config)
                    except Exception as e:
                        print "BeamConfiguration failed for beam %s with error: %s " % (n, str(e))
                pool = Pool(cpu_count() - 1)

                try:
                    try:
                        #Beam1 = Beam(beam_config_list[0])
                        #Beam2 = Beam(beam_config_list[1])
                        self.beam_result = pool.map(Beam, beam_config_list)
                    except:
                        raise Exception("Pathos failed to run beam")

                    """A list of the raw results returned by the multi-processing Pools module"""
                    if self.beam_result is None or len(self.beam_result) == 0:
                        raise Exception("Beams module failed to run beam")
                except Exception:
                    raise Exception("Beams module failed to run beam")

                # Write the beams output file
                if len(self.beam_result) > 0:
                    out_dict = {}

                    with open(inp.beams_out_json, "w") as f:
                        for b in self.beam_result:
                            out_dict[b.name] = self.BeamResult(b.rtang,
                                                               b.beamWidth,
                                                               b.gaussR,
                                                               b.beamE / 1.E3,
                                                               b.beamA,
                                                               b.beamP,
                                                               b.rho.tolist(),
                                                               b.Hofrho.D1(b.rho).tolist(),
                                                               b.Hofrho.D2(b.rho).tolist(),
                                                               b.Hofrho.D3(b.rho).tolist(),
                                                               b.angle_1(b.rho).tolist(),
                                                               b.angle_2(b.rho).tolist(),
                                                               b.angle_3(b.rho).tolist(),
                                                               b.shine)
                        json.dump(out_dict, f, indent=4)

        except:
            print "No Beams JSON file found. Using single-beam data"

    def _plot_single_val(self, title="Title", xLabel=r'$\rho$', yLabel="Y-Axis"):

        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        fig1.set_xlabel(xLabel)
        fig1.set_ylabel(yLabel)
        fig1.set_title(title)
        return fig1

    def plot_S_nbi(self):
        fig = self._plot_single_val(title="NBI Particle Source", yLabel=r'$S_{nbi}[\#/{m s}]$')
        fig.scatter(self.beams_space, self.combined_beam_src_kept.Snbi, color='red')
        plt.show()
        return fig

    def plot_Q_nbi(self):
        fig = self._plot_single_val(title="NBI Energy Source", yLabel=r'$Q_{nbi}[W/m]$')
        fig.scatter(self.beams_space, self.combined_beam_src_kept.Qnbi, color='red')
        plt.show()
        return fig

    def plot_M_nbi(self):
        fig = self._plot_single_val(title="NBI Momentum Source", yLabel=r'$S_{nbi}[M/s}]$')
        fig.scatter(self.beams_space, self.combined_beam_src_kept.Qnbi, color='red')
        plt.show()
        return fig

    def plot_S_dens_nbi(self):
        fig = self._plot_single_val(title="NBI Particle Source Density", yLabel=r'$S_{nbi}[\#/{m^3 s}]$')
        fig.scatter(self.beams_space, self.combined_beam_src_dens_kept.Snbi, color='red')
        plt.show()
        return fig

    def plot_Q_dens_nbi(self):
        fig = self._plot_single_val(title="NBI Energy Source Density", yLabel=r'$Q_{nbi}[W/{m^3}]$')
        fig.scatter(self.beams_space, self.combined_beam_src_dens_kept.Qnbi, color='red')
        plt.show()
        return fig

    def plot_M_dens_nbi(self):
        fig = self._plot_single_val(title="NBI Momentum Source Density", yLabel=r'$S_{nbi}[\frac{M}{s}]$')
        fig.scatter(self.beams_space, self.combined_beam_src_dens_kept.Qnbi, color='red')
        plt.show()
        return fig



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    inp = debugInp()
    core = debugCore()
    beamData = BeamDeposition(inp, core)
    print """Printing NBeams input parameters:
            Core
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            """.format("R0_a", str(core.R0_a),
                       "a", str(core.a),
                       "kappa_axis", str(core.kappa.axis),
                       "kappa_sep", str(core.kappa.sep),
                       "pts.axis.mag[0]", str(core.pts.axis.mag[0]),
                       "shaf_shift", str(core.shaf_shift),
                       )

    print """
            inp           
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            """.format("BT0", str(inp.BT0),
                       "ebeam", str(inp.ebeam),
                       "nbeams_loc", str(inp.nbeams_loc),
                       "pbeam", str(inp.pbeam),
                       "rtang", str(inp.rtang)
                       )

    print """
        D1 info
        {:<10}      {}
        {:<10}      {}
        {:<10}      {}
        """.format("m", str(beamData.beams.D1.m),
                   "rtan", str(beamData.beams.D1.rtan),
                   "zeta", str(beamData.beams.D1.zeta))
    print """
        D2 info
        {:<10}      {}
        {:<10}      {}
        {:<10}      {}
        """.format("m", str(beamData.beams.D2.m),
                   "rtan", str(beamData.beams.D2.rtan),
                   "zeta", str(beamData.beams.D2.zeta))
    print """
        D3 info
        {:<10}      {}
        {:<10}      {}
        {:<10}      {}
        """.format("m", str(beamData.beams.D3.m),
                   "rtan", str(beamData.beams.D3.rtan),
                   "zeta", str(beamData.beams.D3.zeta))

    print """Calculated nbeams values"""
    fig0 = plt.figure(figsize=(12, 8))
    fig0.suptitle(r"Deposition profiles")
    axi = 230
    for a in beamData.debug.keys():
        if type(beamData.debug[a]) is float:
            print "{:<10}      {}".format(str(a), str(beamData.debug[a]))
        elif type(beamData.debug[a]) is np.ndarray:
            if "hofr" in a:
                axi += 1
                ax = fig0.add_subplot(axi)
                ax.set_title(a, fontsize=16)
                ax.set_ylabel("dep", fontsize=16)
                ax.set_xlabel(r"rho", fontsize=16)
                ax.plot(np.linspace(0, 1, 51), beamData.debug[a])
    d1Ptot = UnivariateSpline(core.rho[:, 0], beamData.beams.D1.dPdr.v1D.W).integral(0., 1.)
    d2Ptot = UnivariateSpline(core.rho[:, 0], beamData.beams.D2.dPdr.v1D.W).integral(0., 1.)
    d3Ptot = UnivariateSpline(core.rho[:, 0], beamData.beams.D3.dPdr.v1D.W).integral(0., 1.)
    Ptot = d1Ptot + d2Ptot + d3Ptot
    print "{:<10}      {}".format("D1 Total power", str(d1Ptot))
    print "{:<10}      {}".format("D2 Total power", str(d2Ptot))
    print "{:<10}      {}".format("D3 Total power", str(d3Ptot))
    print "{:<10}      {}".format("Total power from dpdr graph integration", str(Ptot))

    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(r'NBeams Debug Info')
    ax1 = fig.add_subplot(231)
    ax1.set_title(r'ion density', fontsize=16)
    ax1.set_ylabel(r'$n_i[1/{m^3}]$', fontsize=16)
    ax1.set_xlabel(r'rho', fontsize=16)
    ax1.plot(core.rho[:, 0], core.n.i[:, 0])

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'electron density', fontsize=16)
    ax2.set_ylabel(r'$n_e[1/{m^3}]$', fontsize=16)
    ax2.set_xlabel(r'rho', fontsize=16)
    ax2.plot(core.rho[:, 0], core.n.e[:, 0])

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'Ion temp (kev)', fontsize=16)
    ax3.set_ylabel(r'$T_i[kev}$', fontsize=16)
    ax3.set_xlabel(r'rho', fontsize=16)
    ax3.plot(core.rho[:, 0], core.T.i.kev[:, 0])

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'Electron temp (kev)', fontsize=16)
    ax4.set_ylabel(r'$T_e[kev}$', fontsize=16)
    ax4.set_xlabel(r'rho', fontsize=16)
    ax4.plot(core.rho[:, 0], core.T.e.kev[:, 0])

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'Electron temp (kev)', fontsize=16)
    ax5.set_ylabel(r'$T_e[keV]$', fontsize=16)
    ax5.set_xlabel(r'rho', fontsize=16)
    ax5.plot(core.rho[:, 0], core.T.e.kev[:, 0])

    fig2 = plt.figure(figsize=(12, 8))
    fig2.suptitle(r"""$D_1$ data""")

    ax21 = fig2.add_subplot(231)
    ax21.set_title(r'Power density profile', fontsize=16)
    ax21.set_ylabel(r'$\frac{dP}{dr}$', fontsize=16)
    ax21.set_xlabel(r'rho', fontsize=16)
    ax21.plot(core.rho[:, 0], beamData.beams.D1.dPdr.v1D.W)

    ax22 = fig2.add_subplot(232)
    ax22.set_title(r'Power density profile', fontsize=16)
    ax22.set_ylabel(r'$\frac{dP}{dV}$', fontsize=16)
    ax22.set_xlabel(r'rho', fontsize=16)
    ax22.plot(core.rho[:, 0], beamData.beams.D1.dPdV.v1D.W)

    fig3 = plt.figure(figsize=(12, 8))
    fig3.suptitle(r"""$D_2$ data""")

    ax31 = fig3.add_subplot(231)
    ax31.set_title(r'Power density profile', fontsize=16)
    ax31.set_ylabel(r'$\frac{dP}{dr}$', fontsize=16)
    ax31.set_xlabel(r'rho', fontsize=16)
    ax31.plot(core.rho[:, 0], beamData.beams.D2.dPdr.v1D.W)

    ax32 = fig3.add_subplot(232)
    ax32.set_title(r'Power density profile', fontsize=16)
    ax32.set_ylabel(r'$\frac{dP}{dV}$', fontsize=16)
    ax32.set_xlabel(r'rho', fontsize=16)
    ax32.plot(core.rho[:, 0], beamData.beams.D2.dPdV.v1D.W)

    fig4 = plt.figure(figsize=(12, 8))
    fig4.suptitle(r"""$D_3$ data""")

    ax41 = fig4.add_subplot(231)
    ax41.set_title(r'Power density profile', fontsize=16)
    ax41.set_ylabel(r'$\frac{dP}{dr}$', fontsize=16)
    ax41.set_xlabel(r'rho', fontsize=16)
    ax41.plot(core.rho[:, 0], beamData.beams.D3.dPdr.v1D.W)

    ax42 = fig4.add_subplot(232)
    ax42.set_title(r'Power density profile', fontsize=16)
    ax42.set_ylabel(r'$\frac{dP}{dV}$', fontsize=16)
    ax42.set_xlabel(r'rho', fontsize=16)
    ax42.plot(core.rho[:, 0], beamData.beams.D3.dPdV.v1D.W)

    ax43 = fig4.add_subplot(233)
    ax43.set_title(r'Power density profile', fontsize=16)
    ax43.set_ylabel(r'$\frac{dP}{dr}$', fontsize=16)
    ax43.set_xlabel(r'rho', fontsize=16)
    ax43.plot(beamData.debug['rho_nbeams'], beamData.debug['dvol'])

    plt.show()