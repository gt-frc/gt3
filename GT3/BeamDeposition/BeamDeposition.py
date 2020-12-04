#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d
from collections import namedtuple
import os, sys
from .Functions.Beam import Beam
import json
import yaml
from pathos.multiprocessing import ProcessPool as Pool
from pathos.multiprocessing import cpu_count
import GT3.constants as constants
from GT3.utilities.PlotBase import PlotBase
from GT3 import Core
from GT3.Core.Functions.ProfileClasses import OneDProfile

m_d = constants.deuteron_mass
z_d = 1


class BeamDeposition(PlotBase):
    beam_result = None  # type: list
    beams = None  # type: list

    BeamResult = namedtuple('BeamResult', 'rtang beamWidth gaussR beamE beamA beamP rho Hofrho_D1 Hofrho_D2 '
                                          'Hofrho_D3 angle_D1 angle_D2 angle_D3 shine')
    BeamConfiguration = namedtuple('BeamConfiguration', 'name rho r2vol dVdr dVdrho shaf_shift kappa_vals a R0'
                                                        ' Te Ti ne Z_eff rtang beamWidth gaussR beamE beamA '
                                                        'beamP iol coCurr verbose timing iolFlag '
                                                        'beamResult pwrFracOverride')

    BeamSources = namedtuple('BeamSources', 'Snbi Qnbi Mnbi')

    def __init__(self, inp, core, iol=None, reRun=False, pwrFracOverride=None):
        super().__init__()
        self._core = core
        sys.dont_write_bytecode = True
        # If deposition profile is provided as an input, use that.
        # Note, this is designed to read in the dP/dr(r) quantity, not dP/dV(r).
        # When integrated over r, this should equal P_beam, or at least the amount
        # that is deposited in the plasma.

        self.beams_space = np.linspace(0., core.sep_val, 50)
        """The rho space used by the Beams module on [0., 1.] or [0., sep_val] to speed up computation"""

        self.set_plot_rho1d(self.beams_space)

        self.pwrFracOverride = pwrFracOverride
        """Allows us to override the beam power fraction for sensitivity studies"""

        if not reRun:
            # nbi_vals = calc_nbi_vals(inp, core)
            # Look for NBI output file
            try:
                self.load_beams(inp, core, iol)
            except BaseException as e:
                print("Failed to load beams data. Running beams module...")
                self.call_beams(inp, core, iol)
        else:
            self.call_beams(inp, core, iol)

        try:
            # Combine into one beam.
            total_S = self._to_main_grid(np.array([a.part_src.total for a in self.beam_result]).sum(axis=0))
            total_Q = self._to_main_grid(np.array([a.en_src.total for a in self.beam_result]).sum(axis=0))
            total_M = self._to_main_grid(np.array([a.mom_src.total for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_total = self.BeamSources(OneDProfile(core.psi, total_S, core.R, core.Z),
                                                            OneDProfile(core.psi, total_Q, core.R, core.Z),
                                                            OneDProfile(core.psi, total_M, core.R, core.Z))

            lost_S = self._to_main_grid(np.array([a.part_src.lost for a in self.beam_result]).sum(axis=0))
            lost_Q = self._to_main_grid(np.array([a.en_src.lost for a in self.beam_result]).sum(axis=0))
            lost_M = self._to_main_grid(np.array([a.mom_src.lost for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_lost = self.BeamSources(OneDProfile(core.psi, lost_S, core.R, core.Z),
                                                           OneDProfile(core.psi, lost_Q, core.R, core.Z),
                                                           OneDProfile(core.psi, lost_M, core.R, core.Z))

            kept_S = self._to_main_grid(np.array([a.part_src.lost for a in self.beam_result]).sum(axis=0))
            kept_Q = self._to_main_grid(np.array([a.en_src.lost for a in self.beam_result]).sum(axis=0))
            kept_M = self._to_main_grid(np.array([a.mom_src.lost for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_kept = self.BeamSources(OneDProfile(core.psi, kept_S, core.R, core.Z),
                                                           OneDProfile(core.psi, kept_Q, core.R, core.Z),
                                                           OneDProfile(core.psi, kept_M, core.R, core.Z))

            dens_total_S = self._to_main_grid(np.array([a.part_src_dens.total for a in self.beam_result]).sum(axis=0))
            dens_total_Q = self._to_main_grid(np.array([a.en_src_dens.total for a in self.beam_result]).sum(axis=0))
            dens_total_M = self._to_main_grid(np.array([a.mom_src_dens.total for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_dens_total = self.BeamSources(OneDProfile(core.psi, dens_total_S, core.R, core.Z),
                                                                OneDProfile(core.psi, dens_total_Q, core.R, core.Z),
                                                                OneDProfile(core.psi, dens_total_M, core.R, core.Z))

            dens_lost_S = self._to_main_grid(np.array([a.part_src_dens.lost for a in self.beam_result]).sum(axis=0))
            dens_lost_Q = self._to_main_grid(np.array([a.en_src_dens.lost for a in self.beam_result]).sum(axis=0))
            dens_lost_M = self._to_main_grid(np.array([a.mom_src_dens.lost for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_dens_lost = self.BeamSources(OneDProfile(core.psi, dens_lost_S, core.R, core.Z),
                                                                OneDProfile(core.psi, dens_lost_Q, core.R, core.Z),
                                                                OneDProfile(core.psi, dens_lost_M, core.R, core.Z))

            dens_kept_S = self._to_main_grid(np.array([a.part_src_dens.kept for a in self.beam_result]).sum(axis=0))
            dens_kept_Q = self._to_main_grid(np.array([a.en_src_dens.kept for a in self.beam_result]).sum(axis=0))
            dens_kept_M = self._to_main_grid(np.array([a.mom_src_dens.kept for a in self.beam_result]).sum(axis=0))

            self.combined_beam_src_dens_kept = self.BeamSources(OneDProfile(core.psi, dens_kept_S, core.R, core.Z),
                                                                OneDProfile(core.psi, dens_kept_Q, core.R, core.Z),
                                                                OneDProfile(core.psi, dens_kept_M, core.R, core.Z))
        except Exception as e:
            raise Exception("NBI module failed to create main beam profiles with error: %s" % str(e))

    def load_beams(self, inp, core: Core, iol):
        try:
            self.beam_result = []
            with open(os.path.join(os.getcwd(), inp.beams_out_json), "r") as f:
                l = yaml.safe_load(f)
                for a in list(l.keys()):
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
                                                                        UnivariateSpline(core.rho[:, 0], core.T.e.kev.to1D())(
                                                                            self.beams_space),
                                                                        UnivariateSpline(core.rho[:, 0],
                                                                                         core.T.i.kev.to1D())(
                                                                            self.beams_space),
                                                                        UnivariateSpline(core.rho[:, 0], core.n.e.to1D())(
                                                                            self.beams_space),
                                                                        UnivariateSpline(core.rho[:, 0], core.z_eff.to1D())(
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

    def call_beams(self, inp, core: Core, iol):
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
                for n in list(newdata.keys()):
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
                                                             core.R0_a,
                                                             UnivariateSpline(core.rho[:, 0], core.T.e.kev.to1D())(
                                                                 self.beams_space),
                                                             UnivariateSpline(core.rho[:, 0], core.T.i.kev.to1D())(
                                                                 self.beams_space),
                                                             UnivariateSpline(core.rho[:, 0], core.n.e.to1D())(
                                                                 self.beams_space),
                                                             UnivariateSpline(core.rho[:, 0], core.z_eff.to1D())(
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
                        print(("BeamConfiguration failed for beam %s with error: %s " % (n, str(e))))
                pool = Pool(cpu_count() - 1)

                try:
                    try:
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
            print("No Beams JSON file found. Using single-beam data")

    def _to_main_grid(self, val):
        return np.array(interp1d(self.beams_space, val, kind="linear")(self._core.rho[:, 0]))
