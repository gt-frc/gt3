#!/usr/bin/python2.7
from numpy.core.multiarray import ndarray

from GT3.BeamDeposition.Functions.GetMFP import get_mfp, sigeff
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad, fixed_quad, cumtrapz
from math import sqrt, exp
from time import time
import warnings
import numpy as np
from scipy.optimize import root_scalar
from math import pi, isinf
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy import constants
from GT3.Core.Functions.eVConvert import eVConvert

m_d = constants.physical_constants['deuteron mass'][0]

class Beam:

    beamE = None  # type: float
    beamP = None # type: float
    beamA = None # type: float
    beamWidth = None  # type: float
    shaf_shift_0 = None  # type: float
    dVdr = None  # type: UnivariateSpline
    dVdrho = None  # type: UnivariateSpline
    rtang = None  # type: float
    kappa_vals = None  # type: object
    r2vol = None  # type: UnivariateSpline
    timing = None  # type: bool
    rho = None  # type: np.ndarray
    verbose = None  # type: bool
    R0 = None # type: float
    Ti = None # type: np.ndarray
    Te = None # type: np.ndarray
    ne = None # type: np.ndarray
    iolFlag = True # type: bool
    coCurr = True # type: bool

    DepositionProfiles = namedtuple('DepositionProfiles', ['D1', 'D2', 'D3'])
    """A named tuple for deposition profiles for D1, D2 and D3 molecules"""
    EnergySplit = namedtuple('EnergySplit', ['D1', 'D2', 'D3'])
    PowerProfiles = namedtuple('PowerProfiles', ['D1', 'D2', 'D3', 'total'])
    MomentumProfiles = namedtuple('MomentumProfiles', ['D1', 'D2', 'D3', 'total', 'coCurr'])
    IOLSplit = namedtuple('IOLSplit', 'total lost kept')



    #def __init__(self, rho, r2vol, dVdr, rtang, shaf_shift, kappa_vals,
    #             beamWidth, Te, TC, R0, ne, beamE, beamA, beamP, Zeff, a, gaussR, verbose=False, timing=False):

    def __init__(self, config):
        """
        The Beams class for a single neutral beam injecter beam

        This class mimicks NBEAMS and provides the deposition profile H(rho,i) for the three molecular species.

        LIMITATIONS:

        This model assumes a circular beam with gaussian current density J. Beam shinethrough is not currently calculated
        correctly.

        Parameters
        ----------

        :param verbose: Enable verbose output
        :type verbose: bool
        :param timing: Enable timing data
        :type timing: bool
        :param rho: The rho vector
        :type rho:  np.ndarray
        :param r2vol: The r2vol interpolator
        :type r2vol:  UnivariateSpline
        :param dVdr: The dVdr interpolator
        :type dVdr:  UnivariateSpline
        :param rtang: The radius of tangency of this beam
        :type rtang:  float
        :param shaf_shift: The shavronof shift interpolator
        :type shaf_shift:  UnivariateSpline
        :param kappa_vals: The kappa values at the mag axes/sep
        :type kappa_vals:  np.ndarary(2)
        :param beamWidth: The beam width (2x the radius for circular)
        :type beamWidth:  float
        :param Te: The electron temperature np.ndarray
        :type Te:  np.ndarray
        :param TC:  The carbon temperature np.ndarray
        :type TC:  np.ndarray
        :param R0: The major radius at the magnetic axis
        :type R0:  float
        :param ne: The electron density np.ndarray
        :type ne:  np.ndarray
        :param beamE: The beam energy in keV
        :type beamE:  float
        :param beamA: The beam ion mass in u
        :type beamA:  float
        :param beamP: The beam power in MW (as scaled by the duty cycle)
        :type beamP: float
        :param Zeff:  The Z-effective np.array
        :type Zeff:  np.ndarray
        :param a: The plasma radius
        :type a:  float
        :param gaussR: The gaussian radius of hte circular beam
        :type gaussR:  float

        """

        self.config = config
        self.name = config.name
        start = time()
        self.verbose = config.verbose
        self.timing = config.timing
        self.coCurr = config.coCurr
        self.rho = config.rho
        self.iol = config.iol
        self.iolFlag = config.iolFlag

        self.r2vol = config.r2vol
        self.dVdr = config.dVdr
        self.dVdrho = config.dVdrho
        self.rtang = config.rtang
        self.shaf_shift_0 = config.shaf_shift
        self.kappa_vals = config.kappa_vals
        self.beamWidth = config.beamWidth
        self.Te = config.Te
        self.Tc = config.Ti
        self.R0 = config.R0
        self.ne = config.ne
        self.beamE = config.beamE * 1000
        self.beamA = config.beamA
        self.beamP = config.beamP
        self.Zeff = config.Z_eff
        self.a = config.a
        self.gaussR = config.gaussR
        self.quad = "fixed_quad"
        self.prof_root = {}
        self.prof_mfp = {}
        self.prof_D = {}
        self.kappa_vals = self.kappa_vals.axis + (self.kappa_vals.sep - self.kappa_vals.axis) * self.rho**2
        self.shift_vals = self.shaf_shift_0 * (1 - self.rho**2)
        self.shift = UnivariateSpline(self.rho, self.shift_vals)
        self.shift_prime = self.shift.derivative()
        self.kappa = UnivariateSpline(self.rho, self.kappa_vals)
        self.kappa_prime = self.kappa.derivative()
        self.vol = self.r2vol(self.a)
        self.pwrFracOverride = config.pwrFracOverride
        self.sep_val = self.rho[-1]

        self.prof_root['count'] = 0
        self.prof_root['time'] = []
        self.prof_mfp['count'] = 0
        self.prof_mfp['time'] = []
        self.prof_D['count'] = 0
        self.prof_D['time'] = []

        if not config.beamResult:
            self.is_new(config)
        else:
            old_result = config.beamResult
            try:

                self.angle_1 = UnivariateSpline(old_result[6], old_result[10], k=3, s=0)
                self.angle_2 = UnivariateSpline(old_result[6], old_result[11], k=3, s=0)
                self.angle_3 = UnivariateSpline(old_result[6], old_result[12], k=3, s=0)
                self.Hofrho = self.DepositionProfiles(UnivariateSpline(old_result[6], old_result[7], k=3, s=0),
                                                  UnivariateSpline(old_result[6], old_result[8], k=3, s=0),
                                                  UnivariateSpline(old_result[6], old_result[9], k=3, s=0))
                self.shine = old_result[13]
                self.Hofr = self.calc_Hofr(self.Hofrho, self.DepositionProfiles)
                """The normalized H(r) function on [0., a]"""
                if self.pwrFracOverride:
                    self.pwrfrac = self.pwrFracOverride
                    print(("Power fraction overwritten: " + str(self.pwrfrac)))
                else:
                    self.pwrfrac = self.calc_power_frac(self.beamE)
                self.dPdV = self.calc_dPdV(self.PowerProfiles)
                self.energies = self.EnergySplit(eVConvert(self.beamE), eVConvert(self.beamE / 2.),
                                            eVConvert(self.beamE / 3.))

                self.calc_iol(self.iol)
                self.calc_part_sources(self.IOLSplit)
                self.calc_heat_sources(self.IOLSplit)
                self.calc_mom_sources(self.IOLSplit)

            except Exception as e:
                print(("Failed to rebuild beam %s with error: %s" % (self.name, str(e))))


    def is_new(self, config):
        start = time()
        h_res_1 = None
        """The unnormalized H(rho) interpolator for Energy Group 1 on [0., 1.]"""
        h_res_2 = None
        """The unnormalized H(rho) interpolator for Energy Group 2 on [0., 1.]"""
        h_res_3 = None
        """The unnormalized H(rho) interpolator for Energy Group 3 on [0., 1.]"""

        self.angle_1 = None
        """The interpolator for Energy Group 1 for the launch angle on [0., 1.]"""
        self.angle_2 = None
        """The interpolator for Energy Group 2 for the launch angle on [0., 1.]"""
        self.angle_3 = None
        """The interpolator for Energy Group 3 for the launch angle on [0., 1.]"""

        print(("Running %s Beam" % config.name))
        if self.verbose: print("Calculating energy group 1")
        h_res_1, self.angle_1 = self.calcHofRho(self.beamE)

        if self.verbose: print("Calculating energy group 2")
        h_res_2, self.angle_2 = self.calcHofRho(self.beamE / 2.)
        """The H(rho) and launch angle(rho) for Energy Group 2 on [0., 1.]"""

        if self.verbose: print("Calculating energy group 3")
        h_res_3, self.angle_3 = self.calcHofRho(self.beamE / 3.)  # type: (UnivariateSpline, UnivariateSpline)
        """The H(rho) and launch angle(rho) for Energy Group 3 on [0., 1.]"""


        self.Hofrho_unnorm = self.DepositionProfiles(h_res_1,
                                                h_res_2,
                                                h_res_3)

        self.shine = self.calc_shine()
        self.Hofrho = self.normalize_Hofrho(self.DepositionProfiles)
        self.Hofr = self.calc_Hofr(self.Hofrho, self.DepositionProfiles)
        """The normalized H(r) function on [0., a]"""
        if self.pwrFracOverride:
            self.pwrfrac = self.pwrFracOverride
            print(("Power fraction overwritten: " + str(self.pwrfrac)))
        else:
            self.pwrfrac = self.calc_power_frac(self.beamE)
        self.dPdV = self.calc_dPdV(self.PowerProfiles)
        self.energies = self.EnergySplit(eVConvert(self.beamE), eVConvert(self.beamE / 2.), eVConvert(self.beamE / 3.))
        self.calc_iol(self.iol)
        self.calc_part_sources(self.IOLSplit)
        self.calc_heat_sources(self.IOLSplit)
        self.calc_mom_sources(self.IOLSplit)


        end=time()

        if self.timing:
            print(("""
            Timing summary:
            
            Total runtime: %s seconds
            
            Number of rootfinding operations: %s
            Total time for rootfinding: %s seconds
            Number of MFP operations: %s
            Total time of MFP operations: %s seconds
            Number of calculations of D: %s
            Total time of D calculations: %s seconds""" % (end - start,
                                                           self.prof_root['count'],
                                                           np.sum(self.prof_root['time']),
                                                           self.prof_mfp['count'],
                                                           np.sum(self.prof_mfp['time']),
                                                           self.prof_D['count'],
                                                           np.sum(self.prof_D['time']))))


    def calc_iol(self, iol):

        """
        Calculate the Ion Orbit Loss coefficients from the given iol module results. Returns nothing.
        The momentum loss fraction is retuned as the absolute value since there shouldn't realistically be any direction
        change.

        :param iol:
        """
        # TODO: Introduce accurate launch angles
        if iol:
            D1_iol = iol.calc_iol_beams(1, m_d, iol.iol_p, 33,  sqrt(2 * self.beamE * 1.6021E-19 / m_d), -.96, iol.coslist)
            D2_iol = iol.calc_iol_beams(1, m_d, iol.iol_p, 33, sqrt(2 * self.beamE * 1.6021E-19 / (m_d * 2.)), -.96,
                                        iol.coslist)
            D3_iol = iol.calc_iol_beams(1, m_d, iol.iol_p, 33, sqrt(2 * self.beamE * 1.6021E-19 / (3. * m_d)), -.96,
                                        iol.coslist)
            # Since we don't know if the beam meshing will be the same as the iol rho meshing, we map by index

            self.F_orb_D1 = np.zeros(len(self.rho))
            self.M_orb_D1 = np.zeros(len(self.rho))
            self.E_orb_D1 = np.zeros(len(self.rho))
            self.F_orb_D1[np.where(self.rho > iol.rho[:, 0][np.where(D1_iol[0][:, 0] > .1)[0].min()])[0].min():] = 0.5
            self.M_orb_D1[np.where(self.rho > iol.rho[:, 0][np.where(np.abs(D1_iol[1][:, 0]) > .1)[0].min()])[0].min():] = 0.5
            self.E_orb_D1[np.where(self.rho > iol.rho[:, 0][np.where(D1_iol[2][:, 0] > .1)[0].min()])[0].min():] = 0.5

            self.F_orb_D2 = np.zeros(len(self.rho))
            self.M_orb_D2 = np.zeros(len(self.rho))
            self.E_orb_D2 = np.zeros(len(self.rho))
            self.F_orb_D2[np.where(self.rho > iol.rho[:, 0][np.where(D2_iol[0][:, 0] > .1)[0].min()])[0].min():] = 0.5
            self.M_orb_D2[np.where(self.rho > iol.rho[:, 0][np.where(np.abs(D2_iol[1][:, 0]) > .1)[0].min()])[0].min():] = 0.5
            self.E_orb_D2[np.where(self.rho > iol.rho[:, 0][np.where(D2_iol[2][:, 0] > .1)[0].min()])[0].min():] = 0.5

            self.F_orb_D3 = np.zeros(len(self.rho))
            self.M_orb_D3 = np.zeros(len(self.rho))
            self.E_orb_D3 = np.zeros(len(self.rho))
            self.F_orb_D3[np.where(self.rho > iol.rho[:, 0][np.where(D3_iol[0][:,0]>.1)[0].min()])[0].min():] = 0.5
            self.M_orb_D3[np.where(self.rho > iol.rho[:, 0][np.where(np.abs(D3_iol[1][:, 0]) > .1)[0].min()])[0].min():] = 0.5
            self.E_orb_D3[np.where(self.rho > iol.rho[:, 0][np.where(D3_iol[2][:, 0] > .1)[0].min()])[0].min():] = 0.5

        else:
            self.F_orb_D1 = np.zeros(len(self.rho))
            self.F_orb_D2 = np.zeros(len(self.rho))
            self.F_orb_D3 = np.zeros(len(self.rho))

            self.E_orb_D1 = np.zeros(len(self.rho))
            self.E_orb_D2 = np.zeros(len(self.rho))
            self.E_orb_D3 = np.zeros(len(self.rho))

            self.M_orb_D1 = np.zeros(len(self.rho))
            self.M_orb_D2 = np.zeros(len(self.rho))
            self.M_orb_D3 = np.zeros(len(self.rho))


    def calc_part_sources(self, IOLSplit):
        """
        Calculate the particle sources in #/m^3*s as well as dS/d(rho)
        """
        self.part_src_dens_D1 = None  # type: ndarray
        """The particle source density as a function of r in #/(m^3 s) for the D1 beam"""
        self.part_src_dens_D2 = None  # type: ndarray
        """The particle source density as a function of r in #/(m^3 s) for the D2 beam"""
        self.part_src_dens_D3 = None  # type: ndarray
        """The particle source density as a function of r in #/(m^3 s) for the D3 beam"""
        self.part_src_D1 = None  #type: ndarray
        """The particle source as a function of r in #/(s) for the D1 beam"""
        self.part_src_D2 = None  #type: ndarray
        """The particle source as a function of r in #/(s) for the D2 beam"""
        self.part_src_D2 = None  #type: ndarray
        """The particle source as a function of r in #/(s) for the D3 beam"""

        if self.iolFlag:

            self.part_src_dens_D1 = IOLSplit(self.dPdV.D1(self.rho) / self.energies.D1.J,  # total
                                             self.F_orb_D1 * self.dPdV.D1(self.rho) / self.energies.D1.J,  # lost
                                             (1. - self.F_orb_D1) * self.dPdV.D1(self.rho) / self.energies.D1.J)  # kept

            self.part_src_dens_D2 = IOLSplit(self.dPdV.D2(self.rho) / self.energies.D2.J,  # total
                                             self.F_orb_D2 * self.dPdV.D2(self.rho) / self.energies.D2.J,  # lost
                                             (1. - self.F_orb_D2) * self.dPdV.D2(self.rho) / self.energies.D2.J)  # kept

            self.part_src_dens_D3 = IOLSplit(self.dPdV.D3(self.rho) / self.energies.D3.J,  # total
                                        self.F_orb_D3 * self.dPdV.D3(self.rho) / self.energies.D3.J,  # lost
                                        (1. - self.F_orb_D3) * self.dPdV.D3(self.rho) / self.energies.D3.J)  # kept

            self.part_src_D1 = IOLSplit(self.dPdV.D1(self.rho) * self.dVdrho(self.rho) / self.energies.D1.J,  # total
                                        self.F_orb_D1 * self.dPdV.D1(self.rho) * self.dVdrho(self.rho) / self.energies.D1.J,  # lost
                                        (1. - self.F_orb_D1) * self.dPdV.D1(self.rho) * self.dVdrho(self.rho) / self.energies.D1.J)  # kept

            self.part_src_D2= IOLSplit(self.dPdV.D2(self.rho) * self.dVdrho(self.rho) / self.energies.D2.J,  # total
                                        self.F_orb_D2 * self.dPdV.D2(self.rho) * self.dVdrho(self.rho) / self.energies.D2.J,  # lost
                                        (1. - self.F_orb_D2) * self.dPdV.D2(self.rho) * self.dVdrho(self.rho) / self.energies.D2.J)  # kept

            self.part_src_D3 = IOLSplit(self.dPdV.D3(self.rho) * self.dVdrho(self.rho) / self.energies.D3.J,  # total
                                        self.F_orb_D3 * self.dPdV.D3(self.rho) * self.dVdrho(self.rho) / self.energies.D3.J,  # lost
                                        (1. - self.F_orb_D3) * self.dPdV.D3(self.rho) * self.dVdrho(self.rho) / self.energies.D3.J)  # kept


        else:
            self.part_src_dens_D1 = IOLSplit(self.dPdV.D1(self.rho) / self.energies.D1.J,
                                        np.zeros(len(self.rho)),
                                        self.dPdV.D1(self.rho) / self.energies.D1.J)

            self.part_src_dens_D2 = IOLSplit(self.dPdV.D1(self.rho) / self.energies.D2.J,
                                        np.zeros(len(self.rho)),
                                        self.dPdV.D1(self.rho) / self.energies.D2.J)

            self.part_src_dens_D3 = IOLSplit(self.dPdV.D1(self.rho) / self.energies.D3.J,
                                        np.zeros(len(self.rho)),
                                        self.dPdV.D1(self.rho) / self.energies.D3.J)

            self.part_src_D1 = IOLSplit(self.dPdV.D1(self.rho) * self.dVdrho(self.rho) / self.energies.D1.J,  # total
                                        self.dPdV.D1(self.rho) * self.dVdrho(self.rho) / self.energies.D1.J,  # lost
                                        self.dPdV.D1(self.rho) * self.dVdrho(self.rho) / self.energies.D1.J)  # kept

            self.part_src_D2 = IOLSplit(self.dPdV.D2(self.rho) * self.dVdrho(self.rho) / self.energies.D2.J,  # total
                                        self.dPdV.D2(self.rho) * self.dVdrho(self.rho) / self.energies.D2.J,  # lost
                                        self.dPdV.D2(self.rho) * self.dVdrho(self.rho) / self.energies.D2.J)  # kept

            self.part_src_D3 = IOLSplit(self.dPdV.D3(self.rho) * self.dVdrho(self.rho) / self.energies.D3.J,  # total
                                        self.dPdV.D3(self.rho) * self.dVdrho(self.rho) / self.energies.D3.J,  # lost
                                        self.dPdV.D3(self.rho) * self.dVdrho(self.rho) / self.energies.D3.J)  # kept


        self.part_src_dens = IOLSplit(self.part_src_dens_D1.total + self.part_src_dens_D2.total + self.part_src_dens_D3.total,
                                      self.part_src_dens_D1.lost + self.part_src_dens_D2.lost + self.part_src_dens_D3.lost,
                                      self.part_src_dens_D1.kept + self.part_src_dens_D2.kept + self.part_src_dens_D3.kept)
        """The particle source density IOL split as a function of r in #/(m^3 s) """

        self.part_src = IOLSplit(self.part_src_D1.total + self.part_src_D2.total + self.part_src_D3.total,
                                 self.part_src_D1.lost + self.part_src_D2.lost + self.part_src_D3.lost,
                                 self.part_src_D1.kept + self.part_src_D2.kept + self.part_src_D3.kept)
        """The particle source IOL split as a function of r in #/s """

        # Remove any values that are stupidly small but not 0.
        #
        # self.part_src_tot_total[self.part_src_tot_total < 1E8] = 0.0
        # self.part_src_tot_kept[self.part_src_tot_kept < 1E8] = 0.0
        # self.part_src_tot_lost[self.part_src_tot_lost < 1E8] = 0.0
        # self.part_src_tot_of_rho[self.part_src_tot_of_rho < 1E8] = 0.0


    def calc_heat_sources(self, IOLSplit):

        self.en_src_dens_D1 = None  # type: ndarray
        """The energy source density as a function of r in #/(m^3 s) for the D1 beam"""
        self.en_src_dens_D2 = None  # type: ndarray
        """The energy source density as a function of r in #/(m^3 s) for the D2 beam"""
        self.en_src_dens_D3 = None  # type: ndarray
        """The energy source density as a function of r in #/(m^3 s) for the D3 beam"""
        self.en_src_D1 = None  #type: ndarray
        """The energy source as a function of r in #/(s) for the D1 beam"""
        self.en_src_D2 = None  #type: ndarray
        """The energy source as a function of r in #/(s) for the D2 beam"""
        self.en_src_D2 = None  #type: ndarray
        """The energy source as a function of r in #/(s) for the D3 beam"""

        if self.iolFlag:
            self.en_src_dens_D1 = IOLSplit(self.dPdV.D1(self.rho),  # total
                                      self.E_orb_D1 * self.dPdV.D1(self.rho),  # lost
                                      (1. - self.E_orb_D1) * self.dPdV.D1(self.rho))   # kept

            self.en_src_dens_D2 = IOLSplit(self.dPdV.D2(self.rho),  # total
                                      self.E_orb_D2 * self.dPdV.D2(self.rho),  # lost
                                      (1. - self.E_orb_D2) * self.dPdV.D2(self.rho))   # kept

            self.en_src_dens_D3 = IOLSplit(self.dPdV.D3(self.rho),  # total
                                      self.E_orb_D3 * self.dPdV.D3(self.rho),  # lost
                                      (1. - self.E_orb_D3) * self.dPdV.D3(self.rho))   # kept

            self.en_src_D1 = IOLSplit(self.dPdV.D1(self.rho) * self.dVdrho(self.rho),  # total
                                        self.E_orb_D1 * self.dPdV.D1(self.rho) * self.dVdrho(self.rho),  # lost
                                        (1. - self.E_orb_D1) * self.dPdV.D1(self.rho) * self.dVdrho(self.rho))  # kept

            self.en_src_D2 = IOLSplit(self.dPdV.D2(self.rho) * self.dVdrho(self.rho),  # total
                                        self.E_orb_D2 * self.dPdV.D2(self.rho) * self.dVdrho(self.rho),  # lost
                                        (1. - self.E_orb_D2) * self.dPdV.D2(self.rho) * self.dVdrho(self.rho))  # kept

            self.en_src_D3 = IOLSplit(self.dPdV.D3(self.rho) * self.dVdrho(self.rho),  # total
                                        self.E_orb_D3 * self.dPdV.D3(self.rho) * self.dVdrho(self.rho) ,  # lost
                                        (1. - self.E_orb_D3) * self.dPdV.D3(self.rho) * self.dVdrho(self.rho))  # kept

        else:
            self.en_src_dens_D1 = IOLSplit(self.dPdV.D1(self.rho),  # total
                                      np.zeros(len(self.rho)),  # lost
                                      self.dPdV.D1(self.rho))   # kept

            self.en_src_dens_D2 = IOLSplit(self.dPdV.D2(self.rho),  # total
                                      np.zeros(len(self.rho)),  # lost
                                      self.dPdV.D2(self.rho))   # kept

            self.en_src_dens_D3 = IOLSplit(self.dPdV.D3(self.rho),  # total
                                      np.zeros(len(self.rho)),  # lost
                                      self.dPdV.D3(self.rho))   # kept

            self.en_src_D1 = IOLSplit(self.dPdV.D1(self.rho) * self.dVdrho(self.rho),  # total
                                        self.dPdV.D1(self.rho) * self.dVdrho(self.rho),  # lost
                                        self.dPdV.D1(self.rho) * self.dVdrho(self.rho))  # kept

            self.en_src_D2 = IOLSplit(self.dPdV.D2(self.rho) * self.dVdrho(self.rho),  # total
                                        self.dPdV.D2(self.rho) * self.dVdrho(self.rho),  # lost
                                        self.dPdV.D2(self.rho) * self.dVdrho(self.rho))  # kept

            self.en_src_D3 = IOLSplit(self.dPdV.D3(self.rho) * self.dVdrho(self.rho),  # total
                                        self.dPdV.D3(self.rho) * self.dVdrho(self.rho),  # lost
                                        self.dPdV.D3(self.rho) * self.dVdrho(self.rho))  # kept


        self.en_src_dens = IOLSplit(self.en_src_dens_D1.total + self.en_src_dens_D2.total + self.en_src_dens_D3.total,
                                      self.en_src_dens_D1.lost + self.en_src_dens_D2.lost + self.en_src_dens_D3.lost,
                                      self.en_src_dens_D1.kept + self.en_src_dens_D2.kept + self.en_src_dens_D3.kept)
        """The particle source density IOL split as a function of r in J/(m^3 s) """

        self.en_src = IOLSplit(self.en_src_D1.total + self.en_src_D2.total + self.en_src_D3.total,
                                 self.en_src_D1.lost + self.en_src_D2.lost + self.en_src_D3.lost,
                                 self.en_src_D1.kept + self.en_src_D2.kept + self.en_src_D3.kept)
        """The particle source IOL split as a function of r in #/s """

        # Zero out very small values
        #
        # self.en_src_tot_total[self.en_src_tot_total < 1E2] = 0.0
        # self.en_src_tot_kept[self.en_src_tot_kept < 1E2] = 0.0
        # self.en_src_tot_lost[self.en_src_tot_lost < 1E2] = 0.0
        # self.en_src_tot_of_rho[self.en_src_tot_of_rho < 1E2] = 0.0


    def calc_mom_sources(self, IOLSplit):

        if self.iolFlag:
            self.mom_src_dens_D1 = IOLSplit(self.dPdV.D1(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1(self.rho),
                                       self.M_orb_D1 * self.dPdV.D1(self.rho)  * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1(self.rho),
                                       (1. - self.M_orb_D1) * self.dPdV.D1(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1(self.rho))

            self.mom_src_dens_D2 = IOLSplit(self.dPdV.D2(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2(self.rho),
                                       self.M_orb_D2 * self.dPdV.D2(self.rho)  * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2(self.rho),
                                       (1. - self.M_orb_D2) * self.dPdV.D2(self.rho)  * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2(self.rho))

            self.mom_src_dens_D3 = IOLSplit(self.dPdV.D3(self.rho) * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3(self.rho),
                                       self.M_orb_D3 * self.dPdV.D3(self.rho)  * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3(self.rho),
                                       (1. - self.M_orb_D3) * self.dPdV.D3(self.rho)  * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3(self.rho))

            self.mom_src_D1 = IOLSplit(self.dPdV.D1(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1(self.rho),
                                       self.M_orb_D1 * self.dPdV.D1(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1(self.rho),
                                       (1. - self.M_orb_D1) * self.dPdV.D1(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1(self.rho))

            self.mom_src_D2 = IOLSplit(self.dPdV.D2(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2(self.rho),
                                       self.M_orb_D2 * self.dPdV.D2(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2(self.rho),
                                       (1. - self.M_orb_D2) * self.dPdV.D2(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2(self.rho))

            self.mom_src_D3 = IOLSplit(self.dPdV.D3(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3(self.rho),
                                       self.M_orb_D3 * self.dPdV.D3(self.rho)  * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3(self.rho),
                                       (1. - self.M_orb_D3) * self.dPdV.D3(self.rho)  * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3(self.rho))


            if not self.coCurr:
                self.mom_src_D1 = -1. * self.mom_src_D1
                self.mom_src_D2 = -1. * self.mom_src_D2
                self.mom_src_D3 = -1. * self.mom_src_D3
                self.mom_src_D1_of_rho = -1 * self.mom_src_D1_of_rho
                self.mom_src_D2_of_rho = -1 * self.mom_src_D2_of_rho
                self.mom_src_D3_of_rho = -1 * self.mom_src_D3_of_rho

        else:
            self.mom_src_dens_D1 = IOLSplit(self.dPdV.D1(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1,
                                       np.zeros(len(self.rho)),
                                       self.dPdV.D1(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1)

            self.mom_src_dens_D2 = IOLSplit(self.dPdV.D2(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2,
                                       np.zeros(len(self.rho)),
                                       self.dPdV.D2(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2)

            self.mom_src_dens_D3 = IOLSplit(self.dPdV.D3(self.rho) * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3,
                                       np.zeros(len(self.rho)),
                                       self.dPdV.D3 * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3)

            self.mom_src_D1 = IOLSplit(self.dPdV.D1(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1,
                                       np.zeros(len(self.rho)),
                                       self.dPdV.D1(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D1.J) * self.rtang * self.angle_1)

            self.mom_src_D2 = IOLSplit(self.dPdV.D2(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2,
                                       np.zeros(len(self.rho)),
                                       self.dPdV.D2(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D2.J) * self.rtang * self.angle_2)

            self.mom_src_D3 = IOLSplit(self.dPdV.D3(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3,
                                       np.zeros(len(self.rho)),
                                       self.dPdV.D3(self.rho) * self.dVdrho(self.rho) * np.sqrt(0.5 * m_d / self.energies.D3.J) * self.rtang * self.angle_3)
            if not self.coCurr:
                self.mom_src_dens_D1 = -1. * self.mom_src_D1
                self.mom_src_dens_D2 = -1. * self.mom_src_D2
                self.mom_src_dens_D3 = -1. * self.mom_src_D3
                self.mom_src_D1 = -1 * self.mom_src_D1_of_rho
                self.mom_src_D2 = -1 * self.mom_src_D2_of_rho
                self.mom_src_D3 = -1 * self.mom_src_D3_of_rho

        self.mom_src_dens = IOLSplit(self.mom_src_dens_D1.total + self.mom_src_dens_D2.total + self.mom_src_dens_D3.total,
                                      self.mom_src_dens_D1.lost + self.mom_src_dens_D2.lost + self.mom_src_dens_D3.lost,
                                      self.mom_src_dens_D1.kept + self.mom_src_dens_D2.kept + self.mom_src_dens_D3.kept)
        """The momentum source density IOL split as a function of r in #/(m^3 s) """

        self.mom_src = IOLSplit(self.mom_src_D1.total + self.mom_src_D2.total + self.mom_src_D3.total,
                                 self.mom_src_D1.lost + self.mom_src_D2.lost + self.mom_src_D3.lost,
                                 self.mom_src_D1.kept + self.mom_src_D2.kept + self.mom_src_D3.kept)
        """The momentum source IOL split as a function of r in #/s """


    def calc_dPdV(self, Profile):
        """

        :param Profile:
        :return: Univariate splines of D1, D2, D3, and the total
        :rtype: Profile(UnivariateSpline)
        """
        y1 = lambda x: self.Hofrho.D1(x) * self.beamP *self.pwrfrac[0] * 1E6 / self.vol
        y2 = lambda x: self.Hofrho.D2(x) * self.beamP *self.pwrfrac[1] * 1E6 / self.vol
        y3 = lambda x: self.Hofrho.D3(x) * self.beamP *self.pwrfrac[2] * 1E6 / self.vol
        y4 = lambda x: y1(x) + y2(x) + y3(x)

        interp_1 = UnivariateSpline(self.rho, y1(self.rho), k=3, s=2)
        """The dP(rho)/dV"""
        interp_2 = UnivariateSpline(self.rho, y2(self.rho), k=3, s=2)
        interp_3 = UnivariateSpline(self.rho, y3(self.rho), k=3, s=2)
        interp_4 = UnivariateSpline(self.rho, y4(self.rho), k=3, s=2)


        return Profile(interp_1,
                       interp_2,
                       interp_3,
                       interp_4)

    def calc_power_frac(self, beamEeV):
        """
        Calculate the D3D power fraction

        TODO: Accept arbitrary functions into the Beams class to calculate this.

        :param beamE: The beam power
        :type beamE: float
        :return: A list of the D1/D2/D3 power fracs
        :rtype: list
        """
        beamE = beamEeV / 1E3
        pwr_frac = np.zeros(3)
        pwr_frac[0] = (68 + 0.11 * beamE) / 100
        pwr_frac[1] = (-159 + 6.53 * beamE - 0.082 * (beamE ** 2) + 0.00034 * (beamE ** 3)) / 100
        pwr_frac[2] = (191 - 6.64 * beamE + 0.082 * (beamE ** 2) - 0.00034 * (beamE ** 3)) / 100
        return pwr_frac

    def calc_shine(self):
        shine_func_1 = lambda x: self.dVdr(x * self.a) * self.Hofrho_unnorm.D1(x)
        shine_func_2 = lambda x: self.dVdr(x * self.a) * self.Hofrho_unnorm.D2(x)
        shine_func_3 = lambda x: self.dVdr(x * self.a) * self.Hofrho_unnorm.D3(x)

        shine = 1 - (cumtrapz(shine_func_1(self.rho), self.rho)[-1] / self.vol +
                          cumtrapz(shine_func_2(self.rho), self.rho)[-1] / self.vol +
                          cumtrapz(shine_func_3(self.rho), self.rho)[-1] / self.vol)

        if shine <= 0. or shine > 1.:

             warnings.warn("Shinethrough is out of range at %s" % shine)

        return shine

    def normalize_Hofrho(self, DepClass):
        """
        Normalize H(rho) such that integral(H(rho)dV) = Vplasma and provide an interpolator
        Note that to do this integration from 0->1, you must multiply the integrand by dV/d(rho)

        :type DepClass:
        """

        norm_int_D1 = lambda x: self.dVdr(x * self.a) * self.Hofrho_unnorm.D1(x) * self.a
        norm_int_D2 = lambda x: self.dVdr(x * self.a) * self.Hofrho_unnorm.D2(x) * self.a
        norm_int_D3 = lambda x: self.dVdr(x * self.a) * self.Hofrho_unnorm.D3(x) * self.a

        normal_const_D1 = self.vol / cumtrapz(norm_int_D1(self.rho), self.rho)[-1]
        normal_const_D2 = self.vol / cumtrapz(norm_int_D2(self.rho), self.rho)[-1]
        normal_const_D3 = self.vol / cumtrapz(norm_int_D3(self.rho), self.rho)[-1]

        interp_1 = UnivariateSpline(self.rho, self.Hofrho_unnorm.D1(self.rho) * normal_const_D1)
        interp_2 = UnivariateSpline(self.rho, self.Hofrho_unnorm.D2(self.rho) * normal_const_D2)
        interp_3 = UnivariateSpline(self.rho, self.Hofrho_unnorm.D3(self.rho) * normal_const_D3)

        return DepClass(interp_1,
                        interp_2,
                        interp_3)

    def calc_Hofr(self, Hofrho, depClass):
        """
        Calculate H(r) and return the depClass named tuple containing 3 interpolators as functions of r

        :param depClass: The deposition class object to be used
        :param Hofrho: The H(rho) function
        :type Hofrho: DepositionProfiles
        :return: The H(r) np.ndarray
        :rtype: depClass
        """
        y1 = lambda x: Hofrho.D1(x * self.a)
        y2 = lambda x: Hofrho.D2(x * self.a)
        y3 = lambda x: Hofrho.D3(x * self.a)
        interp_space = np.linspace(0., self.a, len(self.rho))

        return depClass(UnivariateSpline(interp_space, y1(self.rho), k=3, s=0),
                             UnivariateSpline(interp_space, y2(self.rho), k=3, s=0),
                             UnivariateSpline(interp_space, y3(self.rho), k=3, s=0))

    def calcHofRho(self, energy):
        """

        Calculate the beam deposition profile

        'quad' - scipy.integrate.quad - Uses Fortran QUADPACK technique
        'fixed_quad' - scipy.integrate.fixed_quad - Fixed-order gaussian quadrature (n=5)
        'trap' - scipy.integrate.cumtrapz - Composite trapazoidal rule integration (n=10)

        Fixed quadrature seems to get the same answer on a few tests in a fraction of the time

        TODO: Run more extensive tests to validate quadrature choice
        :param energy: The energy group (in keV)
        :return: UnivariateSplines of H(rho) and ksi(rho) interpolated on [0., 1.]
        :rtype: list(UnivariateSpline)

        """

        # Circular beam

        result_pos = np.zeros(len(self.rho))

        result_minus = np.zeros(len(self.rho))

        Hresult = np.zeros(len(self.rho))
        angle_result = np.zeros(len(self.rho))
        mfp = get_mfp(self.ne, energy, self.Te, self.beamA, self.Zeff, self.rho)
        smallFlag = False

        for n, val in enumerate(self.rho):
        #for n, val in enumerate([.5]):
            if self.verbose: print(val)
            if float(val)== self.sep_val:
                Hresult[n] = 0.
                angle_result[n] = 0.
                continue
            if val <=0.05:
                smallFlag = True
            else:
                smallFlag = False
            posFlag = True
            Zulimit= min(val * self.a * self.kappa(val), self.beamWidth/2.)
            RLowerLimit = lambda x: self.rtang - sqrt((self.beamWidth/2.)**2 - x**2)
            Rpm = lambda x, n=val: self.get_Rpm(x, n, posFlag)
            hofrpzJLambda = lambda y, x, val=val: self.hofrpzJ(x, y, val, mfp, posFlag, smallFlag)
            angleFunTopLambda = lambda y, x, val=val: self.hofrpzJ(x, y, val, mfp, posFlag, smallFlag) * y / Rpm(x)

            RUpperLimit = lambda x: min(Rpm(x), self.rtang + sqrt((self.beamWidth / 2.0)**2 - x**2))
            if self.verbose: print(("Beginning integration of node n=%s" % n))
            now = time()

            if self.verbose: print("Running positive")
            if smallFlag:
                result_pos[n] = pi * self.kappa(val) * self.beamdblquad(hofrpzJLambda, 0., Zulimit, RLowerLimit, RUpperLimit, val)
            else:
                result_pos[n] = 2.0 * self.beamdblquad(hofrpzJLambda, 0., Zulimit, RLowerLimit, RUpperLimit, val)
            posFlag = False

            if self.verbose: print("Running negative")
            if smallFlag:
                result_minus[n] = pi * self.kappa(val) * self.beamdblquad(hofrpzJLambda, 0., Zulimit, RLowerLimit, RUpperLimit,
                                                                        val)
            else:
                result_minus[n] = 2.0 * self.beamdblquad(hofrpzJLambda, 0.0, Zulimit, RLowerLimit, RUpperLimit, val)
            total = time() - now
            Hresult[n] = result_pos[n] + result_minus[n]

            posFlag = True
            angle_top_pos = self.beamdblquad(angleFunTopLambda, 0., Zulimit, RLowerLimit, RUpperLimit, val)
            angle_bottom_pos = result_pos[n] / 2.
            posFlag = False
            angle_top_neg = self.beamdblquad(angleFunTopLambda, 0., Zulimit, RLowerLimit, RUpperLimit, val)
            angle_bottom_neg = result_minus[n] / 2.

            # TODO: Still don't know which is correct

            angle_result[n] = (angle_top_pos + angle_top_neg)/(angle_bottom_neg + angle_bottom_pos)
            #angle_result[n] = (angle_top_pos / angle_bottom_pos) + (angle_top_neg / angle_bottom_neg)

            if self.verbose: print(("Node n=%s completed in %s s" % (n, total)))
            if self.verbose: print(("Result: %s" % Hresult[n]))
            if self.timing: print(("Root scalar count: %s" % self.prof_root["count"]))
            if self.timing: print(("Root scalar total time: %s" % np.sum(self.prof_root['time'])))
            if self.timing: print(("MFP count: %s" % self.prof_mfp["count"]))
            if self.timing: print(("MFP total time: %s" % np.sum(self.prof_mfp['time'])))
            if self.timing: print(("D calc count: %s" % self.prof_D["count"]))
            if self.timing: print(("D calc total time: %s" % np.sum(self.prof_D['time'])))

        Hofrho=UnivariateSpline(self.rho, Hresult, k=3, s=0)
        """The H(rho) interpolator on [0., 1.]"""
        angle = UnivariateSpline(self.rho, np.nan_to_num(angle_result), k=3, s=0)
        """The launch angle as a function of rho interpolator on [0., 1.]"""
        return Hofrho, angle

    def beamdblquad(self, f, yllimit, yulimit, xllimit, xulimit, *args):
        """
        Custom double integration function for f(Rb, Zb) integration.

        Integration limits are Z -> [yllimit, yulimit], R->[xllimit, xulimit]

        :rtype: float
        :return: Returns the numerical result of this integration
        :param f: The function to be ingrated in R,Z space
        :type f: function
        :param yllimit: The Zb lower limit
        :type yllimit: float
        :param yulimit: The Zb upper limit
        :type yulimit: float
        :param xllimit: The Rb lower limit
        :type xllimit: float
        :param xulimit: The Rb upper limit
        :type xulimit: float
        :param args:
        """
        if yllimit == yulimit:
            return 0
        rho_val = args[0]
        yres = 10
        mesh = np.linspace(yllimit, yulimit, yres)
        results = np.zeros(yres)

        for a in range(0, yres):
            flambda = lambda x, y=mesh[a], val=rho_val: f(x, y, val)
            if xulimit(mesh[a]) < xllimit(mesh[a]):
                results[a] = 0.
            elif (xulimit(mesh[a]) - xllimit(mesh[a])) / xllimit(mesh[a]) <= 10**(-5):
                results[a] = 0.
            elif a == yres-1:
                results[a] = results[a-1]
            else:
                try:
                    results[a] = fixed_quad(flambda, xllimit(mesh[a]), xulimit(mesh[a]), n=4)[0] * (yulimit - yllimit) / yres
                except:
                    print("Something went wrong")
                if self.verbose: print(("Integration result for a = %s: %s" % (a, results[a])))
        return np.sum(results)

    def hofrpzJ(self, Zb, Rb, rho_val, mfp, p, smallFlag):
        """

        :param smallFlag: If rho is small, the Z term will integrate to pi*k/2., which is handle in the main loop
        :param rho_val:  The rho value
        :type rho_val:  float
        :param p: Whether this is for positive (True) or negative (False)
        :type p:  bool
        :return:
        :param Zb: The beam-Z coordinate
        :type Zb: float
        :param Rb:
        :type Rb: THe beam-R coordinate
        :return:
        """
        if self.verbose: print(("Coordinates: (Rb,Z) = (%s, %s)" % (Rb, Zb)))
        Rpm_of_Z = lambda x, n=rho_val: self.get_Rpm(x, n, p)

        RLowerLimit = lambda x: self.rtang - sqrt((self.beamWidth / 2.) ** 2 - x ** 2)

        RUpperLimit = lambda x: min(Rpm_of_Z(x), self.rtang + sqrt((self.beamWidth / 2.0) ** 2 - x ** 2))

        if type(Rb) == np.ndarray:
            result = []
            for a in Rb:
                r = sqrt(Zb ** 2 + (self.rtang - a) ** 2)
                if (self.beamWidth / 2.0) ** 2 - Zb ** 2 <= 0:
                    result.append(0.)
                else:
                    if RUpperLimit(Zb):
                        if RLowerLimit(Zb) > RUpperLimit(Zb):
                            result.append(0.)
                        elif Rpm_of_Z(Zb) <= (self.rtang - sqrt((rho_val*self.a)**2 - ((Zb**2)/(self.kappa(rho_val)**2)))):
                            result.append(0.)
                        elif r > self.beamWidth / 2.:
                            result.append(0.)
                        else:
                        # print "Running for (Zb, Rb, x) = (%s, %s, %s)" % (Zb, Rb, x)
                            result.append(self.hofrpz(a, Zb, rho_val, mfp, p, smallFlag) * self.J(a, Zb))

                    else:
                        result.append(0.)
            return result
        else:
            r = sqrt(Zb ** 2 + (self.rtang - Rb) ** 2)
            if (self.beamWidth / 2.0) ** 2 - Zb ** 2 <= 0:
                return 0
            if RUpperLimit(Zb):
                if r > self.beamWidth / 2.:
                    return 0
                if RLowerLimit(Zb) > RUpperLimit(Zb):
                    return 0
                elif Rpm_of_Z(Zb) <= (self.rtang - sqrt((rho_val*self.a) ** 2 - ((Zb ** 2) / (self.kappa(rho_val) ** 2)))):
                    return 0
                return self.hofrpz(Rb, Zb, rho_val, mfp, p, smallFlag) * self.J(Rb, Zb)
            else:
                return 0.

    def J(self, Rb, Zb):
        """
        The beam current function at a given Rb,Zb coordinate

        :param Rb:
        :type Rb: float
        :param Zb:
        :type Zb: float
        :return: float
        """

        r = sqrt(Zb**2 + (self.rtang - Rb)**2)
        if r > self.beamWidth / 2.:
            return 0.
        else:
            expTermUp = -1. * (r**2 / self.gaussR**2)
            expTermDown = -1. * ((self.beamWidth / 2.)**2) / (self.gaussR**2)
            result = exp(expTermUp)/(pi * self.gaussR**2 * (1. - exp(expTermDown)))
            return result

    def hofrpz(self, Rb, Zb, rho_val, mfp, p, smallFlag):
        """
        The h(rho, Rb, Zb) function.

        :rtype: float
        :param p: Whether this is for the positive or negative half of the integration
        :type p: bool
        :param smallFlag: Flag to indicate that we are at small rho and need to address the singularity
        :type smallFlag: bool
        :return: The numeric value of this function
        :param mfp: The mean free path interpolator for this plasma energy group
        :type mfp: UnivariateSpline
        :param Rb: The beam-R coordinate
        :type float
        :param Zb: The beam-Z coordinate
        :type Zb: float
        :param rho_val: The rho value (0. < rho < 1.)
        :type float

        :return:
        """


        # TODO: This next lien is as in the manual, but NBEAMS was a subtraction
        Rin = self.R0 + self.shift(self.rho[-1]) + sqrt(self.a**2 - Zb**2 / self.kappa(self.a)**2)
        """ The R that the beam entered the plasma at (the plasma edge basically)"""

        lambda_mfp = mfp  # type: UnivariateSpline
        gamma = 1
        r = rho_val * self.a

        if p:
            # This probably should be updated for negative triangularity
            if (r)**2 - Zb**2/(self.kappa(rho_val)**2) >= 0:
                Rplus = self.get_Rpm(Zb, rho_val, True)
                """The largest R of this flux surface"""
                if Rplus:
                    # Is the beam R greater than the minimum R from the tangency radius on the inboard side
                    if Rplus > (self.rtang - sqrt((self.beamWidth / 2.) ** 2 - Zb**2)):
                        D0plus = self.D(Rplus, Rin, Zb, mfp)
                        if Rb >= Rin:
                            D1plus = self.D(Rb, Rplus, Zb, mfp)
                        else:
                            D1plus = 0.
                            gamma = 0.
                        # Updated to use r2vol instead of full plasma volume since that seems to make more sense.
                        # Doing self.rho[-2] beacuse when the separatrix is manually defined, the interpolator seems
                        # to mess up that very last value, so we push inward slightly.
                        #part1plus = ((2 * r * self.r2vol(self.a)) / self.dVdr(r)) * (1 / lambda_mfp(rho_val)) * (Rplus / (sqrt(Rplus ** 2 - self.rtang ** 2)))
                        part1plus = ((2 * r * self.r2vol(self.a * self.rho[-2])) / self.dVdr(r)) * (1 / lambda_mfp(rho_val)) * (Rplus / (sqrt(Rplus ** 2 - Rb ** 2)))
                        if smallFlag:
                            part2plus = 1.
                        else:
                            part2plus = np.nan_to_num((1 + (self.kappa_prime(rho_val) / (r)) * ((Zb ** 2) / self.kappa(rho_val) ** 3)) / (
                                sqrt(r ** 2 - (Zb ** 2) / (self.kappa(rho_val) ** 2))) + (self.shift_prime(rho_val) / (r)))
                        part3plus = exp(-1. * D0plus) + gamma * exp(-1. * (D0plus + 2.0 * D1plus))

                        hplus = part1plus * part2plus * part3plus
                        if self.verbose: print(("Posiitve: Part 1: %s" % part1plus))
                        if self.verbose: print(("Posiitve: Part 2: %s" % part2plus))
                        if self.verbose: print(("Posiitve: Part 3: %s" % part3plus))
                        if self.dVdr(r) <= 0.:
                            warnings.warn("Beams Module - Warning - dV/dr is negative at rho = %s" % rho_val)
                            hplus = 0.
                        if hplus < 0.:
                            raise ValueError("Beams Module - Error - hplus is negative at rho = %s" % rho_val)
                        if isinf(hplus):
                            raise ValueError("Beams Module - Erorr - hplus is infinite at rho = %s" % rho_val)
                        #print "Positive result: %s" % hplus

                        return hplus
                    else:
                        return 0
                else:
                    return 0
            else:
                return 0
        else:
            if (r)**2 - Zb**2/(self.kappa(rho_val)**2) >= 0:
                Rminus = self.get_Rpm(Zb, rho_val, False)
                if Rminus:
                    #if Rminus > (self.rtang - sqrt((self.beamWidth / 2.) ** 2 - Zb**2)):
                    if Rminus > (Rb - sqrt((self.beamWidth / 2.) ** 2 - Zb ** 2)):
                        D0minus = self.D(Rminus, Rin, Zb, mfp)
                        if Rb >= Rin:
                            D1minus = self.D(Rb, Rminus, Zb, mfp)
                        else:
                            D1minus = 0.
                            gamma = 0.
                        #part1minus = ((2 * (r) * self.r2vol(self.a)) / self.dVdr(r)) * (1 / lambda_mfp(rho_val)) * (Rminus / (sqrt(Rminus**2 - self.rtang**2)))
                        part1minus = ((2 * (r) * self.r2vol(self.a)) / self.dVdr(r)) * (1 / lambda_mfp(rho_val)) * (Rminus / (sqrt(Rminus ** 2 - Rb ** 2)))
                        if smallFlag:
                            part2minus = 1.
                        else:
                            part2minus = (1 + (self.kappa_prime(rho_val) / (r)) * ((Zb**2)/self.kappa(rho_val)**3))/(sqrt((r)**2 - (Zb**2)/(self.kappa(rho_val)**2))) + (self.shift_prime(rho_val)/(r))
                        part3minus = exp((-1. * D0minus)) + gamma * exp(-1. * (D0minus + 2.0 * D1minus))

                        hminus = part1minus + part2minus + part3minus
                        if self.dVdr(r) <= 0.:
                            warnings.warn("Beams Module - Warning - dV/dr is negative at rho = %s" % rho_val)
                            hminus = 0.
                        if hminus < 0.:
                            raise ValueError("Beams Module - Error - hminus is negative at rho = %s" % rho_val)
                        if isinf(hminus):
                            raise ValueError("Beams Module - Error - hminus is infinite at rho = %s" % rho_val)
                        if self.verbose: print(("Negative: Part 1: %s" % part1minus))
                        if self.verbose: print(("Negative: Part 2: %s" % part2minus))
                        if self.verbose: print(("Negative: Part 3: %s" % part3minus))
                        if self.verbose: print(("Negative result: %s" % hminus))
                        return hminus
                    else:
                        return 0
                else:
                    return 0
            else:
                return 0

    def get_Rpm(self, Zb, rho_val, p):
        """
        Get R_+- for a given (rho, Zb). This is the major radius of the outer (+) and inner (-) extents of the given
        flux surface a distance Zb above the midplane.

        :param Zb: The Zb coordinate
        :type Zb:  float
        :param rho_val: The rho value
        :type rho_val:  float
        :param p: Is this positive (True) or negative (False)
        :type p:  bool
        :return:
        """
        if p:
            if (rho_val*self.a)**2 - Zb**2/(self.kappa(rho_val)**2) <= 0:
                return False
            else:
                return self.R0 + self.shift(rho_val) + sqrt((rho_val*self.a)**2 - (Zb**2)/(self.kappa(rho_val)**2))
        else:
            if (rho_val*self.a)**2 - (Zb**2)/(self.kappa(rho_val)**2) <= 0:
                return False
            else:
                return self.R0 + self.shift(rho_val) - sqrt((rho_val*self.a)**2 - (Zb**2)/(self.kappa(rho_val)**2))

    def D(self, R1, R2, Zb, mfp):
        """
        This is the D function

        :param R1: The lower bound of the integration
        :type R1: float
        :param R2: The upper bound of the integration
        :type R2: float
        :param Zb: The Zb value that this integration is being performed at
        :type Zb: float
        :param mfp: The mean free path interpolator function passed to the Lambda function
        :type mfp: UnivariateSpline
        :return: The numeric value for this calculation
        :rtype: float
        """
        DfunLambda = lambda x, mfp=mfp: self.Dfun(x, Zb, mfp)

        if self.quad == 'quad':
            result = quad(DfunLambda, R1, R2, epsrel=.05)
            return result[0]

        elif self.quad == 'fixed_quad':
            self.prof_D['count'] += 1
            now = time()
            try:
                result = fixed_quad(DfunLambda, R1, R2, n=2)[0]
            except:
                print("Something went wrong")
            if result < 0.:
                raise ValueError("Attenuation D is negative")
            self.prof_D['time'].append(time() - now)
            return result

        elif self.quad == 'trap':
            try:
                return cumtrapz([DfunLambda(a) for a in np.linspace(R1, R2, 10)], np.linspace(R1, R2, 10))[-1]
            except:
                print("Something went wrong")

    def Dfun(self, x, Zb, mfp):
        """
        The function for use in the integration for calculating D. This takes Rb as a float or array depending on
        the integration scheme.

        :param x: The Rb value for the D integration at this point
        :type x: float
        :param Zb: The Zb value for the D integration at this point
        :type Zb: float
        :param mfp: The mean free path interpolator
        :type mfp: UnivariateSpline
        :return: The result of this function
        """
        if type(x) == np.ndarray:
            return[(a / (sqrt(a**2 - self.rtang**2))) * (1 / self.lam_special(a, Zb, mfp)) for a in x]
        else:
            return (x / (sqrt(x**2 - self.rtang**2))) * (1 / self.lam_special(x, Zb, mfp))

    def lam_special(self, R, Zb, mfp):
        """
        Attempts to find the special lambda for this R, Zb coordinate for this beam energy

        :return: Returns the mean free path in m
        :param mfp: The mean free path interpolator
        :type mfp: UnivariateSpline
        :param R: The major radius coordinate of the beam
        :type R:  float
        :param Zb: The Z coordinate of the beam
        :type Zb:  float
        :return:
        """
        now = time()
        fluxrho_l = lambda x: self.fluxrho(x, R, Zb)
        try:
            result = root_scalar(fluxrho_l, bracket=[0., self.sep_val]).root
        except Exception as e:
            print(("Problem with root calculation: %s" % str(e)))
            result = .9

        self.prof_root['count'] += 1
        self.prof_root['time'].append(time() - now)
        now = time()
        special_mfp = mfp(result)

        self.prof_mfp['count'] += 1
        self.prof_mfp['time'].append(time() - now)
        return special_mfp

    def fluxrho(self, x, R, Zb):
        """
        The F(x) = 0 function for calculating the flux surface that corresponds to the major radius in Eq (5)

        :param x: The rho (in [0., 1.]) value
        :type x:  float
        :param R: The R value from the integration
        :type R:  float
        :param Zb: The Zb value from the integration
        :type Zb:  float
        :return:  The rho value of the flux surface
        """
        return (self.a * x)**2 - (R - (self.R0 + self.shift(x)))**2 - Zb**2 / (self.kappa(x)**2)

    def plot_mfp(self):
        mfp_1 = get_mfp(self.ne, self.beamE, self.Te, self.beamA, self.Zeff, self.rho)
        mfp_2 = get_mfp(self.ne, self.beamE / 2., self.Te, self.beamA, self.Zeff, self.rho)
        mfp_3 = get_mfp(self.ne, self.beamE / 3., self.Te, self.beamA, self.Zeff, self.rho)

        mfp_l_1 = lambda x: mfp_1(x)
        mfp_l_2 = lambda x: mfp_2(x)
        mfp_l_3 = lambda x: mfp_3(x)
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$\lambda [m]$')
        fig1.scatter(self.rho, [mfp_l_1(x) for x in self.rho], color='red')
        fig1.scatter(self.rho, [mfp_l_2(x) for x in self.rho], color='blue')
        fig1.scatter(self.rho, [mfp_l_3(x) for x in self.rho], color='green')
        plt.show()
        return fig1

    def plot_xs(self):
        sig_int_1 = sigeff(self.ne, self.beamE, self.Te, self.beamA, self.Zeff, self.rho)
        sig_int_2 = sigeff(self.ne, self.beamE / 2., self.Te, self.beamA, self.Zeff, self.rho)
        sig_int_3 = sigeff(self.ne, self.beamE / 3., self.Te, self.beamA, self.Zeff, self.rho)
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$\sigma_{eff} [cm^2]$')
        fig1.set_ylim(1E-15, 1E-16)
        fig1.scatter(self.rho, sig_int_1(self.rho), color="red")
        fig1.scatter(self.rho, sig_int_2(self.rho), color="blue")
        fig1.scatter(self.rho, sig_int_3(self.rho), color="green")
        fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        plt.show()
        return fig1

    def plot_HofRho(self):

        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$H(\rho)$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho, self.Hofrho.D1(self.rho), color='red')
        fig1.scatter(self.rho, self.Hofrho.D2(self.rho), color='blue')
        fig1.scatter(self.rho, self.Hofrho.D3(self.rho), color='green')
        plt.show()
        return fig1
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    def plot_Hofr(self):
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        # fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$H(r)$')
        fig1.scatter(self.rho, self.Hofr.D1(self.rho / self.a), color='red')
        fig1.scatter(self.rho, self.Hofr.D2(self.rho / self.a), color='blue')
        fig1.scatter(self.rho, self.Hofr.D3(self.rho / self.a), color='green')
        plt.show()
        return fig1
        # fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    def plot_power(self):

        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$\frac{dP}{dV}_{nb}(\rho) [W/{m^3}]$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho, self.dPdV.D1(self.rho), color='red')
        fig1.scatter(self.rho, self.dPdV.D2(self.rho), color='blue')
        fig1.scatter(self.rho, self.dPdV.D3(self.rho), color='green')
        fig1.scatter(self.rho, self.dPdV.total(self.rho), color='black')
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        plt.show()
        return fig1

    def plot_hofRho(self, rho, energy):
        r = np.linspace(1, 1.3, 50)
        z = np.linspace(0., .1, 50)
        R, Z = np.meshgrid(r, z)
        f = lambda x, y: self.hofrpzJ(x, y, rho, energy,  True) + self.hofrpzJ(x, y, rho, energy, False)
        values = np.zeros((50, 50))
        for i in range(50):
            for j in range(50):
                values[i][j] = f(z[i], r[j])

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(R, Z, values, cmap='viridis', linewidth=0.5)
        plt.show()
        return fig, ax

    def plot_hofrho_unnorm(self):

        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$P_{nb}(\rho)$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho, self.Hofrho_unnorm.D1(self.rho), color='red')
        fig1.scatter(self.rho, self.Hofrho_unnorm.D2(self.rho), color='blue')
        fig1.scatter(self.rho, self.Hofrho_unnorm.D3(self.rho), color='green')
        plt.show()
        return fig1
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    def plot_J(self):
        r = np.linspace(1, 1.3, 50)
        z = np.linspace(0., .1, 50)
        R, Z = np.meshgrid(r, z)
        f = lambda x, y: self.J(x, y)
        values = np.zeros((50, 50))
        for i in range(50):
            for j in range(50):
                values[i][j] = f(r[i], z[j])

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(R, Z, values, cmap='viridis', linewidth=0.5)
        plt.show()
        return fig, ax

    def plot_part_src_kept(self):
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$r$')
        fig1.set_ylabel(r'$S_{nb}(r) [#/m^3 s]$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho * self.a, self.part_src_D1.kept, color='red')
        fig1.scatter(self.rho * self.a, self.part_src_D2.kept, color='blue')
        fig1.scatter(self.rho * self.a, self.part_src_D3.kept, color='green')
        fig1.scatter(self.rho * self.a, self.part_src_total.kept, color='green')
        plt.show()
        return fig1
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    def plot_energy_src_kep(self):
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        # fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$Q_{nb}(r) [W/m^3]$')
        # fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho * self.a, self.en_src_D1.kept, color='red')
        fig1.scatter(self.rho * self.a, self.en_src_D2.kept, color='blue')
        fig1.scatter(self.rho * self.a, self.en_src_D3.kept, color='green')
        fig1.scatter(self.rho * self.a, self.en_src_dens, color='black')
        plt.show()
        return fig1
        # fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    def plot_mom_src_kept(self):
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$r')
        fig1.set_ylabel(r'$M_{nb}(r)$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho * self.a, self.mom_src_D1.kept, color='red')
        fig1.scatter(self.rho * self.a, self.mom_src_D2.kept, color='blue')
        fig1.scatter(self.rho * self.a, self.mom_src_D3.kept, color='green')
        fig1.scatter(self.rho * self.a, self.mom_src_tot_kept, color='black')
        plt.show()
        return fig1
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    def plot_iol(self):
        fig = plt.figure()
        fig1 = fig.add_subplot(311)
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$F_{iol}(\rho)$')
        fig1.scatter(self.rho, self.F_orb_D1, color='red')
        fig1.scatter(self.rho, self.F_orb_D2, color='blue')
        fig1.scatter(self.rho, self.F_orb_D3, color='green')

        fig2 = fig.add_subplot(312)
        fig2.set_xlabel(r'$/rho$')
        fig2.set_ylabel(r'$E_{iol}(\rho)$')
        fig2.scatter(self.rho, self.E_orb_D1, color='red')
        fig2.scatter(self.rho, self.E_orb_D2, color='blue')
        fig2.scatter(self.rho, self.E_orb_D3, color='green')

        fig3 = fig.add_subplot(313)
        fig3.set_xlabel(r'$/rho$')
        fig3.set_ylabel(r'$M_{iol}(\rho)$')
        fig3.scatter(self.rho, self.M_orb_D1, color='red')
        fig3.scatter(self.rho, self.M_orb_D2, color='blue')
        fig3.scatter(self.rho, self.M_orb_D3, color='green')

        return fig

