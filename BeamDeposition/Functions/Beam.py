#!/usr/bin/python2.7


from BeamDeposition.Functions.GetMFP import get_mfp, sigeff
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad, fixed_quad, cumtrapz, dblquad
from math import sqrt, exp
from time import time
import warnings
import exceptions
import logging
import numpy as np
from scipy.optimize import root_scalar
from math import pi, isinf
import matplotlib.pyplot as plt
from collections import namedtuple
from mpl_toolkits import mplot3d
import matplotlib.ticker as mtick




class Beam:
    def __init__(self, rho, r2vol, dVdr, rtang, shaf_shift, kappa_vals,
                 beamHeight, beamWidth, Te, TC, R0, ne, beamE, beamA, beamP, Zeff, a, gaussR, gaussZ, verbose=False, timing=False):
        """
        The Beams class for a single neutral beam injecter beam

        This class mimicks NBEAMS


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
        :param beamHeight: The beam height
        :type beamHeight:  float
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
        :param gaussZ: The gaussian height of the circular beam
        :type gaussZ:  float
        """

        start = time()
        self.verbose = verbose
        self.timing = timing
        self.rho = rho
        self.r2vol = r2vol
        self.dVdr = dVdr
        self.rtang = rtang
        self.shaf_shift_0 = shaf_shift
        self.kappa_vals = kappa_vals
        self.beamHeight = beamHeight
        self.beamWidth = beamWidth
        self.Te = Te
        self.Tc = TC
        self.R0 = R0
        self.ne = ne
        self.beamE = beamE
        self.beamA = beamA
        self.beamP = beamP
        self.Zeff = Zeff
        self.a = a
        self.gaussR = gaussR
        self.gaussZ = gaussZ
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

        self.prof_root['count'] = 0
        self.prof_root['time'] = []
        self.prof_mfp['count'] = 0
        self.prof_mfp['time'] = []
        self.prof_D['count'] = 0
        self.prof_D['time'] = []

        print "Calculating energy group 1"
        h_res_1, self.angle_1 = self.calcHofRho(self.beamE)

        print "Calculating energy group 2"
        h_res_2, self.angle_2 = self.calcHofRho(self.beamE / 2.)

        print "Calculating energy group 3"
        h_res_3, self.angle_3 = self.calcHofRho(self.beamE / 3.)

        DepositionProfiles = namedtuple('DepositionProfiles', ['D1', 'D2', 'D3'])

        self.Hofrho_unnorm = DepositionProfiles(h_res_1,
                                                h_res_2,
                                                h_res_3)

        self.shine = self.calc_shine()

        self.Hofrho = self.normalize_Hofrho(DepositionProfiles)
        self.Hofr = self.calc_Hofr(self.Hofrho, DepositionProfiles)

        self.pwrfrac = self.calc_power_frac(self.beamE)

        PowerProfiles = namedtuple('PowerProfiles', ['D1', 'D2', 'D3', 'total'])  # type: object

        self.power_profile = PowerProfiles(self.beamP * self.pwrfrac[0] * self.Hofrho.D1 / self.vol,
                                           self.beamP * self.pwrfrac[1] * self.Hofrho.D2 / self.vol,
                                           self.beamP * self.pwrfrac[2] * self.Hofrho.D3 / self.vol,
                                           self.beamP * self.pwrfrac[0] * self.Hofrho.D1 / self.vol +
                                           self.beamP * self.pwrfrac[1] * self.Hofrho.D2 / self.vol +
                                           self.beamP * self.pwrfrac[2] * self.Hofrho.D3 / self.vol)

        EnergySplit = namedtuple('EnergySplit', ['D1', 'D2', 'D3'])
        self.energies = EnergySplit(self.beamE, self.beamE / 2., self.beamE / 3.)

        end=time()

        if self.timing:
            print """
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
                                                           np.sum(self.prof_D['time']))


    def calc_power_frac(self, beamE):
        """
        Calculate the D3D power fraction

        TODO: Accept arbitrary functions into the Beams class to calculate this.

        :param beamE: The beam power
        :type beamE: float
        :return:
        """
        pwr_frac = np.zeros(3)
        pwr_frac[0] = (68 + 0.11 * beamE) / 100
        pwr_frac[1] = (-159 + 6.53 * beamE - 0.082 * (beamE ** 2) + 0.00034 * (beamE ** 3)) / 100
        pwr_frac[2] = (191 - 6.64 * beamE + 0.082 * (beamE ** 2) - 0.00034 * (beamE ** 3)) / 100
        return pwr_frac

    def calc_shine(self):
        shine_func_1 = lambda n: self.dVdr(self.rho[n] * self.a) * self.Hofrho_unnorm.D1
        shine_func_2 = lambda n: self.dVdr(self.rho[n] * self.a) * self.Hofrho_unnorm.D2
        shine_func_3 = lambda n: self.dVdr(self.rho[n] * self.a) * self.Hofrho_unnorm.D3

        self.shine = 1 - (cumtrapz(shine_func_1(range(0, len(self.rho))), self.rho)[-1] / self.vol +
                          cumtrapz(shine_func_2(range(0, len(self.rho))), self.rho)[-1] / self.vol +
                          cumtrapz(shine_func_3(range(0, len(self.rho))), self.rho)[-1] / self.vol)

        if self.shine <= 0. or self.shine > 1.:

             warnings.warn("Shinethrough is out of range at %s" % self.shine)

        return self.shine

    def normalize_Hofrho(self, DepClass):
        """
        Normalize H(rho) such that integral(H(rho)dV) = Vpol.
        Note that to do this integration from 0->1, you must multiply the integrand by dV/d(rho)

        :type DepClass:
        """

        norm_int_D1 = lambda n: self.dVdr(self.rho[n] * self.a) * self.Hofrho_unnorm.D1[n] * self.a
        norm_int_D2 = lambda n: self.dVdr(self.rho[n] * self.a) * self.Hofrho_unnorm.D2[n] * self.a
        norm_int_D3 = lambda n: self.dVdr(self.rho[n] * self.a) * self.Hofrho_unnorm.D3[n] * self.a

        normal_const_D1 = self.vol / cumtrapz([norm_int_D1(a) for a in range(0, len(self.rho))], self.rho)[-1]
        normal_const_D2 = self.vol / cumtrapz([norm_int_D2(a) for a in range(0, len(self.rho))], self.rho)[-1]
        normal_const_D3 = self.vol / cumtrapz([norm_int_D3(a) for a in range(0, len(self.rho))], self.rho)[-1]



        return DepClass(self.Hofrho_unnorm.D1 * normal_const_D1,
                        self.Hofrho_unnorm.D2 * normal_const_D2,
                        self.Hofrho_unnorm.D3 * normal_const_D3)

    def calc_Hofr(self, Hofrho, depClass):
        """

        :param Hofrho: The H(rho) function
        :type Hofrho: np.ndarray
        :return: The H(r) Univariate spline
        :rtype: UnivariateSpline
        """
        x = self.rho * self.a
        y1 = Hofrho.D1
        y2 = Hofrho.D2
        y3 = Hofrho.D3
        interp_1 = UnivariateSpline(x, y1, k=3, s=2)
        interp_2 = UnivariateSpline(x, y2, k=3, s=2)
        interp_3 = UnivariateSpline(x, y3, k=3, s=2)
        return depClass(interp_1(np.linspace(0., self.a, len(self.rho))),
                        interp_2(np.linspace(0., self.a, len(self.rho))),
                        interp_3(np.linspace(0., self.a, len(self.rho))))

    def calcHofRho(self, energy):
        """

        Calculate the beam deposition profile

        'quad' - scipy.integrate.quad - Uses Fortran QUADPACK technique
        'fixed_quad' - scipy.integrate.fixed_quad - Fixed-order gaussian quadrature (n=5)
        'trap' - scipy.integrate.cumtrapz - Composite trapazoidal rule integration (n=10)

        Fixed quadrature seems to get the same answer on a few tests in a fraction of the time

        TODO: Run more extensive tests to validate quadrature choice
        :param energy: The energy group (in keV)
        :return: The H(rho) array
        :rtype: np.ndarray

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
            if self.verbose: print val
            if float(val)==1.00:
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
            if self.verbose: print "Beginning integration of node n=%s" % n
            now = time()

            if self.verbose: print "Running positive"
            if smallFlag:
                result_pos[n] = pi * self.kappa(val) * self.beamdblquad(hofrpzJLambda, 0, Zulimit, RLowerLimit, RUpperLimit, val)
            else:
                result_pos[n] = 2.0 * self.beamdblquad(hofrpzJLambda, 0, Zulimit, RLowerLimit, RUpperLimit, val)
            posFlag = False

            if self.verbose: print "Running negative"
            if smallFlag:
                result_minus[n] = pi * self.kappa(val) * self.beamdblquad(hofrpzJLambda, 0, Zulimit, RLowerLimit, RUpperLimit,
                                                                        val)
            else:
                result_minus[n] = 2.0 * self.beamdblquad(hofrpzJLambda, 0.0, Zulimit, RLowerLimit, RUpperLimit, val)
            total = time() - now
            Hresult[n] = result_pos[n] + result_minus[n]

            posFlag = True
            angle_top_pos = self.beamdblquad(angleFunTopLambda, 0, Zulimit, RLowerLimit, RUpperLimit, val)
            angle_bottom_pos = result_pos[n] / 2.
            posFlag = False
            angle_top_neg = self.beamdblquad(angleFunTopLambda, 0, Zulimit, RLowerLimit, RUpperLimit, val)
            angle_bottom_neg = result_minus[n] / 2.
            angle_result[n] = (angle_top_pos + angle_top_neg)/(angle_bottom_neg + angle_bottom_pos)
            #angle_result[n] = (angle_top_pos / angle_bottom_pos) + (angle_top_neg / angle_bottom_neg)

            if self.verbose: print "Node n=%s completed in %s s" % (n, total)
            if self.verbose: print "Result: %s" % Hresult[n]
            if self.timing: print "Root scalar count: %s" % self.prof_root["count"]
            if self.timing: print "Root scalar total time: %s" % np.sum(self.prof_root['time'])
            if self.timing: print "MFP count: %s" % self.prof_mfp["count"]
            if self.timing: print "MFP total time: %s" % np.sum(self.prof_mfp['time'])
            if self.timing: print "D calc count: %s" % self.prof_D["count"]
            if self.timing: print "D calc total time: %s" % np.sum(self.prof_D['time'])


        return Hresult, angle_result

    def beamdblquad(self, f, yllimit, yulimit, xllimit, xulimit, *args):
        if yllimit == yulimit: return 0
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
                results[a] = fixed_quad(flambda, xllimit(mesh[a]), xulimit(mesh[a]), n=4)[0] * (yulimit - yllimit)/ yres
                if self.verbose: print "Integration result for a = %s: %s" % (a, results[a])
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
        if self.verbose: print "Coordinates: (Rb,Z) = (%s, %s)" % (Rb, Zb)
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
                return self.hofrpz(Rb, Zb, rho_val, mfp, p) * self.J(Rb, Zb)
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
        The h(rho, Rb, Zb) function

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
        Rin = self.R0 + sqrt(self.a**2 - Zb**2 / self.kappa(self.a)**2)

        lambda_mfp = mfp  # type: UnivariateSpline
        gamma = 1
        r = rho_val * self.a

        if p:
            if (r)**2 - Zb**2/(self.kappa(rho_val)**2) >= 0:
                Rplus = self.get_Rpm(Zb, rho_val, True)
                if Rplus:
                    #if Rplus > (self.rtang - sqrt((self.beamWidth / 2.) ** 2 - Zb**2)):
                    if Rplus > (Rb - sqrt((self.beamWidth / 2.) ** 2 - Zb ** 2)):
                        D0plus = self.D(Rplus, Rin, Zb, mfp)
                        if Rb >= Rin:
                            D1plus = self.D(Rb, Rplus, Zb, mfp)
                        else:
                            D1plus = 0.
                            gamma = 0.
                        # Updated to use r2vol instead of full plasma volume since that seems to make more sense.
                        part1plus = ((2 * r * self.r2vol(self.a)) / self.dVdr(r)) * (1 / lambda_mfp(rho_val)) * (Rplus / (sqrt(Rplus ** 2 - self.rtang ** 2)))
                        if smallFlag:
                            part2plus = 1.
                        else:
                            part2plus = np.nan_to_num((1 + (self.kappa_prime(rho_val) / (r)) * ((Zb ** 2) / self.kappa(rho_val) ** 3)) / (
                                sqrt(r ** 2 - (Zb ** 2) / (self.kappa(rho_val) ** 2))) + (self.shift_prime(rho_val) / (r)))
                        part3plus = exp(-1. * D0plus) + gamma * exp(-1. * (D0plus + 2.0 * D1plus))

                        hplus = part1plus * part2plus * part3plus
                        if self.verbose: print "Posiitve: Part 1: %s" % part1plus
                        if self.verbose: print "Posiitve: Part 2: %s" % part2plus
                        if self.verbose: print "Posiitve: Part 3: %s" % part3plus
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
                        part1minus = ((2 * (r) * self.r2vol(self.a)) / self.dVdr(r)) * (1 / lambda_mfp(rho_val)) * (Rminus / (sqrt(Rminus**2 - self.rtang**2)))
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
                        if self.verbose: print "Negative: Part 1: %s" % part1minus
                        if self.verbose: print "Negative: Part 2: %s" % part2minus
                        if self.verbose: print "Negative: Part 3: %s" % part3minus
                        if self.verbose: print "Negative result: %s" % hminus
                        return hminus
                    else:
                        return 0
                else:
                    return 0
            else:
                return 0


    def get_Rpm(self, Zb, rho_val, p):
        """
        Get R_+- for a given (rho, Zb)

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

        # type: (float, float, float, UnivariateSpline) -> float

        DfunLambda = lambda x, mfp=mfp: self.Dfun(x, Zb, mfp)

        if self.quad == 'quad':
            result = quad(DfunLambda, R1, R2, epsrel=.05)
            return result[0]

        elif self.quad == 'fixed_quad':
            self.prof_D['count'] += 1
            now = time()
            result = fixed_quad(DfunLambda, R1, R2, n=2)[0]
            if result < 0.:
                raise ValueError("Attenuation D is negative")
            self.prof_D['time'].append(time() - now)
            return result

        elif self.quad == 'trap':
            return cumtrapz([DfunLambda(a) for a in np.linspace(R1, R2, 10)], np.linspace(R1, R2, 10))[-1]

    def Dfun(self, x, Zb, mfp):

        if type(x) == np.ndarray:
            return[(a / (sqrt(a**2 - self.rtang**2))) * (1 / self.lam_special(a, Zb, mfp)) for a in x]
        else:
            return (x / (sqrt(x**2 - self.rtang**2))) * (1 / self.lam_special(x, Zb, mfp))

    def lam_special(self, R, Zb, mfp):
        """
        Attempts to find the special lambda for this R, Zb coordinate for this beam energy

        :param mfp:
        :type mfp: UnivariateSpline
        :param R: The major radius coordinate of the beam
        :type R:  float
        :param Zb: The Z coordinate of the beam
        :type Zb:  float
        :return:
        """
        now = time()
        fluxrho_l = lambda x: self.fluxrho(x, R, Zb)
        result = root_scalar(fluxrho_l, bracket=[0., 1.]).root

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

    def plot_HofRho(self):

        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$H(\rho)$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho, self.Hofrho.D1, color='red')
        fig1.scatter(self.rho, self.Hofrho.D2, color='blue')
        fig1.scatter(self.rho, self.Hofrho.D3, color='green')
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    def plot_power(self):

        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$P_{nb}(\rho)$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho, self.power_profile.D1, color='red')
        fig1.scatter(self.rho, self.power_profile.D2, color='blue')
        fig1.scatter(self.rho, self.power_profile.D3, color='green')
        fig1.scatter(self.rho, self.power_profile.total, color='black')
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))


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

    def plot_hofrho_unnorm(self):

        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #fig1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig1.set_xlabel(r'$/rho$')
        fig1.set_ylabel(r'$P_{nb}(\rho)$')
        #fig1.set_ylim(np.min(y, 0), np.max(y, 0))
        fig1.scatter(self.rho, self.Hofrho_unnorm.D1, color='red')
        fig1.scatter(self.rho, self.Hofrho_unnorm.D2, color='blue')
        fig1.scatter(self.rho, self.Hofrho_unnorm.D3, color='green')
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



