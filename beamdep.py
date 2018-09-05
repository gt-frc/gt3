#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 21:51:19 2018

@author: max
"""
from __future__ import division
import numpy as np
from scipy.interpolate import interp1d
from subprocess import call, Popen, PIPE
import os
import re
import sys
from scipy.interpolate import interp1d, UnivariateSpline
from collections import namedtuple
from scipy.constants import physical_constants
from shapely.geometry import LineString, Point
from math import sqrt
import matplotlib.pyplot as plt

m_d = physical_constants['deuteron mass'][0]
z_d = 1

def calc_pwr_frac(ebeam):
    # calculate pwr fracs. These fits are only useful for DIII-D

    pwr_frac = np.zeros(3)
    pwr_frac[0] = (68 + 0.11 * ebeam) / 100
    pwr_frac[1] = (-159 + 6.53 * ebeam - 0.082 * (ebeam ** 2) + 0.00034 * (ebeam ** 3)) / 100
    pwr_frac[2] = (191 - 6.64 * ebeam + 0.082 * (ebeam ** 2) - 0.00034 * (ebeam ** 3)) / 100

    return pwr_frac


def prep_nbi_infile(inp, core):
    pwr_frac = calc_pwr_frac(inp.ebeam)
    # pwr_frac = [0.7, 0.2, 0.1]

    # f1=open(inp.nbeams_loc+"inbeams.dat", "w")
    if __name__=="__main__":
        f = open("inbeams_test.dat", "w")
    else:
        f = open("inbeams.dat", "w")
    f.write("$nbin\n")
    f.write("nbeams = 1\n")
    f.write("inbfus = 1\n")
    f.write("amb = 2.0\n")
    f.write("zbeam = 1.0\n")
    f.write("ebeam = " + str(inp.ebeam) + "\n")
    f.write("pbeam = " + str(inp.pbeam) + "\n")
    f.write("rtang = " + str(inp.rtang) + "\n")
    f.write("nbshape = 1\n")
    f.write("bwidth = 0.12\n")  # Is this default?
    f.write("bheigh = 0.48\n")
    f.write("bgaussR = 0.066\n")
    f.write("bgaussZ = 0.18\n")
    f.write("bzpos = 0.0\n")
    f.write("nbptype = 1\n")
    f.write("maxiter = 2\n")
    f.write("pwrfrac(1, 1) = " + str(pwr_frac[0]) + "   " + str(pwr_frac[1]) + "   " + str(pwr_frac[2]) + "\n")
    f.write("a = " + str(core.a) + "\n")
    f.write("r0 = " + str(core.R0_a) + "\n")
    f.write("b0 = " + str(inp.BT0) + "\n")
    f.write("n = 51\n")
    f.write("e0 = " + str(core.kappa_vals.axis) + "\n")
    f.write("ea = " + str(core.kappa_vals.sep) + "\n")
    f.write("shft0 = " + str(core.shaf_shift) + "\n")
    f.write("nion = 2\n")
    f.write("aion = 2.0 12.0\n")
    f.write("zion = 1.0 6.0\n")

    rho_nbi = np.linspace(0, 1, 51)

    ni_nbi = interp1d(core.rho[:, 0], core.n.i[:, 0])(rho_nbi)
    ne_nbi = interp1d(core.rho[:, 0], core.n.e[:, 0])(rho_nbi)
    Ti_nbi = interp1d(core.rho[:, 0], core.T.i.kev[:, 0])(rho_nbi)
    Te_nbi = interp1d(core.rho[:, 0], core.T.e.kev[:, 0])(rho_nbi)

    for i, v in enumerate(rho_nbi):
        f.write('ni20(' + str(i + 1) + ', 1) = ' + str(ni_nbi[i] * 1E-20) + '\n')
    for i, v in enumerate(rho_nbi):
        f.write('ne20(' + str(i + 1) + ') = ' + str(ne_nbi[i] * 1E-20) + '\n')
    for i, v in enumerate(rho_nbi):
        f.write('tikev(' + str(i + 1) + ') = ' + str(Ti_nbi[i]) + '\n')
    for i, v in enumerate(rho_nbi):
        f.write('tekev(' + str(i + 1) + ') = ' + str(Te_nbi[i]) + '\n')
    f.write("$end\n")
    f.close()


class read_nbi_outfile:

    def __init__(self, inp, core):
        with open(os.getcwd() + '/outbeams.dat', 'r') as f:
            for count, line in enumerate(f):
                if line.startswith(" Total Absorbed Power"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        P_abs_tot = float(result)
                    except:
                        P_abs_tot = np.NaN

                if line.startswith(" Total Lost Power"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        P_lst_tot = float(result)
                    except:
                        P_lst_tot = np.NaN

                if line.startswith(" Total NB Driven Current"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        I_nbi_tot = float(result)
                    except:
                        I_nbi_tot = np.NaN

                if line.startswith(" Total NBCD Efficiency"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        I_nbi_eff = float(result)
                    except:
                        I_nbi_eff = np.NaN

                if line.startswith(" Total Beam Beta"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        Beta_nbi_tot = float(result)
                    except:
                        Beta_nbi_tot = np.NaN

                if line.startswith(" Taus"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        Taus_nbi = float(result)
                    except:
                        Taus_nbi = np.NaN

                if line.startswith(" Volp"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        Volp_nbi = float(result)
                    except:
                        Volp_nbi = np.NaN

                if line.startswith(" Absorbed Power"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        P_abs_1 = float(result)
                    except:
                        P_abs_1 = np.NaN

                if line.startswith(" Lost Power"):
                    # this will need to be modified for multiple beams

                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        P_lst_1 = float(result)
                    except:
                        P_lst_1 = np.NaN

                if line.startswith(" NB driven current"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        I_nbi_1 = float(result)
                    except:
                        I_nbi_1 = np.NaN

                if line.startswith(" NBCD efficiency"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        I_nbi_eff_1 = float(result)
                    except:
                        I_nbi_eff_1 = np.NaN

                if line.startswith(" NBCD gamma"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        I_nbi_gam_1 = float(result)
                    except:
                        I_nbi_gam_1 = np.NaN

                if line.startswith("    energy group 1"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        st_en1_1 = float(result)
                    except:
                        st_en1_1 = np.NaN

                if line.startswith("    energy group 2"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        st_en2_1 = float(result)
                    except:
                        st_en2_1 = np.NaN

                if line.startswith("    energy group 3"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        st_en3_1 = float(result)
                    except:
                        st_en3_1 = np.NaN

                if line.startswith(" Total Beam-Target Fusion Power"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        fus_pwr_bt = float(result)
                    except:
                        fus_pwr_bt = np.NaN

                if line.startswith(" Total Power to Charged Particles"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        cp_pwr_tot = float(result)
                    except:
                        cp_pwr_tot = np.NaN

                if line.startswith(" Total DT Neutron Rate"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        rate_dt_n = float(result)
                    except:
                        rate_dt_n = np.NaN

                if line.startswith(" Total DD Neutron Rate"):
                    # this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        rate_dd_n = float(result)
                    except:
                        rate_dd_n = np.NaN

        with open(os.getcwd() + '/outbeams.dat', 'r') as f:
            data = f.read().replace('\n', ' ')

        result = re.match(r'.*hofr_3 *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) +Pitch.*', data).group(1)
        array = np.reshape(np.asarray(result.split(), dtype=float), (-1, 4))
        rho_nbeams = array[:, 0]
        hofr_1 = array[:, 1]
        hofr_2 = array[:, 2]
        hofr_3 = array[:, 3]

        result = re.match(r'.*zeta_3 *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) +rho.*', data).group(1)
        array = np.reshape(np.asarray(result.split(), dtype=float), (-1, 4))
        ptch_angl1 = array[:, 1]
        ptch_angl2 = array[:, 2]
        ptch_angl3 = array[:, 3]

        result = re.match(r'.*dA *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) *', data).group(1)
        array = np.reshape(np.asarray(result.split(), dtype=float), (-1, 9))
        jnbtot = array[:, 1]
        pNBe   = array[:, 2]  # in MW/m^3
        pNBi   = array[:, 3]  # in MW/m^3
        nbfast = array[:, 4]
        pressb = array[:, 5]
        pfusb  = array[:, 6]
        dvol   = array[:, 7]
        dA     = array[:, 8]

        #
        # H(r) is a strangely normalized power DENSITY deposition profile
        #
        #   H(r) = (V_tot / P_beam) * dP/dV(r)
        #
        # The two most useful quantities for GT3 are:
        #
        #   dP/dV(r) = H(r) * (P_beam / V_tot)
        #
        # and
        #
        #   dP/dr(r) = dP/dV(r) * dV/dr(r) = H(r) * (P_beam / V_tot) * dV/dr(r)
        #



        pwr_frac = calc_pwr_frac(inp.ebeam)
        self.beam_pwr_1 = P_abs_tot * pwr_frac[0]  # atomic deuterium
        self.beam_pwr_2 = P_abs_tot * pwr_frac[1]  # molecular deuterium D2
        self.beam_pwr_3 = P_abs_tot * pwr_frac[2]  # molecular deuterium D3

        self.beam_en_1 = inp.ebeam
        self.beam_en_2 = inp.ebeam / 2
        self.beam_en_3 = inp.ebeam / 3

        norm_const = 1./UnivariateSpline(rho_nbeams, hofr_1 * dvol / Volp_nbi, k=3, s=0).integral(0, 1.)

        dPdV_1_interp = UnivariateSpline(rho_nbeams,
                                         hofr_1 * self.beam_pwr_1 / Volp_nbi,
                                         k=3,
                                         s=0)

        dPdV_2_interp = UnivariateSpline(rho_nbeams,
                                         hofr_2 * self.beam_pwr_2 / Volp_nbi,
                                         k=3,
                                         s=0)

        dPdV_3_interp = UnivariateSpline(rho_nbeams,
                                         hofr_3 * self.beam_pwr_3 / Volp_nbi,
                                         k=3,
                                         s=0)

        dPdr_1_interp = UnivariateSpline(rho_nbeams,
                                         hofr_1 * self.beam_pwr_1 / Volp_nbi * dvol,
                                         k=3,
                                         s=0)

        dPdr_2_interp = UnivariateSpline(rho_nbeams,
                                         hofr_2 * self.beam_pwr_2 / Volp_nbi * dvol,
                                         k=3,
                                         s=0)

        dPdr_3_interp = UnivariateSpline(rho_nbeams,
                                         hofr_3 * self.beam_pwr_3 / Volp_nbi * dvol,
                                         k=3,
                                         s=0)

        zeta_1_interp = UnivariateSpline(rho_nbeams,
                                         ptch_angl1,
                                         k=3,
                                         s=0)

        zeta_2_interp = UnivariateSpline(rho_nbeams,
                                         ptch_angl2,
                                         k=3,
                                         s=0)

        zeta_3_interp = UnivariateSpline(rho_nbeams,
                                         ptch_angl3,
                                         k=3,
                                         s=0)

        self.dPdV_1 = dPdV_1_interp(core.rho)
        self.dPdV_2 = dPdV_2_interp(core.rho)
        self.dPdV_3 = dPdV_3_interp(core.rho)

        self.dPdr_1 = dPdr_1_interp(core.rho)
        self.dPdr_2 = dPdr_2_interp(core.rho)
        self.dPdr_3 = dPdr_3_interp(core.rho)

        self.dPdV_1_1D = dPdV_1_interp(core.rho[:, 0])
        self.dPdV_2_1D = dPdV_2_interp(core.rho[:, 0])
        self.dPdV_3_1D = dPdV_3_interp(core.rho[:, 0])


        # Introducing a normalization constant that seems to be getting lost in this calculation somewhere
        # TODO: Figure out what the hell is this

        self.dPdr_1_1D = dPdr_1_interp(core.rho[:, 0]) * norm_const
        self.dPdr_2_1D = dPdr_2_interp(core.rho[:, 0]) * norm_const
        self.dPdr_3_1D = dPdr_3_interp(core.rho[:, 0]) * norm_const

        # calculate deposition profile-weighted averaged zeta
        self.zeta_1 = np.sum(zeta_1_interp(np.linspace(0,1,100)) *
                             (1/self.beam_pwr_1) *
                             dPdr_1_interp(np.linspace(0,1,100)))
        self.zeta_2 = np.sum(zeta_1_interp(np.linspace(0,1,100)) *
                             (1/self.beam_pwr_1) *
                             dPdr_1_interp(np.linspace(0,1,100)))
        self.zeta_3 = np.sum(zeta_1_interp(np.linspace(0,1,100)) *
                             (1/self.beam_pwr_1) *
                             dPdr_1_interp(np.linspace(0,1,100)))

        self.debug={"P_abs_tot" : P_abs_tot,
                    "P_lst_tot"  : P_lst_tot,
                    "I_nbi_tot" : I_nbi_tot,
                    "I_nbi_eff" : I_nbi_eff,
                    "Volp_nbi"  : Volp_nbi,
                    "P_abs_1" : P_abs_1,
                    "hofr_1" : hofr_1,
                    "hofr_2" : hofr_2,
                    "hofr_3" : hofr_3,
                    "jnbtot" : jnbtot,
                    "pNBe" : pNBe,
                    "pNBi" : pNBi,
                    "nbfast" : nbfast,
                    "pressb" : pressb,
                    "pfusb" : pfusb,
                    "dvol" : dvol,
                    "dA" : dA}

class calc_nbi_vals:
    def __init__(self, inp, core):
        R0 = core.pts.axis.mag[0]
        R_tan = inp.rtang

        # center of central solenoid is at the origin
        # define vertical line at x=R_tan going from -y to +y, where
        # y = sqrt(R_tan^2 - (R_0(a) + a)^2)

        origin = Point(0,0)
        pt1 = Point(R_tan, -1*sqrt((core.R0_a + core.a)**2 - R_tan**2))
        pt2 = Point(R_tan, sqrt((core.R0_a + core.a)**2 - R_tan**2))

        beamline = LineString((pt1, pt2))
        # interpolate along half of the beam path (the second half will give identical values as the first half)
        # use a hundred points or so, just to get a smooth interpolation
        rho_beam = np.zeros(100)
        zeta_beam = np.zeros(100)
        for i,v in enumerate(np.linspace(0,0.5,100)):
            pt = beamline.interpolate(v, normalized=True)
            R = pt.distance(origin)
            rho_beam[i] = abs(R - R0) / core.a
            zeta_beam[i] = R_tan / R

        rho_zeta = np.column_stack((rho_beam,zeta_beam))
        rho_zeta_sorted = rho_zeta[rho_zeta[:,0].argsort()]

        #
        # A couple things need to be noted here:
        #   1.  If the radius of tangency of a beam is less than the R of the magnetic axis,
        #       then those inner flux surfaces that the beam hits four times will have two
        #       distinct zeta values for the IOL calculation. It hits the inboard side of the
        #       flux surface with a different direction cosine than when it hits the outbard
        #       side of the flux surface. This results in a bifurcated plot of zeta vs rho.
        #       You can see it with the following code:
        #
        #       plt.scatter(rho_zeta_sorted[:,0], rho_zeta_sorted[:,1])
        #
        #   2.  Although in principle, you could include the correct zeta values for each rho
        #       in the beam IOL calculation, the resulting F_orb, M_orb, and E_orb functions are
        #       not going to be terribly different from what they would be if used a single,
        #       appropriately averaged zeta value for the entire beam. This is partially because
        #       all of the zeta values are likely to be fairly high (i.e. ~> 0.8).
        #
        #   3.  If someone down the road wants to treat this more carefully, knock yourself out.
        #       At present, a single, deposition profile-weighted zeta will be calculated and used
        #       for the beam IOL calculation. If you change this, be sure to update these notes.
        #

        def moving_average_sortof(x, y, deltax):
            xmin = np.amin(x)
            xmax = np.amax(x)

            xnew = np.linspace(xmin,xmax,100)
            ynew = np.zeros(100)
            for i,xval in enumerate(xnew):
                ymax = np.amax(y[np.logical_and(x < xval+deltax, x > xval-deltax)])
                ymin = np.amin(y[np.logical_and(x < xval+deltax, x > xval-deltax)])
                ynew[i] = (ymax + ymin)/2
            return xnew, ynew

        rho_av_prof, zeta_av_prof = moving_average_sortof(rho_zeta_sorted[:,0], rho_zeta_sorted[:,1], 0.01)
        dPdr_norm1_interp = UnivariateSpline(inp.dPdr_norm1_data[:,0], inp.dPdr_norm1_data[:,1], k=3, s=0)

        zeta_av_val = np.sum(zeta_av_prof * dPdr_norm1_interp(rho_av_prof))

        # calculat dVdr, which is used in the calculation of dPdV
        dVdr = core.rho2vol.derivative()(core.rho)
        dVdr_1D = core.rho2vol.derivative()(core.rho[:,0])

        #Set final values for beam 1
        self.zeta_1 = zeta_av_val

        # assume dPdr_norm1 given with rho, like most other inputs
        self.dPdr_1 = inp.pbeam * dPdr_norm1_interp(core.rho)
        self.dPdr_1_1D = inp.pbeam * dPdr_norm1_interp(core.rho[:,0])

        self.dPdV_1 = self.dPdr_1 / dVdr
        self.dPdV_1_1D = self.dPdr_1_1D / dVdr_1D

        # TODO: add ability to have more than one beam when specifying deposition profiles as inputs.
        self.zeta_2 = zeta_av_val

        self.dPdr_2 = np.zeros(core.rho.shape)
        self.dPdr_2_1D = np.zeros(core.rho[:, 0].shape)

        self.dPdV_2 = np.zeros(core.rho.shape)
        self.dPdV_2_1D = np.zeros(core.rho[:, 0].shape)

        self.zeta_3 = zeta_av_val

        self.dPdr_3 = np.zeros(core.rho.shape)
        self.dPdr_3_1D = np.zeros(core.rho[:, 0].shape)

        self.dPdV_3 = np.zeros(core.rho.shape)
        self.dPdV_3_1D = np.zeros(core.rho[:, 0].shape)

        pwr_frac = np.array([1, 0, 0])
        self.beam_pwr_1 = inp.pbeam * pwr_frac[0]  # atomic deuterium
        self.beam_pwr_2 = inp.pbeam * pwr_frac[1]  # molecular deuterium D2
        self.beam_pwr_3 = inp.pbeam * pwr_frac[2]  # molecular deuterium D3

        self.beam_en_1 = inp.ebeam
        self.beam_en_2 = inp.ebeam / 2
        self.beam_en_3 = inp.ebeam / 3

class BeamDeposition:
    """"""
    
    def __init__(self, inp, core):
        sys.dont_write_bytecode = True
        # If deposition profile is provided as an input, use that.
        # Note, this is designed to read in the dP/dr(r) quantity, not dP/dV(r).
        # When integrated over r, this should equal P_beam, or at least the amount
        # that is deposited in the plasma.

        try:
            nbi_vals = calc_nbi_vals(inp, core)
            # If no deposition profile is provided, then prepare nbeams input file and run nbeams
        except:

            # set the name of the executable based on the operating system
            if os.name == 'nt':
                nbeams_name = 'nbeams.exe'
            elif os.name == 'posix':
                nbeams_name = 'nbeams'
            else:
                print 'not sure what os you\'re running. If mac, you might need to add some code \
                        to the beamdep module to help it find and run nbeams.'
                sys.exit()

            try:
                #prepare nbeams input file
                prep_nbi_infile(inp, core)

                # call nbeams. Note to those familiar with the old nbeams, I modified
                # the source code to take the input file as a commandline argument. - MH

                # If run as debug script, look at "/inputs" for inbeams_test.dat"
                if __name__ == "__main__":
                    try:
                        # try to find nbeams in the system path
                        p = Popen([nbeams_name, os.path.join(os.getcwd(), 'inputs', 'inbeams_test.dat')], stdin=PIPE, stdout=PIPE).wait()
                    except:
                        try:
                            # otherwise use the location specified in the input file
                            p = Popen([inp.nbeams_loc, os.path.join(os.getcwd(), 'inputs', 'inbeams_test.dat')], stdin=PIPE, stdout=PIPE).wait()
                        except:
                            print 'Unable to find nbeams executable. Stopping.'
                            sys.exit()
                else:
                    try:
                        # try to find nbeams in the system path
                        p = Popen([nbeams_name, os.getcwd()+'/inbeams.dat'], stdin=PIPE, stdout=PIPE).wait()
                    except:
                        try:
                            # otherwise use the location specified in the input file
                            p = Popen([inp.nbeams_loc, os.getcwd()+'/inbeams.dat'], stdin=PIPE, stdout=PIPE).wait()
                        except:
                            print 'Unable to find nbeams executable. Stopping.'
                            sys.exit()
                # instantiate class with nbeams output file information
                nbi_vals = read_nbi_outfile(inp, core)
            except:
                print 'unable to create beam deposition information. Stopping.'
                sys.exit()
        self.debug = nbi_vals.debug
        # create beams object
        self.beams = namedtuple('beam', 'D1 D2 D3')(
            # create D1 beam
            namedtuple('beam_D1', 'z m E P rtan dPdV dPdr zeta')(
                z_d,
                m_d,
                namedtuple('E', 'kev ev J')(
                    nbi_vals.beam_en_1,
                    nbi_vals.beam_en_1 * 1E3,
                    nbi_vals.beam_en_1 * 1E3 * 1.6021E-19
                ),
                namedtuple('P', 'MW W')(
                    nbi_vals.beam_pwr_1,
                    nbi_vals.beam_pwr_1 * 1E6
                ),
                inp.rtang,
                namedtuple('dPdV', 'v1D v2D')(
                    namedtuple('v1D', 'MW W')(
                        nbi_vals.dPdV_1_1D,
                        nbi_vals.dPdV_1_1D * 1E6
                    ),
                    namedtuple('v2D', 'MW W')(
                        nbi_vals.dPdV_1,
                        nbi_vals.dPdV_1 * 1E6
                    )
                ),
                namedtuple('dPdr', 'v1D v2D')(
                    namedtuple('v1D', 'MW W')(
                        nbi_vals.dPdr_1_1D,
                        nbi_vals.dPdr_1_1D * 1E6
                    ),
                    namedtuple('v2D', 'MW W')(
                        nbi_vals.dPdr_1,
                        nbi_vals.dPdr_1 * 1E6
                    )
                ),
                nbi_vals.zeta_1
            ),

            # create D2 beam
            namedtuple('beam_D2', 'z m E P rtan dPdV dPdr zeta')(
                z_d,
                m_d,
                namedtuple('E', 'kev ev J')(
                    nbi_vals.beam_en_2,
                    nbi_vals.beam_en_2 * 1E3,
                    nbi_vals.beam_en_2 * 1E3 * 1.6021E-19
                ),
                namedtuple('P', 'MW W')(
                    nbi_vals.beam_pwr_2,
                    nbi_vals.beam_pwr_2 * 1E6
                ),
                inp.rtang,
                namedtuple('dPdV', 'v1D v2D')(
                    namedtuple('v1D', 'MW W')(
                        nbi_vals.dPdV_2_1D,
                        nbi_vals.dPdV_2_1D * 1E6
                    ),
                    namedtuple('v2D', 'MW W')(
                        nbi_vals.dPdV_2,
                        nbi_vals.dPdV_2 * 1E6
                    )
                ),
                namedtuple('dPdr', 'v1D v2D')(
                    namedtuple('v1D', 'MW W')(
                        nbi_vals.dPdr_2_1D,
                        nbi_vals.dPdr_2_1D * 1E6
                    ),
                    namedtuple('v2D', 'MW W')(
                        nbi_vals.dPdr_2,
                        nbi_vals.dPdr_2 * 1E6
                    )
                ),
                nbi_vals.zeta_2
            ),

            # create D3 beam
            namedtuple('beam_D3', 'z m E P rtan dPdV dPdr zeta')(
                z_d,
                m_d,
                namedtuple('E', 'kev ev J')(
                    nbi_vals.beam_en_3,
                    nbi_vals.beam_en_3 * 1E3,
                    nbi_vals.beam_en_3 * 1E3 * 1.6021E-19
                ),
                namedtuple('P', 'MW W')(
                    nbi_vals.beam_pwr_3,
                    nbi_vals.beam_pwr_3 * 1E6
                ),
                inp.rtang,
                namedtuple('dPdV', 'v1D v2D')(
                    namedtuple('v1D', 'MW W')(
                        nbi_vals.dPdV_3_1D,
                        nbi_vals.dPdV_3_1D * 1E6
                    ),
                    namedtuple('v2D', 'MW W')(
                        nbi_vals.dPdV_3,
                        nbi_vals.dPdV_3 * 1E6
                    )
                ),
                namedtuple('dPdr', 'v1D v2D')(
                    namedtuple('v1D', 'MW W')(
                        nbi_vals.dPdr_3_1D,
                        nbi_vals.dPdr_3_1D * 1E6
                    ),
                    namedtuple('v2D', 'MW W')(
                        nbi_vals.dPdr_3,
                        nbi_vals.dPdr_3 * 1E6
                    )
                ),
                nbi_vals.zeta_3
            )
        )



class debugInp():
    def __init__(self):
        self.ebeam = 76.56296774
        self.pbeam = 4.88
        self.rtang = 1.715
        self.BT0 = 2.0
        self.nbeams_loc='/home/jonathan/Dropbox/GTEDGE/MyPrograms/GTEDGE/lib/beams/NBeamsMDS/NBeams/bin/Release/nbeams'

class debugCore():
    from core import calc_rho1d
    from core import calc_theta1d

    def __init__(self):

        rho1d= np.concatenate((np.linspace(0, 0.8, 10, endpoint=False), np.linspace(0.8, 1, 100, endpoint=True)))
        self.a = 0.5640316391827651
        self.R0_a = 1.691071481388474
        self.kappa = namedtuple('kappa', 'axis sep')(1.2507812426138631,1.6533128554966692)
        self.shaf_shift = 0.04273763321406595
        self.theta, self.rho = np.meshgrid([0., 1.], rho1d)
        loader=np.loadtxt("inputs/nbDebugProfs.dat")
        # This named tuple is generated by:
        #
        # np.savetxt('nbDebugProfs.dat', (core.n.e[:,0] / ((1.+.025*6.0)), core.n.e[:, 0], core.T.i.kev[:, 0], core.T.e.kev[:, 0]))
        #
        # in a module supporting calls to core. 9/4/2018 nbeams debug profiles obtained from shot 118890.1560

        self.n = namedtuple('n','i e')(np.array([loader[0],loader[0]]).transpose(),np.array([loader[1],loader[1]]).transpose())
        self.T = namedtuple('T','i e')(namedtuple('i', 'kev')(np.array([loader[2],loader[2]]).transpose()),
                                       namedtuple('e', 'kev')(np.array([loader[3],loader[3]]).transpose()))


        self.pts = namedtuple('pts','axis')(namedtuple('axis','mag')([1.7151768587049954,None]))
        self.kappa_vals = namedtuple('kappa_vals', 'axis sep')(1.2507812426138631, 1.6533128554966692)

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    inp=debugInp()
    core=debugCore()
    beamData=BeamDeposition(inp,core)
    print """Printing NBeams input parameters:
            Core
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            """.format("R0_a",str(core.R0_a),
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
            """.format("BT0",str(inp.BT0),
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
    fig0 = plt.figure(figsize=(12,8))
    fig0.suptitle(r"Deposition profiles")
    axi=230
    for a in beamData.debug.keys():
        if type(beamData.debug[a]) is float:
            print "{:<10}      {}".format(str(a),str(beamData.debug[a]))
        elif type(beamData.debug[a]) is np.ndarray:
            if "hofr" in a:
                axi+=1
                ax=fig0.add_subplot(axi)
                ax.set_title(a,fontsize=16)
                ax.set_ylabel("dep", fontsize=16)
                ax.set_xlabel(r"rho", fontsize=16)
                ax.plot(np.linspace(0, 1, 51), beamData.debug[a])
    d1Ptot = UnivariateSpline(core.rho[:,0],beamData.beams.D1.dPdr.v1D.W).integral(0.,1.)
    d2Ptot = UnivariateSpline(core.rho[:,0],beamData.beams.D2.dPdr.v1D.W).integral(0.,1.)
    d3Ptot = UnivariateSpline(core.rho[:,0],beamData.beams.D3.dPdr.v1D.W).integral(0.,1.)
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
    ax1.plot(core.rho[:, 0], core.n.i[:,0])

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'electron density', fontsize=16)
    ax2.set_ylabel(r'$n_e[1/{m^3}]$', fontsize=16)
    ax2.set_xlabel(r'rho', fontsize=16)
    ax2.plot(core.rho[:, 0], core.n.e[:,0])

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'Ion temp (kev)', fontsize=16)
    ax3.set_ylabel(r'$T_i[kev}$', fontsize=16)
    ax3.set_xlabel(r'rho', fontsize=16)
    ax3.plot(core.rho[:, 0], core.T.i.kev[:,0])

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'Electron temp (kev)', fontsize=16)
    ax4.set_ylabel(r'$T_e[kev}$', fontsize=16)
    ax4.set_xlabel(r'rho', fontsize=16)
    ax4.plot(core.rho[:, 0], core.T.e.kev[:,0])

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'Electron temp (kev)', fontsize=16)
    ax5.set_ylabel(r'$T_e[keV]$', fontsize=16)
    ax5.set_xlabel(r'rho', fontsize=16)
    ax5.plot(core.rho[:, 0], core.T.e.kev[:,0])

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

    plt.show()