#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 21:51:19 2018

@author: max
"""
from __future__ import division
import numpy as np
from scipy.interpolate import interp1d
from subprocess import call
import os
import re
import sys
from scipy.interpolate import interp1d
from collections import namedtuple
from scipy.constants import physical_constants

m_d = physical_constants['deuteron mass'][0]


def calc_pwr_frac(ebeam):
    # calculate pwr fracs. These fits are only useful for DIII-D

    pwr_frac = np.zeros(3)
    pwr_frac[0] = (68 + 0.11 * ebeam) / 100
    pwr_frac[1] = (-159 + 6.53 * ebeam - 0.082 * (ebeam ** 2) + 0.00034 * (ebeam ** 3)) / 100
    pwr_frac[2] = (191 - 6.64 * ebeam + 0.082 * (ebeam ** 2) - 0.00034 * (ebeam ** 3)) / 100

    return pwr_frac


def prep_nbi_infile(inp, core):
    pwr_frac = calc_pwr_frac(inp.ebeam)

    # f1=open(inp.nbeams_loc+"inbeams.dat", "w")
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


class BeamDeposition:
    """"""
    
    def __init__(self, inp, core):
        sys.dont_write_bytecode = True 
        prep_nbi_infile(inp, core)
        #call nbeams. Note to those familiar with the old nbeams, I modified
        #the source code to take the input file as a commandline argument. - MH
        try:
            call([inp.nbeams_loc+'nbeams', os.getcwd()+'/inbeams.dat'])
        except:
            try:
                call(['nbeams', os.getcwd()+'/inbeams.dat'])
            except:
                print 'Unable to find nbeams executable. Stopping.'
                sys.exit()
        self.read_nbi_outfile(inp, core)
        pass
    

    
    def read_nbi_outfile(self, inp, brnd):
        with open(os.getcwd() + '/outbeams.dat', 'r') as f:
            for count, line in enumerate(f):
                if line.startswith(" Total Absorbed Power"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.P_abs_tot = float(result)
                    except:
                        self.P_abs_tot = np.NaN

                if line.startswith(" Total Lost Power"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.P_lst_tot = float(result)
                    except:
                        self.P_lst_tot = np.NaN

                if line.startswith(" Total NB Driven Current"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.I_nbi_tot = float(result)
                    except:
                        self.I_nbi_tot = np.NaN
                        
                if line.startswith(" Total NBCD Efficiency"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.I_nbi_eff = float(result)
                    except:
                        self.I_nbi_eff = np.NaN
                        

                if line.startswith(" Total Beam Beta"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.Beta_nbi_tot = float(result)
                    except:
                        self.Beta_nbi_tot = np.NaN
                        
                if line.startswith(" Taus"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.Taus_nbi = float(result)
                    except:
                        self.Taus_nbi = np.NaN
                        
                if line.startswith(" Volp"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.Volp_nbi = float(result)
                    except:
                        self.Volp_nbi = np.NaN

                if line.startswith(" Absorbed Power"):
                    #this will need to be modified for multiple beams
                    print 'line = ', line
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.P_abs_1 = float(result)
                    except:
                        self.P_abs_1 = np.NaN
                        
                if line.startswith(" Lost Power"):                   
                    #this will need to be modified for multiple beams
                    
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.P_lst_1 = float(result)
                    except:
                        self.P_lst_1 = np.NaN
                        
                if line.startswith(" NB driven current"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.I_nbi_1 = float(result)
                    except:
                        self.I_nbi_1 = np.NaN
                        
                if line.startswith(" NBCD efficiency"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.I_nbi_eff_1 = float(result)
                    except:
                        self.I_nbi_eff_1 = np.NaN
                        
                if line.startswith(" NBCD gamma"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.I_nbi_gam_1 = float(result)
                    except:
                        self.I_nbi_gam_1 = np.NaN
                        
                if line.startswith("    energy group 1"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.st_en1_1 = float(result)
                    except:
                        self.st_en1_1 = np.NaN
                        
                if line.startswith("    energy group 2"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.st_en2_1 = float(result)
                    except:
                        self.st_en2_1 = np.NaN
                        
                if line.startswith("    energy group 3"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.st_en3_1 = float(result)
                    except:
                        self.st_en3_1 = np.NaN
                        
                if line.startswith(" Total Beam-Target Fusion Power"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.fus_pwr_bt = float(result)
                    except:
                        self.fus_pwr_bt = np.NaN
                        
                if line.startswith(" Total Power to Charged Particles"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.cp_pwr_tot = float(result)
                    except:
                        self.cp_pwr_tot = np.NaN
                        
                if line.startswith(" Total DT Neutron Rate"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.rate_dt_n = float(result)
                    except:
                        self.rate_dt_n = np.NaN
                        
                if line.startswith(" Total DD Neutron Rate"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*', line).group(1)
                    try:
                        self.rate_dd_n = float(result)
                    except:
                        self.rate_dd_n = np.NaN
                        
        with open(os.getcwd() + '/outbeams.dat', 'r') as f:
            data = f.read().replace('\n', ' ')
            #print data
            result = re.match(r'.*hofr_3 *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) +Pitch.*', data).group(1)
            array = np.reshape(np.asarray(result.split(), dtype=float), (-1, 4))
            self.dep_prof1 = interp1d(array[:, 0], array[:, 1])(brnd.rho)
            self.dep_prof2 = interp1d(array[:, 0], array[:, 2])(brnd.rho)
            self.dep_prof3 = interp1d(array[:, 0], array[:, 3])(brnd.rho)
            self.dep_prof1_1D = self.dep_prof1[:, 0]
            self.dep_prof2_1D = self.dep_prof2[:, 0]
            self.dep_prof3_1D = self.dep_prof3[:, 0]

            result = re.match(r'.*zeta_3 *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) +rho.*', data).group(1)
            array = np.reshape(np.asarray(result.split(), dtype=float), (-1, 4))
            self.ptch_angl1 = interp1d(array[:, 0], array[:, 1])(brnd.rho)
            self.ptch_angl2 = interp1d(array[:, 0], array[:, 2])(brnd.rho)
            self.ptch_angl3 = interp1d(array[:, 0], array[:, 3])(brnd.rho)
            self.ptch_angl1_1D = self.ptch_angl1[:, 0]
            self.ptch_angl2_1D = self.ptch_angl2[:, 0]
            self.ptch_angl3_1D = self.ptch_angl3[:, 0]

            result = re.match(r'.*dA *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) *', data).group(1)
            array = np.reshape(np.asarray(result.split(), dtype=float), (-1, 9))
            self.jnbtot = interp1d(array[:, 0], array[:, 1])(brnd.rho)
            self.pNBe = interp1d(array[:, 0], array[:, 2])(brnd.rho)  # in MW/m^3
            self.pNBi = interp1d(array[:, 0], array[:, 3])(brnd.rho)  # in MW/m^3
            self.pNB_tot = self.pNBe + self.pNBi
            self.pNBe_dvol = interp1d(array[:, 0], array[:, 2] * array[:, 7])(brnd.rho)  # in MW
            self.pNBi_dvol = interp1d(array[:, 0], array[:, 3] * array[:, 7])(brnd.rho)  # in MW
            self.nbfast = interp1d(array[:, 0], array[:, 4])(brnd.rho)
            self.pressb = interp1d(array[:, 0], array[:, 5])(brnd.rho)
            self.pfusb = interp1d(array[:, 0], array[:, 6])(brnd.rho)
            self.dvol = interp1d(array[:, 0], array[:, 7])(brnd.rho)
            self.dA = interp1d(array[:, 0], array[:, 8])(brnd.rho)
            self.jnbtot_1D = self.jnbtot[:, 0]
            self.pNBe_1D = self.pNBe[:, 0]
            self.pNBi_1D = self.pNBi[:, 0]
            self.pNB_tot_1D = self.pNB_tot[:, 0]
            self.pNBe_dvol_1D = self.pNBe_dvol[:, 0]
            self.pNBi_dvol_1D = self.pNBi_dvol[:, 0]
            self.nbfast_1D = self.nbfast[:, 0]
            self.pressb_1D = self.pressb[:, 0]
            self.pfusb_1D = self.pfusb[:, 0]
            self.dvol_1D = self.dvol[:, 0]
            self.dA_1D = self.dA[:, 0]

            pwr_frac = calc_pwr_frac(inp.ebeam)
            beam_D_pwr = pwr_frac[0]
            beam_D2_pwr = pwr_frac[0]
            beam_D3_pwr = pwr_frac[0]

            # create beams object
            self.beams = namedtuple('beam', 'D D2 D3')(
                #create D beam
                namedtuple('beam_D', 'z m a num e p rtan dp')(
                    # num = 1 for atomic deuterium - find better way to do this
                    1,
                    m_d,
                    2,
                    1,
                    namedtuple('e', 'kev ev J')(
                        inp.ebeam,
                        inp.ebeam * 1E3,
                        inp.ebeam * 1E3 * 1.6021E-19
                    ),
                    namedtuple('p', 'MW W')(
                        beam_D_pwr,
                        beam_D_pwr * 1E6
                    ),
                    inp.rtang,
                    self.dep_prof1 * self.dvol / self.Volp_nbi
                ),

                # create D2 beam
                namedtuple('beam_D2', 'z m a num e p rtan dp')(
                    # num = 1 for atomic deuterium - find better way to do this
                    1,
                    2 * m_d,
                    4,
                    2,
                    namedtuple('e', 'kev ev J')(
                        inp.ebeam,
                        inp.ebeam * 1E3,
                        inp.ebeam * 1E3 * 1.6021E-19
                    ),
                    namedtuple('p', 'MW W')(
                        beam_D2_pwr,
                        beam_D2_pwr * 1E6
                    ),
                    inp.rtang,
                    self.dep_prof2 * self.dvol / self.Volp_nbi
                ),

                # create D3 beam
                namedtuple('beam_D2', 'z m a num e p rtan dp')(
                    # num = 1 for atomic deuterium - find better way to do this
                    1,
                    3 * m_d,
                    6,
                    3,
                    namedtuple('e', 'kev ev J')(
                        inp.ebeam,
                        inp.ebeam * 1E3,
                        inp.ebeam * 1E3 * 1.6021E-19
                    ),
                    namedtuple('p', 'MW W')(
                        beam_D3_pwr,
                        beam_D3_pwr * 1E6
                    ),
                    inp.rtang,
                    self.dep_prof3 * self.dvol / self.Volp_nbi
                )
            )

            # create beams1D object
            self.beams_1D = namedtuple('beam', 'D D2 D3')(
                #create D beam
                namedtuple('beam_D', 'z m a num e p rtan dp')(
                    # num = 1 for atomic deuterium - find better way to do this
                    1,
                    m_d,
                    2,
                    1,
                    namedtuple('e', 'kev ev J')(
                        inp.ebeam,
                        inp.ebeam * 1E3,
                        inp.ebeam * 1E3 * 1.6021E-19
                    ),
                    namedtuple('p', 'MW W')(
                        beam_D_pwr,
                        beam_D_pwr * 1E6
                    ),
                    inp.rtang,
                    self.dep_prof1_1D * self.dvol_1D / self.Volp_nbi
                ),

                # create D2 beam
                namedtuple('beam_D2', 'z m a num e p rtan dp')(
                    # num = 1 for atomic deuterium - find better way to do this
                    1,
                    2 * m_d,
                    4,
                    2,
                    namedtuple('e', 'kev ev J')(
                        inp.ebeam,
                        inp.ebeam * 1E3,
                        inp.ebeam * 1E3 * 1.6021E-19
                    ),
                    namedtuple('p', 'MW W')(
                        beam_D2_pwr,
                        beam_D2_pwr * 1E6
                    ),
                    inp.rtang,
                    self.dep_prof2_1D * self.dvol_1D / self.Volp_nbi
                ),

                # create D3 beam
                namedtuple('beam_D2', 'z m a num e p rtan dp')(
                    # num = 1 for atomic deuterium - find better way to do this
                    1,
                    3 * m_d,
                    6,
                    3,
                    namedtuple('e', 'kev ev J')(
                        inp.ebeam,
                        inp.ebeam * 1E3,
                        inp.ebeam * 1E3 * 1.6021E-19
                    ),
                    namedtuple('p', 'MW W')(
                        beam_D3_pwr,
                        beam_D3_pwr * 1E6
                    ),
                    inp.rtang,
                    self.dep_prof3_1D * self.dvol_1D / self.Volp_nbi
                )
            )