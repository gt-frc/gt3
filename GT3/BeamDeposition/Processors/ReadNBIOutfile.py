#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import re, os
import numpy as np
from GT3.BeamDeposition.Functions.CalcPwrFrac import calc_pwr_frac
from scipy.interpolate import interp1d, UnivariateSpline


class read_nbi_outfile:

    def __init__(self, inp, core, index = False, indBeam = False):
        if index:
            with open(os.getcwd() + '/outbeams_%s.dat' % str(index), 'r') as f:
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

            with open(os.getcwd() + '/outbeams_%s.dat' % str(index), 'r') as f:
                data = f.read().replace('\n', ' ')
        else:
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



        if index:
            pwr_frac = calc_pwr_frac(indBeam['ebeam'])
        else:
            pwr_frac = calc_pwr_frac(inp.ebeam)
        self.beam_pwr_1 = P_abs_tot * pwr_frac[0]  # atomic deuterium
        self.beam_pwr_2 = P_abs_tot * pwr_frac[1]  # molecular deuterium D2
        self.beam_pwr_3 = P_abs_tot * pwr_frac[2]  # molecular deuterium D3

        if index:
            self.beam_en_1 = indBeam['ebeam']
            self.beam_en_2 = indBeam['ebeam'] / 2
            self.beam_en_3 = indBeam['ebeam'] / 3
        else:
            self.beam_en_1 = inp.ebeam
            self.beam_en_2 = inp.ebeam / 2
            self.beam_en_3 = inp.ebeam / 3



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
                                         (hofr_1 * dvol * self.beam_pwr_1 )/ (.02 * Volp_nbi),
                                         k=3,
                                         s=0)

        dPdr_2_interp = UnivariateSpline(rho_nbeams,
                                         hofr_2  * dvol * self.beam_pwr_2 / (.02 * Volp_nbi),
                                         k=3,
                                         s=0)

        dPdr_3_interp = UnivariateSpline(rho_nbeams,
                                         hofr_3  * dvol * self.beam_pwr_3 / (.02 * Volp_nbi),
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

        self.dPdr_1_1D = dPdr_1_interp(core.rho[:, 0])
        self.dPdr_2_1D = dPdr_2_interp(core.rho[:, 0])
        self.dPdr_3_1D = dPdr_3_interp(core.rho[:, 0])

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
                    "dA" : dA,
                    "rho_nbeams" : rho_nbeams}