#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 13:25:51 2018

@author: max
"""
from __future__ import division
from subprocess import call
from collections import namedtuple
import numpy as np
from scipy.interpolate import interp1d, interp2d, UnivariateSpline, interpn
import os
import re
import sys
import matplotlib.pyplot as plt


class ImpRad:
    """Calculates the impurity charge state densities and resulting impurity emissivities.

    Attributes:
        calc_nz_cs
        Lz
        dLzdT

    Methods:
        __init__
        calc_nz_cs
        calc_Lz
        frac_abun
    """
    def __init__(self, inp1, brnd):
        sys.dont_write_bytecode = True

        # TODO: Convert adpak routines to python binary using f2py so we can eliminate the main.f driver program

        # specify some additional input parameters. These could easily be moved into an input file if necessary.
        inp2_d = {}
        inp2_d['imode'] = 1
        inp2_d['nmin_T'] = -3
        inp2_d['nmax_T'] = 2
        inp2_d['nmin_n'] = 13
        inp2_d['nmax_n'] = 15
        inp2_d['nte'] = 25  # Number of electron temperatures in table
        inp2_d['nne'] = 2  # Number of electron densities in table
        inp2_d['tei'] = np.logspace(inp2_d['nmin_T'], inp2_d['nmax_T'], inp2_d['nte'])  # Array of nte Te values in keV
        inp2_d['anei'] = np.logspace(inp2_d['nmin_n'], inp2_d['nmax_n'], inp2_d['nne'])  # Array of nne ne values in cm^-3
        inp2_d['ncxb'] = 1  # Number of NB components (if 0 --> No NB's_)
        inp2_d['ivunit'] = 2  # Units for NB energy  1 --> cm/s  2 --> keV / amu

        # Some notes about the adpak interface:
        #   - Lz is a function of Te, nn/ne, and Tn
        #   - As a result, you don't need to run multiple ne values, however 2 are necessary for adpak to run smoothly
        #     for reasons that aren't totally clear. We will use ne=1E13 cm^-3 (1E19 m^3)
        #   - To get the values necessary to do a 3-D interpolation, we'll loop over nn/ne and Tn (in kev/amu)
        #   - Lz(nn/ne) doesn't change much for nn/ne > 1, which is a situation that only typically occurs in the PFR
        #       - we will loop over ~10 nn/ne ratios ranging from 0 to 1, i.e. (10^-6, 10^-5, ...etc.)
        #   - for a given nn/ne ratio, increasing Tn tends to INCREASE the effects of the neutrals
        #       - increasing Tn beyond ~100kev/amu has little effect on Lz, even for small nn/ne.
        #       - we will loop over  ~10 Tn values ranging from 0.001 to 100 kev/amu (i.e. 0.001, 0.01, 0.1, ...)

        Tn_num = 8
        nf_num = 8
        Te_num = 100

        Tn_vals = np.logspace(-3, 2, Tn_num)
        nf_vals = np.logspace(-7, 0, nf_num)
        Te_vals = np.logspace(-3, 2, Te_num)

        Lz_complete = np.zeros((Tn_num, nf_num, Te_num))
        dLzdT_complete = np.zeros((Tn_num, nf_num, Te_num))

        for i, nf in enumerate(nf_vals):
            for j, Tn in enumerate(Tn_vals):

                inp2_d['anneut'] = nf * 1*10**inp2_d['nmin_n']
                inp2_d['vneut'] = Tn

                # combine inp1 and inp2 as inp
                try:
                    inp1_d = inp1._asdict()
                except AttributeError:
                    inp1_d = inp1.__dict__
                inp_d = inp2_d.copy()
                inp_d.update(inp1_d)
                inp = namedtuple('inp', inp_d.keys())(*inp_d.values())

                # prepare input file
                f = open('./toadpak', 'w')
                f.write(' &inp')
                f.write('\n' + '  inucz = ' + str(inp.z_imp))
                f.write('\n' + '  zte = ' + str(0))  # this parameter doesn't matter because imode = 1
                f.write('\n' + '  zne = ' + str(0))  # this parameter doesn't matter because imode = 1
                f.write('\n' + '  laden = ' + str(inp.laden))
                f.write('\n' + '  ladtip = ' + str(inp.ladtip))
                f.write('\n' + '  leci = ' + str(inp.leci))
                f.write('\n' + '  ldrmlt = ' + str(inp.ldrmlt))
                f.write('\n' + '  ncxb = ' + str(inp.ncxb))
                f.write('\n' + '  ncxopt = ' + str(inp.ncxopt))
                f.write('\n' + '  ivunit = ' + str(inp.ivunit))
                f.write('\n' + '  anneut = ' + str(inp.anneut))
                f.write('\n' + '  vneut = ' + str(inp.vneut))
                f.write('\n' + '  imode = ' + str(inp.imode))
                f.write('\n' + '  nte = ' + str(inp.nte))
                f.write('\n' + '  nne = ' + str(inp.nne))
                f.write('\n' + '  tei = ' + ' '.join(map(str, inp.tei)))
                f.write('\n' + '  anei = ' + ' '.join(map(str, inp.anei)))
                f.write('\n' + '  nmin = ' + str(inp.nmin_T))
                f.write('\n' + '  nmax = ' + str(inp.nmax_T))
                f.write('\n' + ' $end')
                f.write('\n' + '')
                f.write('\n')
                f.close()

                # call adpak
                try:
                    call([inp.adpak_loc+'adpak', os.getcwd()+'/toadpak'])
                except AttributeError:
                    try:
                        call(['adpak', os.getcwd()+'/toadpak'])
                    except:
                        print 'Unable to locate adpak. Stopping.'
                        sys.exit()

                # read adpak output
                with open(os.getcwd() + '/outblk.dat', 'r') as f:
                    # read entire file into one long string variable
                    data = f.read().replace('\n', ' ').replace('. ', ' ').replace(', ', ' ')

                # the following arrays all have the same dimension despite the fact that alinzr
                # is a function only of Te
                # array[ charge state, density value, temperature value ]
                alinzr = np.zeros((inp.z_imp + 1, inp.nne, inp.nte))
                alradr = np.zeros((inp.z_imp + 1, inp.nne, inp.nte))
                alrecr = np.zeros((inp.z_imp + 1, inp.nne, inp.nte))

                # find altei data
                result = re.match(r'.*?data \(+%s.*?\/(.*?)\/.*' % 'altei', data).group(1)
                altei = np.asarray(result.split(), dtype=float)

                # find alnei data
                result = re.match(r'.*?data \(+%s.*?\/(.*?)\/.*' % 'alnei', data).group(1)
                alnei = np.asarray(result.split(), dtype=float)

                # for each charge state plus the ground state,
                for cs in np.arange(inp.z_imp + 1):
                    # find alinzr data for cs. alinzr is a function of temperature only
                    result = re.match(r'.*?data \(+%s\(\d+ +%s.*?\/(.*?)\/.*' % ('alinzr', cs+1), data).group(1)
                    alinzr[cs, :, :] = np.asarray(result.split(), dtype=float).reshape(-1, inp.nte)

                    # find alradr data for cs. alradr is a function of both temperature and electron density
                    result = re.match(r'.*?data \(+%s\(\d+ +%s.*?\/(.*?)\/.*' % ('alradr', cs+1), data).group(1)
                    alradr[cs, :, :] = np.asarray(result.split(), dtype=float).reshape(-1, inp.nte)

                    # find alrecr data for cs. alrecr is a function of both temperature and electron density
                    result = re.match(r'.*?data \(+%s\(\d+ +%s.*?\/(.*?)\/.*' % ('alrecr', cs+1), data).group(1)
                    alrecr[cs, :, :] = np.asarray(result.split(), dtype=float).reshape(-1, inp.nte)

                # assemble data to pass to Lz and charge state density calculations
                d = {
                    'altei': altei,
                    'alnei': alnei,
                    'alinzr': alinzr,
                    'alradr': alradr,
                    'alrecr': alrecr
                    }

                # d['altei'] = altei
                # d['alnei'] = alnei
                # d['alinzr'] = alinzr
                # d['alradr'] = alradr
                # d['alrecr'] = alrecr

                # convert dictionary to named tuple so keys can be accessed as normal attributes
                data = namedtuple("data", d.keys())(*d.values())

                # calculate impurity charge state densities. Uncomment the following line if desired.
                # self.calc_nz_cs = self.nz_cs(inp, data)

                # create interpolation function for Lz and dLzdT based on coronal equilibrium approximation
                Lz_interp, dLzdT_interp = self.calc_Lz(inp, data)
                Lz_complete[j, i, :] = Lz_interp(Te_vals*1E3*1.6021E-19)
                dLzdT_complete[j, i, :] = dLzdT_interp(Te_vals*1E3*1.6021E-19)

        # Create 3-D interpolation function for Lz_complete and dLzdT_complete arrays

        # Although this module can be run as a standalone code, it's designed for use with GT3. As such it needs
        # to calculate Lz and dLzdT for two different neutral populations, "n_n_slow" and "n_n_thermal." For the
        # situations in which it is being run as a standalone code, it will also look for an array "n_n."
        # If n_n_slow or n_n_thermal are successful, the code will also attempt to update the corresponding Lz
        # and dLzdT values in the brnd instance.

        #first look for brnd.n_n_slow
        try:
            brnd_vals = np.column_stack((np.asarray(brnd.Tn).flatten().T,
                                         (np.asarray(brnd.n_n_slow) / np.asarray(brnd.ne)).flatten().T,
                                         np.asarray(brnd.Te_keV).flatten().T))
        except:
            #brnd.n_n_slow was not found
            pass

        self.Lz = interpn((Tn_vals, nf_vals, Te_vals*1E3*1.6021E-19), Lz_complete, brnd_vals, method='linear')
        self.dLzdT = interpn((Tn_vals, nf_vals, Te_vals*1E3*1.6021E-19), dLzdT_complete, brnd_vals, method='linear')

    @classmethod
    def calc_nz_cs(cls, inp, data):
        """calculates impurity densities by charge state"""

        tei_lin = np.linspace(inp.nmin, inp.nmax, inp.nte)

        nz_cs = []

        # for each charge state,
        for cs in np.arange(inp.z_imp+1):
            # get fractional abundances of each as a function of ne and Te
            frac_abun_interp = interp2d(tei_lin,
                                        inp.anei,
                                        cls.frac_abun(inp, data)[:, :, cs])
            
            # get the fractional abundance for each value in the nz array
            nz_frac_abun = np.zeros(brnd.nz.shape)
            for (i, j), brnd.nz in np.ndenumerate(brnd.nz):
                nz_frac_abun[i, j] = frac_abun_interp(np.log10(brnd.Te_kev[i, j]),
                                                      brnd.ne[i, j])
            
            # multiply the fractional abundance by the overall impurity density to get the
            # density of this charge state
            nz_cs.append(nz_frac_abun * brnd.nz)

        return nz_cs

    @classmethod
    def calc_Lz(cls, inp, data):
        """creates the emissivity function Lz over a range of temperatures"""

        def kev_J(x):
            return x*1E3*1.6021E-19

        tei_lin = np.linspace(inp.nmin_T, inp.nmax_T, inp.nte)
        Lz_tot = np.zeros(inp.nte)

        # for each charge state,
        for cs in np.arange(inp.z_imp + 1):
            # get fractional abundances of each as a function of ne and Te
            frac_abun_interp = interp1d(tei_lin,
                                        cls.frac_abun(inp, data)[0, :, cs])

            Lz_cs_interp = interp1d(tei_lin,
                                    data.alradr[cs, 0, :])

            for i, T in enumerate(tei_lin):
                Lz_tot[i] = Lz_tot[i] + 10**Lz_cs_interp(T) * frac_abun_interp(T)

        # spline fit the logarithm of the emissivity. We'll spline fit the actual values next
        Lz_tot_interp = UnivariateSpline(inp.tei, np.log10(Lz_tot)+100.0, s=0)

        # the above function gives the right values, but we also need the derivatives
        # we will now do a spline fit of the values themselves, rather than their base-10 logarithm
        T_kev_vals = np.logspace(-3, 2, 10000)
        T_J_vals = kev_J(T_kev_vals)
        Lz_vals = 10.0**(Lz_tot_interp(T_kev_vals)-100.0)

        Lz_interp = UnivariateSpline(T_J_vals, Lz_vals, s=0)
        dLzdT_interp = Lz_interp.derivative()

        return Lz_interp, dLzdT_interp

    @staticmethod
    def frac_abun(inp, data):
        """Calculates the fractional abundances of impurity charge states using a coronal equilibrium approximation."""

        def null(M, eps=1e-20):
            u, s, vh = np.linalg.svd(M)
            return vh[np.argmin(np.abs(s))]

        # Create nchrgsr x nchrgsr array
        M = np.zeros((inp.z_imp + 1, inp.z_imp + 1))

        frac_abun = np.zeros((inp.nne, inp.nte, inp.z_imp + 1))

        for i, dens in enumerate(inp.anei):
            for j, temp in enumerate(inp.tei):
                for (cs, col), value in np.ndenumerate(M):
                    if cs == 0:  # do the first row differently
                        M[cs, cs] = -(10 ** data.alinzr[cs, i, j] + 10 ** data.alrecr[cs, i, j])
                        M[cs, cs + 1] = 10 ** data.alrecr[cs + 1, i, j]
                    elif cs == inp.z_imp:  # do the last row differently
                        M[cs, cs - 1] = 10 ** data.alinzr[cs - 1, i, j]
                        M[cs, cs] = -(10 ** data.alinzr[cs, i, j] + 10 ** data.alrecr[cs, i, j])
                    else:
                        M[cs, cs - 1] = 10 ** data.alinzr[cs - 1, i, j]
                        M[cs, cs] = -(10 ** data.alinzr[cs, i, j] + 10 ** data.alrecr[cs, i, j])
                        M[cs, cs + 1] = 10 ** data.alrecr[cs + 1, i, j]

                # the following line may help convergence. Uncomment it if necessary. Everything
                # gets normalized anyway.
                # M = M*1E10

                result = np.abs(null(M).T)  # solve(M, b)
                frac_abun[i, j, :] = result / np.sum(result)

        return frac_abun


if __name__ == '__main__':
    """Example script."""
    # input parameters. Note that some parameters are hard coded in ImpRad.__init__
    z_imp = 6       # Atomic number, Z, of desired element
    laden = 1       # Flag for excitation energy and ion. potential calc. If 0 --> use Mayer formalism. If 1 --> use More formalism
    ladtip = 1      # If 1 --> use tabulated ionization potentials
    leci = 1        # Flag for calculation of ionization rates: 1: XSNQ, 2: Belfast group (for H through O) 3: Younger (Scandium and Fe)
    ldrmlt = 2      # Dielectronic multiplier. 0: Use CDNN, CDNM arrays as given# 1 --> Set CDNN, CDNM equal to 1, 2 --> Use Y. Hahn factors (dut) to set up CDNN, CDNM
    ncxopt = 1      # Selects cross sections to be used   1 --> OSAS   2 --> GJ    3 --> OSCT

    # convert input parameters to dictionary. Will make a namedtuple after nn and Tn are added to the dictionary.
    d_inp = {}
    d_inp['z_imp'] = z_imp
    d_inp['laden'] = laden
    d_inp['ladtip'] = ladtip
    d_inp['leci'] = leci
    d_inp['ldrmlt'] = ldrmlt
    d_inp['ncxopt'] = ncxopt

    # specify density and temperature parameters at which you want to
    # evaluate Lz, dLzdT, etc. and convert to dict and namedtuple
    Te_kev = np.logspace(-3, 2, 100)

    d_brnd = {}
    d_brnd['Te_kev'] = Te_kev
    brnd = namedtuple('brnd', d_brnd.keys())(*d_brnd.values())

    #prepare figure if necessary
    Lz_fig = plt.figure(figsize=(6, 4))
    ax1 = Lz_fig.add_subplot(1, 1, 1)
    ax1.set_title('Impurity Radiation Emissivity as a Function of Temperature')
    ax1.set_xlabel(r'Electron Temperature ($keV$)')
    ax1.set_ylabel(r'Carbon Radiation Emissivity ($W*m^3$)')

    # run adpak for each neutral density over a range.
    # there might be a way to pass adpak an array of neutral densities, but I can't make it work.
    #for i, (nn, Tn) in enumerate(zip(anneut, vneut)):

    # finish the d_inp dictionary with the nn and Tn and convert to namedtuple
    # This allows keys to be accessed as normal attributes

    inp = namedtuple('inp', d_inp.keys())(*d_inp.values())

    # instantiate ImpRad
    inst = ImpRad(inp, brnd)

    # plot
    #ax1.loglog(brnd.Te_kev, inst.Lz)

    # clean up and show plot
    plt.tight_layout()
    plt.show()
