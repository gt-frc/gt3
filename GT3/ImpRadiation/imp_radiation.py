#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 13:25:51 2018

@author: max
"""

from subprocess import call
from collections import namedtuple
import numpy as np
from scipy.interpolate import interp1d, interp2d, UnivariateSpline, LinearNDInterpolator
import os
import re
import sys
import pickle
from GT3 import Core


class ImpRad:
    """Calculates the impurity charge state densities and resulting impurity emissivities.

    Attributes:
        Lz
        dLzdT

    Methods:
        __init__
        calc_Lz
        frac_abun
    """

    def __init__(self, z=None, core=None):
        """

        :param z:
        :param core: Core
        """
        sys.dont_write_bytecode = True

        # list of impurity names and pickled impurity interpolation object filenames, if they exist. Will be used later.
        imp_names = {}
        imp_names[2] = 'Helium'
        imp_names[4] = 'Beryllium'
        imp_names[6] = 'Carbon'
        imp_names[8] = 'Oxygen'
        imp_names[10] = 'Neon'
        imp_names[18] = 'Argon'
        imp_names[36] = 'Krypton'
        imp_names[54] = 'Xenon'
        imp_names[74] = 'Tungsten'

        if (core is not None) and (z is not None):
            print("""Both core and z were specified. z will be ignored and core Lz parameters will be updated.'
                  'If you would like a specific z only, then instantiate ImpRad as a standalone instance and'
                  'pass z, but not a core instance.""")
        if core is not None:
            update_dict = {4: core.Lz.update_Be,
                           6: core.Lz.update_C,
                           74: core.Lz.update_W,
                           10: core.Lz.update_Ne,
                           18: core.Lz.update_Ar,
                           36: core.Lz.update_Kr}
            for z in update_dict.keys():  # list of all z elements in core. Update this list as necessary.
                try:  # before running adpak, try to find a pickled interpolator somewhere in the main directory.
                    Lz, dLzdT = self.find_interp(z, imp_names)
                except:
                    print('Pickled interpolater not found for {}. Running adpak.'.format(imp_names[z]))
                    Lz, dLzdT = self.run_adpak(z, imp_names)
                update_dict[z](core.n, core.T, Lz, dLzdT)
            cool_rate = core.n.e * core.n.C * np.nan_to_num(core.Lz.C.s) + core.n.e * core.n.C * np.nan_to_num(
                core.Lz.C.t)
            core.cool_rate.update(cool_rate)



        elif z is not None:
            try:  # before running adpak, try to find a pickled interpolator somewhere in the main directory.
                print('trying to find interpolation object')
                self.Lz, self.dLzdT = self.find_interp(z, imp_names)
            except:
                print('didn\'t work')
                print('Pickled interpolater not found for {}. Running adpak.'.format(imp_names[z]))
                self.Lz, self.dLzdT = self.run_adpak(z, imp_names)
        else:
            raise Exception("Neither core nor z were specified. I can\'t read minds. Stopping.")

        # TODO: Convert adpak routines to python binary using f2py so we can eliminate the main.f driver program

    def find_interp(self, z, imp_names):
        # create filename
        Lz_pkl_file = imp_names[z] + '_Lz.pkl'
        dLzdT_pkl_file = imp_names[z] + '_dLzdT.pkl'

        # search for pkl_file
        outfile_found = 0
        for root, subdirs, files in os.walk(os.getcwd()):
            for filename in files:
                if filename == Lz_pkl_file:
                    Lz_found = True
                    os.path.join(root, filename)
                    pkl_file_loc = os.path.join(root, filename)
                    pickle_in = open(pkl_file_loc, "rb")
                    Lz_interp = pickle.load(pickle_in, encoding='latin1')
                    pickle_in.close()
                if filename == dLzdT_pkl_file:
                    dLzdT_found = True
                    os.path.join(root, filename)
                    pkl_file_loc = os.path.join(root, filename)
                    pickle_in = open(pkl_file_loc, "rb")
                    dLzdT_interp = pickle.load(pickle_in, encoding='latin1')
                    pickle_in.close()
        if Lz_found and dLzdT_found:
            return Lz_interp, dLzdT_interp
        else:
            # pickle file was not found. Bummer.
            raise ValueError('One or more of the pickle files was not found for {}.'.format(imp_names[z]))

    def run_adpak(self, z, imp_names):

        # specify input parameters. These could easily be moved into an input file if necessary.
        inp = {}
        inp['z_imp'] = z  # impurity z
        inp['laden'] = laden  # Flag for excitation energy and ion. potential calc. If 0 --> use Mayer formalism. If 1 --> use More formalism
        inp['ladtip'] = ladtip  # If 1 --> use tabulated ionization potentials
        inp[
            'leci'] = leci  # Flag for calculation of ionization rates: 1: XSNQ, 2: Belfast group (for H through O) 3: Younger (Scandium and Fe)
        inp[
            'ldrmlt'] = ldrmlt  # Dielectronic multiplier. 0: Use CDNN, CDNM arrays as given# 1 --> Set CDNN, CDNM equal to 1, 2 --> Use Y. Hahn factors (dut) to set up CDNN, CDNM
        inp['ncxopt'] = ncxopt  # Selects cross sections to be used   1 --> OSAS   2 --> GJ    3 --> OSCT
        inp['imode'] = 1
        inp['nmin_T'] = -3
        inp['nmax_T'] = 2
        inp['nmin_n'] = 13
        inp['nmax_n'] = 15
        inp['nte'] = 25  # Number of electron temperatures in table
        inp['nne'] = 2  # Number of electron densities in table
        inp['tei'] = np.logspace(inp['nmin_T'], inp['nmax_T'], inp['nte'])  # Array of nte Te values in keV
        inp['anei'] = np.logspace(inp['nmin_n'], inp['nmax_n'], inp['nne'])  # Array of nne ne values in cm^-3
        inp['ncxb'] = 1  # Number of NB components (if 0 --> No NB's_)
        inp['ivunit'] = 2  # Units for NB energy  1 --> cm/s  2 --> keV / amu

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

                inp['anneut'] = nf * 1 * 10 ** 13  # inp['nmin_n']
                inp['vneut'] = Tn

                inp2 = namedtuple('inp2', list(inp.keys()))(*list(inp.values()))

                # prepare input file
                f = open('./toadpak', 'w')
                f.write(' &inp')
                f.write('\n' + '  inucz = ' + str(inp2.z_imp))
                f.write('\n' + '  zte = ' + str(0))  # this parameter doesn't matter because imode = 1
                f.write('\n' + '  zne = ' + str(0))  # this parameter doesn't matter because imode = 1
                f.write('\n' + '  laden = ' + str(inp2.laden))
                f.write('\n' + '  ladtip = ' + str(inp2.ladtip))
                f.write('\n' + '  leci = ' + str(inp2.leci))
                f.write('\n' + '  ldrmlt = ' + str(inp2.ldrmlt))
                f.write('\n' + '  ncxb = ' + str(inp2.ncxb))
                f.write('\n' + '  ncxopt = ' + str(inp2.ncxopt))
                f.write('\n' + '  ivunit = ' + str(inp2.ivunit))
                f.write('\n' + '  anneut = ' + str(inp2.anneut))
                f.write('\n' + '  vneut = ' + str(inp2.vneut))
                f.write('\n' + '  imode = ' + str(inp2.imode))
                f.write('\n' + '  nte = ' + str(inp2.nte))
                f.write('\n' + '  nne = ' + str(inp2.nne))
                # f.write('\n' + '  tei = ' + ' '.join(map(str, inp2.tei)))
                # f.write('\n' + '  tei = ' +  np.array2string(inp2.tei, precision=3, separator=' ',suppress_small = False))
                print(np.array2string(inp2.tei, formatter={'float_kind': lambda x: "%.5f" % x})[1:-1])
                f.write(
                    '\n' + '  tei = ' + np.array2string(inp2.tei, formatter={'float_kind': lambda x: "%.5f" % x})[1:-1])
                f.write('\n' + '  anei = ' + ' '.join(map(str, inp2.anei)))
                f.write('\n' + '  nmin = ' + str(inp2.nmin_T))
                f.write('\n' + '  nmax = ' + str(inp2.nmax_T))
                f.write('\n' + ' $end')
                f.write('\n' + '')
                f.write('\n')
                f.close()
                # sys.exit()
                # call adpak
                try:
                    call([inp2.adpak_loc + 'adpak', os.getcwd() + '/toadpak'])
                except AttributeError:
                    try:
                        call(['adpak', os.getcwd() + '/toadpak'])
                    except:
                        raise RuntimeError("Unable to locate adpak.")

                # read adpak output
                with open(os.getcwd() + '/outblk.dat', 'r') as f:
                    # read entire file into one long string variable
                    data = f.read().replace('\n', ' ').replace('. ', ' ').replace(', ', ' ')

                # the following arrays all have the same dimension despite the fact that alinzr
                # is a function only of Te
                # array[ charge state, density value, temperature value ]
                alinzr = np.zeros((inp2.z_imp + 1, inp2.nne, inp2.nte))
                alradr = np.zeros((inp2.z_imp + 1, inp2.nne, inp2.nte))
                alrecr = np.zeros((inp2.z_imp + 1, inp2.nne, inp2.nte))

                # find altei data
                result = re.match(r'.*?data \(+%s.*?\/(.*?)\/.*' % 'altei', data).group(1)
                altei = np.asarray(result.split(), dtype=float)

                # find alnei data
                result = re.match(r'.*?data \(+%s.*?\/(.*?)\/.*' % 'alnei', data).group(1)
                alnei = np.asarray(result.split(), dtype=float)

                # for each charge state plus the ground state,
                for cs in np.arange(inp2.z_imp + 1):
                    # find alinzr data for cs. alinzr is a function of temperature only
                    result = re.match(r'.*?data \(+%s\(\d+ +%s.*?\/(.*?)\/.*' % ('alinzr', cs + 1), data).group(1)
                    alinzr[cs, :, :] = np.asarray(result.split(), dtype=float).reshape(-1, inp2.nte)

                    # find alradr data for cs. alradr is a function of both temperature and electron density
                    result = re.match(r'.*?data \(+%s\(\d+ +%s.*?\/(.*?)\/.*' % ('alradr', cs + 1), data).group(1)
                    alradr[cs, :, :] = np.asarray(result.split(), dtype=float).reshape(-1, inp2.nte)

                    # find alrecr data for cs. alrecr is a function of both temperature and electron density
                    result = re.match(r'.*?data \(+%s\(\d+ +%s.*?\/(.*?)\/.*' % ('alrecr', cs + 1), data).group(1)
                    alrecr[cs, :, :] = np.asarray(result.split(), dtype=float).reshape(-1, inp2.nte)

                # assemble data to pass to Lz and charge state density calculations
                data = namedtuple("data", 'altei alnei alinzr alradr alrecr')(altei, alnei, alinzr, alradr, alrecr)

                # calculate impurity charge state densities. Uncomment the following line if desired.
                # self.calc_nz_cs = self.nz_cs(inp, data)

                # create interpolation function for Lz and dLzdT based on coronal equilibrium approximation
                Lz_interp, dLzdT_interp = self.calc_Lz(inp2, data)

                Lz_complete[i, j, :] = Lz_interp(Te_vals * 1E3 * 1.6021E-19)
                dLzdT_complete[i, j, :] = dLzdT_interp(Te_vals * 1E3 * 1.6021E-19)

        Tn_mesh, nf_mesh, Te_mesh = np.meshgrid(Tn_vals, nf_vals, Te_vals)

        points = np.column_stack((np.log10(Tn_mesh).flatten(),
                                  np.log10(nf_mesh).flatten(),
                                  np.log10(Te_mesh).flatten()))
        Lz = LinearNDInterpolator(points, Lz_complete.flatten(), rescale=False)
        dLzdT = LinearNDInterpolator(points, dLzdT_complete.flatten(), rescale=False)

        # pickle the interpolation object to save time in the future
        filename = os.getcwd() + '/Lz_interpolators/' + imp_names[z] + '_Lz.pkl'
        outfile = open(filename, 'wb')
        pickle.dump(Lz, outfile)
        outfile.close()

        filename = os.getcwd() + '/Lz_interpolators/' + imp_names[z] + '_dLzdT.pkl'
        outfile = open(filename, 'wb')
        pickle.dump(dLzdT, outfile)
        outfile.close()

        # cleanup
        files = ['toadpak', 'outplt.txt', 'outblk.dat', 'adfits.txt', 'outadpk.txt']
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass

        return Lz, dLzdT

    @classmethod
    def calc_nz_cs(cls, inp, data):
        """calculates impurity densities by charge state"""

        tei_lin = np.linspace(inp.nmin, inp.nmax, inp.nte)

        nz_cs = []

        # for each charge state,
        for cs in np.arange(inp.z_imp + 1):
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
            return x * 1E3 * 1.6021E-19

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
                Lz_tot[i] = Lz_tot[i] + 10 ** Lz_cs_interp(T) * frac_abun_interp(T)

        # spline fit the logarithm of the emissivity. We'll spline fit the actual values next
        Lz_tot_interp = UnivariateSpline(inp.tei, np.log10(Lz_tot) + 100.0, s=0)

        # the above function gives the right values, but we also need the derivatives
        # we will now do a spline fit of the values themselves, rather than their base-10 logarithm
        T_kev_vals = np.logspace(-3, 2, 10000)
        T_J_vals = kev_J(T_kev_vals)
        Lz_vals = 10.0 ** (Lz_tot_interp(T_kev_vals) - 100.0)

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

        # for i, dens in enumerate(inp.anei):
        #     for j, temp in enumerate(inp.tei):
        #         for (cs, col), value in np.ndenumerate(M):
        #             if cs == 0:  # do the first row differently
        #                 M[cs, cs] = -(10**data.alinzr[cs, i, j] + 10**data.alrecr[cs, i, j])
        #                 M[cs, cs + 1] = 10 ** data.alrecr[cs + 1, i, j]
        #             elif cs == inp.z_imp:  # do the last row differently
        #                 M[cs, cs - 1] = 10 ** data.alinzr[cs - 1, i, j]
        #                 M[cs, cs] = -(10 ** data.alinzr[cs, i, j] + 10 ** data.alrecr[cs, i, j])
        #             else:
        #                 M[cs, cs - 1] = 10 ** data.alinzr[cs - 1, i, j]
        #                 M[cs, cs] = -(10 ** data.alinzr[cs, i, j] + 10 ** data.alrecr[cs, i, j])
        #                 M[cs, cs + 1] = 10 ** data.alrecr[cs + 1, i, j]

        for i, dens in enumerate(inp.anei):
            for j, temp in enumerate(inp.tei):

                for (cs, col), value in np.ndenumerate(M):
                    if cs == 0:  # do the first row differently
                        M[cs, cs] = -(10 ** data.alinzr[cs, i, j])
                        M[cs, cs + 1] = 10 ** data.alrecr[cs + 1, i, j]
                    elif cs == inp.z_imp:  # do the last row differently
                        M[cs, cs - 1] = 10 ** data.alinzr[cs - 1, i, j]
                        M[cs, cs] = -(10 ** data.alrecr[cs, i, j])
                    else:
                        M[cs, cs - 1] = 10 ** data.alinzr[cs - 1, i, j]
                        M[cs, cs] = -(10 ** data.alinzr[cs, i, j] + 10 ** data.alrecr[cs, i, j])
                        M[cs, cs + 1] = 10 ** data.alrecr[cs + 1, i, j]
                # the following line may help convergence. Uncomment it if necessary. Everything
                # gets normalized anyway.
                # M = M*1E20

                # print M
                # sys.exit()

                result = np.abs(null(M).T)  # solve(M, b)
                # print 'result = ',result
                result = np.where(result < 0, -result, result)
                frac_abun[i, j, :] = result / np.sum(result)

        return frac_abun


if __name__ == '__main__':
    """Example script."""

    # input parameters. Note that some parameters are hard coded in ImpRad.__init__
    z_imp = 18  # Atomic number, Z, of desired element
    laden = 1  # Flag for excitation energy and ion. potential calc. If 0 --> use Mayer formalism. If 1 --> use More formalism
    ladtip = 1  # If 1 --> use tabulated ionization potentials
    leci = 1  # Flag for calculation of ionization rates: 1: XSNQ, 2: Belfast group (for H through O) 3: Younger (Scandium and Fe)
    ldrmlt = 2  # Dielectronic multiplier. 0: Use CDNN, CDNM arrays as given# 1 --> Set CDNN, CDNM equal to 1, 2 --> Use Y. Hahn factors (dut) to set up CDNN, CDNM
    ncxopt = 1  # Selects cross sections to be used   1 --> OSAS   2 --> GJ    3 --> OSCT

    # Helium
    He_2 = ImpRad(z=2)

    # Beryllium
    Be_4 = ImpRad(z=4)

    # Carbon
    C_6 = ImpRad(z=6)

    # Oxygen
    O_8 = ImpRad(z=8)

    # #Neon
    Ne_10 = ImpRad(z=10)

    # #Argon
    Ar_18 = ImpRad(z=18)

    # #Krypton
    Kr_36 = ImpRad(z=36)

    # Xenon
    Xe_54 = ImpRad(z=54)

    # Tungsten
    W_74 = ImpRad(z=74)


