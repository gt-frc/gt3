#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
"""
import sys, os
from neutpy_prep import Neutrals
from iol import IOL
from sol import Sol
from read_infile import ReadInfile
from imp_rad import ImpRad
from core import Core
from beamdep import BeamDeposition
from dens_lim import DensityLimit
from marfe import Marfe
from radial_transport import RadialTransport
from collections import namedtuple

class gt3:
    """GT3 calculates various tokamak-related quantities

    Methods:
        allthethings
        justneutrals
        plotstuff

    Attributes:

    External Dependencies:
        Triangle:       Used to create the triangular mesh for neutpy. Source
                        code and documentation can be found at
                        https://www.cs.cmu.edu/~quake/triangle.html. Can be
                        installed from Ubuntu repositories as well.
        Neutpy:
        nbeams:
        adpack:
    """
    def __init__(self, shotlabel=None, mode=None, iolFlag=True, neutFlag = True, debugRT= False):
        sys.dont_write_bytecode = True 
        # Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel
        self.inp = ReadInfile(self.shotlabel)
        self.core = Core(self.inp)
        self.iolFlag = iolFlag
        self.neutFlag = neutFlag
        self.debugRT = debugRT

        if mode == 'coreonly':
            pass

        if mode == 'coreandsol':
            self.sol = Sol(self.inp, self.core)
        elif mode == 'thermaliol':
            self.iol = IOL(self.inp, self.core)
        elif mode == 'fulliol':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
        elif mode == 'imp':
            self.imp = ImpRad(core=self.core)
        elif mode == 'ntrls':
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'ntrlsandiol':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'nbi':
            self.nbi = BeamDeposition(self.inp, self.core)
        elif mode == 'marfe_denlim':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.core, self.nbi)
            self.mar = Marfe(core=self.core)
        elif mode == 'marfe':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.mar = Marfe(core=self.core)
        elif mode == 'allthethings':
            self.nbi = BeamDeposition(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(core=self.core)
            self.dl = DensityLimit(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = Marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            self.sol = Sol(self.inp, self.core)
            self.iol = IOL(self.inp, self.core)
            self.nbi = BeamDeposition(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(z=None, core=self.core)
            self.rtrans = RadialTransport(self.inp, self.core, self.iol, self.nbi,self.iolFlag, self.neutFlag, debugFlag = self.debugRT)

def getNum(prompt, valType):
    """

    Function to obtain valid inputs from keyboard

    :param prompt: Text prompt to user
    :param valType: Type of value expected from user
    :return: Your mom
    """

    if valType=='f':
        while True:
            try:
                return float(input(prompt))
            except:
                print("Invalid input \n")
    elif valType=='i':
        while True:
            try:
                x = input(prompt)
                if isinstance(x, (int, long)):
                    return(int(x))
                else:
                    print "Invalid input \n"
                    continue
            except:
                print "Invalid input \n"

def getVals(s, t, f):
    """

    This function interactively obtains the 0D parameters for your input file to be created. The commented-out section
    is data for a shot just for debugging purposes

    :param s: shot id
    :param t: time id
    :param f: file to be generated
    :return: Your mom
    """

    data={}
    """ Function that gets 0D plasma values from user input """

    print """Enter 0D plasma parameters for shot %s.%s
    
             File will be generated at inputs/%s""" % (str(s), str(t), str(f))

    # data={'eq1' : 1.6E-19,
    #       'eq2' : 9.6E-19,
    #       'xmas1' : 3.35E-27,
    #       'xmas2' : 2.01E-26,
    #       'ephia' : .04,
    #       'xk' : 1.6E-19,
    #       'delma' : 0.005,
    #       'aminor' : 0.598,
    #       'bphi' : -2.04,
    #       'rmajor' : 1.7,
    #       'kappaup': 1.82,
    #       'kappalo' : 1.82,
    #       'triup': .237,
    #       'trilo' : 0.,
    #       'thetapts' : 30,
    #       'rhopts_edge' : 100,
    #       'rhopts_core' : 10,
    #       'xptR' : 1.48,
    #       'xptZ' : -1.24,
    #       'jknot' : 2850000,
    #       'plasmaCur' : 1.38,
    #       'q95' : 3.57,
    #       'ebeam' : 77.48,
    #       'pbeam' : 4.6,
    #       'rtang' : 1.09,
    #       'bknot' : 2.0,
    #       'pwrfrac1' : .76,
    #       'pwrfrac2' : .13,
    #       'pwrfrac3' : .11,
    #       'epsknot' : 1.265,
    #       'epssep' : 1.82,
    #       'shftknot' : 0.033
    #       }

    data['eq1'] = getNum("Enter charge of main ion: ", "f")
    data['eq2'] = getNum("Enter charge of main impurity species: ", "f")
    data['xmas1'] = getNum("Enter main ion mass: ", "f")
    data['xmas2'] = getNum("Enter main impurity mass: ", "f")
    data['ephia'] = getNum("Enter the toroidal electric field: ", "f")
    data['xk'] = 1.6E-19
    data['delma'] = .005
    data['aminor'] = getNum("Enter the plasma radius (a minor): ","f")
    data['bphi'] = getNum("Enter the toroidal plasma field strength in T (abs mag): ", "f")
    data['rmajor'] = getNum("Enter R major: ", "f")
    data['kappaup'] = getNum("Enter upper elongation (kappa): ", "f")
    data['kappalo'] = getNum("Enter lower elongation (kappa): ", "f")
    data['triup'] = getNum("Enter upper triangularity (delta): ", "f")
    data['trilo'] = getNum("Enter lower triangularity (delta): ", "f")
    data['thetapts'] = getNum("Enter approximate number of theta points (typically 30): ", 'i')
    data['rhopts_edge'] = getNum("Enter rho points in the edge (typically 100): ", 'i')
    data['rhopts_core'] = getNum("Enter rho points in the core (typically 10): ", 'i')
    data['xptR'] = getNum("Enter the X-point R coordinate: ", 'f')
    data['xptZ'] = getNum("Enter the X-point Z coordinate: ", 'f')
    data['jknot'] = getNum("Enter the r=0 plasma current density (in A/m^3): ", 'f')
    data['plasmaCur'] = getNum("Enter the plasma current (in MA): ", 'f')
    data['q95'] = getNum("Enter q95: ", 'f')
    data['ebeam'] = getNum("Enter beam ion energy in eV: ", 'f')
    data['pbeam'] = getNum("Enter beam power in MWYup: ", 'f')
    data['rtang'] = getNum("Enter radius of tangency in cm: ", 'f')
    data['bknot'] = abs(data['bphi'])
    data['pwrfrac1'] = getNum("Enter fraction of beam power to D1: ", 'f')
    data['pwrfrac2'] = getNum("Enter fraction of beam power to D2: ", 'f')
    data['pwrfrac3'] = getNum("Enter fraction of beam power to D3: ", 'f')
    data['epsknot'] = getNum("Enter epsilon at the plasma center: ", 'f')
    data['epssep'] = getNum("Enter epsilon at separatrix: ", 'f')
    data['shftknot'] = getNum("Enter the Shavranof shift at the plasma center: ", 'f')

    return namedtuple('data', sorted(data))(**data)

def writeFile(s, t, fpath, data, reNeut=False):

    """
    This function writes an input file and saves it in the inputs/shot_timeid folder.

    :param s: shot id
    :param t: time id
    :param fpath: This is the file that is to be written. See the folder structure in getVals() for file structure
    :param data: These data are from the interactive getVals() code that lets you interactively give parameters for
                 setting up input files
    :param reNeut: This flag tells GT3 whether a new run of neutpy needs to be performed. Concerns ntrl_switch in gt3
    :return: Your mom
    """
    with open(fpath, "w") as f:
        f.write("d3d_iter = 1 \n")
        f.write("exp_inp = 1 \n")
        f.write("\n")
        f.write("#1D PROFILE INPUT FILES \n")
        f.write("ne_file         = %s_%s/gt3_%s_%s_ne.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("nD_file         = %s_%s/gt3_%s_%s_nD.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("Te_file         = %s_%s/gt3_%s_%s_Te.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("Ti_file         = %s_%s/gt3_%s_%s_Ti.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("er_file         = %s_%s/gt3_%s_%s_er.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("fz1_file        = %s_%s/gt3_%s_%s_fz1.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("fracz_file      = %s_%s/gt3_%s_%s_fracz.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("exlti_file      = %s_%s/gt3_%s_%s_exlti.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("exlte_file      = %s_%s/gt3_%s_%s_exlte.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("exlni_file      = %s_%s/gt3_%s_%s_exlni.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("vpolC_file      = %s_%s/gt3_%s_%s_vpolC.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("vtorC_file      = %s_%s/gt3_%s_%s_vtorC.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("vpolD_file      = %s_%s/gt3_%s_%s_vpolD.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("vtorD_file      = %s_%s/gt3_%s_%s_vtorD.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("q_file          = %s_%s/gt3_%s_%s_q.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("zbar2_file      = %s_%s/gt3_%s_%s_zbar2.dat    \n" % (str(s), str(t), str(s), str(t)))
        f.write("\n")
        f.write("2D QUANTITIES INPUT FILES \n")
        f.write("bpol_file       = %s_%s/gt3_%s_%s_bpol.dat  \n" % (str(s), str(t), str(s), str(t)))
        f.write("btor_file       = %s_%s/gt3_%s_%s_btor.dat  \n" % (str(s), str(t), str(s), str(t)))
        f.write("psirz_file      = %s_%s/gt3_%s_%s_psirz.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("\n")
        f.write("2D LINE INPUT FILES \n")
        f.write("sep_file        = %s_%s/gt3_%s_%s_sep.dat \n" % (str(s), str(t), str(s), str(t)))
        f.write("wall_file       = %s_%s/gt3_diiid_wall.dat \n" % (str(s), str(t)))
        f.write("\n")

        f.write("# CONSTANTS \n")
        f.write("eq1 = %s \n" % str(data.eq1))
        f.write("eq2 = %s \n" % str(data.eq2))
        f.write("xmas1 = %s \n" % str(data.xmas1))
        f.write("xmas2 = %s \n" % str(data.xmas2))
        f.write("ephia = %s \n" % str(data.ephia))
        f.write("xk = %s \n" % str(data.xk))
        f.write("delma = %s \n" % str(data.delma))
        f.write("xnuioni = 0.0 \n")
        f.write("xnuati = 0.0 \n")
        f.write("\n")

        f.write("#NEUTRAL BEAM DEPOSITION \n")
        f.write(
            "nbeams_loc = /home/jonathan/Dropbox/GTEDGE/MyPrograms/GTEDGE/lib/beams/NBeamsMDS/NBeams/bin/Release/nbeams \n")
        f.write("ebeam = %s \n" % str(data.ebeam))
        f.write("abeam = 2 \n")
        f.write("alphain = .6475 \n")
        f.write("pbeam = %s \n" % str(data.pbeam))
        f.write("rtang = %s \n" % str(data.rtang))
        f.write("bknot = %s \n" % str(data.bknot))
        f.write("pwrfrac1 = %s \n" % str(data.pwrfrac1))
        f.write("pwrfrac2 = %s \n" % str(data.pwrfrac2))
        f.write("pwrfrac3 = %s \n" % str(data.pwrfrac3))
        f.write("epsknot = %s \n" % str(data.epsknot))
        f.write("eps_sep = %s \n" % str(data.epssep))
        f.write("shftknot = %s \n" % str(data.shftknot))
        f.write("\n")

        f.write("#IMPURITY CALCULATION \n")
        f.write("adpak_loc = /home/jonathan/Dropbox/GTEDGE/MyPrograms/GTEDGE/MaxPlasma/adpak/ \n")
        f.write("# ADPAK PARAMETERS \n")
        f.write("z_imp = 6 \n")
        f.write("laden = 1 \n")
        f.write("ladtip = 1 \n")
        f.write("leci = 1 \n")
        f.write("ldrmlt = 2 \n")
        f.write("ncxopt = 1 \n")
        f.write("anneut = 1.0E9 \n")
        f.write("vneut = 0.001 \n")
        f.write("\n")

        f.write("#GENERAL GEOMETRY \n")
        f.write("a = %s \n" % str(data.aminor))
        f.write("BT0 = %s \n" % str(data.bphi * -1.))
        f.write("R0_a = %s \n" % str(data.rmajor))
        f.write("Z0 = 0.0 \n")
        f.write("kappa_up = %s \n" % str(data.kappaup))
        f.write("kappa_lo = %s \n" % str(data.kappalo))
        f.write("s_k_up = 0.0 \n")
        f.write("s_k_lo = 0.0 \n")
        f.write("tri_up = %s \n" % str(data.triup))
        f.write("tri_lo = %s \n" % str(data.trilo))
        f.write("thetapts_approx = %s \n" % str(int(data.thetapts)))
        f.write("rhopts = %s \n" % str(201))
        f.write("edge_rho = 0.8 \n")
        try:
            f.write("rhopts_edge = %s \n" % str(data.rhopts_edge))
        except:
            f.write("rhopts_edge = %s \n" % str(100))
        try:
            f.write("rhopts_core = %s \n" % str(data.rhopts_core))
        except:
            f.write("rhopts_core = %s \n" % str(10))
        try:
            f.write("thetapts_approx = %s \n" % str(data.thetapts))
        except:
            f.write("thetapts_approx = %s \n" % str(30))
        f.write("\n")
        f.write("#NEUTRALS CALCULATION \n")
        if reNeut == True:
            f.write("ntrl_switch = 2 \n")
        elif reNeut == False:
            f.write("ntrl_switch = 1 \n")
        else:
            raise Exception("reNeut not defined")
        f.write("edge_rho_ntrl = 0.8 \n")
        f.write("neut_outfile =  inputs/%s_%s/gt3_%s_%s_neut.dat \n" % (str(s), str(t), str(s), str(t)))
        #        f.write("neut_outfile = gt3_%s_%s_neut.dat \n" % (str(s),str(t)))
        f.write("rhopts_edge_ntrl = %s \n" % str(10))
        f.write("ntrl_thetapts = 33 \n")
        f.write("\n")

        f.write("#XMILLER PARAMETERS \n")
        f.write("xmil = 1 \n")
        f.write("xpt_R = %s \n" % str(data.xptR))
        f.write("xpt_Z = %s \n" % str(data.xptZ))
        f.write("\n")

        f.write("#BACKGROUND DENSITIES AND TEMPERATURES (IF NOT READING FROM INPUT FILE) \n")
        f.write("#ni0 = 3.629E19 \n")
        f.write("#ni9 = 1.523E19\n")
        f.write("#ni_sep = 0.3E19\n")
        f.write("#ni_dp = 1E17\n")
        f.write("#nu_ni = 3.0\n")
        f.write("#ne0 = 3.629E19\n")
        f.write("#ne9 = 1.523E19\n")
        f.write("#ne_sep = 0.3E19\n")
        f.write("#ne_dp = 1E17\n")
        f.write("#nu_ne = 2.5\n")
        f.write("#Ti0 = 35\n")
        f.write("#Ti9 = 6\n")
        f.write("#Ti_sep = 0.6\n")
        f.write("#Ti_dp = 0.03\n")
        f.write("#nu_Ti = 3.5\n")
        f.write("#Te0 = 36\n")
        f.write("#Te9 = 6\n")
        f.write("#Te_sep = 0.6\n")
        f.write("#Te_dp = 0.01\n")
        f.write("#nu_Te = 3.5 \n")
        f.write("\n")

        f.write("#CURRENT-RELATED PARAMETERS \n")
        f.write("j0 = %s \n" % str(data.jknot))
        f.write("j_sep = 0 \n")
        f.write("nu_j = 1. \n")
        f.write("IP = %s \n" % str(data.plasmaCur))
        f.write("q95 = %s \n" % str(data.q95))
        f.write("\n")
        f.write("rmeshnum_p = 5 \n")
        f.write("xtheta1 = 3.6 \n")
        f.write("xtheta2 = 2.7 \n")
        f.write("xtheta3 = 0.1 \n")
        f.write("xtheta4 = -1.0 \n")
        f.write("\n")
        f.write("#ION ORBIT LOSS CALCULATION \n")
        f.write("numcos = 8 \n")
        f.write("R_loss = 0.5 \n")
        f.write("\n")

        f.write("pfr_ni_val = 1.0E14 \n")
        f.write("pfr_ne_val = 1.0E14 \n")
        f.write("pfr_Ti_val = 0.002  \n")
        f.write("pfr_Te_val = 0.002  \n")

        f.write("############################################################################### \n")
        f.write("# CONFIG - YOU CAN PROBABLY LEAVE THESE ALONE IF YOU DON'T KNOW WHAT THEY ARE \n")
        f.write("############################################################################### \n")
        f.write("verbose = 1 \n")

        f.write("corelines_begin = 0.75 \n")
        f.write("num_corelines = 10 \n")

        f.write("sollines_psi_max = 1.07 \n")
        f.write("num_sollines = 6 \n")

        f.write("xi_sep_pts = 50 \n")
        f.write("ib_trim_off = 0.1 \n")
        f.write("ob_trim_off = 0.1 \n")

        f.write("xi_ib_pts = 10 \n")
        f.write("xi_ob_pts = 10 \n")

        f.write("core_pol_pts = 30 \n")
        f.write("core_thetapts_ntrl = 50 \n")

        f.write("#rhopts_ntrl = 100 \n")
        f.write("edge_rho_ntrl = 0.8 \n")
        f.write("rhopts_edge_ntrl = 5 \n")
        f.write("rhopts_core_ntrl = 10 \n")
        f.write("ib_div_pol_pts = 7 \n")
        f.write("ob_div_pol_pts = 7 \n")
        f.write("wall_ni_min = 1.0E15 \n")
        f.write("wall_ne_min = 1.0E15 \n")
        f.write("wall_Ti_min = %s \n" % str(0.02 * 1.0E3 * 1.6021E-19))
        f.write("wall_Te_min = %s \n" % str(0.02 * 1.0E3 * 1.6021E-19))
        f.write("core_thetapts_ntrl = 30 \n")
        f.write("ib_thetapts_ntrl = 10 \n")
        f.write("ob_thetapts_ntrl = 10 \n")
        f.write("tri_min_angle = 20 \n")
        f.write("tri_min_area = 0.01 \n")

        f.close()


def gt3Prep(s, t, r, m, quiet = False, genFiles=False):

    """
    This function allows you to prep GT3 in an interactive manner.

    Your file structure should look something like this if you want to use this template, or you could create your own

├── inputs
│   ├── 118888_1525
│   ├── 118888_1570
│   │   ├── gt3_118888_1570_bpol.dat
│   │   ├── gt3_118888_1570_btor.dat
│   │   ├── gt3_118888_1570_er.dat
│   │   ├── gt3_118888_1570_exlni.dat
│   │   ├── gt3_118888_1570_exlte.dat
│   │   ├── gt3_118888_1570_exlti.dat
│   │   ├── gt3_118888_1570_fracz.dat
│   │   ├── gt3_118888_1570_fz1.dat
│   │   ├── gt3_118888_1570_nD.dat
│   │   ├── gt3_118888_1570_ne.dat
│   │   ├── gt3_118888_1570_neut.dat
│   │   ├── gt3_118888_1570_psirz.dat
│   │   ├── gt3_118888_1570_q.dat
│   │   ├── gt3_118888_1570_R.dat
│   │   ├── gt3_118888_1570_Te.dat
│   │   ├── gt3_118888_1570_Ti.dat
│   │   ├── gt3_118888_1570_vpolC.dat
│   │   ├── gt3_118888_1570_vpolD.dat
│   │   ├── gt3_118888_1570_vtorC.dat
│   │   ├── gt3_118888_1570_vtorD.dat
│   │   ├── gt3_118888_1570_zbar2.dat
│   │   ├── gt3_118888_1570_Z.dat
│   │   └── gt3_diiid_wall.dat
│   ├── 144977_3000
│   ├── 164436_3720
│   ├── 164988_1915
│   ├── 175826_2010
│   ├── 92980_3000
│   ├── 92980_3600
│   ├── togt3_d3d_118888_1570
│   ├── togt3_d3d_118890_1560
│   ├── togt3_d3d_164436_3720



    :param s: shot id
    :param t: time id
    :param r: run id (unsure if necessary)
    :param m: filename for input file to be created
    :param quiet: If True, does not provide interactive requests. Will use input file that it gt3 finds. If no data file,
                  system will halt
    :param genFiles: Flag to tell code in the future to transform DIII-D data into format we want for gt3
    :return: Your mom
    """

    varDict = {}

    #    direct="MaxPlasma\\inputs\%s\%s" % (str(s),str(t))
    direct = os.path.join("inputs", "%s_%s") % (str(s), str(t))
    if not os.path.exists(direct):
        print """Directory for shot %s.%s does not exist at inputs/%s_%s
                 Directory used: %s
                 Shutting down now.""" % (str(s),str(t),str(s),str(t), str(os.getcwd()+direct))
        raise Exception
    fileName = m
    fpath = os.path.join("inputs", fileName)
    if not os.path.exists(os.path.join("inputs", fileName)):
        if quiet:
            print """Silent mode active: No input file found.
            
                     Shutting down now."""
            raise Exception
        data=getVals(s, t, fileName)
        writeFile(s, t, fpath, data, reNeut=False)
    else:
        print """Input file for shot %s.%s exists""" % (str(s), str(t))
        while True:
            if quiet:
                print """Silent mode active: Using current input file for shot %s.%s""" % (str(s),str(t))
                break
            ch = raw_input("Use current input file? (Y/N) ")
            if (str(ch) == 'Y' or str(ch) == 'y'): break
            elif (str(ch) == 'N' or str(ch) == 'n'):
                os.remove(os.path.join("inputs", fileName))
                data=getVals(s, t, fileName)
                writeFile(s, t, fpath, data, reNeut=False)
                break
            else:
                print "Invalid selection \n"




        #####################################################################################################
        #
        #   This section will probably be usable for when we do not have data ready for gt3
        #
        #   TODO: Implement this somehow
        #
        #####################################################################################################

    if genFiles:
        fileList = {
            'data.xne': ('gt3_%s_%s_ne.dat' % (str(s), str(t)), data.xne),
            'data.xni': ('gt3_%s_%s_ni.dat' % (str(s), str(t)), data.xni),
            'data.xte': ('gt3_%s_%s_Te.dat' % (str(s), str(t)), data.xte),
            'data.xti': ('gt3_%s_%s_Ti.dat' % (str(s), str(t)), data.xti),
            'data.xer': ('gt3_%s_%s_er.dat' % (str(s), str(t)), data.xer),
            'data.fz1': ('gt3_%s_%s_fz1.dat' % (str(s), str(t)), data.fz1),
            'data.fracz': ('gt3_%s_%s_fracz.dat' % (str(s), str(t)), data.fracz),
            'data.exlti': ('gt3_%s_%s_exlti.dat' % (str(s), str(t)), data.exlti),
            'data.exlte': ('gt3_%s_%s_exlte.dat' % (str(s), str(t)), data.exlte),
            'data.exlni': ('gt3_%s_%s_exlni.dat' % (str(s), str(t)), data.exlni),
            'data.vpolC': ('gt3_%s_%s_vpolC.dat' % (str(s), str(t)), data.vpolC),
            'data.vtorC': ('gt3_%s_%s_vtorC.dat' % (str(s), str(t)), data.vtorC),
            'data.q': ('gt3_%s_%s_q.dat' % (str(s), str(t)), data.q),
            'data.zbar2': ('gt3_%s_%s_zbar2.dat   ' % (str(s), str(t)), data.zbar2)}
        try:
            os.chdir("inputs")
        except:
            raise Exception("gt3 input directory not found")

        #####################################################################################################
        #
        #   This section will probably be usable for when we do not have data ready for gt3
        #
        #   TODO: Implement this somehow
        #
        #####################################################################################################

        unitCon = {'data.xti': 1E-3,
                   'data.xte': 1E-3,
                   'data.xer': 1E-3}
        for var in fileList.keys():
            if var in unitCon.keys():
                with open(fileList[var][0], "w") as openFile:
                    for a, b in zip(data.rhor, fileList[var][1]):
                        openFile.write("%s    %s \n" % (str(a), str(b * unitCon[var])))
            else:
                with open(fileList[var][0], "w") as openFile:
                    for a, b in zip(data.rhor, fileList[var][1]):
                        openFile.write("%s    %s \n" % (str(a), str(b)))

        if hasattr(data, 'vpolD'):
            with open('gt3_%s_%s_vpolD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, data.vpolD):
                    openFile.write("%s   %s \n" % (a, b))
        else:
            with open('gt3_%s_%s_vpolD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, [0.] * len(data.rhor)):
                    openFile.write("%s   %s \n" % (a, b))

        if hasattr(data, 'vtorD'):
            with open('gt3_%s_%s_vtorD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, data.vtorD):
                    openFile.write("%s   %s \n" % (a, b))
        else:
            with open('gt3_%s_%s_vtorD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, [0.] * len(data.rhor)):
                    openFile.write("%s   %s \n" % (a, b))
