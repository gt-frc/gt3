#!/usr/bin/env python2
# -*- coding: utf-8 -*-

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
        f.write("wall_file       = %s_%s/gt3_diiid_wall.dat \n" % (str(s), str(t)))
        f.write("\n")

        f.write("# CONSTANTS \n")
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

        f.write("#GENERAL GEOMETRY \n")
        f.write("a = %s \n" % str(data.aminor))
        f.write("BT0 = %s \n" % str(data.bphi * -1.))
        f.write("Z0 = 0.0 \n")
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
        if reNeut:
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
        f.write("IP = %s \n" % str(data.plasmaCur))
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

        f.write("sollines_psi_max = 1.07 \n")
        f.write("num_sollines = 6 \n")

        f.write("xi_ib_pts = 10 \n")
        f.write("xi_ob_pts = 10 \n")

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

        f.close()