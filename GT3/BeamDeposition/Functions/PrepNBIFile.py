#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from GT3.BeamDeposition.Functions.CalcPwrFrac import calc_pwr_frac
from scipy.interpolate import interp1d
import numpy as np
import os

def prep_nbi_infile(inp, core, index = False, indBeam = False):
    i = index
    if indBeam: pwr_frac = calc_pwr_frac(indBeam['ebeam'])
    else: pwr_frac = calc_pwr_frac(inp.ebeam)

    # pwr_frac = [0.7, 0.2, 0.1]

    # f1=open(inp.nbeams_loc+"inbeams.dat", "w")
    if __name__=="__main__":
        f = open("inputs/inbeams_test.dat", "w")
    else:
        if indBeam: f=open(os.path.dirname(__file__) + "/../inbeams_%s.dat" % str(i), "w")
        else: f = open(os.path.dirname(__file__) + "/../inbeams.dat", "w")
    f.write("$nbin\n")
    f.write("nbeams = 1\n")
    f.write("inbfus = 1\n")
    f.write("amb = 2.0\n")
    f.write("zbeam = 1.0\n")

    if indBeam: f.write("ebeam = " + str(indBeam["ebeam"]) + "\n")
    else: f.write("ebeam = " + str(inp.ebeam) + "\n")

    if indBeam: f.write("pbeam = " + str(indBeam["pbeam"]) + "\n")
    else: f.write("pbeam = " + str(inp.pbeam) + "\n")

    if indBeam: f.write("pbeam = " + str(indBeam["rtang"]) + "\n")
    else: f.write("rtang = " + str(inp.rtang) + "\n")

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