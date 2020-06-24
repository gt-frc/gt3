#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from Functions.PrepNBIFile import prep_nbi_infile
from Processors.ReadNBIOutfile import read_nbi_outfile
from CalcNBIValues import calc_nbi_vals
from Test.DebugCore import debugCore
from Test.DebugInp import debugInp
import numpy as np
from scipy.interpolate import UnivariateSpline
from subprocess import Popen, PIPE
from collections import namedtuple
from scipy.constants import physical_constants
import os, sys
from Functions.Beam import Beam

m_d = physical_constants['deuteron mass'][0]
z_d = 1

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
            # Call signature: def calcHofRho(rho, rho2vol, dVdrho, rtang, shaf_shift, kappa_vals, beamHeight, beamWidth, Te, TC, R0, ne, beamE, Zeff, a, gaussR, gaussZ)
            inp.beamWidth = 0.11
            inp.beamHeight = 0.28
            inp.beamE = 77
            inp.beamA = 2.
            inp.beamP = 0.9
            inp.bgaussZ = .188
            inp.bgaussR = .066
            inp.rtang = 1.2


            nbeams_space = np.linspace(0., 1.0, 50)
            self.Beam1 = Beam(nbeams_space, core.r2vol, core.dVdr, inp.rtang, core.shaf_shift, core.kappa_vals,
                                     inp.beamHeight,
                                     inp.beamWidth,
                                     UnivariateSpline(core.rho[:, 0], core.T_fsa.e.kev)(nbeams_space),
                                     UnivariateSpline(core.rho[:, 0], core.T_fsa.i.kev)(nbeams_space),
                                     core.R0_a,
                                     UnivariateSpline(core.rho[:, 0], core.n_fsa.e)(nbeams_space),
                                     inp.beamE,
                                     inp.beamA,
                                     inp.beamP,
                                     UnivariateSpline(core.rho[:, 0], core.z_eff_fsa)(nbeams_space),
                                     core.a,
                                     inp.bgaussR,
                                     inp.bgaussZ,
                                     verbose=True,
                                     timing=False)

            self.Beam2 = Beam(nbeams_space, core.r2vol, core.dVdr, inp.rtang, core.shaf_shift, core.kappa_vals,
                                     inp.beamHeight,
                                     inp.beamWidth,
                                     UnivariateSpline(core.rho[:, 0], core.T_fsa.e.kev)(nbeams_space),
                                     UnivariateSpline(core.rho[:, 0], core.T_fsa.i.kev)(nbeams_space),
                                     core.R0_a,
                                     UnivariateSpline(core.rho[:, 0], core.n_fsa.e)(nbeams_space),
                                     inp.beamE,
                                     inp.beamA,
                                     inp.beamP,
                                     UnivariateSpline(core.rho[:, 0], core.z_eff_fsa)(nbeams_space),
                                     core.a,
                                     inp.bgaussR,
                                     inp.bgaussZ,
                                     verbose=True,
                                     timing=False)


        try:
            if nbi_vals: self.debug = nbi_vals.debug
        except:
            self.debug = nbi_vals_list

        try:
            nbi_vals_list
            # create multi-beam nbeams object
            self.beams = namedtuple('beam', 'D1 D2 D3')(
                # create D1 beam
                namedtuple('beam_D1', 'z m E P rtang dPdV dPdr zeta coI')(
                    z_d,
                    m_d,
                    namedtuple('E', 'kev ev J')(
                        [nbi_vals_list[i].beam_en_1 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_en_1 * 1E3 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_en_1 * 1E3 * 1.6021E-19 for i in range(len(nbi_vals_list))]
                    ),
                    namedtuple('P', 'MW W')(
                        [nbi_vals_list[i].beam_pwr_1 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_pwr_1 * 1E6 for i in range(len(nbi_vals_list))]
                    ),
                    [data[i]['rtang'] for i in range(len(nbi_vals_list))],
                    namedtuple('dPdV', 'v1D v2D')(
                        namedtuple('v1D', 'MW W')(
                            [nbi_vals_list[i].dPdV_1_1D for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdV_1_1D * 1E6 for i in range(len(nbi_vals_list))]
                        ),
                        namedtuple('v2D', 'MW W')(
                            [nbi_vals_list[i].dPdV_1 for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdV_1 * 1E6 for i in range(len(nbi_vals_list))]
                        )
                    ),
                    namedtuple('dPdr', 'v1D v2D')(
                        namedtuple('v1D', 'MW W')(
                            [nbi_vals_list[i].dPdr_1_1D for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdr_1_1D * 1E6 for i in range(len(nbi_vals_list))]
                        ),
                        namedtuple('v2D', 'MW W')(
                            [nbi_vals_list[i].dPdr_1 for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdr_1 * 1E6 for i in range(len(nbi_vals_list))]
                        )
                    ),
                    [nbi_vals_list[i].zeta_1 for i in range(len(nbi_vals_list))],
                    [data[i]['coI'] for i in range(len(data))]
                ),

                # create D2 beam
                namedtuple('beam_D2', 'z m E P rtang dPdV dPdr zeta coI')(
                    z_d,
                    m_d,
                    namedtuple('E', 'kev ev J')(
                        [nbi_vals_list[i].beam_en_2 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_en_2 * 1E3 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_en_2 * 1E3 * 1.6021E-19 for i in range(len(nbi_vals_list))]
                    ),
                    namedtuple('P', 'MW W')(
                        [nbi_vals_list[i].beam_pwr_2 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_pwr_2 * 1E6 for i in range(len(nbi_vals_list))]
                    ),
                    [data[i]['rtang'] for i in range(len(nbi_vals_list))],
                    namedtuple('dPdV', 'v1D v2D')(
                        namedtuple('v1D', 'MW W')(
                            [nbi_vals_list[i].dPdV_2_1D for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdV_2_1D * 1E6 for i in range(len(nbi_vals_list))]
                        ),
                        namedtuple('v2D', 'MW W')(
                            [nbi_vals_list[i].dPdV_2 for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdV_2 * 1E6 for i in range(len(nbi_vals_list))]
                        )
                    ),
                    namedtuple('dPdr', 'v1D v2D')(
                        namedtuple('v1D', 'MW W')(
                            [nbi_vals_list[i].dPdr_2_1D for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdr_2_1D * 1E6 for i in range(len(nbi_vals_list))]
                        ),
                        namedtuple('v2D', 'MW W')(
                            [nbi_vals_list[i].dPdr_2 for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdr_2 * 1E6 for i in range(len(nbi_vals_list))]
                        )
                    ),
                    [nbi_vals_list[i].zeta_2 for i in range(len(nbi_vals_list))],
                    [data[i]['coI'] for i in range(len(data))]
                ),

                # create D3 beam
                namedtuple('beam_D3', 'z m E P rtang dPdV dPdr zeta coI')(
                    z_d,
                    m_d,
                    namedtuple('E', 'kev ev J')(
                        [nbi_vals_list[i].beam_en_3 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_en_3 * 1E3 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_en_3 * 1E3 * 1.6021E-19 for i in range(len(nbi_vals_list))]
                    ),
                    namedtuple('P', 'MW W')(
                        [nbi_vals_list[i].beam_pwr_3 for i in range(len(nbi_vals_list))],
                        [nbi_vals_list[i].beam_pwr_3 * 1E6 for i in range(len(nbi_vals_list))]
                    ),
                    [data[i]['rtang'] for i in range(len(data))],
                    namedtuple('dPdV', 'v1D v2D')(
                        namedtuple('v1D', 'MW W')(
                            [nbi_vals_list[i].dPdV_3_1D for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdV_3_1D * 1E6 for i in range(len(nbi_vals_list))]
                        ),
                        namedtuple('v2D', 'MW W')(
                            [nbi_vals_list[i].dPdV_3 for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdV_3 * 1E6 for i in range(len(nbi_vals_list))]
                        )
                    ),
                    namedtuple('dPdr', 'v1D v2D')(
                        namedtuple('v1D', 'MW W')(
                            [nbi_vals_list[i].dPdr_3_1D for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdr_3_1D * 1E6 for i in range(len(nbi_vals_list))]
                        ),
                        namedtuple('v2D', 'MW W')(
                            [nbi_vals_list[i].dPdr_3 for i in range(len(nbi_vals_list))],
                            [nbi_vals_list[i].dPdr_3 * 1E6 for i in range(len(nbi_vals_list))]
                        )
                    ),
                    [nbi_vals_list[i].zeta_3 for i in range(len(nbi_vals_list))],
                    [data[i]['coI'] for i in range(len(data))]
                )
            )
        except:
            # create beams object
            self.beams = namedtuple('beam', 'D1 D2 D3')(
                # create D1 beam
                namedtuple('beam_D1', 'z m E P rtang dPdV dPdr zeta')(
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
                namedtuple('beam_D2', 'z m E P rtang dPdV dPdr zeta')(
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
                namedtuple('beam_D3', 'z m E P rtang dPdV dPdr zeta')(
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

    def print_something(self):
        self.beams.D1.E

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

    ax43 = fig4.add_subplot(233)
    ax43.set_title(r'Power density profile', fontsize=16)
    ax43.set_ylabel(r'$\frac{dP}{dr}$', fontsize=16)
    ax43.set_xlabel(r'rho', fontsize=16)
    ax43.plot(beamData.debug['rho_nbeams'], beamData.debug['dvol'])

    plt.show()