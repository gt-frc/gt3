#!/usr/bin/python2.7

import numpy as np
import os
import texttable
import sys
from shapely.geometry import LineString
import matplotlib.pyplot as plt


def getGT3Test(iolFlag=True, verbose=False, neutFlag=True, mode="coreonly"):
    try:
        from GT3 import gt3
    except ImportError:
        ImportError("Unable to load GT3 module")
    return gt3(preparedInput=TestClass(), iolFlag=iolFlag, verbose=verbose, neutFlag=neutFlag, mode=mode)

class TestClass:

    def __init__(self, shotid=164436, timeid=3740):
        self.shotid = shotid
        self.timeid = timeid
        print("This is the GT3 Test Base Class")
        print("This class uses DIII-D shot %s.%s for debugging purposes" % (self.shotid, self.timeid))
        print("Use the print_summary method to see a summary of this shot")

        """
        Generate various 0D parameters for the plasma from shot 164436.3740
        """
        self.iolFlag = True
        self.BT0 = 2.004
        self.wall_ni_min = 1.0E15
        self.wall_ne_min = 1.0E15
        self.wall_Ti_min = 0.005
        self.wall_Te_min = 0.005
        self.pfr_ni_val = 1.0E14
        self.pfr_ne_val = 1.0E14
        self.pfr_Ti_val = 0.002
        self.pfr_Te_val = 0.002
        self.ebeam = 64.72
        self.abeam = 2
        self.alphain = 0.6475
        self.pbeam = 0.8
        self.rtang = 1.146
        self.R_loss = .5
        self.sep_val = 1.0


        """
        Define the meshing parameters
        """
        self.thetapts_approx = 50
        self.edge_rho = 0.7
        self.rhopts_edge = 100
        self.rhopts_core = 10
        self.core_thetapts_ntrl = 50
        self.edge_rho_ntrl = 0.8
        self.rhopts_edge_ntrl = 5
        self.rhopts_core_ntrl = 10
        self.sollines_psi_max = 1.07
        self.num_sollines = 3
        self.xi_ib_pts = 10
        self.xi_ob_pts = 10
        self.numcos = 20
        self.sep_val = 1.0

        """
        Define the meshing parameters
        """

        self.er_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_er.dat' % (self.shotid, self.timeid), comments='#')
        self.ne_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_ne.dat' % (self.shotid, self.timeid), comments='#')
        self.nD_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_nd.dat' % (self.shotid, self.timeid), comments='#')
        self.nC_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_nc.dat' % (self.shotid, self.timeid), comments='#')
        self.Te_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_Te.dat' % (self.shotid, self.timeid), comments='#')
        self.Ti_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_Ti.dat' % (self.shotid, self.timeid), comments='#')
        self.vpolC_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_vpolC.dat' % (self.shotid, self.timeid), comments='#')
        self.vpolD_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_vpolD.dat' % (self.shotid, self.timeid), comments='#')
        self.vtorD_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_vtorD.dat' % (self.shotid, self.timeid), comments='#')
        self.vtorC_data = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_vtorC.dat' % (self.shotid, self.timeid), comments='#')
        self.psirz_exp = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_psirz.txt' % (self.shotid, self.timeid), comments='#')
        self.wall_exp = np.genfromtxt(os.path.dirname(__file__) + '/TestBaseProfiles/gt3_diiid_wall.dat', comments='#')
        self.wall_line = LineString(self.wall_exp)
        self.nbeams_loc = os.path.dirname(__file__) + '/../nbeams/bin/Release/nbeams'
        self.beams_json = os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_beams.json' % (self.shotid, self.timeid)
        self.beams_out_json = os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_outbeams.json' % (self.shotid, self.timeid)
        self.neutfile_loc = os.path.dirname(__file__) + '/TestBaseProfiles/gt3_%s_%s_outneuts.json' % (self.shotid, self.timeid)

    def print_summary(self):
        print("0D Parameters \n")
        table1 = texttable.Texttable()
        table1.set_cols_align(["l", "l", "m", "m"])
        table1.set_cols_valign(["m", "m", "m", "m"])
        table1.add_rows([["Variable name", "Description", "Value", "Units"],
                        ["shotid", "DIII-D Shot ID", self.shotid, "n/a"],
                        ["timeid", "DIII-D Time ID", self.timeid, "n/a"],
                        ["BT0", "Toroidal magnetic field at rho=0", self.BT0, "T"],
                        ["wall_ni_min", "Minimum wall ion density", self.wall_ni_min, "#/m^3"],
                        ["wall_ne_min", "Minimum wall electron density", self.wall_ne_min, "#/m^3"],
                        ["wall_Ti_min", "Minimum wall ion temperature", self.wall_Ti_min, "keV"],
                        ["wall_Te_min", "Minimum wall electron temperature", self.wall_Te_min, "keV"],
                        ["pfr_ni_val", "Ion density in the private flux region", self.pfr_ni_val, "#/m^3"],
                        ["pfr_ne_val", "Electron density in the private flux region", self.pfr_ne_val, "#/m^3"],
                        ["pfr_ti_val", "Ion temperature in the private flux region", self.pfr_Ti_val, "keV"],
                        ["pfr_te_val", "Electron temperature in the private flux region", self.pfr_Te_val, "keV"],
                        ["ebeam", "Power-averaged NBI beam energy", self.ebeam, "keV"],
                        ["abeam", "Some value for NBeams",  self.abeam, "1"],
                        ["alphain", "Power-averaged launch angle (maybe)", self.alphain, "??"],
                        ["pbeam", "Duty-cycle-averaged beam power", self.pbeam, "MW"],
                        ["rtang", "Radius of tangency of beam ions", self.rtang, "m"],
                        ["R_loss", "Fraction of IOL particles that do not return on loss orbits", self.R_loss, "1"]
                        ])

        print(table1.draw() + "\n")

        print("1D Parameters \n")
        table2 = texttable.Texttable()
        table2.set_cols_align(["l", "l", "c"])
        table2.set_cols_valign(["m", "m", "m"])
        table2.add_rows([["Variable name", "Description", "Units"],
                         ["Er", "Radial electric field", "V/m"],
                         ["ne", "Electron density", "#/m^3"],
                         ["nD", "Main ion density", "#/m^3"],
                         ["nC", "Carbon density", "#/m^3"],
                         ["Te", "Electron temperature", "keV"],
                         ["Ti", "Ion temperature", "keV"],
                         ["VpolC", "Poloidal carbon velocity", "m/s"],
                         ["VpolD", "Poloidal main ion velocity", "m/s"],
                         ["VtorC", "Toroidal carbon velocity", "m/s"],
                         ["VtorD", "Toroidal main ion velocity", "m/s"],
                         ["wall_line", "LineString object for the machine wall", "n/a"]
                         ])
        print(table2.draw() + "\n")

        print("Computation-related Parameters \n")
        table3 = texttable.Texttable()
        table3.set_cols_align(["l", "l", "m"])
        table3.set_cols_valign(["m", "m", "m"])
        table3.add_rows([["Variable name", "Description", "Value"],
                         ["thetapts_approx", "Approximate number of theta points (may be modified in situ)", self.thetapts_approx],
                         ["edge_rho", "Definition of the plasma edge", self.edge_rho],
                         ["rhopts_edge", "Number of rho values in the edge", self.rhopts_edge],
                         ["rhopts_core", "Number of rho values in the core", self.rhopts_core],
                         ["core_thetapts_ntrl", "Number of theta values in the neutrals core calculation", self.core_thetapts_ntrl],
                         ["edge_rho_ntrl", "Definition of the edge in the neutrals calculation", self.edge_rho_ntrl],
                         ["rhopts_edge_ntrl", "Number of rho values in the edge in the neutrals calculation", self.rhopts_edge_ntrl],
                         ["rhopts_core_ntrl", "Number of rho values in the core in the neutrals calculation", self.rhopts_core_ntrl],
                         ["sollines_psi_max", "Maximum psi value after the scrape-off layer", self.sollines_psi_max],
                         ["num_sollines", "Number of lines after the scrape-off layer", self.num_sollines],
                         ["numcos", "Number of cosine angles used in IOL calculation", self.numcos]
                         ])
        print(table3.draw())

        print("Nbeams location: %s" % self.nbeams_loc)

if __name__ == "__main__":

    plasma = getGT3Test()
    plasma.run_NBI()
    print("DONE")
