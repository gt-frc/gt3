#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:10:25 2017

@author: max
"""
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from GT3.utilities.PlotBase import PlotBase
import configparser

from shapely.geometry import LineString


class ReadInfile:
    """Reads main GT3 input file.

    Methods:
        read_vars
        read_exp
        showparams

    Attributes:
        a                (float)    tokamak minor radius (m)
        BT0          (float)    toroidal field strength at mag. axis (T)
        R0_a             (float)    tokamak major radius (m)
        Z0               (float)    vertical height of the magnetic axis (m)
        xpt_R            (float)
        xpt_Z            (float)
        xpt
        thetapts_approx  (int)
        thetapts
        rmeshnum_p       (int)
        rpts             (int)
        ni0              (float)
        ni9              (float)
        ni_sep           (float)
        nu_ni            (float)
        ne0              (float)
        ne9              (float)
        ne_sep           (float)
        nu_ne            (float)
        Ti0              (float)
        Ti9              (float)
        Ti_sep           (float)
        nu_Ti            (float)
        Te0              (float)
        Te9              (float)
        Te_sep           (float)
        nu_Te            (float)
        j0               (float)
        wallfile         (str)
    """

    def __init__(self, infile, **kwargs):
        """
        Initializes the read_infile class with the filename for the to_gt3 file.
        :param infile:
        :type infile: str
        """
        super().__init__()
        sys.dont_write_bytecode = True
        self.read_vars(infile, **kwargs)

    def _float_loader(self, parser, section, name, **kwargs):

        try:
            result = parser.getfloat(section, name)
            try:
                if kwargs.get("inp_override").get(name).get("multi"):
                    result = result * kwargs.get("inp_override").get(name).get("multi")
                    print(name + " is being scaled")
                if kwargs.get("inp_override").get(name).get("add"):
                    result = result + kwargs.get("inp_override").get(name).get("add")
                    print(name + " is being translated")
                if kwargs.get("inp_override").get(name).get("replace"):
                    result = kwargs.get("inp_override").get(name).get("replace")
                    print(name + " is being replaced")
                return result
            except AttributeError:
                return result
        except configparser.NoOptionError as e:
            print("%s not found" % name)
            return e
        except IOError as e:
            print("%s not found" % name)
            return e

    def _profile_loader(self, parser, section, name, **kwargs):
        """ A helper function for loading profiles.

        :param parser: The ConfigParser
        :type parser: ConfigParser.RawConfigParser
        :param section: The section name
        :type section: str
        :param name: The variable name
        :type name: str
        :return:
        """

        try:
            filename = parser.get(section, name)
            filepath = os.path.join(os.getcwd(), filename)
            result = np.genfromtxt(filepath, comments='#')
            try:
                if kwargs.get("inp_override").get(name):
                    temp = np.array((result[:, 0], result[:, 1]))
                    if kwargs.get("inp_override").get(name).get("multi") is not None:
                        multi = kwargs.get("inp_override").get(name).get("multi")
                        print(name+" is being scaled")
                        temp = np.array((result[:, 0], result[:, 1] * multi)).T
                    if kwargs.get("inp_override").get(name).get("add") is not None:
                        add = kwargs.get("inp_override").get(name).get("add")
                        if float(add) == 0.0:
                            return result
                        print(name + " is being translated")
                        temp = np.array((result[:, 0], result[:, 1] + add)).T
                    return temp
                else:
                    return result
            except AttributeError:
                return result
        except configparser.NoOptionError as e:
            print("%s not found" % name)
            return e
        except IOError as e:
            print("%s not found" % name)
            return e

    def read_vars(self, infile, *args, **kwargs):
        """
        Reads variables given infile filename
        :param infile:
        :type infile: str
        """

        config = configparser.RawConfigParser()
        config.read(infile)

        # Grid construction parameters

        self.rhopts = config.getint('Mesh', 'rhopts')
        self.edge_rho = config.getfloat('Mesh', 'edge_rho')
        self.rhopts_edge = config.getint('Mesh', 'rhopts_edge')
        self.rhopts_core = config.getint('Mesh', 'rhopts_core')
        self.thetapts_approx = config.getint('Mesh', 'thetapts_approx')
        try:
            self.Er_scale = config.getfloat('Mesh', 'Er_scale')
        except configparser.NoOptionError:
            self.Er_scale = 1.
        try:
            self.psi_scale = config.getfloat('Mesh', 'psi_scale')
        except configparser.NoOptionError:
            self.psi_scale = 1.
        self.sollines_psi_max = config.getfloat('Mesh', 'sollines_psi_max')
        self.num_sollines = config.getint('Mesh', 'num_sollines')
        self.xi_ib_pts = config.getint('Mesh', 'xi_ib_pts')
        self.xi_ob_pts = config.getint('Mesh', 'xi_ob_pts')
        self.numcos = config.getint('Mesh', 'numcos')


        # Plasma
        self.BT0 = self._float_loader(config, 'Plasma', 'Bt0', **kwargs)
        self.pfr_ni_val = config.getfloat('Plasma', 'pfr_ni_val')
        self.pfr_ne_val = config.getfloat('Plasma', 'pfr_ne_val')
        self.pfr_Ti_val = config.getfloat('Plasma', 'pfr_Ti_val')
        self.pfr_Te_val = config.getfloat('Plasma', 'pfr_Te_val')
        self.R_loss = config.getfloat('Plasma', 'R_loss')
        try:
            self.sep_val = config.getfloat('Plasma', 'sep_val')
        except configparser.NoOptionError:
            self.sep_val = 1.0


        # Neutrals Output File
        self.neutfile_loc = config.get('1DProfiles', 'neutfile_loc')

        # NBeams multi-beam location
        self.beams_json = config.get('1DProfiles', 'beams_json')
        self.beams_out_json = config.get('1DProfiles', 'beams_out_json')

        self.er_data = self._profile_loader(config, '1DProfiles', 'er_file', **kwargs)
        self.jr_data = self._profile_loader(config, '1DProfiles', 'jr_file', **kwargs)
        self.ne_data = self._profile_loader(config, '1DProfiles', 'ne_file', **kwargs)
        self.nD_data = self._profile_loader(config, '1DProfiles', 'nD_file', **kwargs)
        self.nT_data = self._profile_loader(config, '1DProfiles', 'nT_file', **kwargs)
        self.nW_data = self._profile_loader(config, '1DProfiles', 'nW_file', **kwargs)
        self.nBe_data = self._profile_loader(config, '1DProfiles', 'nBe_file', **kwargs)
        self.na_data = self._profile_loader(config, '1DProfiles', 'na_file', **kwargs)
        self.nC_data = self._profile_loader(config, '1DProfiles', 'nC_file', **kwargs)
        self.Te_data = self._profile_loader(config, '1DProfiles', 'Te_file', **kwargs)
        self.Ti_data = self._profile_loader(config, '1DProfiles', 'Ti_file', **kwargs)
        self.TC_data = self._profile_loader(config, '1DProfiles', 'TC_file', **kwargs)
        self.frac_C_data = self._profile_loader(config, '1DProfiles', 'frac_C_file', **kwargs)

        self.vpolC_data = self._profile_loader(config, '1DProfiles', 'vpolC_file', **kwargs)
        self.vtorC_data = self._profile_loader(config, '1DProfiles', 'vtorC_file', **kwargs)
        self.vpolD_data = self._profile_loader(config, '1DProfiles', 'vpolD_file', **kwargs)
        self.vtorD_data = self._profile_loader(config, '1DProfiles', 'vtorD_file', **kwargs)

        if isinstance(self.vpolD_data, (IOError, configparser.NoOptionError)):
            self.omegDpol_data = self._profile_loader(config, '1DProfiles', 'omegDpol_file', **kwargs)
        if isinstance(self.vtorD_data, (IOError, configparser.NoOptionError)):
            self.omegDtor_data = self._profile_loader(config, '1DProfiles', 'omegDtor_file', **kwargs)
        if isinstance(self.vpolC_data, (IOError, configparser.NoOptionError)):
            self.omegCpol_data = self._profile_loader(config, '1DProfiles', 'omegCpol_file', **kwargs)
        if isinstance(self.vtorC_data, (IOError, configparser.NoOptionError)):
            self.omegCtor_data = self._profile_loader(config, '1DProfiles', 'omegCtor_file', **kwargs)

        self.psirz_exp = self._profile_loader(config, '2DProfiles', 'psirz_file')
        self.wall_exp = self._profile_loader(config, 'Wall', 'wall_file')

        self.wall_line = LineString(self.wall_exp)


    def showparams(self):
        """
        Show the parameters of this shot
        """
        raise NotImplementedError
        print('**PARAMETERS FOR SHOT \'{}\'.'.format(self.shotlabel))
        for key in vars(self).items():
            if key[0][0] != '_' and key[0] != 'line' and key[0] != 'infile' and key[0] != 'variable' and key[
                0] != 'value':
                print(('{} = {}'.format(key[0], key[1])))
        print('**END OF PARAMETERS**')

    def _plot_base(self, xLabel=r'$\rho$', yLabel="Value", title="Title", edge=False):

        plot = plt.figure()
        fig = plot.add_subplot(111)
        fig.set_xlabel(xLabel, fontsize=20)
        fig.set_ylabel(yLabel, fontsize=20)
        fig.set_title(title)
        if edge:
            fig.set_xlim(0.85, 1.0)
        plt.show()
        return fig

    def plot_raw_Ti(self, edge=False):
        fig = self._plot_base(self.Ti_data.T[0], edge=edge)
        fig.scatter(self.Ti_data.T[0], self.Ti_data.T[1], color="black", s=8)
        return fig

    def plot_raw_Te(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.Te_data.T[0], self.Te_data.T[1], color="black", s=8)
        return fig

    def plot_raw_Er(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.er_data.T[0], self.er_data.T[1], color="black", s=8)
        return fig

    def plot_raw_ni(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.nD_data.T[0], self.nD_data.T[1], color="black", s=8)
        return fig

    def plot_raw_ne(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.ne_data.T[0], self.ne_data.T[1], color="black", s=8)
        return fig

    def plot_raw_vpol_C(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.vpolC_data.T[0], self.vpolC_data.T[1], color="black", s=8)
        return fig

    def plot_raw_vpol_D(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.vpolD_data.T[0], self.vpolD_data.T[1], color="black", s=8)
        return fig

    def plot_raw_vtor_C(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.vtorC_data.T[0], self.vtorC_data.T[1], color="black", s=8)
        return fig

    def plot_raw_vtor_D(self, edge=False):
        fig = self._plot_base(edge=edge)
        fig.scatter(self.vtorD_data.T[0], self.vtorD_data.T[1], color="black", s=8)
        return fig
