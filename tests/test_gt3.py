#!/usr/bin/python

import unittest
import GT3
from tests.ShotBase import *
from matplotlib.axes._axes import Axes
from matplotlib.pyplot import Figure
from GT3.Core.Functions.ProfileClasses import TwoDProfile, OneDProfile
from shapely.geometry import Point
import numpy as np


class CommonFunctions(object):
    """Tests to see if a shot has the expected attributes typical for a fully run shot."""

    def test_gt3_has_core(cls):
        cls.assertTrue(hasattr(cls.plasma, "core"))

    def test_gt3_has_iol(cls):
        cls.assertTrue(hasattr(cls.plasma, "iol"))

    def test_gt3_has_nbi(cls):
        cls.assertTrue(hasattr(cls.plasma, "nbi"))

    def test_gt3_has_rtrans(cls):
        cls.assertTrue(hasattr(cls.plasma, "rtrans"))


class SingleNullRun(SingleLowerNullTest, CommonFunctions):

    @classmethod
    def setUpClass(cls):
        super(SingleNullRun, cls).setUpClass()
        cls.plasma.run_radial_transport()


class DoubleNullRun(DoubleNullTest, CommonFunctions):

    @classmethod
    def setUpClass(cls):
        super(DoubleNullRun, cls).setUpClass()
        cls.plasma.run_radial_transport()


class NegativeTriangularityRun(NegativeTriangularityTest, CommonFunctions):

    @classmethod
    def setUpClass(cls):
        super(NegativeTriangularityRun, cls).setUpClass()
        cls.plasma.run_radial_transport()


class RunModificationTest(SingleLowerNullTest):

    # def test_sol_exists(self):
    #     self.plasma.run_SOL()
    #     self.assertTrue(hasattr(self.plasma, "sol"))
    #     self.assertIsInstance(self.plasma.sol, GT3.Sol)

    def test_iol_exists(self):
        self.plasma.run_IOL()
        self.assertTrue(hasattr(self.plasma, "iol"))
        self.assertIsInstance(self.plasma.iol, GT3.IOL)

    def test_nbi_exists(self):
        self.plasma.run_NBI()
        self.assertTrue(hasattr(self.plasma, "nbi"))
        self.assertIsInstance(self.plasma.nbi, GT3.BeamDeposition)

    def test_rtrans_exists(self):
        self.plasma.run_radial_transport()
        self.assertTrue(hasattr(self.plasma, "rtrans"))
        self.assertIsInstance(self.plasma.rtrans, GT3.RadialTransport)

    def test_rtrans_modifiers(self):
        try:
            del self.plasma.rtrans
        except:
            pass
        kwargs = {
            'rtrans_override': {
                'splines': {
                    'T_i': True,
                    'T_e': True,
                    'n_i': True,
                    'n_e': True,
                }
            }
        }
        self.plasma.run_radial_transport(**kwargs)
        self.assertTrue(self.plasma.rtrans._T.i.kev._data_overwritten)
        self.assertTrue(self.plasma.rtrans._T.e.kev._data_overwritten)
        self.assertTrue(self.plasma.rtrans._T.i.ev._data_overwritten)
        self.assertTrue(self.plasma.rtrans._T.e.ev._data_overwritten)
        self.assertTrue(self.plasma.rtrans._n.i._data_overwritten)
        self.assertTrue(self.plasma.rtrans._n.e._data_overwritten)
        self.plasma.rtrans.set_chi_asymmetry(0.05, 0.05)


class PlotCoreTest(DoubleNullTest):

    @classmethod
    def setUpClass(cls):
        super(PlotCoreTest, cls).setUpClass()
        import matplotlib.pyplot as plt

        cls.plt = plt
        cls.plasma.run_radial_transport()
        cls.plt.ion()

    def plot_tester(self, plotter, edge=False):
        import inspect
        args = inspect.getfullargspec(plotter)
        if 'logPlot' in args and 'edge' in args:
            fig = plotter(logPlot=True, edge=True)
        elif 'logPlot' in args:
            fig = plotter(logPlot=True)
        elif 'edge' in args:
            fig = plotter(edge=True)
        else:
            fig = plotter()
        self.assertIsInstance(fig, (Figure, Axes))
        self.plt.close(fig.get_figure())

    def test_plot_core(self):
        """
        Plot all plots in the Core module
        """

        plot_vars = [self.plasma.core.n.i.fsa.plot,
                     self.plasma.core.n.e.fsa.plot,
                     self.plasma.core.n.n.s.plot2D,
                     self.plasma.core.n.n.t.plot2D,
                     self.plasma.core.n.n.tot.plot2D,
                     self.plasma.core.T.i.ev.plot2D,
                     self.plasma.core.T.i.J.plot2D,
                     self.plasma.core.T.i.kev.plot2D,
                     self.plasma.core.T.e.ev.plot2D,
                     self.plasma.core.T.i.ev.L.plot2D,
                     self.plasma.core.T.e.J.L.plot2D,
                     self.plasma.core.n.i.L.plot2D,
                     self.plasma.core.n.n.s.L.plot2D,
                     self.plasma.core.n.n.tot.L.plot2D,
                     self.plasma.core.v.D.pol.plot2D,
                     self.plasma.core.v.C.tor.plot2D]

        for v in plot_vars:
            self.plot_tester(v)

    def test_plot_beams(self):
        """
        Plot all plots in the NBI module
        """
        plot_vars = [self.plasma.nbi.combined_beam_src_dens_lost.Snbi.plot,
                     self.plasma.nbi.combined_beam_src_dens_lost.Qnbi.plot,
                     self.plasma.nbi.combined_beam_src_dens_lost.Mnbi.plot,
                     self.plasma.nbi.combined_beam_src_kept.Snbi.plot,
                     self.plasma.nbi.combined_beam_src_kept.Qnbi.plot,
                     self.plasma.nbi.combined_beam_src_kept.Mnbi.plot]
        for v in plot_vars:
            self.plot_tester(v)

    def test_plot_rtrans(self):
        """
        Plot all plots in the Radial Transport module
        """
        plot_vars = [self.plasma.rtrans.gamma.D.diff.plot,
                     self.plasma.rtrans.gamma.D.int.plot,
                     self.plasma.rtrans.gamma.e.int.plot,
                     self.plasma.rtrans.gamma.C.int.plot,
                     self.plasma.rtrans.gamma.plot_D,
                     self.plasma.rtrans.gamma.plot_C,
                     self.plasma.rtrans.plot_Q_sources,
                     self.plasma.rtrans.plot_S_sources,
                     self.plasma.rtrans.plot_chi_terms]
        for v in plot_vars:
            self.plot_tester(v)

    @classmethod
    def tearDownClass(cls):
        cls.plt.clf()
        cls.plt.close()


class PlotIOLTest(DoubleNullTest):

    @classmethod
    def setUpClass(cls):
        super(PlotIOLTest, cls).setUpClass()
        import matplotlib.pyplot as plt
        cls.plt = plt
        cls.plasma.run_IOL()

    def test_plot_iol_F_i(self):
        self.plasma.iol.plot_F_i(edge=True)
        self.assertIsInstance(self.plasma.iol.plot_F_i(), Axes)


class GT3TestClassTest(unittest.TestCase, CommonFunctions):
    @classmethod
    def setUpClass(cls):
        super(GT3TestClassTest, cls).setUpClass()
        from GT3 import gt3
        from GT3.TestBase.testbase import TestClass
        cls.plasma = gt3(preparedInput=TestClass())
        cls.plasma.run_radial_transport()
        TestClass().print_summary()


class ProfileClassTest(SingleLowerNullTest):
    plasma = None

    @classmethod
    def setUpClass(cls):
        import matplotlib.pyplot as plt
        super(ProfileClassTest, cls).setUpClass()
        cls.plt = plt
        cls.plasma.run_radial_transport()
        cls.plt.ion()
        cls.data = cls.plasma.core.T.i.kev * 0.1
        cls.raw_data = cls.data[:,0]
        cls.rho = cls.plasma.core.T.i.kev.rho[:,0]
        cls.docs = "A test TwoDProfile"
        cls.units = "/m^3"
        cls.title = "Test TwoDProfile"
        cls.xLabel = r"$/rho$"
        cls.yLabel = r"Test"

        cls.dataProfile = cls._get_default_twod_profile(cls)


    def test_TwoDProfile(self):
        self.assertEquals(self.dataProfile.__doc__, self.docs)

        self.assertIsInstance(self.dataProfile.fsa, OneDProfile)
        self.assertIsInstance(self.dataProfile.fsa, OneDProfile)
        self.assertIsInstance(self.dataProfile.L, TwoDProfile)
        self.assertIsInstance(self.dataProfile.plot_raw_data(), Axes)
        self.assertIsInstance(self.dataProfile.plot_raw_data(fsa=False), Axes)
        self.assertIsInstance(self.dataProfile.plot_fsa(), Axes)
        self.assertIsInstance(self.dataProfile.getGrid(), np.ndarray)

    def test_TwoDProfile_errors(self):
        self.assertRaises(IndexError, self.dataProfile.set_raw_data, "")
        del self.dataProfile._raw_data
        del self.dataProfile._raw_rho
        self.assertWarns(Warning, self.dataProfile.plot_raw_data)

        self.assertRaises(TypeError, self._create_2d_type_error)
        self.assertWarns(Warning, self.dataProfile.set_wall, Point((0., 0.)))
        del self.dataProfile.wall
        self.assertWarns(Warning, self.dataProfile.plot2D)
        self.dataProfile.set_wall(self.plasma.core.wall_line)

    def _create_2d_type_error(self):
        return TwoDProfile(self.plasma.core.psi,
                                      [3., 4., 5.],
                                      self.plasma.core.R,
                                      self.plasma.core.Z)

    def test_TwoDProfile_mods(self):
        self.dataProfile_mod = self._get_default_twod_profile()
        self.dataProfile_mod.set_to_zeros()
        self.assertFalse(self.dataProfile_mod.isNonZero())

    def _get_default_twod_profile(self):
        return TwoDProfile(self.plasma.core.psi,
                           self.data,
                           self.plasma.core.R,
                           self.plasma.core.Z,
                           wall=self.plasma.core.wall_line,
                           docs=self.docs,
                           units=self.units,
                           plotTitle=self.title,
                           xLabel=self.xLabel,
                           yLabel=self.yLabel,
                           raw=(self.rho, self.raw_data))

    @classmethod
    def tearDownClass(cls):
        pass



if __name__ == '__main__':
    unittest.main()
