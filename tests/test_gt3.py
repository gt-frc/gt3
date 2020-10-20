#!/usr/bin/python

import unittest
import GT3
from tests.ShotBase import *
from matplotlib.axes._axes import Axes

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

    def test_sol_exists(self):
        self.plasma.run_SOL()
        self.assertTrue(hasattr(self.plasma, "sol"))
        self.assertIsInstance(self.plasma.sol, GT3.Sol)

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
        self.assertIsInstance(fig, Axes)
        self.plt.close(fig.get_figure())

    def test_plot_core(self):
        """
        Plot all plots in the Core module
        """
        plot_vars = [self.plasma.core.plot_t,
                     self.plasma.core.plot_n,
                     self.plasma.core.plot_ni,
                     self.plasma.core.plot_ntrl_s,
                     self.plasma.core.plot_izn_rate_total,
                     self.plasma.core.plot_izn_rate_t,
                     self.plasma.core.plot_ntrl_total,
                     self.plasma.core.plot_ntrl_s,
                     self.plasma.core.plot_ntrl_t,
                     self.plasma.core.plot_er,
                     self.plasma.core.plot_L_n_e,
                     self.plasma.core.plot_L_n_i,
                     self.plasma.core.plot_L_t_e,
                     self.plasma.core.plot_L_t_i,
                     self.plasma.core.plot_pressure_C,
                     self.plasma.core.plot_pressure_D,
                     self.plasma.core.plot_pressure_e,
                     self.plasma.core.plot_vpol_C,
                     self.plasma.core.plot_vpol_D,
                     self.plasma.core.plot_vtor_C,
                     self.plasma.core.plot_vtor_D]
        for v in plot_vars:
            self.plot_tester(v)



    def test_plot_beams(self):
        """
        Plot all plots in the NBI module
        """
        plot_vars = [self.plasma.nbi.plot_S_nbi,
                     self.plasma.nbi.plot_S_dens_nbi,
                     self.plasma.nbi.plot_M_nbi,
                     self.plasma.nbi.plot_M_dens_nbi,
                     self.plasma.nbi.plot_Q_nbi,
                     self.plasma.nbi.plot_Q_dens_nbi]
        for v in plot_vars:
            self.plot_tester(v)

    def test_plot_rtrans(self):
        """
        Plot all plots in the Radial Transport module
        """
        plot_vars = [self.plasma.rtrans.plot_n,
                     self.plasma.rtrans.plot_ni,
                     self.plasma.rtrans.plot_D,
                     self.plasma.rtrans.plot_gamma_diff,
                     self.plasma.rtrans.plot_T,
                     self.plasma.rtrans.plot_gamma,
                     self.plasma.rtrans.plot_Chi_i_comp,
                     self.plasma.rtrans.plot_Q_sources,
                     self.plasma.rtrans.plot_S_sources,
                     self.plasma.rtrans.plot_chi_terms,]
        self.plot_tester(self.plasma.rtrans.plot_gamma_diff)


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

if __name__ == '__main__':
    unittest.main()
