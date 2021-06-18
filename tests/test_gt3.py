#!/usr/bin/python

import unittest
import GT3
from tests.ShotBase import *
from matplotlib.axes._axes import Axes
from matplotlib.pyplot import Figure


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


if __name__ == '__main__':
    unittest.main()
