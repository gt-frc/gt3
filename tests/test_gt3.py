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


    def test_plot_core(self):
        fig = self.plasma.core.plot_t(edge=True)
        self.assertIsInstance(fig, Axes)

    def test_plot_core_n(self):
        self.assertIsInstance(self.plasma.core.plot_n(edge=True), Axes)

    def plot_core_n_s(self):
        self.assertIsInstance(self.plasma.core.plot_ntrl_s(edge=True), Axes)


    def test_plot_beams(self):
        self.plasma.nbi.plot_S_dens_nbi()
        self.plasma.nbi.plot_S_nbi()

    def test_plot_rtrans_gamma(self):
        self.plasma.rtrans.plot_gamma_diff()

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
