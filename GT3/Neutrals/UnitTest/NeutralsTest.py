#!/usr/bin/python

from GT3.TestBase import TestClass
import unittest
from GT3 import gt3
import neutpy

class NeutralsUnitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testPlasma = gt3(preparedInput=TestClass())

class NeutralsTestCase(NeutralsUnitTest):
    @classmethod
    def setUpClass(cls):
        super(NeutralsTestCase, cls).setUpClass()
    def test_get_plasma(self):
        self.assertIsInstance(self.testPlasma, gt3)
    def test_gt3_has_objs(self):
        self.assertTrue(hasattr(self.testPlasma, "core"))
        self.assertTrue(hasattr(self.testPlasma, "inp"))
    def test_neutpy_has_gt3_method(self):
        self.assertTrue(hasattr(neutpy.neutrals, "from_gt3"))
    def test_neutpy_run(self):
        self.npi = neutpy.neutrals()
        self.npi.from_gt3(self.testPlasma.core, self.testPlasma.inp)
    def test_neutpy_has_nn(self):
        self.assertTrue(hasattr(self.npi.nn, "s"))
        self.assertTrue(hasattr(self.npi.nn, "t"))
        self.assertTrue(hasattr(self.npi.nn, "tot"))
    def test_neutpy_has_izn(self):
        self.assertTrue(hasattr(self.npi.izn_rate, "s"))
        self.assertTrue(hasattr(self.npi.izn_rate, "t"))
        self.assertTrue(hasattr(self.npi.izn_rate, "tot"))

if __name__ == '__main__':
    unittest.main()