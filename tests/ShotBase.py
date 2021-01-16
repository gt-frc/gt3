#!/usr/bin/python

import unittest



class SingleLowerNullTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from GT3 import gt3
        cls.plasma = gt3(inputFile="tests/TestProfiles/unit_test_d3d_144977_3000")


class DoubleNullTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from GT3 import gt3
        cls.plasma = gt3(inputFile="tests/TestProfiles/unit_test_d3d_175826_2010")


class NegativeTriangularityTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from GT3 import gt3
        cls.plasma = gt3(inputFile="tests/TestProfiles/unit_test_d3d_170672_1900")