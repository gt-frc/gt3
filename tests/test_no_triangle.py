#!/usr/bin/python

import unittest
from tests.ShotBase import *

"""
Unit tests to be run on a system without triangle to ensure GT3 doesn't crash.
"""

class NoTriangleTest(SingleLowerNullTest):

    @classmethod
    def setUpClass(cls):
        super(NoTriangleTest, cls).setUpClass()
        from shutil import which
        if which("triangle"):
            raise Exception("Triangle has been found on this system. This test will not run")
        cls.plasma.run_neutrals()

    def test_no_npi(self):
        self.assertFalse(hasattr(self.plasma.ntrl, "npi"))


