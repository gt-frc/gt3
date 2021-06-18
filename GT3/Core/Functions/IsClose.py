#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# isclose is included in python3.5+, so you can delete this if the code ever gets ported into python3.5+
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)