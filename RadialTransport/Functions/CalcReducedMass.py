#!/usr/bin/env python2
# -*- coding: utf-8 -*-

def calc_reduced_mass(m1, m2):
    return m1 * (1 + (m1 / m2))  # * 1E3  # TODO: not sure why this 1E3 is here
