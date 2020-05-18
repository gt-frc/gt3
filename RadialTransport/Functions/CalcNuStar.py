#!/usr/bin/env python2
# -*- coding: utf-8 -*-


def calc_nustar(nu_c, q95, R0_a, vpol):
    nustar = nu_c * abs(q95) * R0_a / vpol
    return nustar