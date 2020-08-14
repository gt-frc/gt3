#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import e


def calc_vsep(z, m, p):
    """Calculates V_sep"""
    a = (np.abs(p.B1 / p.B0) * p.f0 / p.f1 * p.zeta0) ** 2 - 1 + (1 - p.zeta0 ** 2) * np.abs(p.B1 / p.B0)
    b = 2 * z * e * (p.psi0 - p.psi1) / (p.R1 * m * p.f1) * np.abs(p.B1 / p.B0) * p.f0 / p.f1 * p.zeta0
    c = (z * e * (p.psi0 - p.psi1) / (p.R1 * m * p.f1)) ** 2 - 2 * z * e * (p.phi0 - p.phi1) / m

    v_sep_1 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    v_sep_2 = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

    v_sep = np.zeros(p.r0.shape)
    v_sep = np.where(
        np.logical_and(
            np.logical_or(v_sep_1 <= 0, np.isnan(v_sep_1)),
            np.logical_or(v_sep_2 <= 0, np.isnan(v_sep_2))
        ),
        np.nan, v_sep)
    v_sep = np.where(
        np.logical_or(
            np.logical_and(v_sep_1 > 0, np.logical_or(v_sep_2 <= 0, np.isnan(v_sep_2))),
            np.logical_and(np.logical_and(v_sep_1 > 0, v_sep_2 > 0), v_sep_1 < v_sep_2)
        )
        , v_sep_1, v_sep)
    v_sep = np.where(
        np.logical_or(
            np.logical_and(v_sep_2 > 0, np.logical_or(v_sep_1 <= 0, np.isnan(v_sep_1))),
            np.logical_and(np.logical_and(v_sep_1 > 0, v_sep_2 > 0), v_sep_2 <= v_sep_1)
        )
        , v_sep_2, v_sep)

    v_sep_min = np.nanmin(np.nanmin(v_sep, axis=0), axis=2).T
    v_sep_min[-1] = 0
    return v_sep, v_sep_min
