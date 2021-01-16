#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from deprecation import deprecated


@deprecated(deprecated_in="0.0.3", removed_in="0.0.4", details="NBI-related calculations are now performed in the NBI module")
def calc_torque(beam, fforb, index=False):
    """ Calculates torque from a neutral beam (or beam component)

    torque = F * r_tan = (P/v) * r_tan = (P/sqrt(2E/m)) * r_tan = P * sqrt(m/(2E)) * r_tan

    :param fforb:
    :param index:
    :param beam: beam object with attributes z, m, a, en, pwr, rtan
    :return: torque
    """
    if index is not False:

        power = beam.P.W[index]
        energy = beam.E.J[index]
        mass = beam.m
        rtan = beam.rtang[index]
        torque = power * np.sqrt(0.5 * mass / energy) * rtan * (1.0 - fforb)  # Piper Changes: Included fast ion losses.
        return torque
    else:
        power = beam.P.W
        energy = beam.E.J
        mass = beam.m
        rtan = beam.rtang

        torque = power * np.sqrt(0.5 * mass / energy) * rtan * (1.0-fforb)  # Piper Changes: Included fast ion losses.
        return torque