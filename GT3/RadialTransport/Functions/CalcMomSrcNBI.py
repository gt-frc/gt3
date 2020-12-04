#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from GT3.RadialTransport.Functions.CalcTorque import calc_torque
from math import cos
from deprecation import deprecated

@deprecated(deprecated_in="0.0.3", removed_in="0.0.4", details="NBI-related calculations are now performed in the NBI module")
def calc_mom_src_nbi(beam,  n, z_eff, R0_a, fforb):
    """Calculates toroidal momentum input from a neutral beam (or beam component)

    :param fforb:
    :param beam:
    :param n:
    :param z_eff:
    :param rpts:
    :param R0_a:
    :return:
    """

    def calc_atten(beamE, n, z_eff, rpts=len(n.i)):  # TODO: Figure out what this is. has to do with beam attenuation
        delma = 1 / (rpts - 1)
        alphain = 0.6475  # TODO: no idea what this is
        xlam = 5.5E17 * (beamE / 1) / (2.0 * 0.5 * n.i * z_eff ** 0.5)
        atten = 1 - np.exp(-delma / cos(alphain) / xlam)
        return atten

#
#   Note: Legacy GTEDGE had the beam energy in keV I THINK. Using Joules just gives 0 torque.
#
    try:
        beam.E.J[1] #checks to see if multiple beams
        atten, unatten, torque, tor_mom = [], [], [], []
        for i in range(len(beam.E.J)):
            atten.append(calc_atten(beam.E.kev[i] , n, z_eff))
            unatten.append(1-atten[i])
            torque.append(calc_torque(beam, fforb, index=i))
            if beam.coI[i]:
                tor_mom.append(torque[i]*unatten*atten/R0_a)
            else:
                tor_mom.append(-1 * torque[i] * atten * unatten / R0_a)

        return np.asarray([sum(x) for x in zip(*tor_mom)][0])
    except:
        atten = calc_atten(beam.E.kev, n, z_eff)
        unatten = 1-atten  # might need to flip atten...
        torque = calc_torque(beam, fforb)

        tor_mom = unatten*atten*torque/R0_a
    return tor_mom