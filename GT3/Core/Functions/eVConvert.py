#!/usr/bin/python
from collections import namedtuple

def eVConvert(e):
    """
    Given a value in eV, returns an EnergyConversion named tuple containing the value in eV, keV, and J.

    :param e:
    :type e: float
    :return:
    """
    EnergyConversions = namedtuple('EnergyConversion', 'eV keV J')
    return EnergyConversions(e, e * 1E-3, e * 1.0618E-19)