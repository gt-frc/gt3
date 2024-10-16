#!/usr/bin/python


from numpy import log as ln
from scipy.interpolate import UnivariateSpline
import numpy as np

def get_mfp(ne, beamE, Te, beamA, Zeff, rho):
    """
    Returns the mean free path as a Univariate Spline

    :return: A univariate spline of the MFP function
    :param x: The rho value (in [0., 1]
    :type x: float
    :param ne: The elecron density ndarray
    :type ne:  np.ndarray
    :param beamE: The beam energy
    :type beamE:  float
    :param Te: The electron temperature array
    :type Te:  np.ndarray
    :param beamA: The beam ion mass in u
    :type beamA:  float
    :param Zeff:  The Z-effective array
    :type Zeff:  np.ndarray
    :return:
    :param rho: The rho array
    :type rho:  np.ndarray

    :return:
    """

    ne_spline = UnivariateSpline(rho, ne / 1E6)
    sig_eff_spline = sigeff(ne, beamE, Te, beamA, Zeff, rho)
    result = np.array([(1. / (ne_spline(x) * sig_eff_spline(x))) * 1E-2 for x in rho])

    if not (result > 0).all():
        raise ValueError("Negative MFP encountered")
    else:
        return UnivariateSpline(rho, result, k=3, s=2)

def sigeff(ne, beamE, Te, beamA, Zeff, rho):

    """
    Calculates effective cross sections for a single-impurity plasma. Returns a Univariate spline.

    Reference:
    R.K. Janev, C.D. Boley, D.E. Post. Penetration of energetic neutral beams into fusion plasmas.
    NUCLEAR  FUSION,  Vol.29,  No. 12  (1989)

    :param kwargs:
    :return: A Univariate spline
    :param rho: The rho array
    :type rho: np.array
    :param E: The beam energy
    :type float
    :param ne: The electron density
    :type ne: np.array
    :param Te: The electron temperature
    :type Te: np.array
    :param Zeff: The effective Z
    :type Zeff: np.array
    :param nz: The impurity density
    :type nz: np.array
    """


    expTerm = np.exp(S1(ne, beamE, Te, beamA, rho)) / (beamE / beamA)
    rest = (1. + (Zeff - 1.) * SC(ne, beamE, Te, beamA)) * 1E-16
    fun = expTerm * rest
    return UnivariateSpline(rho, fun)


def S1(ne, E, Te, beamA, rho):
    """
    Calculates S1 coefficient for hydrogenic plasma impurities

    :param ne: Electron density (in m^-3)
    :type ne: np.ndarray
    :param E: Beam ion temperature  (in eV)
    :type E: np.ndarray
    :param Te: Plasma electron temperature  (in keV)
    :type Te: np.ndarray
    :param beamA: Atomic mass of beam ions (in u)
    :type beamA: float
    :return result: The rho-dependent S1
    :type result: np.array
    """

    A1 = np.zeros((3, 3, 2))
    A1[0, 0, 0] = 3.95E+00
    A1[0, 0, 1] = 1.60E-02
    A1[0, 1, 0] = -3.84E-02
    A1[0, 1, 1] = -5.98E-03
    A1[0, 2, 0] = -3.10E-03
    A1[0, 2, 1] = -1.09E-03
    A1[1, 0, 0] = 3.67E-01
    A1[1, 0, 1] = -2.15E-02
    A1[1, 1, 0] = 3.07E-02
    A1[1, 1, 1] = 1.78E-03
    A1[1, 2, 0] = 3.16E-03
    A1[1, 2, 1] = 3.47E-04
    A1[2, 0, 0] = -9.95E-03
    A1[2, 0, 1] = 6.19E-04
    A1[2, 1, 0] = -2.36E-03
    A1[2, 1, 1] = -1.67E-04
    A1[2, 2, 0] = -1.31E-04
    A1[2, 2, 1] = -2.28E-05

    result = np.zeros(len(rho))
    for n, val in enumerate(rho):
        for i in [0,1]:
            for j in [0, 1, 2]:
                for k in [0, 1]:
                    result[n] += A1[i, j, k] * (ln(E * 1.E-3/beamA)**(i)) * ( ( ln((ne[n] / 1E6) / 1E13) ) )**(j) * (ln(Te[n]))**(k)
    return result

def SC(ne, E, Te, beamA):
    """
    Calculates S2 coefficient for Carbon impurities

    :param ne: Electron density (in m^-3)
    :type ne: np.ndarray
    :param E: Beam ion temperature  (in eV)
    :type E: np.ndarray
    :param Te: Electron temperature  (in keV)
    :type Te: np.ndarray
    :param beamA: Atomic mass of beam ions (in u)
    :type beamA: float
    :return result: The rho-dependent SC
    :type result: np.array
    :return:
    """
    # TODO: Get improved coefficients since Mandrekas used improved coefficients in S1

    BC = np.zeros((4, 2, 2))
    # BC[0, 0, 0] = -1.49E0
    # BC[0, 0, 1] = -1.54E-2
    # BC[0, 1, 0] = -1.19E-1
    # BC[0, 1, 1] = -1.50E-2
    # BC[1, 0, 0] = 5.18E-1
    # BC[1, 0, 1] = 7.18E-3
    # BC[1, 1, 0] = 2.92E-2
    # BC[1, 1, 1] = 3.66E-3
    # BC[2, 0, 0] = -3.36E-2
    # BC[2, 0, 1] = 3.41E-4
    # BC[2, 1, 0] = -1.79E-3
    # BC[2, 1, 1] = -2.14E-4

    # Updated coefficients from S Suzuki et al 1998 Plasma Phys. Control. Fusion40 209

    # BC[0, 0, 0] = -1.54E0
    # BC[0, 0, 1] = -8.68E-2
    # BC[0, 1, 0] = -1.83E-1
    # BC[0, 1, 1] = -1.53E-2
    # BC[1, 0, 0] = 5.67E-1
    # BC[1, 0, 1] = 3.13E-3
    # BC[1, 1, 0] = 4.60E-2
    # BC[1, 1, 1] = 4.21E-3
    # BC[2, 0, 0] = -3.86E-2
    # BC[2, 0, 1] = -1.60E-3
    # BC[2, 1, 0] = -2.68E-3
    # BC[2, 1, 1] = -2.41E-4

    # Updated coefficients via NBEAMS from Mandrekas/Boyl

    BC[0, 0, 0] = -1.89E-01
    BC[0, 0, 1] = -3.22E-02
    BC[0, 1, 0] = 5.43E-02
    BC[0, 1, 1] = 3.83E-03
    BC[1, 0, 0] = -1.41E-03
    BC[1, 0, 1] = 8.98E-03
    BC[1, 1, 0] = -3.34E-02
    BC[1, 1, 1] = -1.97E-03
    BC[2, 0, 0] = 3.10E-02
    BC[2, 0, 1] = 7.43E-04
    BC[2, 1, 0] = 5.08E-03
    BC[2, 1, 1] = 1.96E-04
    BC[3, 0, 0] = -2.54E-03
    BC[3, 0, 1] = -4.21E-05
    BC[3, 1, 0] = -2.20E-04
    BC[3, 1, 1] = 8.32E-07


    result = np.zeros(len(ne))
    for rho in range(0, len(ne)):
        for i in [0, 1, 2, 3]:
            for j in [0, 1]:
                for k in [0, 1]:
                    result[rho] += BC[i, j, k] * (ln(E * 1.E-3/beamA)**(i)) * (ln((ne[rho] / 1E6) / 1E13))**(j)\
                                   * (ln(Te[rho]))**(k)
    return result