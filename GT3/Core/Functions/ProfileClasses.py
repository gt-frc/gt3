#!/usr/bin/python

from GT3.utilities.PlotBase import PlotBase
from scipy.interpolate import UnivariateSpline, interp1d, griddata
from GT3.Core.Functions.CalcFSA import calc_fsa
from GT3.Core.Functions.CalcGrad import calc_grad
from GT3.Core.Functions.CalcRho2PsiInterp import calc_rho2psi_interp
from shapely.geometry import LineString
from numpy import ndarray, array
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
import GT3.constants as constants
from GT3 import Core
from typing import Union

e = constants.elementary_charge

SlowFastSplit = namedtuple("SlowFastSplit", "s t tot")
TemperatureSplit = namedtuple("TemperatureSplit", "kev ev J")

class BaseMath:
    """
    The BaseMath class provides overrides for how mathematical operations interact with the OneDProfile and
    TwoDProfile classes. This allows instances to, in general, be added, subtracted, multiplied, etc. against other
    instances and against things like int, float, np.ndarray, etc.

    OneDProfile and TwoDProfile each have val attributes that provide the actual arrays and thus is used extensively
    here.
    """
    def __init__(self):
        # Set array priority so that these results attempt to be followed first.
        self.__array_priority__ = 1000.

    def __add__(self, other):

        if isinstance(other, TwoDProfile) or isinstance(other, OneDProfile):
            return self.val + other.val
        else:
            return self.val + other

    def __rtruediv__(self, other):
        div = other / self.val
        return div

    def __pow__(self, power, modulo=None):
        return self.val**power

    def __mul__(self, other):
        if isinstance(other, TwoDProfile) or isinstance(other, OneDProfile) or isinstance(other, np.ndarray):
            if isinstance(other, np.ndarray):
                return self.val * other
            else:
                return self.val * other.val
        else:
            return self.val * other

    def __eq__(self, other):
        return np.equal(self.val, other)

    def __ne__(self, other):
        return np.not_equal(self.val, other)

    def __sub__(self, other):
        return self.val - other

    # Rsub is used when for a - b when a is not a One/TwoDProfile. Python reflects the operators as well.
    # Thus, the order is reversed.
    def __rsub__(self, other):
        return other - self.val

    def __rmul__(self, other):
        if isinstance(other, TwoDProfile) or isinstance(other, OneDProfile):
            return other.val * self.val
        else:
            return other * self.val

    def __truediv__(self, other):
        if isinstance(other, TwoDProfile) or isinstance(other, OneDProfile):
            return self.val / other.val
        else:
            return self.val / other

    def __neg__(self):
        return -1. * self.val

    def __getitem__(self, key):
        return self.val[key]

    def __len__(self):
        return len(self.val)


class OneDProfile(PlotBase, BaseMath):
    """
    This is the OneDProfile class for 1-D variables. It inherits many plotting functions and the BaseMath class for
    mathematical operations.

     R, Z, and Wall are taken in because we allow TemperatureProfiles, DensityProfiles, etc. to be utilized by
     both 1D and 2D profile classes. Instead of adding logic to TemperatureProfiles, etc. to figure out if OneDProfile
     or TwoDProfile is being used, we just accept R,Z, wall at the start and ultimately just not use them.

     TODO: Probably do the TemperatureProfiles, etc. logic instead
     """

    def __init__(self, psi, val, R, Z, docs="", units="", plotTitle="", xLabel=r"$\rho$", yLabel="", raw=None, *args, **kwargs):
        """

        :param psi: This is the Psi class that will be used to
        :param val:
        :param R:
        :param Z:
        :param docs:
        :param units:
        :param plotTitle:
        :param xLabel:
        :param yLabel:
        :param wall:
        :param raw:
        """
        super(OneDProfile, self).__init__()
        super(BaseMath, self).__init__()
        __array_priority__ = 1000
        self.__array_priority__ = 1000

        # If there are no values given, we set this to 0s
        if not np.any(val):
            val = np.zeros(psi.rho[:, 0].shape)
        self.val = val
        self._R = R
        self._Z = Z
        self.a = psi.a
        rho = psi.rho
        if type(rho) is list:
            self._rho1D = array(rho)
        elif type(rho) is ndarray:
            if len(rho.shape) == 2:
                self._rho1D = rho[:, 0]
                self._rho2D = rho
            elif len(rho.shape) == 1:
                self._rho1D = rho
        else:
            raise TypeError("Rho needs to be a 1D np.array, 2D np.array, or list ")

        if type(val) is list:
            self._rho1D = array(val)
        elif type(val) is ndarray:
            if len(val.shape) == 2:
                self.val = val[:, 0]
            elif len(rho.shape) == 1:
                self.val = val
        elif type(val) is OneDProfile:
            self.val = val.val
        else:
            raise TypeError("Value needs to be a 1D np.array, 2D np.array, or list ")

        self._docs = docs
        self._psi = psi
        self.units = units
        self.xLabel, self.yLabel, self.plotTitle = xLabel, yLabel, plotTitle
        self.set_plot_rho1d(self._rho1D)
        self._Spline = None
        if raw:
            try:
                self._raw_rho = raw[0]
                self._raw_data = raw[1]
            except TypeError as err:
                raise(str(err))

        #del self.plot_contours_with_wall
        #del self.plot_with_wall

    def __doc__(self):
        return self._docs

    def plot(self, edge=True, color="red", **kwargs):
        return self._plot_base(self.val, xLabel=self.xLabel, yLabel=self.yLabel,
                               title=self.plotTitle, color=color, edge=edge, **kwargs)

    @property
    def L(self):
        L = -1.0 * self.val / self.derivative()
        self._L = OneDProfile(self._psi, L, self._R, self._Z)
        return self._L

    @L.getter
    def L(self):
        try:
            return self._L
        except:
            L = -1.0 * self.val / self.derivative()
            self._L = OneDProfile(self._psi, L, self._R, self._Z)
            return self._L

    def derivative(self):
        if hasattr(self, "_derivative"):
            return self._derivative
        else:
            derivative = self.Spline.derivative()(self._rho1D * self.a)
            self._derivative = derivative
            return self._derivative


    @property
    def Spline(self, k=3, s=0):
        self._Spline = UnivariateSpline(self._rho1D * self.a, self.val, k=k, s=s)
        return self._Spline

    @Spline.getter
    def Spline(self):
        if self._Spline:
            return self._Spline
        else:
            self._Spline = UnivariateSpline(self._rho1D * self.a, self.val, k=3, s=0)
            return self._Spline

    @Spline.setter
    def Spline(self, val):
        """Input a (k,s) tuple and set the spline fit"""
        if type(val) is not tuple:
            print("Value is not tuple. Spline not modified.")
        else:
            self._Spline = UnivariateSpline(self._rho1D * self.a, self.val, k=val[0], s=val[1])

    @property
    def OneDInterp(self):
        self._1Dinterp = interp1d(self._rho1D * self.a, self.val, kind='cubic')
        return

    @OneDInterp.getter
    def OneDInterp(self):
        return self._1Dinterp

    @OneDInterp.setter
    def OneDInterp(self, kind):
        self._1Dinterp = interp1d(self._rho1D * self.a, self.val, kind=kind)

    def isNonZero(self):
        return not np.all(self.val == 0)

    def set_raw_data(self, raw):
        try:
            self._raw_rho = raw[0]
            self._raw_data = raw[1]
        except TypeError as err:
            raise (str(err))
        return self._raw_rho, self._raw_data

    def plot_raw_data(self, fsa=True, edge=False, **kwargs):
        """
        Plots the experimental data vs. the flux-surface-averaged or [:, 0] data for comparison.
        Set fsa to False to plot the [:, 0] data instead.
        """
        try:
            self._raw_data
        except AttributeError as err:
            print("""
            Raw data not found. If this is a temperature profile, use the kev version, i.e., T.x.kev.plot_raw_data().
            Otherwise, please verify that raw data are attached to this instance.
            """)
            raise err
        if fsa:
            fig = self._plot_base(self.fsa, xLabel=self.xLabel, yLabel=self.yLabel, edge=edge, line=True)
            fig.scatter(self._raw_rho, self._raw_data, marker="x", s=self._markerSize, **kwargs)
        else:
            fig = self._plot_base(self.to1D(), xLabel=self.xLabel, yLabel=self.yLabel, edge=edge, line=True)
            fig.scatter(self._raw_rho, self._raw_data, marker="x", s=self._markerSize, **kwargs)
        return fig


class TwoDProfile(PlotBase, BaseMath):

    def __init__(self, psi, val, R, Z, wall=None, docs="", units="", plotTitle="", xLabel=r"$\rho$", yLabel="", raw=None):
        super(TwoDProfile, self).__init__()
        self.rho = psi.rho
        self._psi = psi
        self.a = psi.a
        if not np.any(val):
            val = np.zeros(psi.rho.shape)
        if type(val) is ndarray:
            if len(val.shape) == 2:
                self.val = val
            else:
                raise TypeError("Value needs to be 2D np.array")
        else:
            raise TypeError("Value needs to be 2D np.array")

        self.val, self.R, self.Z, self.wall = val, R, Z, wall
        self._docs = docs
        self.units = units
        self.plotTitle, self.xLabel, self.yLabel = plotTitle, xLabel, yLabel

        self.set_plot_rho1d(self.rho[: ,0])
        if self.wall:
            self.set_plot_wall(wall)
        self.set_plot_RZ(self.R[:, 0], self.Z[:, 0])
        self.fsa = OneDProfile(self._psi, calc_fsa(self.val, self.R, self.Z), self.R, self.Z)

        if raw:
            try:
                self._raw_rho = raw[0]
                self._raw_data = raw[1]
            except TypeError as err:
                raise(str(err))

    def __doc__(self):
        return self._docs


    @property
    def L(self):
        """
        Returns the gradient scale length of the the value. Currently, we get the FSA and broadcast, but this
        is not accurate and needs to be done properly.

        @TODO: Correctly implement
        :rtype: TwoDProfile
        """
        L = np.broadcast_to(self.fsa.L.val, (self.rho.shape[1], len(self.fsa.L.val))).T
        self.L = TwoDProfile(self._psi, L, self.R, self.Z, wall=self.wall)
        # L = -1.0 * self.val / self.derivative()
        # self._L = TwoDProfile(self._psi, L, self.R, self.Z, wall=self.wall)
        return self._L

    @L.getter
    def L(self):
        try:
            return self._L
        except:
            L = np.broadcast_to(self.fsa.L.val, (self.rho.shape[1], len(self.fsa.L.val))).T
            self.L = TwoDProfile(self._psi, L, self.R, self.Z, wall=self.wall)
            # L = -1.0 * self.val / self.derivative()
            # self._L = TwoDProfile(self._psi, L, self.R, self.Z, wall=self.wall)
            return self._L

    def set_raw_data(self, raw):
        try:
            self._raw_rho = raw[0]
            self._raw_data = raw[1]
        except TypeError as err:
            raise (str(err))
        return self._raw_rho, self._raw_data

    def plot_raw_data(self, fsa=True, edge=False, **kwargs):
        """
        Plots the experimental data vs. the flux-surface-averaged or [:, 0] data for comparison.
        Set fsa to False to plot the [:, 0] data instead.
        """
        try:
            self._raw_data
        except AttributeError as err:
            print("""
            Raw data not found. If this is a temperature profile, use the kev version, i.e., T.x.kev.plot_raw_data().
            Otherwise, please verify that raw data are attached to this instance.
            """)
            raise err
        if fsa:
            fig = self._plot_base(self.fsa, xLabel=self.xLabel, yLabel=self.yLabel, edge=edge, line=True)
            fig.scatter(self._raw_rho, self._raw_data, marker="x", s=self._markerSize, **kwargs)
        else:
            fig = self._plot_base(self.to1D(), xLabel=self.xLabel, yLabel=self.yLabel, edge=edge, line=True)
            fig.scatter(self._raw_rho, self._raw_data, marker="x", s=self._markerSize, **kwargs)
        return fig

    def set_wall(self, wall):
        if type(wall) is LineString:
            self.wall = wall
            self.set_plot_wall(wall)
        else:
            print("Wall must be a LineString. Wall not updated.")

    def plot2D(self, res=50):
        if self.wall:
            return self.plot_contours_with_wall(self.val, res=res)
        else:
            print("Wall not defined. ")

    def plot_fsa(self, color='red', edge=True, **kwargs):
        """Returns the flux-surface-averaged values plotted"""
        return self._plot_base(self.fsa, xLabel=self.xLabel, yLabel=self.yLabel, color=color, edge=edge, **kwargs)

    def to1D(self, l=0):
        """Returns the [:, l] array on np.array's that are broadcast originally. Normally, l=0"""
        return self.val[:, l]

    def getGrid(self):
        """Returns the 2D grid"""
        return self.val

    def Spline(self, vals, k=3, s=0):
        return UnivariateSpline(self.rho, self.val, k=k, s=s)(vals)

    def TwoDInterp(self, vals, method='cubic'):
        return griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
                        self.val.flatten(),
                        vals,
                        method=method)

    def derivative(self):
        if hasattr(self, "_derivative"):
            return self._derivative
        else:
            derivative = calc_grad(self.val, self._psi, self.R, self.Z)
            self._derivative = derivative
            return self._derivative

    def isNonZero(self):
        return not np.all(self.val == 0)

    def update(self, val):
        self.val = val
        self.fsa = OneDProfile(self._psi, calc_fsa(self.val, self.R, self.Z), self.R, self.Z)

    def update_from_1D(self, val):
        if hasattr(self, "_derivative"):
            del self._derivative
        self.val = np.broadcast_to(val, (self.rho.shape[1], len(val))).T
        self.fsa = OneDProfile(self._psi, calc_fsa(self.val, self.R, self.Z), self.R, self.Z)

    @L.setter
    def L(self, value):
        self._L = value


class TemperatureProfiles(PlotBase):

    def __init__(self, psi, R, Z, wall=None, ProfileType=TwoDProfile, *args, **kwargs):
        super(TemperatureProfiles, self).__init__()
        self._profileType = ProfileType
        if 'i' in kwargs:
            self.i = self._builder(psi, kwargs, R, Z, wall, "i", plotTitle="Ion Temperature")
            self.D = self.i
        if 'D' in kwargs:
            self.D = self._builder(psi, kwargs, R, Z, wall, "D", plotTitle="Deuterium Temperature")
        if 'T' in kwargs:
            self.T = self._builder(psi, kwargs, R, Z, wall, "T", plotTitle="Tritium Temperature")
        if 'e' in kwargs:
            self.e = self._builder(psi, kwargs, R, Z, wall, "e", plotTitle="Electron Temperature")
        if 'C' in kwargs:
            self.C = self._builder(psi, kwargs, R, Z, wall, "C", plotTitle="Carbon Temperature")
        if 'ns' in kwargs and 'nt' in kwargs:
            self.n = OneDNeutralsProfiles(self._builder(psi, kwargs, R, Z, wall, "ns", plotTitle="Slow Neutrals Temperature"),
                                          self._builder(psi, kwargs, R, Z, wall, "nt", plotTitle="Hot Neutrals Temperature"),
                                          None)
        #  Other ions

        if "W" in kwargs:
            self.W = self._builder(psi, kwargs, R, Z, wall, "W", plotTitle="Tungsten Temperature")
        if "Be" in kwargs:
            self.Be = self._builder(psi, kwargs, R, Z, wall, "Be", plotTitle="Berylium Temperature")
        if "Ne" in kwargs:
            self.Ne = self._builder(psi, kwargs, R, Z, wall, "Ne", plotTitle="Neon Temperature")
        if "Ar" in kwargs:
            self.Ar = self._builder(psi, kwargs, R, Z, wall, "Ar", plotTitle="Argon Temperature")
        if "Kr" in kwargs:
            self.Kr = self._builder(psi, kwargs, R, Z, wall, "Kr", plotTitle="Krypton Temperature")
        if "alpha" in kwargs:
            self.alpha = self._builder(psi, kwargs, R, Z, wall, "alpha", plotTitle="Alphas Temperature")

        self.set_plot_RZ(R, Z)
        self.set_plot_rho1d(psi.rho[:,0])

        # If raw data are supplied, add it to the density data
        if 'raw' in kwargs:
            raw = kwargs.get('raw')  # type: dict
            if 'i' in raw:
                self.i.kev.set_raw_data(raw.get('i'))
            if 'D' in raw:
                self.D.kev.set_raw_data(raw.get('D'))
            if 'T' in raw:
                self.T.kev.set_raw_data(raw.get('T'))
            if 'e' in raw:
                self.e.kev.set_raw_data(raw.get('e'))
            if 'C' in raw:
                self.C.kev.set_raw_data(raw.get('C'))
            if 'W' in raw:
                self.W.kev.set_raw_data(raw.get('W'))
            if 'Be' in raw:
                self.Be.kev.set_raw_data(raw.get('Be'))
            if 'Ne' in raw:
                self.Ne.kev.set_raw_data(raw.get('Ne'))
            if 'Ar' in raw:
                self.Ar.kev.set_raw_data(raw.get('Ar'))
            if 'Kr' in raw:
                self.Kr.kev.set_raw_data(raw.get('Kr'))
            if 'alpha' in raw:
                self.alpha.kev.set_raw_data(raw.get('alpha'))

    def _builder(self, psi, args, R, Z, wall, s, plotTitle=""):
        return TemperatureSplit(self._profileType(psi, args.get(s), R, Z, units="keV", wall=wall, plotTitle=plotTitle),
                                self._profileType(psi, args.get(s) * 1E3, R, Z, units="eV", wall=wall, plotTitle=plotTitle),
                                self._profileType(psi, args.get(s) * 1E3 * e, R, Z, units="J", wall=wall, plotTitle=plotTitle))


    def plot_contours_with_wall(self, obj, res=50):
        pass

    def plot_with_wall(self, obj=None):
        pass

    def plot(self, e=True, C=True, edge=True, **kwargs):
        fig = self._plot_base(self.i.kev.fsa, yLabel=r"$T[keV]$", edge=edge, title="Temperature Profiles", **kwargs)
        legend = [r"$T_i$"]
        if e:
            if self.e:
                fig.scatter(self._plot_rho1d, self.e.kev.fsa, color="blue", s=self._markerSize, **kwargs)
                legend.append(r"$T_e$")
        if C:
            if self.C:
                fig.scatter(self._plot_rho1d, self.C.kev.fsa, color="green", s=self._markerSize, **kwargs)
                legend.append(r"$T_C$")
        if len(legend) >1:
            fig.legend(legend)
        plt.show()
        return fig


class DensityProfiles(PlotBase):

    def __init__(self, psi, R, Z, wall=None, ProfileType=TwoDProfile, *args, **kwargs):
        super(DensityProfiles, self).__init__()
        self.wall = wall
        self._psi = psi
        if 'i' in kwargs:
            self.i = ProfileType(psi, kwargs.get("i"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Ion Density")
            self.D = self.i
        if "D" in kwargs:
            self.D = ProfileType(psi, kwargs.get("D"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Deuterium Density")
        if "T" in kwargs:
            self.T = ProfileType(psi, kwargs.get("T"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Tritium Density")
        if 'e' in kwargs:
            self.e = ProfileType(psi, kwargs.get("e"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Electron Density")
        if 'C' in kwargs:
            self.C = ProfileType(psi, kwargs.get("C"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Carbon Density")
        if 'ns' in kwargs and 'nt' in kwargs:
            self.n = OneDNeutralsProfiles(ProfileType(psi, kwargs.get("ns"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Ion Density"),
                                          ProfileType(psi, kwargs.get("nt"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Ion Density"),
                                          ProfileType(psi, kwargs.get("ns") + kwargs.get("nt"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Ion Density"))
        #  Other ions

        if "W" in kwargs:
            self.W = ProfileType(psi, kwargs.get("W"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Tungsten Density")
        if "Be" in kwargs:
            self.Be = ProfileType(psi, kwargs.get("Be"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Beryleum Density")
        if "Ne" in kwargs:
            self.Ne = ProfileType(psi, kwargs.get("Ne"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Neon Density")
        if "Ar" in kwargs:
            self.Ar = ProfileType(psi, kwargs.get("Ar"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Argon Density")
        if "Kr" in kwargs:
            self.Kr = ProfileType(psi, kwargs.get("Kr"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Krypton Density")
        if "alpha" in kwargs:
            self.alpha = ProfileType(psi, kwargs.get("alpha"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Alphas Density")

        # If deuterium and tritium is given, add their densities together
        # We add a try in case a deuterium density is not given period.
        try:
            self.T
            if self.D and self.T.isNonZero():
                self.D = ProfileType(psi, kwargs.get("D") + kwargs.get("T"), R, Z, units=r"$m^-3$", wall=wall, plotTitle="Ion Density")
        except:
            pass
        self.set_plot_RZ(R, Z)
        self.set_plot_rho1d(psi.rho[:,0])

        # If raw data are supplied, add it to the density data
        if 'raw' in kwargs:
            raw = kwargs.get('raw')  # type: dict
            if 'i' in raw:
                self.i.set_raw_data(raw.get('i'))
            if 'D' in raw:
                self.D.set_raw_data(raw.get('D'))
            if 'T' in raw:
                self.T.set_raw_data(raw.get('T'))
            if 'e' in raw:
                self.e.set_raw_data(raw.get('e'))
            if 'C' in raw:
                self.C.set_raw_data(raw.get('C'))
            if 'W' in raw:
                self.W.set_raw_data(raw.get('W'))
            if 'Be' in raw:
                self.Be.set_raw_data(raw.get('Be'))
            if 'Ne' in raw:
                self.Ne.set_raw_data(raw.get('Ne'))
            if 'Ar' in raw:
                self.Ar.set_raw_data(raw.get('Ar'))
            if 'Kr' in raw:
                self.Kr.set_raw_data(raw.get('Kr'))
            if 'alpha' in raw:
                self.alpha.set_raw_data(raw.get('alpha'))

    def update_neutrals(self, ns, nt):
        self.n = OneDNeutralsProfiles(TwoDProfile(self._psi, ns, self._R, self._Z, units=r"$m^-3$", wall=self.wall),
                                      TwoDProfile(self._psi, nt, self._R, self._Z, units=r"$m^-3$", wall=self.wall),
                                      TwoDProfile(self._psi, ns + nt, self._R, self._Z, units=r"$m^-3$",
                                                  wall=self.wall))


class LzTwoDProfile(TwoDProfile):
    def __init__(self, core: Core, val, slow: bool):
        super().__init__(core.psi, val, core.R, core.Z)
        self.slow = slow
        self._core = core
        self.val = TwoDProfile(core.psi, np.zeros(core.psi.rho.shape), core.R, core.Z)
        self.ddT = TwoDProfile(core.psi, np.zeros(core.psi.rho.shape), core.R, core.Z)

    def update_Lz(self, Lz, T: TemperatureProfiles, n: DensityProfiles):
        if self.slow:
            self.update(Lz(np.log10(T.n.s.kev.val),
                          np.log10(n.n.s.val / n.e.val),
                          np.log10(T.e.kev.val)))

        else:
            self.update(Lz(np.log10(T.n.t.kev.val),
                        np.log10(n.n.t / n.e.val),
                        np.log10(T.e.kev.val)))
        try:
            del self._derivative
        except AttributeError:
            pass
        return self.val

    def update_ddT(self, LzddT, T: TemperatureProfiles, n: DensityProfiles):
        if self.slow:
            self.update(LzddT(np.log10(T.n.s.kev.val),
                             np.log10(n.n.s.val / n.e.val),
                             np.log10(T.e.kev.val)))
        else:
            self.update(LzddT(np.log10(T.n.t.kev.val),
                             np.log10(n.n.t.val / n.e.val),
                             np.log10(T.e.kev.val)))
        return self.ddT


class Psi(PlotBase):
    psi = None
    psi_norm = None

    def __init__(self, pts, psi_data, sep_val, rho, a):
        super().__init__()
        self.psi_data = psi_data
        self.rho = rho
        self.a = a
        interp_fns = calc_rho2psi_interp(pts, psi_data, sep_val)
        self.rho2psi = interp_fns[0]  # type: interp1d
        self.rho2psinorm = interp_fns[1]  # type: interp1d
        self.psi2rho = interp_fns[2]  # type: interp1d
        self.psi2psinorm = interp_fns[3]  # type: interp1d
        self.psinorm2rho = interp_fns[4]  # type: interp1d
        self.psinorm2psi = interp_fns[5]  # type: interp1d
        self.psi = self.rho2psi(self.rho)
        self.psi_norm = self.rho2psinorm(self.rho)


    def plot_exp_psi(self, res=50):
        try:
            ax = self._plot_with_wall()
        except:
            return
        try:
            ax.contour(self.R, self.Z, self.psi.psi, res)
            return ax
        except NameError:
            print("Psi not defined")
            pass


class ImpurityProfiles(PlotBase, BaseMath):
    def __init__(self, core: Core, *args, **kwargs):
        super().__init__()
        self.wall = core.wall_line
        self._psi = core.psi
        LZSplit = namedtuple("LZData", "s t")
        self.C = LZSplit(LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=True),
                         LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=False))

        self.Be = LZSplit(LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=True),
                         LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=False))

        self.W = LZSplit(LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=True),
                         LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=False))

        self.Ne = LZSplit(LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=True),
                         LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=False))

        self.Ar = LZSplit(LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=True),
                         LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=False))

        self.Kr = LZSplit(LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=True),
                         LzTwoDProfile(core, np.zeros(core.psi.rho.shape), slow=False))


    def update_C(self, n: DensityProfiles, T: TemperatureProfiles, Lz, ddT):
        self.C.s.update_Lz(Lz, T, n)
        self.C.s.update_ddT(ddT, T, n)

    def update_Be(self, n: DensityProfiles, T: TemperatureProfiles, Lz, ddT):
        self.Be.s.update_Lz(Lz, T, n)
        self.Be.s.update_ddT(ddT, T, n)

    def update_W(self, n: DensityProfiles, T: TemperatureProfiles, Lz, ddT):
        self.W.s.update_Lz(Lz, T, n)
        self.W.s.update_ddT(ddT, T, n)

    def update_Ne(self, n: DensityProfiles, T: TemperatureProfiles, Lz, ddT):
        self.Ne.s.update_Lz(Lz, T, n)
        self.Ne.s.update_ddT(ddT, T, n)

    def update_Ar(self, n: DensityProfiles, T: TemperatureProfiles, Lz, ddT):
        self.Ar.s.update_Lz(Lz, T, n)
        self.Ar.s.update_ddT(ddT, T, n)

    def update_Kr(self, n: DensityProfiles, T: TemperatureProfiles, Lz, ddT):
        self.Kr.s.update_Lz(Lz, T, n)
        self.Kr.s.update_ddT(ddT, T, n)


    def plot_contours_with_wall(self, obj, res=50):
        pass

    def plot_with_wall(self, obj=None):
        pass


    def plot(self, e=True, C=True, edge=True, **kwargs):
        fig = self._plot_base(self.i.fsa, yLabel="r$n[#/m^{-3}]$", edge=edge, **kwargs)
        legend = [r"$n_i$"]
        if e:
            if self.e:
                fig.scatter(self._plot_rho1d, self.e.fsa, color="blue", s=self._markerSize)
                legend.append(r"$n_e$")
        if C:
            if self.C:
                fig.scatter(self._plot_rho1d, self.C.fsa, color="green", s=self._markerSize)
                legend.append(r"$n_C$")
        if len(legend) >1:
            fig.legend(legend)
        plt.show()
        return fig


class PressureProfiles(PlotBase):

    def __init__(self, psi, R, Z, wall=None, ProfileType=TwoDProfile, *args, **kwargs):
        super(PressureProfiles, self).__init__()

        if 'i' in kwargs:
            self.i = ProfileType(psi, kwargs.get("i"), R, Z, units="Pa", wall=wall, plotTitle="Ion Pressure")
            self.D = self.i
        if 'e' in kwargs:
            self.e = ProfileType(psi, kwargs.get("e"), R, Z, units="Pa", wall=wall, plotTitle="Electron Pressure")
        if 'C' in kwargs:
            self.C = ProfileType(psi, kwargs.get("C"), R, Z, units="Pa", wall=wall, plotTitle="Carbon Pressure")

        self.set_plot_RZ(R, Z)
        self.set_plot_rho1d(psi.rho[:,0])

    def plot_contours_with_wall(self, obj, res=50):
        pass

    def plot_with_wall(self, obj=None):
        pass

    def plot1D(self, e=True, C=True, edge=True, **kwargs):
        fig = self._plot_base(self.i.fsa, yLabel="r$P[Pa]$", edge=edge, **kwargs)
        legend = [r"$p_i$"]
        if e:
            if self.e:
                fig.scatter(self._plot_rho1d, self.e.fsa, color="blue", s=self._markerSize)
                legend.append(r"$p_e$")
        if C:
            if self.C:
                fig.scatter(self._plot_rho1d, self.C.fsa, color="green", s=self._markerSize)
                legend.append(r"$p_C$")
        if len(legend) >1:
            fig.legend(legend)
        plt.show()
        return fig


class NeutralsProfiles(PlotBase, BaseMath):

    def __init__(self,
                 s : Union[OneDProfile, TwoDProfile, TemperatureSplit],
                 t : Union[OneDProfile, TwoDProfile, TemperatureSplit],
                 tot : Union[OneDProfile, TwoDProfile, TemperatureSplit, None]
                 ):
        """
        :param s: The slow neutron profile
        :param t: The thermal neutron profile
        :param tot: the total neutron profile
        :type s: Union[OneDProfile, TwoDProfile]
        :type t: Union[OneDProfile, TwoDProfile]
        :type tot: Union[OneDProfile, TwoDProfile]

        """
        super(NeutralsProfiles, self).__init__()
        self.s = s
        self.t = t
        self.tot = tot

    def set_plot_wall(self, wall):
        self.s.set_plot_wall(wall)
        self.t.set_plot_wall(wall)
        self.tot.set_plot_wall(wall)
        pass

    def set_plot_rho1d(self, rho1d):
        self.s.set_plot_rho1d(rho1d)
        self.t.set_plot_rho1d(rho1d)
        self.tot.set_plot_rho1d(rho1d)


class OneDNeutralsProfiles(NeutralsProfiles):

    def __init__(self, s, t, tot):
        super(OneDNeutralsProfiles, self).__init__(s, t, tot)

    # def plot_contours_with_wall(self, obj, res=50):
    #     pass

    def plot_add_scatter(self, fig, val, color="blue"):
        pass

class VectorialProfiles:
    tor_D = None
    tor_C = None
    tor_e = None
    pol_D = None
    pol_C = None
    pol_e = None
    def __init__(self, psi, R, Z, wall=None, ProfileType=TwoDProfile, *args, **kwargs):
        """
        The container for the velocity profiles. The velocities should be given by kwarg as follows:

        tor_D = Toroidal deuterium profile
        tor_C = Toroidal carbon profile
        tor_e = Toroidal electron profile
        pol_D = Poloidal deuterium profile
        pol_C = Poloidal carbon profile
        pol_e = Poloidal electron profile


        @TODO: Expand to other particles?

        :param psi:
        :param R:
        :param Z:
        :param wall:
        :param ProfileType:
        :param args:
        :param kwargs:
        """
        super().__init__()
        self._psi = psi
        self._R = R
        self._Z = Z
        self._wall = wall
        self._profileType = ProfileType
        self._set_vals(**kwargs)


    def update_D(self, pol=False, tor=False):
        if not np.any(pol) and not np.any(tor):
            return self
        if np.any(pol) and np.any(tor):
            self._set_vals(tor_D=tor, pol_D=pol)
            return self
        if not np.any(pol) and np.any(tor):
            self._set_vals(tor_D=tor, pol_D=self.i.pol.val)
            return self
        if np.any(pol) and not np.any(tor):
            self._set_vals(tor_D=self.i.tor.val, pol_D=pol)
            return self

    def update_e(self, pol=False, tor=False):
        if not np.any(pol) and not np.any(tor):
            return self
        if np.any(pol) and np.any(tor):
            self._set_vals(tor_e=tor, pol_e=pol)
            return self
        if not np.any(pol) and np.any(tor):
            self._set_vals(tor_e=tor, pol_e=self.e.pol.val)
            return self
        if np.any(pol) and not np.any(tor):
            self._set_vals(tor_e=self.e.tor.val, pol_e=pol)
            return self

    def update_C(self, pol=False, tor=False):
        if not np.any(pol) and not np.any(tor):
            return self
        if np.any(pol) and np.any(tor):
            self._set_vals(tor_C=tor, pol_C=pol)
            return self
        if not np.any(pol) and np.any(tor):
            self._set_vals(tor_C=tor, pol_C=self.C.pol.val)
            return self
        if np.any(pol) and not np.any(tor):
            self._set_vals(tor_C=self.C.tor.val, pol_C=pol)
            return self

    def _set_vals(self, **kwargs):
        PoloidalToroidalSplit = namedtuple("PoloidalToroidalSplit", "pol tor tot")

        if 'tor_D' in kwargs or 'pol_D' in kwargs:
            if 'tor_D' in kwargs and 'pol_D' in kwargs:
                self.i = PoloidalToroidalSplit(self._profileType(self._psi, kwargs.get('pol_D'), self._R, self._Z,
                                                                 units=r"m/s", wall=self._wall),
                                               self._profileType(self._psi, kwargs.get('tor_D'), self._R, self._Z,
                                                                 units=r"m/s", wall=self._wall),
                                               self._profileType(self._psi,
                                                                 np.sqrt(kwargs.get('tor_D')**2 + kwargs.get('pol_D')**2),
                                                                 self._R, self._Z, units=r"m/s", wall=self._wall))
            self.D = self.i

        if 'tor_C' in kwargs or 'pol_C' in kwargs:
            self.C = PoloidalToroidalSplit(self._profileType(self._psi, kwargs.get('pol_C'), self._R, self._Z,
                                                             units=r"m/s", wall=self._wall),
                                           self._profileType(self._psi, kwargs.get('tor_C'), self._R, self._Z,
                                                             units=r"m/s", wall=self._wall),
                                           self._profileType(self._psi, np.sqrt(kwargs.get('tor_C')**2 + kwargs.get('pol_C')**2),
                                                             self._R, self._Z, units=r"m/s", wall=self._wall))

        if 'tor_e' in kwargs or 'pol_e' in kwargs:
            self.e = PoloidalToroidalSplit(self._profileType(self._psi, kwargs.get('pol_e'), self._R, self._Z,
                                                             units=r"m/s", wall=self._wall),
                                           self._profileType(self._psi, kwargs.get('tor_e'), self._R, self._Z,
                                                             units=r"m/s", wall=self._wall),
                                           self._profileType(self._psi, np.sqrt(kwargs.get('tor_e')**2 + kwargs.get('pol_e')**2),
                                                             self._R, self._Z, units=r"m/s", wall=self._wall))

class VectorialBase(PlotBase, BaseMath):

    def __init__(self, tor, pol, psi, R, Z, wall, profileType=TwoDProfile, units=""):
        super().__init__()
        self._profileType = profileType
        self._R = R
        self._Z = Z
        self._wall = wall
        self._psi = psi

        self.pol = profileType(psi, pol, R, Z, units=units, wall=wall)
        self.tor = profileType(psi, tor, R, Z, units=units, wall=wall)
        self.tot = profileType(psi, np.sqrt(pol**2 + tor**2), R, Z, units=units, wall=wall)

class Flux(PlotBase, BaseMath):

    def __init__(self, core: Core, label="", **kwargs):
        super().__init__()
        DiffIntSplit = namedtuple("DiffIntSplit", "diff int")
        self.label = label
        self.set_plot_rho1d(core.rho[:, 0])
        D_int = kwargs.get("D_int")
        D_diff = kwargs.get("D_diff")
        C_int = kwargs.get("C_int")
        C_diff = kwargs.get("C_diff")
        e_int = kwargs.get("e_int")
        e_diff = kwargs.get("e_diff")
        self._core = core

        self.D = DiffIntSplit(OneDProfile(core.psi, D_diff, core.R, core.Z),
                              OneDProfile(core.psi, D_int, core.R, core.Z))
        self.C = DiffIntSplit(OneDProfile(core.psi, C_diff, core.R, core.Z),
                              OneDProfile(core.psi, C_int, core.R, core.Z))
        self.e = DiffIntSplit(OneDProfile(core.psi, e_diff, core.R, core.Z),
                              OneDProfile(core.psi, e_int, core.R, core.Z))

    def plot_D(self, **kwargs):
        fig = self._plot_base(self.D.diff.val, yLabel=self.label, edge=True, title="", **kwargs)
        fig.scatter(self._plot_rho1d, self.D.int.val, color="green")
        legend = [r"integral",
                  r"differential"]
        fig.legend(legend)
        return fig

    def plot_C(self, **kwargs):
        fig = self._plot_base(self.C.diff.val, yLabel=r"$\Gamma$", edge=True, title="", **kwargs)
        fig.scatter(self._plot_rho1d, self.C.int.val, color="green")
        legend = [r"$\Gamma^{int}$",
                  r"$\Gamma^{diff}$"]
        fig.legend(legend)
        return fig

    def plot_e(self, **kwargs):
        fig = self._plot_base(self.e.diff.val, yLabel=self.label, edge=True, title="", **kwargs)
        fig.scatter(self._plot_rho1d, self.e.int.val, color="green")
        legend = [r"integral",
                  r"differential"]
        fig.legend(legend)
        return fig





