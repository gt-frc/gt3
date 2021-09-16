#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp2d, Rbf
from shapely.geometry import LineString, Point, MultiPoint
from matplotlib import Path
from warnings import warn
from GT3.utilities.GT3Figure import GT3FigureSinglePlot

from GT3.utilities import MARKERSIZE_SMALL, MARKERSIZE_MEDIUM, MARKERSIZE_LARGE, PLOTCOLORS, PLOTMARKERS

class PlotBase:

    def __init__(self):
        super(PlotBase, self).__init__()
        self._markerSize = MARKERSIZE_SMALL
        self._markers = False
        self._defColor = PLOTCOLORS[0]

    def set_marker_size(self, s):
        self._markerSize = s

    def set_default_color(self, c):
        self._defColor = c

    def set_plot_rho1d(self, rho1d):
        self._plot_rho1d = rho1d
        return self

    def set_plot_RZ(self, R, Z):
        self._R = R
        self._Z = Z

    def set_plot_wall(self, wall):
        self._wall_line = wall

    def _plot_base(self, val, title="Title", color=None, edge=False, show=True,
                   line=False, **kwargs):

        if not color:
            color = self._defColor

        fig = GT3FigureSinglePlot()
        #plot = plt.figure()
        #fig = plot.add_subplot(111)


        if kwargs.get("multiplier"):
            if isinstance(kwargs.get("multiplier"), (float, int)):
                val = val * kwargs.get("multiplier")
            else:
                warn("Multiplier kwarg is not a float or int", UserWarning)
        if kwargs.get("logPlot"):
            fig.toggle_semilog()
        if line:
            fig.add_line(self._plot_rho1d, val)
        else:
            if kwargs.get("legend"):
                legend = kwargs.get("legend")
            else:
                legend = ""
            fig.add_scatter(self._plot_rho1d, val, legend=legend)
        fig.set_marker_size(self._markerSize)
        if kwargs.get("xLabel"):
            fig.set_xLabel(kwargs.get("xLabel"))
            fig.set_xLabel_fontsize(30)
        else:
            fig.set_xLabel(r'$\rho$')
            fig.set_xLabel_fontsize(30)

        if kwargs.get("yLabel"):
            fig.set_yLabel(kwargs.get("yLabel"))
            fig.set_yLabel_fontsize(30)
        else:
            fig.set_yLabel_fontsize(30)

        fig.set_xticks_fontsize(20)
        fig.set_yticks_fontsize(20)
        fig.set_title(title)
        if edge:
            fig.set_xlim(0.85, 1.0)
        if show:
            plt.show()
        return fig

    def _plot_without_wall(self):

        """
        Generates a Matplotlib plot without the wall pre-plotted.

        :return: An Axis object with the wall line plotted
        """

        fig_width = 6.0
        # Check to see if self.Z/R have been defind yet to generate a figure height, as they aren't immediately
        # calculated
        try:
            fig_height = (np.amax(self._Z) - np.amin(self._Z)) / (np.amax(self._R) - np.amin(self._R)) * fig_width
        except:
            fig_height = 9.0

        fig = plt.figure(figsize=(0.975 * fig_width, fig_height))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.axis('equal')
        return ax1

    def _plot_with_wall(self):
        """
        Generates a Matplotlib plot with the wall pre-plotted.

        :return: An Axis object with the wall line plotted
        """

        if not hasattr(self, "_wall_line") and not hasattr(self, "wall_line"):
            print("Wall Linestring has not been instantiated yet.")
            return self._plot_without_wall()

        if not hasattr(self, "_Z") and not hasattr(self, "Z"):
            print("Z coordinates have not been instantiated yet.")
            return

        if not hasattr(self, "_R") and not hasattr(self, "R"):
            print("R coordinates have not been instantiated yet.")
            return

        if not hasattr(self, "_wall_line"):
            self._wall_line = self.wall_line

        if not hasattr(self, "_Z"):
            self._Z = self.Z

        if not hasattr(self, "_R"):
            self._R = self.R

        fig_width = 6.0
        # Check to see if self.Z/R have been defind yet to generate a figure height, as they aren't immediately
        # calculated
        try:
            fig_height = (np.amax(self._Z) - np.amin(self._Z)) / (np.amax(self._R) - np.amin(self._R)) * fig_width
        except:
            fig_height = 9.0

        fig = plt.figure(figsize=(0.975*fig_width, fig_height))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.axis('equal')
        ax1.plot(np.asarray(self._wall_line.xy)[0], np.asarray(self._wall_line.xy)[1], color='black', lw=1.5)
        return ax1

    def _shapely_obj_plot_handler(self, obj, ax, color=None):
        if not color:
            color = self._defColor
        if isinstance(obj, Point):
            ax.scatter(obj.x, obj.y, color=color, marker='o', s=self._markerSize)
            return ax
        if isinstance(obj, MultiPoint):
            for p in obj:
                ax.scatter(p.x, p.y, color=color, marker='o', s=self._markerSize)
            return ax
        if isinstance(obj, Path):
            ls = LineString(obj.vertices)
            ax.plot(*ls.xy, color=color, marker='o', s=self._markerSize)
            return ax
        if isinstance(obj, LineString):
            ax.plot(*obj.xy, color=color, marker='o', markersize=self._markerSize)
            return ax

    def _unknown_data_plot_helper(self, obj, ax):
        # Is it a shapely or similar object?
        if isinstance(obj, (Point, MultiPoint, Path, LineString)):
            return self._shapely_obj_plot_handler(obj, ax)

        # Does it have the vertices property to convert to a LineString?
        try:
            ls = LineString(obj.vertices)
            return self._shapely_obj_plot_handler(ls, ax)
        except:
            pass

        # Does it have an xy property?

        try:
            ax.plot(*obj.xy, s=self._markerSize)
            return ax
        except:
            pass

        # Can I turn it into a point?
        try:
            p = Point(obj)
            return self._shapely_obj_plot_handler(p, ax)
        except:
            pass

        # Can I turn it into a LineString directly?
        try:
            ls = LineString(obj)
            return self._shapely_obj_plot_handler(ls, ax)
        except:
            pass

        # No dice. Raising error.
        raise

    def plot_with_wall(self, obj=None):
        if not hasattr(self, "_wall_line"):
            print("Wall Linestring has not been instantiated yet. Plotting without Wall")
            self.plot_with_wall(obj)
            return
        ax = self._plot_with_wall()

        try:
            return self._unknown_data_plot_helper(obj, ax)
        except:
            pass

        # Is it an iterable?
        try:
            obj.__iter__
            for p in obj:
                try:
                    ax = self._unknown_data_plot_helper(p, ax)
                except:
                    raise
            return ax
        except:
            pass
    def plot_without_wall(self, obj=None):
        ax = self._plot_without_wall()
        try:
            return self._unknown_data_plot_helper(obj, ax)
        except:
            pass

        # Is it an iterable?
        try:
            obj.__iter__
            for p in obj:
                try:
                    ax = self._unknown_data_plot_helper(p, ax)
                except:
                    raise
            return ax
        except:
            pass
    def plot_contours_with_wall(self, obj, res=50):
        if not hasattr(self, "_wall_line"):
            print("Wall Linestring has not been instantiated yet.")
            return
        ax = self._plot_with_wall()
        try:
            ax.contour(self.R, self.Z, obj, levels=res)
            return ax
        except:
            print("Could not plot contours")
            return

    def plot_add_scatter(self, fig, val, color=None):
        if not color:
            color = self._defColor
        fig.scatter(self._plot_rho1d, val, color=color, s=self._markerSize)
        return fig

class PlotBaseWithHeatMap(PlotBase):
    def __init__(self):
        super(PlotBaseWithHeatMap, self).__init__()

    def _plot_HM(self, x, y, z, aspect=1, cmap=plt.cm.rainbow, *args, **kwargs):
        # x, y, z = self._remove_duplicates_HM(x,y,z)
        xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
        xi, yi = np.meshgrid(xi, yi)
        # interp = interp2d(x, y, z)
        if kwargs.get("no_RBF"):
            interp = interp2d(x, y, z)
        else:
            try:
                interp = Rbf(x, y, z)
            except np.linalg.LinAlgError:
                interp = interp2d(x, y, z)

        zi = interp(xi[0], yi[:, 0])
        if kwargs.get("logScale"):
            zi[zi == 0] = 0.00001
            zi = np.log10(z)
        fig, ax = plt.subplots(figsize=(6, 6))
        hm = ax.imshow(zi, interpolation='nearest', cmap=cmap, extent=[x.min(), x.max(), y.max(), y.min()])
        ax.set_aspect(aspect)
        return fig, hm

    def plot_HM(self, **kwargs):
        self._plot_HM_with_wall(self.R, self.Z, self.val, **kwargs)

    def _plot_HM_with_wall(self, R, Z, val, *args, **kwargs):

        fig, ax = self._plot_HM(R, Z, val, **kwargs)
        fig.colorbar(ax)
        ax.axes.plot(np.asarray(self._wall_line.xy).T[:, 0], np.asarray(self._wall_line.xy).T[:, 1], color='black', lw=1.5)
        return ax