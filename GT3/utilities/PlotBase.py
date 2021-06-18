#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString, Point, MultiPoint
from matplotlib import Path

MARKERSIZE = 10

class PlotBase:

    def __init__(self):
        super(PlotBase, self).__init__()
        self._markerSize = MARKERSIZE
        self._markers = False

    def set_marker_size(self, s):
        self._markerSize = s

    def set_plot_rho1d(self, rho1d):
        self._plot_rho1d = rho1d
        return self

    def set_plot_RZ(self, R, Z):
        self._R = R
        self._Z = Z

    def set_plot_wall(self, wall):
        self._wall_line = wall

    def _plot_base(self, val, xLabel=r'$\rho$', yLabel="Value", title="Title", color='red', edge=False, show=True,
                   line=False, **kwargs):
        plot = plt.figure()
        fig = plot.add_subplot(111)
        fig.set_xlabel(xLabel, fontsize=30)
        fig.set_ylabel(yLabel, fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        fig.set_title(title)
        if kwargs.get("logPlot"):
            fig.set_yscale("log")
            val = np.abs(val)
        if edge:
            fig.set_xlim(0.85, 1.0)
        if line:
            fig.plot(self._plot_rho1d, val, color=color)
        else:
            fig.scatter(self._plot_rho1d, val, color=color, s=self._markerSize)
        if show:
            plt.show()
        return fig

    def _plot_with_wall(self):
        """
        Generates a Matplotlib plot with the wall pre-plotted.

        :return: An Axis object with the wall line plotted
        """

        if not hasattr(self, "_wall_line") and not hasattr(self, "wall_line"):
            print("Wall Linestring has not been instantiated yet.")
            return

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

    def _shapely_obj_plot_hanlder(self, obj, ax):
        if isinstance(obj, Point):
            ax.scatter(obj.x, obj.y, color='red', marker='o', s=self._markerSize)
            return ax
        if isinstance(obj, MultiPoint):
            for p in obj:
                ax.scatter(p.x, p.y, color='red', marker='o', s=self._markerSize)
            return ax
        if isinstance(obj, Path):
            ls = LineString(obj.vertices)
            ax.plot(*ls.xy, color='red', marker='o', s=self._markerSize)
            return ax
        if isinstance(obj, LineString):
            ax.plot(*obj.xy, color='red', marker='o', s=self._markerSize)
            return ax

    def _unknown_data_plot_helper(self, obj, ax):
        # Is it a shapely or similar object?
        if isinstance(obj, (Point, MultiPoint, Path, LineString)):
            return self._shapely_obj_plot_hanlder(obj, ax)

        # Does it have the vertices property to convert to a LineString?
        try:
            ls = LineString(obj.vertices)
            return self._shapely_obj_plot_hanlder(ls, ax)
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
            return self._shapely_obj_plot_hanlder(p, ax)
        except:
            pass

        # Can I turn it into a LineString directly?
        try:
            ls = LineString(obj)
            return self._shapely_obj_plot_hanlder(ls, ax)
        except:
            pass

        # No dice. Raising error.
        raise

    def plot_with_wall(self, obj=None):
        if not hasattr(self, "_wall_line"):
            print("Wall Linestring has not been instantiated yet.")
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

        print("Could not plot given data")

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

    def plot_add_scatter(self, fig, val, color="blue"):
        fig.scatter(self._plot_rho1d, val, color=color, s=self._markerSize)
        return fig