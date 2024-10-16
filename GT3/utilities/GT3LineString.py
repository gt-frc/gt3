#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString, Point, MultiPoint
from GT3.utilities.PlotBase import PlotBase

MARKERSIZE = 10


class GT3LineString(LineString, PlotBase):

    def __init__(self, coordinates=None, wall=None):
        super().__init__(coordinates)
        self._markerSize = MARKERSIZE
        self._markers = False
        self._defColor = 'red'
        self.ls = LineString(coordinates)
        if wall:
            self.set_plot_wall(wall)

    def plot_add_scatter(self, fig, val, color=None):
        return

    def plot_contours_with_wall(self, obj, res=50):
        return

    def plot_ls_with_wall(self):
        self.plot_with_wall(self.ls)

    def plot_ls_without_wall(self):
        self.plot_without_wall(self.ls)
