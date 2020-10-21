#/usr/bin/python

from matplotlib import pyplot as plt
from numpy import genfromtxt, savetxt
import numpy as np
import os
from scipy.interpolate import LSQUnivariateSpline

class dataSmoother:

    def __init__(self, file, knots=[]):
        self.filename = file
        with open(file, "r") as f:
            inData = genfromtxt(f, float, comments='#').T
            self.rho = inData[0]
            self.orig_data = inData[1]

        self.spline = LSQUnivariateSpline(self.rho, self.orig_data, knots, k=5)

    def plot_orig(self):
        fig = self._plot_base()
        fig.scatter(self.rho, self.orig_data, color='red')
        plt.show()
        return fig

    def plot(self):
        fig = self._plot_base()

        fig.scatter(self.rho, self.orig_data, color='red')
        fig.plot(self.rho, self.spline(self.rho), color='black')
        fig.scatter(self.rho, self.spline(self.rho), color='green', marker="x")

        fig.legend(["Spline", "Original", "Cleaned"])

        plt.show()
        return fig

    def save(self):
        answer = ""
        while answer not in ["y", "n"]:
            answer = input("OK to continue [Y/N]? ").lower()
        print("Backing up old file")
        os.rename(self.filename, self.filename + '.bak')
        contents = np.array([self.rho, self.spline(self.rho)]).T
        savetxt(self.filename + "_clean", contents, delimiter='    ')
    def get_spline(self):
        return self.spline

    def set_knots(self, knots):
        self.spline = LSQUnivariateSpline(self.rho, self.orig_data, knots, k=5)
        print("knots set")

    def _plot_base(self):
        plot = plt.figure()
        fig = plot.add_subplot(111)
        fig.set_xlabel(r'$\rho', fontsize=30)
        fig.set_ylabel(r'Data', fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        return fig

    def plot_derivative(self):
        fig = self._plot_base()
        fig.scatter(self.rho, self.spline.derivative()(self.rho), color="blue")
        return fig
