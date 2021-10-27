#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from GT3.utilities import MARKERSIZE_SMALL, MARKERSIZE_MEDIUM, MARKERSIZE_LARGE, PLOTCOLORS, PLOTMARKERS

class GT3FigureSinglePlot:

    def __init__(self, *args, **kwargs):
        self._markerSize = MARKERSIZE_SMALL
        self._markers = False
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.legend = []
        self._sciNotation = False
        self._dpi = 450
        self.filename = "gt3figure.png"
        self._legendLoc = 0
        self._xLabelFontSize = 30
        self._yLabelFontSize = 30
        self._xMin = 0.
        self._xMax = 1.
        self._yMin = 0.
        self._yMax = 1.
        self._legendFontSize = 10
        self._legendScale = 1.0
        self._plotSemilog = False
        self._showTitle = False
        self._maskZeros = False
        self.imgSaveDir = None

        if kwargs.get("title"):
            self.title = kwargs.get("title")
            self.ax.set_title(self.title)
        else:
            self.title = ""

        if kwargs.get("yLabel"):
            self.yLabel = kwargs.get("yLabel")
            self.ax.set_ylabel(self.yLabel)
        else:
            self.yLabel = ""

        if kwargs.get("xLabel"):
            self.xLabel = kwargs.get("xLabel")
            self.ax.set_xlabel(self.xLabel)
        else:
            self.xLabel = ""

        if kwargs.get("markers"):
            self._attachMarkers = True
        else:
            self._attachMarkers = False


        self.showLegend = False
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        self.scatter_data = []
        self.line_data = []


    def set_title(self, title):
        self.title = title
        self._replot()
        return self

    def set_xLabel(self, xLabel):
        self.xLabel = xLabel
        self.ax.set_xlabel(xLabel, fontsize=self._xLabelFontSize)
        return self

    def get_xLabel(self):
        return self.ax.get_xlabel(), self._xLabelFontSize

    def set_xLabel_fontsize(self, s):
        self._xLabelFontSize = s
        self.ax.set_xlabel(self.get_xLabel()[0], fontsize=self._xLabelFontSize)
        return self

    def set_yLabel(self, yLabel):
        self.yLabel = yLabel
        self.ax.set_ylabel(yLabel, fontsize=self._yLabelFontSize)
        return self

    def get_yLabel(self):
        return self.ax.get_ylabel(), self._yLabelFontSize

    def set_yLabel_fontsize(self, s):
        self._yLabelFontSize = s
        self.ax.set_ylabel(self.get_yLabel()[0], fontsize=self._yLabelFontSize)
        return self

    def set_xlim(self, xmin, xmax, replot=True):
        self._xMin = xmin
        self._xMax = xmax
        if replot:
            self._replot()
        return self

    def set_ylim(self, ymin, ymax, replot=True):
        self._yMin = ymin
        self._yMax = ymax
        if replot:
            self._replot()
        return self

    def add_line(self, x, y, *args, **kwargs):
        if kwargs.get("legend"):
            self.legend.append(kwargs.get("legend"))
            self.ax.legend(self.legend)
            self.showLegend = True
        else:
            self.legend.append(None)
        self.line_data.append([x, y])
        self._replot(axisStick=True)
        self._xMin, self._xMax = self.ax.get_xlim()
        self._yMin, self._yMax = self.ax.get_ylim()
        return self

    def add_scatter(self, x, y, *args, **kwargs):

        if kwargs.get("legend"):
            self.legend.append(kwargs.get("legend"))
            self.ax.legend(self.legend)
            self.showLegend = True
        else:
            self.legend.append(None)

        #self.ax.scatter(x, y, color=color, s=self._markerSize)
        self.scatter_data.append([x, y])

        self._replot(axisStick=True)
        #self._xMin, self._xMax = self.ax.get_xlim()
        self._yMin, self._yMax = self.ax.get_ylim()
        return self

    def set_marker_size(self, s):

        self._markerSize = s
        self._replot()
        return self

    def set_number_xticks(self, num):
        self.ax.locator_params(axis="x", nbins=num)
        return self

    def set_number_yticks(self, num):
        self.ax.locator_params(axis="y", nbins=num)
        return self

    def set_xticks_fontsize(self, size):
        self.ax.tick_params(axis='x', labelsize=size)
        return self

    def set_yticks_fontsize(self, size):
        self.ax.tick_params(axis='y', labelsize=size)
        return self

    def set_legend_fontsize(self, size):
        self._legendFontSize = size
        self._replot()
        return self

    def set_legend_scale(self, scale):
        self._legendScale = scale
        self._replot()
        return self

    def toggle_legend(self):
        self.showLegend = not self.showLegend
        self._replot()

    def toggle_markers(self):
        self._markers = not self._markers
        self._replot()

    def toggle_sci_notation(self):
        self._sciNotation = not self._sciNotation
        self._replot()
        return self

    def toggle_semilog(self):
        self._plotSemilog = not self._plotSemilog
        self._replot()
        return self

    def toggle_title(self):
        self._showTitle = not self._showTitle
        self._replot()
        return self

    def toggle_mask_zeros(self):

        self._maskZeros = not self._maskZeros
        self._replot()
        return self

    def thesis_figure(self):

        self.set_xlim(0.85, 1.0)
        self.toggle_markers()
        self.set_xLabel(r"Normalized minor radius, $\rho$")
        self.set_legend_scale(1.5)
        self.set_legend_fontsize(20)
        self.set_yticks_fontsize(26)
        self.set_xticks_fontsize(26)
        self.set_marker_size(80)
        self.set_number_xticks(4)
        self.set_number_yticks(4)
        self.fig.tight_layout()
        return self

    def set_fig_dpi(self, d=450):
        self._dpi = d

    def set_save_name(self, s):
        self.filename = s

    def set_save_dir(self, s):
        self.imgSaveDir = s

    def save_gt3_fig(self):
        if self.imgSaveDir:
            import os
            self.fig.savefig(os.path.join(self.imgSaveDir, self.filename), dpi=self._dpi, pad_inches=0)
        else:
            self.fig.savefig(self.filename, dpi=self._dpi, pad_inches=0)
        return self

    def _replot(self, *args, **kwargs):
        self.ax.cla()
        if self._plotSemilog:
            self.ax.set_yscale("log")
        for n, data in enumerate(self.scatter_data):
            x = data[0]
            y = data[1]
            if self._plotSemilog:
                pass
                #y = np.log10(y)
            if self._maskZeros:
                try:
                    y[np.abs(y) < 1E-10] = np.nan
                except TypeError:
                    try:
                        temp_y = y.val
                        temp_y[np.abs(temp_y) < 1E-10] = np.nan
                        y = temp_y
                    except:
                        raise
            if self._markers:
                self.ax.scatter(x, y, color=PLOTCOLORS[n], s=self._markerSize, marker=PLOTMARKERS[n])
            else:
                self.ax.scatter(x, y, color=PLOTCOLORS[n], s=self._markerSize)

        for n, data in enumerate(self.line_data):
            x = data[0]
            y = data[1]
            if self._plotSemilog:
                #y = np.log10(y)
                pass
            if self._markers:
                self.ax.plot(x, y, color=PLOTCOLORS[n])
            else:
                self.ax.plot(x, y, color=PLOTCOLORS[n])

        if self.showLegend:
            self.ax.legend(self.legend, prop={'size': self._legendFontSize},
                           markerscale=self._legendScale,
                           loc=self._legendLoc)
        if not kwargs.get("axisStick"):
            self.ax.set_xlim(self._xMin, self._xMax)
            self.ax.set_ylim(self._yMin, self._yMax)
        if self._showTitle:
            self.ax.set_title(self.title)
        else:
            self.ax.set_title("")
        self.ax.set_xlabel(self.xLabel)
        self.ax.set_ylabel(self.yLabel)
        if self._sciNotation:
            self.ax = format_exponent(self.ax)
        return self

    def set_legend_location(self, loc):
        self._legendLoc = loc
        self._replot()

class GT3FigureSBSPlot:
    @staticmethod
    def _dict_gen(keys, val):
        result = {}
        for key in keys:
            result[key] = val
        return result

    def __init__(self, keys=None, *args, **kwargs):
        if keys is None:
            self._keys = ['L', 'R']
        else:
            if type(keys) == list:
                self._numPlots = len(keys)
                self._keys = keys
            elif type(keys) == int:
                self._numPlots = keys
                self._keys = range(0, len(keys))
            else:
                raise ValueError("keys needs to be an integer or list")

        self.fig = plt.figure()
        self.fig.tight_layout(pad=3.0)
        self.ax = {}
        for n, key in enumerate(self._keys):
            self.ax[key] = self.fig.add_subplot(1, self._numPlots, n+1)

        self._markerSize = self._dict_gen(self._keys, MARKERSIZE_LARGE)
        self._markers = self._dict_gen(self._keys, False)
        self.legend = self._dict_gen(self._keys, [])
        self._legendLoc = self._dict_gen(self._keys, 0)
        self._sciNotation = self._dict_gen(self._keys, False)
        self._xLabelFontSize = self._dict_gen(self._keys, 30)
        self._yLabelFontSize = self._dict_gen(self._keys, 30)
        self._xtickNum = self._dict_gen(self._keys, 4)
        self._ytickNum = self._dict_gen(self._keys, 4)
        self._xMin = self._dict_gen(self._keys, 0.)
        self._xMax = self._dict_gen(self._keys, 1.)
        self._yMin = self._dict_gen(self._keys, 0.)
        self._yMax = self._dict_gen(self._keys, 1.)
        self._legendFontSize = self._dict_gen(self._keys, 10)
        self._legendScale = self._dict_gen(self._keys, 1.0)
        self._plotSemilog = self._dict_gen(self._keys, False)
        self._showTitle = self._dict_gen(self._keys, False)
        self.title = self._dict_gen(self._keys, "")
        self.xLabel = self._dict_gen(self._keys, "")
        self.yLabel = self._dict_gen(self._keys, "")
        self.fig_title = ""
        self._maskZeros = self._dict_gen(self._keys, False)

        self.showLegend = self._dict_gen(self._keys, False)
        for key in self._keys:
            self.ax[key].tick_params(axis="x", labelsize=16)
            self.ax[key].tick_params(axis="y", labelsize=16)


        self.scatter_data = self._dict_gen(self._keys, [])
        self.line_data = self._dict_gen(self._keys, [])

    def set_fig_title(self, title):
        self.fig_title = title
        self._replot()
        return self

    def set_title(self, title, key=None):
        if key is None:
            for key in self._keys:
                self.title[key] = title
        else:
            self.title[key] = title
        self._replot()
        return self

    def set_xLabel(self, xLabel, key=None):
        if key is None:
            for key in self._keys:
                self.xLabel[key] = xLabel
                self.ax[key].set_xlabel(xLabel, fontsize=self._xLabelFontSize[key])
        else:
            self.xLabel[key] = xLabel
            self.ax[key].set_xlabel(xLabel, fontsize=self._xLabelFontSize[key])
        return self

    def get_xLabel(self, key=None):
        if key is None:
            results = []
            for key in self._keys:
                results.append([self.ax[key].get_xlabel(), self._xLabelFontSize[key]])
            return results
        else:
            return [self.ax[key].get_xlabel(), self._xLabelFontSize[key]]

    def set_xLabel_fontsize(self, s, key=None):
        if key is None:
            for key in self._keys:
                self._xLabelFontSize[key] = s
                self.ax[key].set_xlabel(self.get_xLabel(key=key)[0], fontsize=self._xLabelFontSize[key])
        else:
            self._xLabelFontSize[key] = s
            self.ax[key].set_xlabel(self.get_xLabel(key=key)[0], fontsize=self._xLabelFontSize[key])
        return self

    def set_yLabel(self, yLabel, key=None):
        if key is None:
            for key in self._keys:
                self.yLabel[key] = yLabel
                self.ax[key].set_ylabel(yLabel, fontsize=self._yLabelFontSize[key])
        else:
            self.yLabel[key] = yLabel
            self.ax[key].set_ylabel(yLabel, fontsize=self._yLabelFontSize[key])
        return self

    def get_yLabel(self, key=None):
        if key is None:
            results = []
            for key in self._keys:
                results.append([self.ax[key].get_ylabel(), self._yLabelFontSize[key]])
            return results
        else:
            return [self.ax[key].get_ylabel(), self._yLabelFontSize[key]]

    def set_yLabel_fontsize(self, s, key=None):
        if key is None:
            for key in self._keys:
                self._yLabelFontSize[key] = s
                self.ax[key].set_ylabel(self.get_yLabel(key=key)[0], fontsize=self._yLabelFontSize[key])
        else:
            self._yLabelFontSize[key] = s
            self.ax[key].set_ylabel(self.get_yLabel(key=key)[0], fontsize=self._yLabelFontSize[key])
        return self

    def set_xlim(self, xmin, xmax, key=None, replot=True):
        if key is None:
            for key in self._keys:
                self._xMin[key] = xmin
                self._xMax[key] = xmax
        else:
            self._xMin[key] = xmin
            self._xMax[key] = xmax
        if replot:
            self._replot()
        return self

    def set_ylim(self, ymin, ymax, key=None, replot=True):
        if key is None:
            for key in self._keys:
                self._yMin[key] = ymin
                self._yMax[key] = ymax
        else:
            self._yMin[key] = ymin
            self._yMax[key] = ymax
        if replot:
            self._replot()
        return self

    def set_legend_location(self, loc, key=None):
        if key is None:
            for key in self._keys:
                self._legendLoc[key] = loc
        else:
            self._legendLoc[key] = loc
        self._replot()
        return self

    def add_line(self, x, y, key=None, *args, **kwargs):
        if key is None:
            key = self._keys[0]
        if kwargs.get("legend"):
            self.legend[key].append(kwargs.get("legend"))
            self.ax[key].legend(self.legend[key])
            self.showLegend[key] = True
        else:
            self.legend[key].append(None)
        self.line_data[key].append([x, y])
        self._replot(axisStick=True)
        self._xMin[key], self._xMax[key] = self.ax[key].get_xlim()
        self._yMin[key], self._yMax[key] = self.ax[key].get_ylim()
        return self

    def add_scatter(self, x, y, key=None, *args, **kwargs):
        if key is None:
            key = self._keys[0]
        if kwargs.get("legend"):
            if len(self.legend[key]) == 0:
                self.legend[key] = [kwargs.get("legend")]
            else:
                self.legend[key].append(kwargs.get("legend"))
            self.ax[key].legend(self.legend[key])
            self.showLegend[key] = True
        else:
            self.legend[key].append(None)

        #self.ax.scatter(x, y, color=color, s=self._markerSize)
        if len(self.scatter_data[key]) == 0:
            self.scatter_data[key] = [[x,y]]
        else:
            self.scatter_data[key].append([x,y])

        self._replot(axisStick=True)
        #self._xMin, self._xMax = self.ax.get_xlim()
        self._yMin[key], self._yMax[key] = self.ax[key].get_ylim()
        return self

    def set_marker_size(self, s, key=None):
        if key is None:
            for key in self._keys:
                self._markerSize[key] = s
        else:
            self._markerSize[key] = s
        self._replot()
        return self

    def set_number_xticks(self, num, key=None):
        if key is None:
            for key in self._keys:
                self._xtickNum[key] = num
        else:
            self._xtickNum[key] = num
        self._replot()
        return self

    def set_number_yticks(self, num, key=None):
        if key is None:
            for key in self._keys:
                self._ytickNum[key] = num
        else:
            self._ytickNum[key] = num
        self._replot()
        return self
    def set_xticks_fontsize(self, size, key=None):
        if key is None:
            for key in self._keys:
                self.ax[key].tick_params(axis='x', labelsize=size)
        else:
            self.ax[key].tick_params(axis='x', labelsize=size)
        return self

    def set_yticks_fontsize(self, size, key=None):
        if key is None:
            for key in self._keys:
                self.ax[key].tick_params(axis='y', labelsize=size)
        else:
            self.ax[key].tick_params(axis='y', labelsize=size)
        return self

    def set_legend_fontsize(self, size, key=None):
        if key is None:
            for key in self._keys:
                self._legendFontSize[key] = size
        else:
            self._legendFontSize[key] = size
        self._replot()
        return self

    def set_legend_scale(self, scale, key=None):
        if key is None:
            for key in self._keys:
                self._legendScale[key] = scale
        else:
            self._legendScale[key] = scale
        self._replot()
        return self

    def toggle_legend(self, key=None):
        if key is None:
            for key in self._keys:
                self.showLegend[key] = not self.showLegend[key]
        else:
            self.showLegend[key] = not self.showLegend[key]
        self._replot()

    def toggle_markers(self, key=None):
        if key is None:
            for key in self._keys:
                self._markers[key] = not self._markers[key]
        else:
            self._markers[key] = not self._markers[key]
        self._replot()

    def toggle_sci_notation(self, key=None):
        if key is None:
            for key in self._keys:
                self._sciNotation[key] = not self._sciNotation[key]
        else:
            self._sciNotation[key] = not self._sciNotation[key]
        self._replot()
        return self

    def toggle_semilog(self, key=None):
        if key is None:
            for key in self._keys:
                self._plotSemilog[key] = not self._plotSemilog[key]
        else:
            self._plotSemilog[key] = not self._plotSemilog[key]
        self._replot()
        return self

    def toggle_title(self, key=None):
        if key is None:
            for key in self._keys:
                self._showTitle[key] = not self._showTitle[key]
        else:
            self._showTitle[key] = not self._showTitle[key]
        self._replot()
        return self

    def toggle_mask_zeros(self, key=None):
        if key is None:
            for key in self._keys:
                self._maskZeros[key] = not self._maskZeros[key]
        else:
            self._maskZeros[key] = not self._maskZeros[key]
        self._replot()
        return self

    def _replot(self, *args, **kwargs):
        for key in self._keys:
            self.ax[key].cla()
            if self._plotSemilog[key]:
                self.ax[key].set_yscale("log")
            for n, data in enumerate(self.scatter_data[key]):
                x = data[0]
                y = data[1]
                if self._plotSemilog[key]:
                    pass
                    #y = np.log10(y)
                if self._maskZeros[key]:
                    try:
                        y[np.abs(y) < 1E-10] = np.nan
                    except TypeError:
                        try:
                            temp_y = y.val
                            temp_y[np.abs(temp_y) < 1E-10] = np.nan
                            y = temp_y
                        except:
                            raise
                if self._markers[key]:
                    self.ax[key].scatter(x, y, color=PLOTCOLORS[n], s=self._markerSize[key], marker=PLOTMARKERS[n])
                else:
                    self.ax[key].scatter(x, y, color=PLOTCOLORS[n], s=self._markerSize[key])

            for n, data in enumerate(self.line_data[key]):
                x = data[0]
                y = data[1]
                if self._plotSemilog[key]:
                    #y = np.log10(y)
                    pass
                if self._markers[key]:
                    self.ax[key].plot(x, y, color=PLOTCOLORS[n])
                else:
                    self.ax[key].plot(x, y, color=PLOTCOLORS[n])

            if self.showLegend[key]:
                self.ax[key].legend(self.legend[key],
                                    prop={'size': self._legendFontSize[key]},
                                    markerscale=self._legendScale[key],
                                    loc=self._legendLoc[key])
            if not kwargs.get("axisStick"):
                self.ax[key].set_xlim(self._xMin[key], self._xMax[key])
                self.ax[key].set_ylim(self._yMin[key], self._yMax[key])
            if self._showTitle[key]:
                self.ax[key].set_title(self.title[key])
            else:
                self.ax[key].set_title("")
            self.ax[key].set_xlabel(self.xLabel[key])
            self.ax[key].set_ylabel(self.yLabel[key])
            self.ax[key].locator_params(axis="x", nbins=self._xtickNum[key])
            if not self._plotSemilog[key]:
                self.ax[key].locator_params(axis="y", nbins=self._ytickNum[key])
            if self._sciNotation[key]:
                self.ax[key] = format_exponent(self.ax[key])
        return self

def format_exponent(ax, axis='y'):

    # Change the ticklabel format to scientific format
    ax.ticklabel_format(axis=axis, style='sci', scilimits=(-2, 2))

    # Get the appropriate axis
    if axis == 'y':
        ax_axis = ax.yaxis
        x_pos = 0.0
        y_pos = 1.0
        horizontalalignment='left'
        verticalalignment='bottom'
    else:
        ax_axis = ax.xaxis
        x_pos = 1.0
        y_pos = -0.05
        horizontalalignment='right'
        verticalalignment='top'

    # Run plt.tight_layout() because otherwise the offset text doesn't update
    plt.tight_layout()
    ##### THIS IS A BUG
    ##### Well, at least it's sub-optimal because you might not
    ##### want to use tight_layout(). If anyone has a better way of
    ##### ensuring the offset text is updated appropriately
    ##### please comment!

    # Get the offset value
    offset = ax_axis.get_offset_text().get_text()

    if len(offset) > 0:
        # Get that exponent value and change it into latex format
        minus_sign = u'\u2212'
        expo = np.float(offset.replace(minus_sign, '-').split('e')[-1])
        offset_text = r'x$\mathregular{10^{%d}}$' %expo

        # Turn off the offset text that's calculated automatically
        ax_axis.offsetText.set_visible(False)

        # Add in a text box at the top of the y axis
        ax.text(x_pos, y_pos, offset_text, transform=ax.transAxes,
               horizontalalignment=horizontalalignment,
               verticalalignment=verticalalignment,
                fontsize=24)
    return ax

if __name__ == "__main__":
    gtf = GT3FigureSinglePlot()
    x = np.linspace(0., 1., 100)
    linPlot = np.linspace(0., 3., 100)
    sinPlot = np.sin(np.pi * 2 * x)
    cosPlot = np.cos(np.pi * 2 * x)
    gtf.add_scatter(x, linPlot, legend="Linpl")
    gtf.add_scatter(x, sinPlot, legend="SinPlot")
    gtf.add_scatter(x, cosPlot)

    gtf2 = GT3FigureSinglePlot()
    linPlot = np.linspace(1., 3., 100)
    linPlotHuge = np.linspace(.01, 1E5, 100)
    gtf2.add_scatter(x, linPlot, legend="Linplot")
    gtf2.add_scatter(x, linPlotHuge, legend="Wooot")


