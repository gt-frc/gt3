#/usr/bin/python

from GT3 import gt3
import matplotlib.pyplot as plt


def _plot_base(rho, val, xLabel=r'$\rho$', yLabel="Value", title="Title", color='red', edge=False):

    plot = plt.figure()
    fig = plot.add_subplot(111)
    fig.set_xlabel(xLabel, fontsize=30)
    fig.set_ylabel(yLabel, fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    fig.set_title(title)
    if edge:
        fig.set_xlim(0.85, 1.0)
    fig.scatter(rho, val, color=color, s=8)
    plt.show()
    return fig


class IOLComparison:

    def __init__(self, filename):

        self.plasmaWithIOL = gt3(inputFile=filename)
        self.plasmaWithIOL.run_radial_transport()

        self.plasmaWithoutIOL = gt3(inputFile=filename)
        self.plasmaWithoutIOL.disable_IOL()
        self.plasmaWithoutIOL.run_radial_transport()

    def comp_gamma_diff(self, edge=True):
        fig = _plot_base(self.plasmaWithIOL.rtrans.rhor,
                         self.plasmaWithIOL.rtrans.gamma_diff_D,
                         xLabel=r'$\rho$',
                         yLabel=r"$\Gamma_{i,r}[\#/m^2 s]$",
                         title="Radial particle flux w/ & w/out IOL",
                         color="red", edge=edge)
        fig.scatter(self.plasmaWithoutIOL.rtrans.rhor,
                    self.plasmaWithoutIOL.rtrans.gamma_diff_D,
                    color="blue", s=8)
        fig.legend(["IOL Corrected", "Uncorrected"])
        return fig

    def comp_D(self, edge=True):
        fig = _plot_base(self.plasmaWithIOL.rtrans.rhor,
                         self.plasmaWithIOL.rtrans.D_i,
                         xLabel=r'$\rho$',
                         yLabel=r"$D{i,r}[m/s^2]$",
                         title="Radial Diffusion Coefficient w/ & w/out IOL",
                         color="red", edge=edge)
        fig.scatter(self.plasmaWithoutIOL.rtrans.rhor,
                    self.plasmaWithoutIOL.rtrans.D_i,
                    color="blue", s=8)
        fig.legend(["IOL Corrected", "Uncorrected"])
        return fig
