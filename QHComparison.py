#!/usr/bin/python

from GT3 import gt3
import matplotlib.pyplot as plt
from GT3.Psi import UPPER_XPT, LOWER_XPT
import numpy as np

kwargs = {
    'rtrans_override': {
        'splines': {
            'T_i': False
        }
    },
    'psi_args': {
        'debug': True,
        'xpt_select': LOWER_XPT
    }
}

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



plasma163477 = gt3(inputFile="inputs/togt3_d3d_163477_1800")
plasma163477.run_radial_transport()
plasma163518 = gt3(inputFile="inputs/togt3_d3d_163518_2350", **kwargs)
plasma163518.run_radial_transport()
plasma163477.rtrans.set_marker_size(100)
plasma144977 = gt3(inputFile="inputs/togt3_d3d_144977_3000")
plasma144977.run_radial_transport()
fig=plasma163477.rtrans._plot_base(plasma163477.core.E_r.fsa.val,
                                   yLabel=r'$E_r \left[\frac{kV}{m}\right]$', title="", edge=False)
fig.scatter(plasma163518.rtrans.rhor, plasma163518.core.E_r.fsa.val, color="blue", s=100, marker="x")
fig.legend(["163477",
            "163518"], prop={'size': 20}, markerscale=1.5)
ax = fig.axes
# axins = zoomed_inset_axes(ax, 10, loc="upper right")  # type: Axes
# axins.scatter(plasmaRMP.rtrans.rhor[85:90], plasmaRMP.rtrans.chi.i.chi4[85:90], color='red', s=200)
# axins.scatter(plasmaRMP.rtrans.rhor[85:90], plasma163518.rtrans.chi.i.chi4[85:90], color='blue', s=200, marker="x")
# axins.set_xlim(plasmaRMP.rtrans.rhor[86], plasmaRMP.rtrans.rhor[89])
# axins.set_ylim(2.915, 2.975)
#plt.xticks(visible=False)
#plt.yticks(visible=False)
# mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")

figIOL = plasma163477.iol.plot_all_i(noTitle=True, edge=False, noColor=True)
figIOL.axes[0].scatter(plasma163518.rtrans.rhor, plasma163518.iol.forb_d_therm_1D, color="black", marker="x")
figIOL.axes[1].scatter(plasma163518.rtrans.rhor, plasma163518.iol.eorb_d_therm_1D, color="black", marker="x")
figIOL.axes[2].scatter(plasma163518.rtrans.rhor, plasma163518.iol.morb_d_therm_1D, color="black", marker="x")

figIOL.axes[0].legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)
figIOL.axes[1].legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)
figIOL.axes[2].legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)

figB = plasma163477.rtrans._plot_base(plasma163477.core.B.tot.fsa.val,
                                      yLabel=r'$B^{tot} \left[T\right]$', title="", edge=False)
figB.scatter(plasma163518.rtrans.rhor, plasma163518.core.B.tot.fsa.val, color="blue", s=100, marker="x")
figB.legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)


# figni = plasma163477.rtrans._plot_base(plasma163477.core.n.D.fsa.val,
#                                       yLabel=r'$n_{j} \left[\frac{\#}{m^2}\right]$', title="", edge=False)
# figni.scatter(plasma163518.rtrans.rhor, plasma163518.core.n.D.fsa.val, color="blue", s=100, marker="x")
# figni.legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)

figTi = plasma163477.rtrans._plot_base(plasma163477.core.T.i.kev.fsa.val,
                                      yLabel=r'$T_{j} \left[keV\right]$', title="", edge=False)
figTi.scatter(plasma163518.rtrans.rhor, plasma163518.core.T.i.kev.fsa.val, color="blue", s=100, marker="x")
figTi.legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)
figTi.figure.canvas.set_window_title("Ion Temperature")

figTe = plasma163477.rtrans._plot_base(plasma163477.core.T.e.kev.fsa.val,
                                      yLabel=r'$T_{e} \left[keV\right]$', title="", edge=False)
figTe.scatter(plasma163518.rtrans.rhor, plasma163518.core.T.e.kev.fsa.val, color="blue", s=100, marker="x")
figTe.legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)
figTe.figure.canvas.set_window_title("Electron Temperature")

figvtorD = plasma163477.rtrans._plot_base(plasma163477.core.v.D.tor.fsa.val,
                                      yLabel=r'$V_{j,\phi} \left[\frac{m}{s}\right]$', title="", edge=False)
figvtorD.scatter(plasma163518.rtrans.rhor, plasma163518.core.v.D.tor.fsa.val, color="blue", s=100, marker="x")
figvtorD.legend(["163477", "163518"], prop={'size': 20}, markerscale=1.5)
figvtorD.figure.canvas.set_window_title("Toroidal Velocity")

figvmin = plasma163477.rtrans._plot_base(plasma163477.iol.vsep_min_d.val,
                                      yLabel=r'$v_{min,sep}$', title="", edge=True)
figvmin.scatter(plasma163518.rtrans.rhor, plasma163518.iol.vsep_min_d.val, color="blue", s=100, marker="o")
figvmin.scatter(plasma163477.rtrans.rhor, plasma163477.core.v.D.tor.fsa.val, color="red", s=100, marker="x")
figvmin.scatter(plasma163518.rtrans.rhor, plasma163518.core.v.D.tor.fsa.val, color="blue", s=100, marker="x")
figvmin.legend(["163477", "163518", "163477_vmin", "163518_vmin"], prop={'size': 20}, markerscale=1.5)
figvmin.figure.canvas.set_window_title("VMin Compare")

figeps = plasma163477.rtrans._plot_base(plasma163477.iol.vsep_min_d.val,
                                      yLabel=r'$v_{min,sep}$', title="", edge=True)
figeps.scatter(plasma163518.rtrans.rhor, plasma163518.iol.vsep_min_d.val, color="blue", s=100, marker="o")
figeps.scatter(plasma163477.rtrans.rhor, plasma163477.core.v.D.tor.fsa.val, color="red", s=100, marker="x")
figeps.scatter(plasma163518.rtrans.rhor, plasma163518.core.v.D.tor.fsa.val, color="blue", s=100, marker="x")
figeps.legend(["163477", "163518", "163477_vmin", "163518_vmin"], prop={'size': 20}, markerscale=1.5)
figvmin.figure.canvas.set_window_title("VMin Compare")

plot = plt.figure()
plot.canvas.set_window_title("COOL")
figStacked1 = plot.add_subplot(511)
figStacked1.set_ylabel(r"$T_{e} \left[keV\right]$")
figStacked1.scatter(plasma163518.rtrans.rhor, plasma163477.core.T.e.kev.fsa.val, color="blue", s=100, marker="x")
figStacked1.scatter(plasma163518.rtrans.rhor, plasma163518.core.T.e.kev.fsa.val, color="red", s=100, marker="x")
figStacked1.legend(["163477", "163518"], prop={'size': 12}, markerscale=1)

figStacked2 = plot.add_subplot(512)
figStacked2.set_ylabel(r"$T_{i} \left[keV\right]$")
figStacked2.scatter(plasma163518.rtrans.rhor, plasma163477.core.T.i.kev.fsa.val, color="blue", s=100, marker="x")
figStacked2.scatter(plasma163518.rtrans.rhor, plasma163518.core.T.i.kev.fsa.val, color="red", s=100, marker="x")
figStacked2.legend(["163477", "163518"], prop={'size': 12}, markerscale=1)

figStacked3 = plot.add_subplot(513)
figStacked3.set_ylabel(r"$n_e \left[\frac{\#}{m^3}\right]$")
figStacked3.scatter(plasma163518.rtrans.rhor, plasma163477.core.n.e.fsa.val, color="blue", s=100, marker="x")
figStacked3.scatter(plasma163518.rtrans.rhor, plasma163518.core.n.e.fsa.val, color="red", s=100, marker="x")
figStacked3.legend(["163477", "163518"], prop={'size': 12}, markerscale=1)

figStacked4 = plot.add_subplot(514)
figStacked4.set_ylabel(r"$n_i \left[\frac{\#}{m^3}\right]$")
figStacked4.scatter(plasma163518.rtrans.rhor, plasma163477.core.n.i.fsa.val, color="blue", s=100, marker="x")
figStacked4.scatter(plasma163518.rtrans.rhor, plasma163518.core.n.i.fsa.val, color="red", s=100, marker="x")
figStacked4.legend(["163477", "163518"], prop={'size': 12}, markerscale=1)

figStacked5 = plot.add_subplot(515)
figStacked5.set_ylabel(r"$B^{tot}[T]$")
figStacked5.set_xlabel(r"$\rho$")
figStacked5.scatter(plasma163518.rtrans.rhor, plasma163477.core.B.tot.fsa.val, color="blue", s=100, marker="x")
figStacked5.scatter(plasma163518.rtrans.rhor, plasma163518.core.B.tot.fsa.val, color="red", s=100, marker="x")
figStacked5.legend(["163477", "163518"], prop={'size': 12}, markerscale=1)

plt.draw()
plt.show()


