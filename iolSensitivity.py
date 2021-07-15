#!/usr/bin/python

from GT3 import gt3
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import matplotlib.pyplot as plt
from matplotlib.pyplot import Axes
from GT3.utilities.PlotBase import MARKERSIZE

filename = "inputs/togt3_d3d_170672_1900"

plasmaIOL = gt3(inputFile=filename)
plasmaIOL.override_ntrl_cpus(6)
plasmaIOL.run_neutrals()
plasmaIOL.run_radial_transport()
plasmaNoIOL = gt3(inputFile=filename)
plasmaNoIOL.disable_IOL()
plasmaIOL.rtrans.set_marker_size(150)
fig=plasmaIOL.rtrans._plot_base(plasmaIOL.rtrans.chi.i.chi4,
                                yLabel=r'$\chi_{r}$', title="", edge=True)
fig.scatter(plasmaNoIOL.rtrans.rhor, plasmaNoIOL.rtrans.chi.i.chi4, color="blue", s=200, marker="x")
fig.set_ylim(2.0, 5.0)
fig.legend(["w/ IOL",
            "w/out IOL"], prop={'size': 20}, markerscale=1.5)
ax = fig.axes
# axins = zoomed_inset_axes(ax, 10, loc="upper right")  # type: Axes
# axins.scatter(plasmaIOL.rtrans.rhor[85:90], plasmaIOL.rtrans.chi.i.chi4[85:90], color='red', s=200)
# axins.scatter(plasmaIOL.rtrans.rhor[85:90], plasmaNoIOL.rtrans.chi.i.chi4[85:90], color='blue', s=200, marker="x")
# axins.set_xlim(plasmaIOL.rtrans.rhor[86], plasmaIOL.rtrans.rhor[89])
# axins.set_ylim(2.915, 2.975)
#plt.xticks(visible=False)
#plt.yticks(visible=False)
# mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")
plt.draw()
plt.show()

fig = plasmaIOL.rtrans.gamma.D.diff.plot()
fig.scatter(plasmaNoIOL.rtrans.rhor, plasmaNoIOL.rtrans.gamma.D.diff.val, color="blue", s=MARKERSIZE)

figD = plasmaIOL.rtrans._plot_base(plasmaIOL.rtrans.D_i,
                                   yLabel=r'$D_r$', title="D_r", edge=True)
figD.set_ylim(0.05, 0.1)
figD.scatter(plasmaNoIOL.rtrans.rhor, plasmaNoIOL.rtrans.D_i, color="blue", s=200, marker="x")
figD.legend(["w/ Neutrals",
            "w/out Neutrals"], prop={'size': 20}, markerscale=1.5)

