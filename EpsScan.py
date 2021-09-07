#!/usr/bin/python

from GT3 import gt3
import matplotlib.pyplot as plt
from GT3.Psi import UPPER_XPT, LOWER_XPT
import numpy as np

#import matplotlib
#matplotlib.use('Agg')

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

shot = "inputs/togt3_d3d_163518_2350"
############################################################
#
# Sensitivity study of the Radial Electric Field.
#
#############################################################

Ersensitivity = [(0.1, None),
                 (0.01, None),
                 (0.5, None),
                 (1.0, None),
                 (5.0, None),
                 (10.0, None),
                 (50.0, None),
                 (-1.0, None)]

plasma163477 = gt3(inputFile="inputs/togt3_d3d_163477_1800")
plasma163518 = gt3(inputFile="inputs/togt3_d3d_163518_2350", **kwargs)
plasma144977 = gt3(inputFile="inputs/togt3_d3d_144977_3000")
plasma123302 = gt3(inputFile="inputs/togt3_d3d_123302_2810")
plasma123301 = gt3(inputFile="inputs/togt3_d3d_123301_2800")
plasma163518.run_IOL()
plasma163477.run_IOL()
plasma144977.run_IOL()
plasma123302.run_IOL()
plasma123301.run_IOL()

plasma144977.iol.set_marker_size(50)
figeps=plasma144977.iol._plot_base(plasma144977.iol.eps_d.val,
                                     yLabel=r'$E_r \left[\epsilon_{min}\right]$', title=r"\epsilon comparison", edge=False)
figeps.scatter(plasma144977.iol.rho[:, 0], plasma123301.iol.eps_d.val, color="purple", s=50, marker="o")
figeps.scatter(plasma144977.iol.rho[:, 0], plasma123302.iol.eps_d.val, color="green", s=50, marker="o")
figeps.scatter(plasma144977.iol.rho[:, 0], plasma163518.iol.eps_d.val, color="blue", s=50, marker="o")
figeps.scatter(plasma144977.iol.rho[:, 0], plasma163477.iol.eps_d.val, color="orange", s=50, marker="o")


figeps.legend(["144977", "123301", "123302", "163518", "163477"], prop={'size': 20}, markerscale=1.5)
figeps.axhline(y=1.0, color='black')
figeps.figure.canvas.set_window_title("epsilon Comparison")

fig=plasma144977.iol._plot_base(plasma144977.core.E_pot.fsa.val,
                                     yLabel=r'$E_r \left[\epsilon_{min}\right]$', title=r"\epsilon comparison", edge=False)
fig.scatter(plasma144977.iol.rho[:, 0], plasma123301.core.E_pot.fsa.val, color="purple", s=50, marker="o")
fig.scatter(plasma144977.iol.rho[:, 0], plasma123302.core.E_pot.fsa.val, color="green", s=50, marker="o")
fig.scatter(plasma144977.iol.rho[:, 0], plasma163518.core.E_pot.fsa.val, color="blue", s=50, marker="o")
fig.scatter(plasma144977.iol.rho[:, 0], plasma163477.core.E_pot.fsa.val, color="orange", s=50, marker="o")

fig.legend(["144977", "123301", "123302", "163518", "163477"], prop={'size': 20}, markerscale=1.5)
fig.figure.canvas.set_window_title("E_pot Comparison")
