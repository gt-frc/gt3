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
    'inp_override': {
      'er_file': 1.0,
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
plasma163518 = gt3(inputFile="inputs/togt3_d3d_163518_2350", **kwargs)
plasma144977 = gt3(inputFile="inputs/togt3_d3d_144977_3000")
plasma163518.run_radial_transport()
plasma163477.run_radial_transport()
plasma144977.run_radial_transport()
plasma144977.rtrans.set_marker_size(100)
plasma163477.rtrans.set_marker_size(100)
plasma163518.rtrans.set_marker_size(100)

figER=plasma163477.rtrans._plot_base(plasma163477.core.E_r.fsa.val,
                                     yLabel=r'$E_r \left[\frac{kV}{m}\right]$', title="", edge=False)
figER.scatter(plasma163518.rtrans.rhor, plasma163518.core.E_r.fsa.val, color="blue", s=100, marker="o")
figER.scatter(plasma144977.rtrans.rhor, plasma144977.core.E_r.fsa.val, color="green", s=100, marker="o")
figER.legend(["163477", "163518", "144977"], prop={'size': 20}, markerscale=1.5)
figER.figure.canvas.set_window_title("Electric Field")

figPhi=plasma163477.rtrans._plot_base(plasma163477.core.E_pot.fsa.val,
                                     yLabel=r'$\phi_r \left[kV\right]$', title="", edge=False)
figPhi.scatter(plasma163518.rtrans.rhor, plasma163518.core.E_pot.fsa.val, color="blue", s=100, marker="o")
figPhi.scatter(plasma144977.rtrans.rhor, plasma144977.core.E_pot.fsa.val, color="green", s=100, marker="o")
figPhi.legend(["163477", "163518", "144977"], prop={'size': 20}, markerscale=1.5)
figPhi.figure.canvas.set_window_title("Electric Potential")

figPsi=plasma163477.rtrans._plot_base(plasma163477.core.psi.psi[:,-1],
                                     yLabel=r'$\psi$', title="", edge=False)
figPsi.scatter(plasma163518.rtrans.rhor, plasma163518.core.psi.psi[:,-1], color="blue", s=100, marker="o")
figPsi.scatter(plasma144977.rtrans.rhor, plasma144977.core.psi.psi[:,-1], color="green", s=100, marker="o")
figPsi.legend(["163477", "163518", "144977"], prop={'size': 20}, markerscale=1.5)
figPsi.figure.canvas.set_window_title("Psi")

plt.draw()
plt.show()