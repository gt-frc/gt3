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

for n, val in enumerate(Ersensitivity):
    kwargs['inp_override']['er_file'] = val[0]
    plasma = gt3(inputFile=shot, **kwargs)
    plasma.run_IOL()
    Ersensitivity[n] = (val, plasma.iol.vsep_min_d.val)


############################################################
#
# Sensitivity study of the Toroidal magnetic field.
#
#############################################################

BTsensitivity = [(0.8, None),
                 (0.9, None),
                 (1.0, None),
                 (1.1, None),
                 (1.2, None),
                 (-1.0, None)]
kwargs['inp_override']['er_file'] = None
for n, val in enumerate(BTsensitivity):
    kwargs['inp_override']['Bt0'] = val[0]
    plasma = gt3(inputFile=shot, **kwargs)
    plasma.run_IOL()
    BTsensitivity[n] = (val, plasma.iol.vsep_min_d.val)

############################################################
#
# Sensitivity study of the Ion Temperature
#
#############################################################

TIsensitivity = [(0.8, None),
                 (0.9, None),
                 (1.0, None),
                 (1.1, None),
                 (1.5, None),
                 (2.0, None),
                 (3.0, None)]
kwargs['inp_override']['Bt0'] = None
for n, val in enumerate(TIsensitivity):
    kwargs['inp_override']['Ti_file'] = val[0]
    plasma = gt3(inputFile=shot, **kwargs)
    plasma.run_IOL()
    TIsensitivity[n] = (val, plasma.iol.vsep_min_d.val)


############################################################
#
# Print results
#
#############################################################

plasma.iol.set_marker_size(50)

# Print Bt

figBt = plasma.iol._plot_base(BTsensitivity[0][1],
                            yLabel=r'$v_{sep-min} \left[\frac{m}{s}\right]$', title=r"$B_{\phi}$ Sensitivity", edge=False)
for n, val in enumerate(BTsensitivity):
    if n == 0:
        continue
    figBt.scatter(plasma.iol.rho[:, 0], val[1], s=50)
figBt.legend([a[0] for a, b in BTsensitivity])

# Print Er

figEr = plasma.iol._plot_base(Ersensitivity[0][1],
                            yLabel=r'$v_{sep-min} \left[\frac{m}{s}\right]$', title=r'$E_r$ Sensitivity', edge=False)
for n, val in enumerate(Ersensitivity):
    if n == 0:
        continue
    figEr.scatter(plasma.iol.rho[:, 0], val[1], s=50)
figEr.legend([a[0] for a, b in Ersensitivity])

# Print TI

figTI = plasma.iol._plot_base(TIsensitivity[0][1],
                            yLabel=r'$v_{sep-min} \left[\frac{m}{s}\right]$', title=r'$T_i$ Sensitivity', edge=False)
for n, val in enumerate(TIsensitivity):
    if n == 0:
        continue
    figTI.scatter(plasma.iol.rho[:, 0], val[1], s=50)
figTI.legend([a[0] for a, b in TIsensitivity])