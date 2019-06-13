from __future__ import print_function
from __future__ import division

import numpy as np

from PlotUtils import *
from CutsComb import *


##
##  LOAD DATA
##

# Single-parameter result
infile = "combination_1.npz"
npzfile = np.load(infile)

alpha_total_1 = npzfile['alpha_total']
sigma_stat_1 = npzfile['sigma_stat']
sigma_total_1 = npzfile['sigma_total']


# Twelve-parameter result
infile = "combination_12.npz"
npzfile = np.load(infile)

alpha_total_12 = npzfile['alpha_total']
sigma_stat_12 = npzfile['sigma_stat']
sigma_total_12 = npzfile['sigma_total']



##
##  LABELS
##

xlabels = [r'$4\mu$', r'$2\mu 2\mathrm{e}$', r'$4\mathrm{e}$'] * 4
xpoints = np.array([1, 2, 3,    5, 6, 7,    9, 10, 11,  13, 14, 15])
xpoints_fill = np.arange(-1, xpoints[-1] + 2)

yearpoints = np.array([2, 6, 10, 14])
yearlabels = ["2012", "2016", "2017", "2018"]



##
##  PLOT
##

fig = plt.figure(figsize = (8,4))
ax = plt.axes()

p_total_12 = ax.errorbar(   xpoints,    alpha_total_12,     yerr = sigma_total_12,
        linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth12,
        marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize4l,
        markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
        )

#p_stat_12 = ax.errorbar(   xpoints,  alpha_total_12,    yerr = sigma_stat_12,
#        linewidth = 0,  ecolor = lOrange,       elinewidth = lErrorLineWidth12,
#        marker = 'o',   capsize = lCapSize12,   markersize = lMarkerSize4l,
#        markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
#        )

p_total_12 = ax.fill_between( xpoints_fill,
        alpha_total_1 + sigma_total_1,          alpha_total_1 - sigma_total_1,  
        facecolor = lBlue,                      alpha = 0.4
        )

p_stat_12 = ax.fill_between( xpoints_fill,
        alpha_total_1 + sigma_stat_1,           alpha_total_1 - sigma_stat_1,  
        facecolor = lBlue,                      alpha = 0.4
        )

p_result_12 = ax.axhline(  alpha_total_1,       color = lMarkerColor,
        linewidth = lErrorLineWidth12)


# Fix axes
ax.set_ylim(0.5, 1.5)
ax.set_ylabel(r'$\hat{\alpha}$')

ax.set_xlim(xpoints[0] - 1, xpoints[-1] + 1)
plt.xticks(xpoints, xlabels)

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(yearpoints)
ax2.set_xticklabels(yearlabels)
ax2.tick_params(length=0)


# CMS text
ax.text(0.025,  0.95,   "CMS",
        size = "xx-large",  weight = "bold",
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)
ax.text(0.025,  0.875,  "Work in Progress",
        size = "x-large",   style = "italic",
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)

fig.savefig("combination_1.pdf")
