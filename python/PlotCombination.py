from __future__ import print_function
from __future__ import division

import numpy as np

from PlotUtils import *


##
##  LOAD DATA
##

# Single-parameter result
infile = "combination_1.npz"
npzfile = np.load(infile)

alpha_total_1 = npzfile['alpha_total']
sigma_stat_1 = npzfile['sigma_stat']
sigma_total_1 = npzfile['sigma_total']


# Three-parameter result
infile = "combination_3.npz"
npzfile = np.load(infile)

alpha_total_3 = npzfile['alpha_total']
sigma_stat_3 = npzfile['sigma_stat']
sigma_total_3 = npzfile['sigma_total']


# Four-parameter result
infile = "combination_4.npz"
npzfile = np.load(infile)

alpha_total_4 = npzfile['alpha_total']
sigma_stat_4 = npzfile['sigma_stat']
sigma_total_4 = npzfile['sigma_total']


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

# Single parameter

fig = plt.figure(figsize = (8,3))
ax = plt.axes()

p_total_1 = ax.fill_between( xpoints_fill,
        alpha_total_1 + sigma_total_1,          alpha_total_1 - sigma_total_1,  
        facecolor = lRed,                       alpha = 0.6
        )

p_stat_1 = ax.fill_between( xpoints_fill,
        alpha_total_1 + sigma_stat_1,           alpha_total_1 - sigma_stat_1,  
        hatch = '///',                         alpha = 0
        )

p_result_1 = ax.hlines(  alpha_total_1,     0, 16,
        color = lMarkerColor,   linewidth = lErrorLineWidth12)

p_total_12 = ax.errorbar(   xpoints,    alpha_total_12,     yerr = sigma_total_12,
        linewidth = 0,      ecolor = '#C0C0C0',       elinewidth = 4 * lErrorLineWidth12,
        marker = None,      capsize = 0
        )

ax.errorbar(   xpoints,  alpha_total_12,        yerr = sigma_stat_12,
        linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth12,
        marker = 'o',   capsize = 0,            markersize = lMarkerSize4l,
        markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
        )


# Fix axes
ax.set_ylim(0.5, 1.5)
ax.set_ylabel('Signal strength')

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
ax.text(0.025,  0.85,  "Work in Progress",
        size = "x-large",   style = "italic",
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)

fig.savefig("combination_1.pdf")

fig.clf()




# Single parameter

fig = plt.figure(figsize = (8,3))
ax = plt.axes()

for i in range(len(alpha_total_4)):
    p_total_4 = ax.fill_between( [xpoints[i*3] - 0.5, xpoints[i*3 + 2] + 0.5],
        alpha_total_4[i] + sigma_total_4[i],    alpha_total_4[i] - sigma_total_4[i],  
        facecolor = lRed,                       alpha = 0.6
        )
    p_stat_4 = ax.fill_between( [xpoints[i*3] - 0.5, xpoints[i*3 + 2] + 0.5],
        alpha_total_4[i] + sigma_stat_4[i],     alpha_total_4[i] - sigma_stat_4[i],  
        hatch = '///',                         alpha = 0
        )
    p_result_4 = ax.hlines(  alpha_total_4[i],      xpoints[i*3] - 0.5, xpoints[i*3 + 2] + 0.5,
        color = lMarkerColor, linewidth = lErrorLineWidth12)

p_total_12 = ax.errorbar(   xpoints,    alpha_total_12,     yerr = sigma_total_12,
        linewidth = 0,      ecolor = '#C0C0C0',       elinewidth = 4 * lErrorLineWidth12,
        marker = None,      capsize = 0
        )

ax.errorbar(   xpoints,  alpha_total_12,        yerr = sigma_stat_12,
        linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth12,
        marker = 'o',   capsize = 0,            markersize = lMarkerSize4l,
        markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
        )


# Fix axes
ax.set_ylim(0.5, 1.5)
ax.set_ylabel('Signal strength')

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
ax.text(0.025,  0.85,  "Work in Progress",
        size = "x-large",   style = "italic",
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)

fig.savefig("combination_4.pdf")

fig.clf()



##
##  THREE PARAMETER
##

xlabels = ["2012", "2016", "2017", "2018"] * 3
xpoints = np.array([1, 2, 3, 4,     6, 7, 8, 9,     11, 12, 13, 14])
xpoints_fill = np.arange(-1, xpoints[-1] + 2)

yearlabels = [r'$4\mu$', r'$2\mu 2\mathrm{e}$', r'$4\mathrm{e}$']
yearpoints = np.array([2.5, 7.5, 12.5])


alpha_total_12 = np.reshape(alpha_total_12, (4, 3)).transpose().flatten()
sigma_total_12 = np.reshape(sigma_total_12, (4, 3)).transpose().flatten()
sigma_stat_12 = np.reshape(sigma_stat_12, (4, 3)).transpose().flatten()


fig = plt.figure(figsize = (8,3))
ax = plt.axes()

for i in range(len(alpha_total_3)):
    p_total_3 = ax.fill_between( [xpoints[i*4] - 0.5, xpoints[i*4 + 3] + 0.5],
        alpha_total_3[i] + sigma_total_3[i],    alpha_total_3[i] - sigma_total_3[i],  
        facecolor = lRed,                       alpha = 0.6
        )
    p_stat_3 = ax.fill_between( [xpoints[i*4] - 0.5, xpoints[i*4 + 3] + 0.5],
        alpha_total_3[i] + sigma_stat_3[i],     alpha_total_3[i] - sigma_stat_3[i],  
        hatch = '///',                         alpha = 0
        )
    p_result_3 = ax.hlines(  alpha_total_3[i],      xpoints[i*4] - 0.5, xpoints[i*4 + 3] + 0.5,
        color = lMarkerColor, linewidth = lErrorLineWidth12)

p_total_12 = ax.errorbar(   xpoints,    alpha_total_12,     yerr = sigma_total_12,
        linewidth = 0,      ecolor = '#C0C0C0',       elinewidth = 4 * lErrorLineWidth12,
        marker = None,      capsize = 0
        )

ax.errorbar(   xpoints,  alpha_total_12,        yerr = sigma_stat_12,
        linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth12,
        marker = 'o',   capsize = 0,            markersize = lMarkerSize4l,
        markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
        )


# Fix axes
ax.set_ylim(0.5, 1.5)
ax.set_ylabel('Signal strength')

ax.set_xlim(xpoints[0] - 1, xpoints[-1] + 1)
plt.xticks(xpoints, xlabels)

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(yearpoints)
ax2.set_xticklabels(yearlabels)
ax2.tick_params(length=0)

#for i in range(len(yearlabels)):
#    ax.text(yearpoints[i], 0.05, yearlabels[i], verticalalignment = 'top', 
#            size = "xx-large", transform = ax.transAxes, usetex = True)


# CMS text
ax.text(0.025,  0.95,   "CMS",
        size = "xx-large",  weight = "bold",
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)
ax.text(0.025,  0.85,  "Work in Progress",
        size = "x-large",   style = "italic",
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)

fig.savefig("combination_3.pdf")

fig.clf()
