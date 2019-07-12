from __future__ import print_function
from __future__ import division

import numpy as np
from matplotlib.tri import Triangulation, UniformTriRefiner
from scipy.interpolate import griddata

from PlotUtils import *



##
##  LOAD DATA
##

selection = ["4l", "4m", "4e"]
models = ["VectorU", "ScalarU"]

mU, g, lhs, rhs = {}, {}, {}, {}

for mod in models:
    mU[mod], g[mod], lhs[mod] = {}, {}, {}

    for sel in selection:
        if sel in ["4m", "4e"] and mod == "ScalarU":
            continue

        infile = "exclusion_" + mod + "_" + sel + ".npz"
        npzfile = np.load(infile)
        mU[mod][sel], g[mod][sel], lhs[mod][sel] = npzfile["mU"], npzfile["g"], npzfile["lhs"]
        rhs[sel] = npzfile["rhs"]



##
##  PLOT
##

fig = plt.figure()
ax = plt.axes()

colors = {"4l":lPurple, "4m":lBlue, "4e":lRed}
linestyles = {"VectorU":"solid", "ScalarU":"dashed"}


p = {}
for mod in models:
    p[mod] = {}

    for sel in selection:
        if sel in ["4m", "4e"] and mod == "ScalarU":
            continue

        if False:
            tri = Triangulation(mU[mod][sel], g[mod][sel])
            refiner = UniformTriRefiner(tri)
            tri_refi, lhs_refi = refiner.refine_field(lhs[mod][sel], subdiv=2)
            p[mod][sel] = ax.tricontour(tri_refi, lhs_refi, levels=[np.squeeze(rhs[sel])],
                    colors=[colors[sel]], linestyles=[linestyles[mod]])

        else:
            p[mod][sel] = ax.tricontour(mU[mod][sel], g[mod][sel], lhs[mod][sel],
                    levels=[np.squeeze(rhs[sel])], colors=[colors[sel]], linestyles=[linestyles[mod]])

#       x = np.linspace(0, 90, 10000)
#       y = np.linspace(0, 0.5, 10000)
#       xv, yv = np.meshgrid(x, y)
#       lhs_grid = griddata((mU[mod][sel], g[mod][sel]), lhs[mod][sel], (xv, yv), method='cubic')
#       p[mod][sel] = ax.contour(xv, yv, lhs_grid,
#               levels=[np.squeeze(rhs[sel])], colors=[colors[sel]], linestyles=[linestyles[mod]])


# Dummy plots for legend
d = {}
for mod in models:
    d[mod] = ax.axhline(-5, linestyle=linestyles[mod], color='k')
for sel in selection:
    d[sel] = ax.axhline(-5, linestyle='solid', color=colors[sel])
d["blank"] = ax.axhline(-5, color='w')


# Axes and labels
ax.set_ylabel(r"$g_{\ell}$")
ax.set_ylim(0, 0.5)
ax.yaxis.set_ticks(np.arange(0, 0.55, step = 0.05))
ax.yaxis.set_ticks(np.arange(0, 0.51, step = 0.01), minor=True)

ax.set_xlabel(r"$m_{\mathrm{U}}$ (GeV)")
ax.set_xlim(0, 90)
ax.xaxis.set_ticks(np.arange(0, 100, step = 10))
ax.xaxis.set_ticks(np.arange(0, 92, step = 2), minor=True)

plt.grid()


# Text and legend
handles = [d["4e"], d["4m"], d["4l"], d["blank"], d["ScalarU"], d["VectorU"]]
labels = [r"$g_\ell = g_\mathrm{e}$; $g_\mu = 0$", r"$g_\ell = g_\mu$; $g_\mathrm{e} = 0$",
        r"$g_\ell = g_\mathrm{e} = g_\mu$", "", "Scalar U", "Vector U"]

ax.legend(handles, labels, loc = 'center left',
        frameon = True, facecolor='w', edgecolor='w', framealpha=1)

ax.text(0.05,  0.95,   "CMS",
        size = "xx-large",  weight = "bold", backgroundcolor='w',
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)
ax.text(0.05,  0.89,    "Work in Progress",
        size = "x-large",   style = "italic", backgroundcolor='w',
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)


# Save
fig_name = "Exclusion.pdf"
fig.savefig(fig_name)

print("Wrote plot to", fig_name)



##
##  LOG-LOG VERSION
##

plt.yscale('log')
ax.set_ylim(0.007, 0.5)

plt.xscale('log')
ax.set_xlim(3, 100)
plt.grid(which='minor')


# Save
fig_name = "Exclusion_log.pdf"
fig.savefig(fig_name)
#fig.savefig(fig_name, transparent=True)

print("Wrote plot to", fig_name)

