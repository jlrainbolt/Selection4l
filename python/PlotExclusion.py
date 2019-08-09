from __future__ import print_function
from __future__ import division

import numpy as np
from scipy.interpolate import SmoothBivariateSpline

from PlotUtils import *
import matplotlib.patches as pat



##
##  LOAD DATA
##

selection = ["4l", "4m", "4e"]
models = ["VectorU", "ScalarU"]

c_up = {"ScalarU":{"4l":1.10, "4m":1.10, "4e":1.22}, "VectorU":{"4l":5.00, "4m":5.00, "4e":5.00}}
c_dn = {"ScalarU":{"4l":0.94, "4m":0.93, "4e":0.94}, "VectorU":{"4l":0.65, "4m":0.50, "4e":0.50}}

mU, g, lhs, rhs = {}, {}, {}, {}

for mod in models:
    mU[mod], g[mod], lhs[mod] = {}, {}, {}

    for sel in selection:
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


p, x, y, z, m = {}, {}, {}, {}, {}

for mod in models:
    p[mod], x[mod], y[mod], z[mod], m[mod] = {}, {}, {}, {}, {}

    for sel in selection:
        x_, y_, z_ = mU[mod][sel], g[mod][sel], lhs[mod][sel]
        c = np.squeeze(rhs[sel])

        f = SmoothBivariateSpline(np.log10(x_), np.log10(y_), z_, kx=3, ky=2)
        h = SmoothBivariateSpline(np.log10(x_), np.log10(y_), z_, kx=1, ky=1)

        xv = np.linspace(4, 90, 500)
        yv = 5 * np.logspace(-3, -1, 500)
        zf = f(np.log10(xv), np.log10(yv)).T
        zh = h(np.log10(xv), np.log10(yv)).T

        xv, yv = np.meshgrid(xv, yv)
        mask = np.logical_or(np.less(zh, c * c_dn[mod][sel]), np.greater(zh, c * c_up[mod][sel]))
        if mod == "VectorU":
            m2 = np.logical_and(np.less(zh, c * 0.9), np.greater(xv, 70))
            mask = np.logical_or(mask, m2)
        zv = np.ma.array(zf, mask=mask)

        p[mod][sel] = ax.contour(xv, yv, zv,
                levels=[c], colors=[colors[sel]], linestyles=[linestyles[mod]])

        x[mod][sel], y[mod][sel], z[mod][sel], m[mod][sel] = xv, yv, zf, mask
        



# Dummy plots for legend
d = {}
for mod in models:
    d[mod] = ax.axhline(-5, linestyle=linestyles[mod], color='k')
for sel in selection:
    d[sel] = ax.axhline(-5, linestyle='solid', color=colors[sel])
d["blank"] = ax.axhline(-5, color='w')


# Axes and labels
ax.set_axisbelow(True)

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

r = pat.Rectangle((4.25, 0.055), 23, 0.4, color='w')
ax.add_patch(r)

ax.legend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.025, 0.85),
        frameon = True, facecolor='w', edgecolor='w', framealpha=1)

ax.text(0.05,  0.95,    "CMS",
        size = "xx-large",  weight = "bold", backgroundcolor='w',
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)
ax.text(0.20,  0.945,   "Work in Progress",
        size = "x-large",   style = "italic", backgroundcolor='w',
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)


# Save
fig_name = "Exclusion.pdf"
#fig.savefig(fig_name)

#print("Wrote plot to", fig_name)



##
##  LOG-LOG VERSION
##

plt.yscale('log')
ax.set_ylim(0.007, 0.5)

plt.xscale('log')
ax.set_xlim(4, 100)
plt.grid(which='minor')


# Save
#fig_name = "Exclusion_log.pdf"
fig.savefig(fig_name)

print("Wrote plot to", fig_name)



##
##  SAVE
##

outfile = "exclusion_4m.npz"
np.savez(outfile, x=x["VectorU"]["4m"], y=y["VectorU"]["4m"], z=z["VectorU"]["4m"],
        m=m["VectorU"]["4m"])

print("Wrote arrays to", outfile)
