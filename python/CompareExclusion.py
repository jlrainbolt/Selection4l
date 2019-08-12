from __future__ import print_function
from __future__ import division

import numpy as np
from scipy.interpolate import SmoothBivariateSpline
from ROOT import TFile, TCanvas, TGraph, Double

from PlotUtils import *
import matplotlib.patches as pat



##
##  LOAD DATA
##

selection = ["smp", "exo"]
types = ["exp", "obs"]

infile = "curves_VectorU_4m.npz"
npzfile = np.load(infile)
x, y, z, m = npzfile["x"], npzfile["y"], npzfile["z"], npzfile["m"]
print("Got arrays from", infile)

infile = "exclusion_VectorU_4m.npz"
npzfile = np.load(infile)
c = np.squeeze(npzfile["rhs"])
print("Got data from", infile)



##
##  SMP-19-007
##

fig = plt.figure()
ax = plt.axes()

colors = {"smp":lBlue, "exo":lOrange}
linestyles = {"obs":"solid", "exp":"dashed"}
smp = {}

z = np.ma.array(z, mask=m)
smp["exp"] = ax.contour(x, y, z, levels=[c], colors=[colors["smp"]], linestyles=[linestyles["exp"]])



##
##  EXO-18-008
##

# Get limits from file
exo_p = {}
rootfile = TFile.Open("../data/limits_EXO-18-008.root")
canvas = rootfile.Get("c")
exo_p["exp"] = canvas.GetListOfPrimitives().At(4)
exo_p["obs"] = canvas.GetListOfPrimitives().At(6)
rootfile.Close()

# Read TGraphs into numpy arrays
exo_x, exo_y = {}, {}
for t in types:
    N = exo_p[t].GetN()
    exo_x[t], exo_y[t] = np.zeros(N), np.zeros(N)

    x_, y_ = Double(0), Double(0)
    for i in range(N):
        exo_p[t].GetPoint(i, x_, y_)
        exo_x[t][i], exo_y[t][i] = x_, y_

# Now make the plots
exo = {}
for t in types:
    exo[t] = ax.plot(exo_x[t], exo_y[t], linestyle=linestyles[t], color=colors["exo"])



# Dummy plots for legend
d = {}
for t in types:
    d[t] = ax.axhline(-5, linestyle=linestyles[t], color='k')
for sel in selection:
    d[sel] = ax.axhline(-5, linestyle='solid', color=colors[sel])
d["blank"] = ax.axhline(-5, color='w')


# Axes and labels
ax.set_axisbelow(True)

ax.set_ylabel(r"$g_{\mu}$")
plt.yscale('log')
ax.set_ylim(0.002, 0.5)

ax.set_xlabel(r"$m_{\mathrm{U}}$ (GeV)")
plt.xscale('log')
ax.set_xlim(4, 100)

plt.grid(which='both')

r = pat.Rectangle((4.25, 0.035), 23, 0.42, color='w')
ax.add_patch(r)


# Text and legend
handles = [d["smp"], d["exo"], d["blank"], d["exp"], d["obs"]]
labels = [r"SMP-19-007", r"EXO-18-008", "", "Expected", "Observed"]

ax.legend(handles, labels, loc = 'upper left', bbox_to_anchor = (0.025, 0.85),
        frameon = True, facecolor='w', edgecolor='w', framealpha=1)

ax.text(0.05,  0.95,    "CMS",
        size = "xx-large",  weight = "bold", backgroundcolor='w',
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)
ax.text(0.20,  0.945,   "Work in Progress",
        size = "x-large",   style = "italic", backgroundcolor='w',
        verticalalignment = 'top', transform = ax.transAxes, usetex = False)



fig_name = "Exclusion_EXO.pdf"
fig.savefig(fig_name)

print("Wrote plot to", fig_name)
