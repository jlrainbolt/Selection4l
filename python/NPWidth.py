from __future__ import print_function
from __future__ import division


import sys
import numpy as np
from os import listdir
from scipy.stats import norm
from scipy.interpolate import SmoothBivariateSpline
from ROOT import TFile, TTree, TParameter

from PlotUtils import *
import matplotlib.patches as pat

from Cuts2018 import EOS_PATH



##
##  OPTIONS
##

model = "VectorU"
selection = "4m"


##
##  COMBINATION RESULTS
##

# Single-parameter result
infile = "combination_1.npz"
npzfile = np.load(infile)
delta_total_1 = npzfile['delta_total']
bf_pred = npzfile['bf_pred']


# Three-parameter result
infile = "combination_3.npz"
npzfile = np.load(infile)
delta_total_3 = npzfile['delta_total']

if selection == "4l":
    delta = delta_total_1[0]
    pred = np.sum(bf_pred)
elif selection == "4m":
    delta = delta_total_3[0]
    pred = bf_pred[0]
elif selection == "4e":
    delta = delta_total_3[2]
    pred = bf_pred[2]



##
##  LIMITS
##

inpath = EOS_PATH + "/Generated/"

# Standard model
sm_name = "StandardModel"
sm_file = TFile.Open(inpath + sm_name + "_" + selection + ".root")
tree = sm_file.Get(selection + "_" + sm_name) 

sm_acc = tree.GetEntries("isFiducial") / tree.GetEntries()
sm_xsec = round(sm_file.Get("xsec").GetVal(), 8)
#sm_xsec = 0.18786147 * 1.195 / 4.7
sm_file.Close()

print("SM acceptance:", sm_acc)
print("SM cross section:", sm_xsec)

rhs = (1 + delta * norm.ppf(0.95)) * sm_acc * sm_xsec

print("95% CL UL:", (1 + delta * norm.ppf(0.95)) * pred, "* 10^-6")


# BSM models
inpath = inpath + "Width/"

wU_nom, wU_all = [], []
acc_nom, acc_all = [], []

print("\nLooping over", model, "files...")
for filename in listdir(inpath.replace("root://cmseos.fnal.gov/", "/eos/uscms")):

    rootfile = TFile.Open(inpath + filename)
    tree = rootfile.Get(selection + "_" + model)

    if tree.GetEntries() < 1:
        continue

    acc = tree.GetEntries("isFiducial") / tree.GetEntries()

    mU = round(rootfile.Get("mU").GetVal(), 1)
    wU = round(rootfile.Get("wU").GetVal(), 8)
    ge = round(rootfile.Get("ge").GetVal(), 6)
    gmu = round(rootfile.Get("gmu").GetVal(), 6)
    xsec = round(rootfile.Get("xsec").GetVal(), 8)

    rootfile.Close()

    if mU != 10:
        continue

    if model == "VectorU":
        ge /= 2
        gmu /= 2

    if ge == 0:
        if selection != "4m":
            continue
        g = gmu

    elif gmu == 0:
        if selection != "4e":
            continue
        g = ge

    elif gmu == ge:
        if selection != "4l":
            continue
        g = ge

    else:
        print("ge != gmu")

    print("mU (GeV):", mU, "\twU (GeV):", wU, "\tg:", g, "\txsec:", xsec, "\tacc:", acc)

    lhs = acc * xsec

    if (mU == 10 and 1e-04 < wU < 1e-3) or (mU == 50 and 5e-2 < wU < 5e-1):
        the_acc = rhs / xsec
        wU_nom.append(wU)
        acc_nom.append(acc * 100)
    else:
        wU_all.append(wU)
        acc_all.append(acc * 100)

mU = 10



##
##  PLOT
##

print("Exclusion for nominal acceptance:", the_acc)

fig = plt.figure(figsize = (8,3))
ax = plt.axes()

#p_nom = ax.plot(wU_nom, acc_nom, 'r*', markersize=9)
p_nom = ax.plot(wU_nom, acc_nom, 'bo')
p_all = ax.plot(wU_all, acc_all, 'ro')

ax.axhline(sm_acc * 100, color='k', linestyle='--')


# Axis labels and limits
ax.set_axisbelow(True)

ax.set_ylabel(r"$A^\mathrm{NP}$ (\%)")
ax.set_xlabel(r"$\Gamma_{\mathrm{U}}$ (GeV)")
plt.xscale('log')

if mU == 10:
#ax.set_ylim(15, 40)
#ax.set_xlim(5e-07, 2e2)
    ax.set_ylim(10, 50)
    ax.set_xlim(5e-07, 2e2)

if mU == 50:
    ax.set_ylim(10, 50)
    ax.set_xlim(2.5e-06, 1e03)



# Legend

r = pat.Rectangle((3.25, 0.065), 11, 0.38, color='w')
ax.add_patch(r)

ax.legend([r"Automatic $\Gamma_\mathrm{U}$",# for $m_\mathrm{U} = 10$ GeV; $g_\mu = 0.01$; $g_\mathrm{e} = 0$",
            r"Manual $\Gamma_\mathrm{U}$",# with $m_\mathrm{U} = 10$ GeV; $g_\mu = 0.01$; $g_\mathrm{e} = 0$"],
            r"SM acceptance"],
        loc = 'upper left', fontsize='x-large', #bbox_to_anchor = (0, 0.6),
        numpoints = 1, frameon = True, facecolor='w', edgecolor='w', framealpha=1, ncol=2)
#plt.grid(which='both')
'''
modelTeX = {"ScalarU":"Scalar U", "VectorU":"Vector U"}
selTeX = {"4l":r"$g_\ell = g_\mathrm{e} = g_\mu$", "4m":r"$g_\mu = 0.01$; $g_\mathrm{e} = 0$",
        "4e":r"$g_\ell = g_\mathrm{e}$; $g_\mu = 0$"}

ax.text(0.04,  0.9,    modelTeX[model] + r"\\ \\" + selTeX[selection] + r"\\ \\ $m_\mathrm{U} = 10$ GeV",
        size = "xx-large",
        weight = "bold", verticalalignment = 'top', transform = ax.transAxes)
'''


# Save
fig.tight_layout()
fig_name = "exclusion_width_M-{:g}.pdf".format(mU)
fig.savefig(fig_name)

print("Wrote plot to", fig_name) 



##
##  SAVE
##
'''
x_arr, y_arr, z_arr = np.array(x), np.array(y), np.array(z)
rhs_arr = np.array(rhs)

outfile = "exclusion_" + model +  "_" + selection + ".npz"
np.savez(outfile, mU=x_arr, g=y_arr, lhs=z_arr, rhs=rhs_arr)

print("\nWrote arrays to", outfile)
'''
