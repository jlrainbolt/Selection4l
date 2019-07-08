from __future__ import print_function
from __future__ import division


import sys
import numpy as np
from os import listdir
from scipy.stats import norm
from ROOT import TFile, TTree, TParameter

from PlotUtils import *

from Cuts2018 import EOS_PATH



##
##  OPTIONS
##

model = sys.argv[1]
selection = sys.argv[2]


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
inpath = inpath + model + "/"

mU_exc, g_exc = [], []
mU_all, g_all = [], []

x, y, z = [], [], []

print("\nLooping over", model, "files...")
for filename in listdir(inpath.replace("root://cmseos.fnal.gov/", "/eos/uscms")):

    rootfile = TFile.Open(inpath + filename)
    tree = rootfile.Get(selection + "_" + model)

    if tree.GetEntries() < 1:
        continue

    acc = tree.GetEntries("isFiducial") / tree.GetEntries()

    mU = round(rootfile.Get("mU").GetVal(), 1)
    ge = round(rootfile.Get("ge").GetVal(), 6)
    gmu = round(rootfile.Get("gmu").GetVal(), 6)
    xsec = round(rootfile.Get("xsec").GetVal(), 8)

    rootfile.Close()

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

    print("mU:",  mU, "\tg:", g, "\txsec:", xsec, "\tacc:", acc)

    lhs = acc * xsec

    if lhs > rhs:
        mU_exc.append(mU)
        g_exc.append(g)
    else:
        mU_all.append(mU)
        g_all.append(g)

    x.append(mU)
    y.append(g)
    z.append(lhs)



##
##  PLOT
##

fig = plt.figure()
ax = plt.axes()

p_exc = ax.plot(mU_exc, g_exc, 'ro')
p_all = ax.plot(mU_all, g_all, 'bo')


# Axis labels and limits
plt.yscale('log')
ax.set_ylabel(r"$g_{\ell}$")
ax.set_ylim(0.005, 0.5)

plt.xscale('log')
ax.set_xlabel(r"$m_{\mathrm{U}}$")
ax.set_xlim(3, 100)


# Legend
ax.legend(["Excluded", "Allowed"], loc = 'upper left', numpoints = 1, frameon = False)
plt.grid(which='both')


# Save
fig_name = "exclusion_" + model + "_" + selection + ".pdf"
fig.savefig(fig_name)

print("Wrote plot to", fig_name) 



##
##  SAVE
##

x_arr, y_arr, z_arr = np.array(x), np.array(y), np.array(z)
rhs_arr = np.array(rhs)

outfile = "exclusion_" + model +  "_" + selection + ".npz"
np.savez(outfile, mU=x_arr, g=y_arr, lhs=z_arr, rhs=rhs_arr)

print("\nWrote arrays to", outfile)
