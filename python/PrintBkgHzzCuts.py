from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

from Cuts2018 import *
#from Cuts2017 import *
#from Cuts2016 import *

selection   = ["4l", "4m", "2m2e", "4e"]

cut = "(nLooseLeptons == 0)"

#cut2 = "(((fabs(z2p4.M() - 91.2) < fabs(z1p4.M() - 91.2)) * (z2p4.M() > 40) * (z2p4.M() < 120) * (z1p4.M() > 12) * (z1p4.M() < 120)) || ((fabs(z2p4.M() - 91.2) > fabs(z1p4.M() - 91.2)) * (z1p4.M() > 40) * (z1p4.M() < 120) * (z2p4.M() > 12) * (z2p4.M() < 120))) * (zzp4.M() > 70)"
#cutstr = "HZZCuts"

#cut2 = "(((fabs(z2p4.M() - 91.2) < fabs(z1p4.M() - 91.2)) * (z2p4.M() > 40) * (z2p4.M() < 120) * (z1p4.M() > 12) * (z1p4.M() < 120)) || ((fabs(z2p4.M() - 91.2) > fabs(z1p4.M() - 91.2)) * (z1p4.M() > 40) * (z1p4.M() < 120) * (z2p4.M() > 12) * (z2p4.M() < 120))) * (zzp4.M() > 118) * (zzp4.M() < 130)"
#cutstr = "HZZCuts_118-130"

cut2 = "(((fabs(z2p4.M() - 91.2) < fabs(z1p4.M() - 91.2)) * (z2p4.M() > 40) * (z2p4.M() < 120) * (z1p4.M() > 12) * (z1p4.M() < 120)) || ((fabs(z2p4.M() - 91.2) > fabs(z1p4.M() - 91.2)) * (z1p4.M() > 40) * (z1p4.M() < 120) * (z2p4.M() > 12) * (z2p4.M() < 120))) * (zzp4.M() > 80) * (zzp4.M() < 100)"
cutstr = "HZZCuts_80-100"


##
##  SAMPLE INFO
##

selTeX      = { "4l":r"\fL",    "4m":r"\fM",    "4e":r"\fE",    "2m2e":r"\tMtE" }
T = np.dtype([(sel, 'f4') for sel in selection])


##
##  DATA
##

inPath = EOS_PATH + "/Extended/" + YEAR_STR + "_v1/"
prefix = "background"

# Muon file
muName = prefix + "_" + MU_SUFF + ".root"
muFile = TFile.Open(inPath + muName)
print("Opened", inPath + muName)

# Electron file
elName = prefix + "_" + EL_SUFF + ".root"
elFile = TFile.Open(inPath + elName)
print("Opened", inPath + elName)

# Get yields
data = np.zeros(1, dtype=T)

for sel in selection:
    muTree = muFile.Get(sel + "_" + MU_SUFF)
    elTree = elFile.Get(sel + "_" + EL_SUFF)
    data[sel] = muTree.GetEntries(cut + " * " + cut2) + elTree.GetEntries(cut + " * " + cut2)

muFile.Close()
elFile.Close()



##
##  MONTE CARLO
##

mc_arr, mc_unc_arr = np.zeros(N_MC, dtype=T), np.zeros(N_MC, dtype=T)
mc, mc_unc = {}, {}
row = 0
weight = "weight"

# Loop over all samples
for suff in MC_SUFF:
    inName = prefix + "_" + suff + ".root"
    inFile = TFile.Open(inPath + inName)
    print("Opened", inPath + inName)

    # Get histograms
    for sel in selection:
        sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]

        hist = TH1D("hist", "", 1, 0, 2)
        hist.Sumw2()

        tree = inFile.Get(sel + "_" + suff)
        tree.Draw("1>>hist", cut + " * " + weight + " * " + cut2, "goff")
        mc_arr[row][sel] = sf * hist.Integral()
        tree.Draw("1>>hist", cut + " * " + weight + " * " + weight + " * " + cut2, "goff")            

        mc_unc_arr[row][sel] = sf * np.sqrt(hist.GetBinContent(1))

        mc_arr[row][sel] = np.abs(mc_arr[row][sel])

        hist.Delete()

    mc[suff] = mc_arr[row]
    mc_unc[suff] = mc_unc_arr[row]
    row = row + 1
    inFile.Close()



##
##  ADD SAMPLES
##

# Take average of ttbar (inclusive) and tt_2l2nu
if (YEAR_STR != "2012"):
    for sel in selection:
        mc['ttbar'][sel] = mc['ttbar'][sel] + mc['tt_2l2nu'][sel]
        mc_unc['ttbar'][sel] = np.sqrt(mc_unc['ttbar'][sel] ** 2 + mc_unc['tt_2l2nu'][sel] ** 2)
        mc['tt_2l2nu'][sel] = 0
        mc_unc['tt_2l2nu'][sel] = 0

# Get total expected and background events
exp, exp_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
npt, npt_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bg, bg_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
diff, diff_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
sf, sf_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    exp[sel]        = np.sum(mc_arr[sel])
    exp_unc[sel]    = np.sqrt(np.sum(mc_unc_arr[sel] ** 2))

    npt[sel]        = mc['zjets_m-50'][sel] + mc['ttbar'][sel]
    npt_unc[sel]    = mc_unc['zjets_m-50'][sel] ** 2 + mc_unc['ttbar'][sel] ** 2

    bg[sel]         = exp[sel] - npt[sel]
    bg_unc[sel]     = np.sum(mc_unc_arr[sel] ** 2) - npt_unc[sel]

    diff[sel]       = data[sel] - bg[sel]
    diff_unc[sel]   = bg_unc[sel] + data[sel]

    if diff[sel] < 0:
        diff_unc[sel] = diff_unc[sel] - diff[sel]
        diff[sel] = 0

    if npt[sel] == 0 or diff[sel] == 0:
        sf[sel]     = 0
        sf_unc[sel] = 0
    else:
        sf[sel]     = diff[sel] / npt[sel]
        sf_unc[sel] = np.abs(sf[sel]) * np.sqrt(1 / diff[sel] + npt_unc[sel] / npt[sel] ** 2)

    bg_unc[sel]     = np.sqrt(bg_unc[sel])
    npt_unc[sel]    = np.sqrt(npt_unc[sel])
    diff_unc[sel]   = np.sqrt(diff_unc[sel])



##
##  WRITE TEX FILES
##

prefix = "Background" + YEAR_STR + cutstr

fileName = prefix + ".tex"

fmt = '%.2f'

f = open(fileName, "w")

f.write(r"\begin{tabular}{lll r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Type}")
for sel in selection:
    f.write(r" & \multicolumn{2}{l}{$N_{" + selTeX[sel] + r"}$}")
f.write(r" \\" + "\n")

f.write(r"\midrule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Observed}")
for sel in selection:
    f.write(r" & \multicolumn{2}{l}{" + '%.0f' % np.squeeze(data[sel]) + "}")
f.write(r" \\" + "\n")

f.write(r"\addlinespace\addlinespace" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Expected}")
for sel in selection:
    f.write(r" & " + fmt % np.squeeze(exp[sel]) + r" & " + fmt % np.squeeze(exp_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t&\t" + r"\multicolumn{2}{l}{Nonprompt MC}")
for sel in selection:
    f.write(r" & " + fmt % np.squeeze(npt[sel]) + r" & " + fmt % np.squeeze(npt_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

for suff in ["zjets_m-50", "ttbar"]:
    f.write("\t&\t&\t" + MC_TEX[suff])
    for sel in selection:
        f.write(r" & " + fmt % np.squeeze(mc[suff][sel]) + r" & "
                + fmt % np.squeeze(mc_unc[suff][sel]))
    f.write(r" \\" + "\n")

f.write(r"\addlinespace\addlinespace" + "\n")

f.write("\t&\t" + r"\multicolumn{2}{l}{Prompt MC}")
for sel in selection:
    f.write(r" & " + fmt % np.squeeze(bg[sel]) + r" & " + fmt % np.squeeze(bg_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

for suff in MC_SUFF:
    if suff in ["zjets_m-50", "tt_2l2nu", "ttbar"]:
        continue
    else:
        f.write("\t&\t&\t" + MC_TEX[suff])
        for sel in selection:
#           f.write(r" & " + fmt % np.squeeze(np.abs(mc[suff][sel])) + r" & "
            f.write(r" & " + fmt % np.squeeze(mc[suff][sel]) + r" & "
                    + fmt % np.squeeze(mc_unc[suff][sel]))
        f.write(r" \\" + "\n")

f.write(r"\midrule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Observed $-$ Prompt MC}")
for sel in selection:
    f.write(r" & " + fmt % np.squeeze(diff[sel]) + r" & " + fmt % np.squeeze(diff_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{$N^\text{obs}_\text{npr} / N^\text{exp}_\text{npr}$}")
for sel in selection:
    f.write(r" & " + fmt % np.squeeze(sf[sel]) + r" & " + fmt % np.squeeze(sf_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)



##
##  SAVE
##

print("")
outfile = "nonprompt" + YEAR_STR + cutstr + ".npz"
np.savez(outfile, npt=diff, npt_unc=diff_unc)
print("Wrote nonprompt arrays to", outfile)
