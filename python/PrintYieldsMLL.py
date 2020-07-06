from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

#from Cuts2018 import *
#from Cuts2017 import *
from Cuts2016 import *



##
##  SAMPLE INFO
##

#selection   = ["4l", "4m", "2m2e", "4e"]
selection   = ["mumu", "ee", "4l", "4m", "2m2e", "4e"]
selTeX      = {"mumu":r"\MM", "ee":r"\EE", "4l":r"\fL", "4m":r"\fM", "2m2e":r"\tMtE", "4e":r"\fE"}
T = np.dtype([(sel, 'f4') for sel in selection])

MC_SUFF = MC_SUFF_MLL



##
##  DATA
##

inPath = EOS_PATH + "/Extended/" + YEAR_STR + "_mll/"
prefix = "selected"

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
    data[sel] = muTree.GetEntries() + elTree.GetEntries()

muFile.Close()
elFile.Close()



##
##  MONTE CARLO
##

mc_arr, mc_stat_arr = np.zeros(N_MC, dtype=T), np.zeros(N_MC, dtype=T)
mc_sys_arr, mc_unc_arr = np.zeros(N_MC, dtype=T), np.zeros(N_MC, dtype=T)
sig, sig_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
sig_stat, sig_sys = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
mc, mc_stat, mc_sys, mc_unc = {}, {}, {}, {}
row = 0

# Loop over all samples
for suff in MC_SUFF:
    inName = prefix + "_" + suff + ".root"
    inFile = TFile.Open(inPath + inName)
    print("Opened", inPath + inName)

    # Get histograms
    for sel in selection:
        if suff in ["zjets_m-50", "ttbar", "tt_2l2nu"] and sel in ["4l", "4m", "2m2e", "4e"]:
            continue
        
        sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]

        cut = ""
        weight = "weight"

        hist = TH1D("hist", "", 1, 0, 2)
        tree = inFile.Get(sel + "_" + suff)

        
        # Get signal
        if (suff == "zz_4l_m-1" and sel in ["4l", "4m", "2m2e", "4e"]) or (suff == "zjets_m-50" and sel in ["mumu", "ee"]):
            tree.Draw("1>>hist", "!hasTauDecay * " + weight, "goff")
            sig[sel] = sf * hist.Integral()
            tree.Draw("1>>hist", "!hasTauDecay * " + weight + " * " + weight, "goff")
            sig_stat[sel] = sf * np.sqrt(hist.Integral())
            sig_sys[sel] = MC_UNC[suff] * sig[sel]
            sig_unc[sel] = np.sqrt(sig_stat[sel] ** 2 + sig_sys[sel] ** 2)
            cut = "hasTauDecay"
            weight = cut + " * " + weight

        tree.Draw("1>>hist", weight, "goff")

        mc_arr[row][sel] = sf * hist.Integral()

        tree.Draw("1>>hist", weight + " * " + weight, "goff")
        mc_stat_arr[row][sel] = sf * np.sqrt(hist.Integral())
        mc_sys_arr[row][sel] = MC_UNC[suff] * mc_arr[row][sel]
        mc_unc_arr[row][sel] = np.sqrt(mc_stat_arr[row][sel] ** 2 + mc_sys_arr[row][sel] ** 2)

        mc_arr[row][sel] = np.abs(mc_arr[row][sel])

        hist.Delete()

    mc[suff] = mc_arr[row]
    mc_stat[suff] = mc_stat_arr[row]
    mc_sys[suff] = mc_sys_arr[row]
    mc_unc[suff] = mc_unc_arr[row]
    row = row + 1
    inFile.Close()



##
##  ADD SAMPLES
##

# Take average of ttbar (inclusive) and tt_2l2nu

if YEAR_STR != "2012":
    for sel in ["mumu", "ee"]:
        mc['ttbar'][sel] = mc['ttbar'][sel] + mc['tt_2l2nu'][sel]
        mc_stat['ttbar'][sel] = np.sqrt(mc_stat['ttbar'][sel] ** 2 + mc_stat['tt_2l2nu'][sel] ** 2)
        mc_sys['ttbar'][sel] = mc['ttbar'][sel] * MC_UNC['ttbar']
        mc_unc['ttbar'][sel] = np.sqrt(mc_stat['ttbar'][sel] ** 2 + mc_sys['ttbar'][sel] ** 2)

        mc['tt_2l2nu'][sel] = 0
        mc_stat['tt_2l2nu'][sel] = 0
        mc_sys['tt_2l2nu'][sel] = 0
        mc_unc['tt_2l2nu'][sel] = 0



# Get nonprompt background
infile = "nonprompt" + YEAR_STR + "MLL.npz"
npzfile = np.load(infile)
npt_, npt_stat_ = npzfile['npt'], npzfile['npt_unc']
npt, npt_stat, npt_sys = np.zeros(1, dtype=T), np.zeros(1, dtype=T), np.zeros(1, dtype=T)
npt_unc = np.zeros(1, dtype=T)

for sel in ["4l", "4m", "2m2e", "4e"]:
    npt[sel] = npt_[sel]
    npt_stat[sel] = npt_stat_[sel]
    npt_sys[sel] = DELTA_LAMBDA * npt[sel]
    npt_unc[sel] = np.sqrt(npt_stat[sel] **2 + npt_sys[sel] ** 2)


# Get total expected and background events
exp, exp_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bg, bg_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bg_stat, bg_sys = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
pur, pur_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
sf, sf_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    exp[sel]        = np.sum(mc_arr[sel]) + sig[sel] + npt[sel]
    exp_unc[sel]    = np.sqrt(np.sum(mc_unc_arr[sel] ** 2) + sig_unc[sel] ** 2 + npt_unc[sel] ** 2)
    bg[sel]         = np.sum(mc_arr[sel]) + npt[sel]
    bg_stat[sel]    = np.sqrt(np.sum(mc_stat_arr[sel] ** 2) + npt_stat[sel] ** 2)
    bg_sys[sel]     = np.sqrt(np.sum(mc_sys_arr[sel] ** 2) + npt_sys[sel] ** 2)
    bg_unc[sel]     = np.sqrt(bg_stat[sel] ** 2 + bg_sys[sel] ** 2)
    pur[sel]        = sig[sel] / exp[sel] * 100
    pur_unc[sel]    = pur[sel] * np.sqrt(sig_unc[sel] / sig[sel] ** 2 + exp_unc[sel] / exp[sel] ** 2)
    sf[sel]         = data[sel] / exp[sel] * 100
    sf_unc[sel]     = sf[sel] * np.sqrt(1 / data[sel] + exp_unc[sel] / exp[sel] ** 2)



##
##  WRITE TEX FILES
##

fileName = "Yield" + YEAR_STR + "MLL.tex"
f = open(fileName, "w")

f.write(r"\begin{tabular}{lll r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Type}")
for sel in selection:
    f.write(r" & \multicolumn{2}{l}{$N_{" + selTeX[sel] + r"}$}")
f.write(r" \\" + "\n")

f.write(r"\midrule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Observed}")
for sel in selection:
    f.write(r" & \multicolumn{2}{l}{" + '%i' % np.squeeze(data[sel]) + "}")
f.write(r" \\" + "\n")

f.write(r"\addlinespace\addlinespace" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Expected}")
for sel in selection:
    if sel in ["mumu", "ee"]:
        fmt = '%i'
    else:
        fmt = '%.2f'
    f.write(" & " + fmt % np.squeeze(exp[sel]) + " & " + fmt % np.squeeze(exp_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t&\t" + r"\multicolumn{2}{l}{Signal}")
for sel in selection:
    if sel in ["mumu", "ee"]:
        fmt = '%i'
    else:
        fmt = '%.2f'
    f.write(" & " + fmt % np.squeeze(sig[sel]) + " & " + fmt % np.squeeze(sig_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t&\t" + r"\multicolumn{2}{l}{Background}")
for sel in selection:
    if sel in ["mumu", "ee"]:
        fmt = '%i'
    else:
        fmt = '%.2f'
    f.write(" & " + fmt % np.squeeze(bg[sel]) + " & " + fmt % np.squeeze(bg_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t&\t&\t" + r"Nonprompt & \multicolumn{4}{c}{}")
for sel in ["4l", "4m", "2m2e", "4e"]:
    if sel in ["mumu", "ee"]:
        fmt = '%i'
    else:
        fmt = '%.2f'
    f.write(" & " + fmt % np.squeeze(npt[sel]) + " & " + fmt % np.squeeze(npt_unc[sel]))
f.write(r" \\" + "\n")

for suff in MC_SUFF:
    if suff == "tt_2l2nu":
        continue
    if suff == "zz_4l_m-1" and sel in ["4l", "4m", "2m2e", "4e"]:
        f.write("\t&\t&\t" + r"$\ZZtofL$")
    elif suff == "zjets_m-50" and sel in ["mumu", "ee"]:
        f.write("\t&\t&\t" + r"$\ZtoTT$")
    else:
        f.write("\t&\t&\t" + MC_TEX[suff])

    for sel in selection:
        if sel in ["mumu", "ee"]:
            fmt = '%i'
        else:
            fmt = '%.2f'
        if suff in ["zjets_m-50", "ttbar"] and sel in ["4l", "4m", "2m2e", "4e"]:
            continue
        else:
            f.write(" & " + fmt % np.squeeze(mc[suff][sel]) + " & "
                    + fmt % np.squeeze(mc_unc[suff][sel]))
    f.write(r" \\" + "\n")

f.write(r"\midrule" + "\n")

fmt = '%.2f'

f.write("\t" + r"\multicolumn{3}{l}{Signal purity (\%)}")
for sel in selection:
    f.write(" & " + fmt % np.squeeze(pur[sel]) + " & " + fmt % np.squeeze(pur_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{$N^\text{obs} / N^\text{exp}$ (\%)}")
for sel in selection:
    f.write(" & " + fmt % np.squeeze(sf[sel]) + " & " + fmt % np.squeeze(sf_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)



##
##  SAVE
##

print("")
outfile = "yields" + YEAR_STR + "MLL.npz"
np.savez(outfile, data=data, exp=exp, exp_unc=exp_unc, sig=sig, sig_unc=sig_unc, sig_stat=sig_stat, sig_sys=sig_sys, 
        bg=bg, bg_unc=bg_unc, bg_stat=bg_stat, bg_sys=bg_sys, mc=mc_arr, mc_stat=mc_stat_arr,
        mc_sys=mc_sys_arr, mc_unc=mc_unc_arr, npt=npt, npt_stat=npt_stat, npt_sys=npt_sys, npt_unc=npt_unc)

print("Wrote arrays to", outfile)
