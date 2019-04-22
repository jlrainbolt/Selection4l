from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

#from Cuts2017 import *
#from Cuts2016 import *
from Cuts2012 import *

signalOnly = True
tightOnly = True

if signalOnly:
    weight = "(isSameSign) && (!isDiffFlavor)"
    selection   = ["4l", "4m", "2m2e", "2e2m", "4e"]
else:
    weight = "(isSameSign || isDiffFlavor)"
    selection   = ["4l", "4m", "2m2e", "2e2m", "4e", "3m1e", "1m3e"]

if tightOnly:
    weight = "(nLooseLeptons == 0) && " + weight
else:
    weight = "(nLooseLeptons <= 1) && " + weight


##
##  SAMPLE INFO
##

selTeX      = { "4l":r"4\Pell", "4m":r"4\PGm",  "4e":r"4\Pe",   "2m2e":r"2\PGm 2\Pe",
                "3m1e":r"3\PGm\Pe",             "1m3e":r"\PGm 3\Pe"
                }
selDef      = { "4l":"4L",  "4m":"4M",  "4e":"4E",  "2m2e":"2M2E",  "3m1e":"3M1E",  "1m3e":"1M3E"}
T = np.dtype([(sel, 'f4') for sel in selection])


##
##  DATA
##

inPath = EOS_PATH + "/Selected/" + YEAR_STR + "/"
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
    if sel == "4l":
        continue
    elif sel in ["mumu", "4m", "2m2e", "3m1e"]:
        tree = muFile.Get(sel + "_" + MU_SUFF)
    elif sel in ["ee", "4e", "2e2m", "1m3e"]:
        tree = elFile.Get(sel + "_" + EL_SUFF)

    data[sel] = tree.GetEntries(weight)

muFile.Close()
elFile.Close()



##
##  MONTE CARLO
##

mc_arr, mc_unc_arr = np.zeros(N_MC, dtype=T), np.zeros(N_MC, dtype=T)
mc, mc_unc = {}, {}
row = 0

# Loop over all samples
for suff in MC_SUFF:
    inName = prefix + "_" + suff + ".root"
    inFile = TFile.Open(inPath + inName)
    print("Opened", inPath + inName)

    # Get histograms
    for sel in selection:
        if sel == "4l":
            continue
        elif sel in ["mumu", "4m", "2m2e", "3m1e"]:
            lumi = MUON_TRIG_LUMI
        elif sel in ["ee", "4e", "2e2m", "1m3e"]:
            lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF
        
        sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

        hist = TH1D("hist", "", 1, 0, 2)
        hist.Sumw2()

        tree = inFile.Get(sel + "_" + suff)
        if suff == "zjets_m-50":
            tree.Draw("1>>hist", weight + " * weight/trigWeight", "goff")
        else:
            tree.Draw("1>>hist", weight + " * weight/trigWeight/qtWeight", "goff")

        mc_arr[row][sel] = sf * hist.GetBinContent(1);
        mc_unc_arr[row][sel] = sf * hist.GetBinError(1);

        hist.Delete()

    mc[suff] = mc_arr[row]
    mc_unc[suff] = mc_unc_arr[row]
    row = row + 1
    inFile.Close()



##
##  ADD CHANNELS
##

for sample in [data, mc_arr]:
    sample['2m2e']  = sample['2m2e'] + sample['2e2m']
    sample['4l']    = sample['4m'] + sample['2m2e'] + sample['4e']
    if not signalOnly:
        sample['4l']    = sample['4l'] + sample['3m1e'] + sample['1m3e']
    sample['2e2m']  = 0

# Handle uncertainty
mc_unc_arr['2m2e']  = np.sqrt(mc_unc_arr['2m2e'] ** 2 + mc_unc_arr['2e2m'] ** 2)
mc_unc_arr['4l']    = mc_unc_arr['4m'] ** 2 + mc_unc_arr['2m2e'] ** 2 + mc_unc_arr['4e'] ** 2
if not signalOnly:
    mc_unc_arr['4l'] = mc_unc_arr['4l'] + mc_unc_arr['3m1e'] ** 2 + mc_unc_arr['1m3e'] ** 2
mc_unc_arr['4l'] = np.sqrt(mc_unc_arr['4l'])
mc_unc_arr['2e2m']  = 0



##
##  ADD SAMPLES
##

# Take average of ttbar (inclusive) and tt_2l2nu
if (YEAR_STR != "2012"):
    for sel in selection:
        mc['ttbar'][sel] = (mc['ttbar'][sel] + mc['tt_2l2nu'][sel])
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

    if (npt[sel] == 0):
        sf[sel]     = 0
        sf_unc[sel] = 0
    else:
        sf[sel]     = diff[sel] / npt[sel]
        sf_unc[sel] = np.abs(sf[sel]) * np.sqrt(diff_unc[sel] / diff[sel] ** 2 + npt_unc[sel] / npt[sel] ** 2)

    bg_unc[sel]     = np.sqrt(bg_unc[sel])
    npt_unc[sel]    = np.sqrt(npt_unc[sel])
    diff_unc[sel]   = np.sqrt(diff_unc[sel])



##
##  WRITE TEX FILES
##

if signalOnly:
    selection   = ["4l", "4m", "2m2e", "4e"]
    fileName = "Background" + YEAR_STR + ".tex"
else:
    selection   = ["4l", "4m", "3m1e", "2m2e", "1m3e", "4e"]    # removed 2e2m
    fileName = "Background" + YEAR_STR + "_all.tex"
if not tightOnly:
    fileName = "Background" + YEAR_STR + "_loose.tex"
fmt = '%.2f'

f = open(fileName, "w")

if signalOnly:
    f.write(r"\begin{tabular}{lll r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r}" + "\n")
else:
    f.write(r"\begin{tabular}{lll r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Type}")
for sel in selection:
    f.write(r" & \multicolumn{2}{l}{$N_{" + selTeX[sel] + r"}$ (events)}")
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

f.write(r"\addlinespace\addlinespace" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Exp. $\mapsto$ Obs.~nonprompt}")
for sel in selection:
    f.write(r" & " + fmt % np.squeeze(sf[sel]) + r" & " + fmt % np.squeeze(sf_unc[sel]))
f.write(r" \\" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
