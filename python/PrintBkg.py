from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

#from Cuts2017 import *
from Cuts2016 import *



##
##  SAMPLE INFO
##

selection   = ["4l", "4m", "2m2e", "2e2m", "4e", "3m1e", "1m3e"]
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

    data[sel] = tree.GetEntries("(isSameSign || isDiffFlavor)")

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
        tree = inFile.Get(sel + "_" + suff)
        if suff == "zjets_m-50":
            tree.Draw("1>>hist", "(isSameSign || isDiffFlavor) * weight/trigWeight", "goff")
        else:
            tree.Draw("1>>hist", "(isSameSign || isDiffFlavor) * weight/trigWeight/qtWeight", "goff")

        mc_arr[row][sel] = sf * hist.Integral()
        mc_unc_arr[row][sel] = sf * np.sqrt(tree.GetEntries())

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
    sample['4l']    = sample['4m'] + sample['2m2e'] + sample['4e'] + sample['3m1e'] + sample['1m3e']
    sample['2e2m']  = 0

# Handle uncertainty
mc_unc_arr['2m2e']  = np.sqrt(mc_unc_arr['2m2e'] ** 2 + mc_unc_arr['2e2m'] ** 2)
mc_unc_arr['4l']    = np.sqrt(mc_unc_arr['4m']**2 + mc_unc_arr['2m2e']**2 + mc_unc_arr['4e']**2
                                + mc_unc_arr['3m1e']**2 + mc_unc_arr['1m3e']**2)
mc_unc_arr['2e2m']  = 0



##
##  ADD SAMPLES
##

# Get events from signal process
sig, sig_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in ["4l", "4m", "2m2e", "4e", "3m1e", "1m3e"]:
    sig[sel]        = mc['zz_4l'][sel]
    sig_unc[sel]    = mc_unc['zz_4l'][sel]

# Get total expected and background events
exp, exp_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bg, bg_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in selection:
    exp[sel]        = np.sum(mc_arr[sel])
    exp_unc[sel]    = np.sqrt(np.sum(mc_unc_arr[sel] ** 2))
    bg[sel]         = exp[sel] - sig[sel]
    bg_unc[sel]     = np.sqrt(np.sum(mc_unc_arr[sel] ** 2) - sig_unc[sel] ** 2)



##
##  WRITE TEX FILES
##

selection   = ["4l", "4m", "3m1e", "2m2e", "1m3e", "4e"]    # removed 2e2m
fmt = '%.2f'

fileName = "Background" + YEAR_STR + ".tex"
f = open(fileName, "w")


f.write(r"\begin{tabular}{lll r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Type}")
for sel in selection:
    f.write(r" & \multicolumn{2}{l}{$N_{" + selTeX[sel] + r"}$ (events)}")
f.write(r" \\" + "\n")

f.write(r"\midrule" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Observed}")
for sel in selection:
    f.write(r" & \multicolumn{2}{l}{\num{" + '%.0f' % np.squeeze(data[sel]) + "}}")
f.write(r" \\" + "\n")

f.write(r"\addlinespace\addlinespace" + "\n")

f.write("\t" + r"\multicolumn{3}{l}{Expected}")
for sel in selection:
    f.write(r" & \num{" + fmt % np.squeeze(exp[sel]) + r"} & \num{"
            + fmt % np.squeeze(exp_unc[sel]) + "}")
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t&\t" + r"\multicolumn{2}{l}{Signal}")
for sel in selection:
    f.write(r" & \num{" + fmt % np.squeeze(sig[sel]) + r"} & \num{"
            + fmt % np.squeeze(sig_unc[sel]) + "}")
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

f.write("\t&\t" + r"\multicolumn{2}{l}{Background}")
for sel in selection:
    f.write(r" & \num{" + fmt % np.squeeze(bg[sel]) + r"} & \num{"
            + fmt % np.squeeze(bg_unc[sel]) + "}")
f.write(r" \\" + "\n")

f.write(r"\addlinespace" + "\n")

for suff in MC_SUFF:
    if suff == "zz_4l":
        continue
    else:
        f.write("\t&\t&\t" + MC_TEX[suff])
        for sel in selection:
            f.write(r" & \num{" + fmt % np.squeeze(mc[suff][sel]) + r"} & \num{"
                    + fmt % np.squeeze(mc_unc[suff][sel]) + "}")
        f.write(r" \\" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
