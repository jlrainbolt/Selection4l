from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

#from Cuts2018 import *
#from Cuts2017 import *
from Cuts2016 import *
#from Cuts2012 import *



##
##  SAMPLE INFO
##

selection   = ["mumu", "ee", "4l", "4m", "2m2e", "2e2m", "4e"]
selTeX      = { "mumu":r"\MM",  "ee":r"\EE",
                "4l":r"4\Pell", "4m":r"4\PGm",  "4e":r"4\Pe",   "2m2e":r"2\PGm 2\Pe"
                }
selDef      = { "mumu":"MM",    "ee":"EE",  "4l":"4L",  "4m":"4M",  "4e":"4E",  "2m2e":"2M2E"   }
T = np.dtype([(sel, 'f4') for sel in selection])

lumiUp = False
lumiDown = False


##
##  DATA
##

inPath = EOS_PATH + "/Selected/" + YEAR_STR + "/"
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
    if sel == "4l":
        continue
    elif sel in ["mumu", "4m", "2m2e"]:
        tree = muFile.Get(sel + "_" + MU_SUFF)
    elif sel in ["ee", "4e", "2e2m"]:
        tree = elFile.Get(sel + "_" + EL_SUFF)

    data[sel] = tree.GetEntries()

muFile.Close()
elFile.Close()



##
##  MONTE CARLO
##

mc_arr, mc_unc_arr = np.zeros(N_MC, dtype=T), np.zeros(N_MC, dtype=T)
sig, sig_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
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
        elif sel in ["mumu", "4m", "2m2e"]:
            lumi = MUON_TRIG_LUMI
        elif sel in ["ee", "4e", "2e2m"]:
            lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

        if lumiUp:
            lumi = lumi * (1 + LUMI_UNC)
        elif lumiDown:
            lumi = lumi * (1 - LUMI_UNC)

        if suff in ["zjets_m-50", "ttbar", "tt_2l2nu"] and sel in ["4m", "2m2e", "2e2m", "4e"]:
            continue
        
        sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

        if suff == "zjets_m-50":
            weight = "weight/trigWeight"
        else:
            weight = "weight/trigWeight/qtWeight"

        cut = ""

        hist = TH1D("hist", "", 1, 0, 2)
        tree = inFile.Get(sel + "_" + suff)

        
        # Get signal
        if (suff == "zz_4l" and sel in ["4m", "2m2e", "2e2m", "4e"]) or (suff == "zjets_m-50" and sel in ["mumu", "ee"]):
            tree.Draw("1>>hist", "!hasTauDecay * " + weight, "goff")
            sig[sel] = sf * hist.Integral()
            sig_unc[sel] = sf * np.sqrt(tree.GetEntries("!hasTauDecay"))
            cut = "hasTauDecay"
            weight = cut + " * " + weight

            
        tree.Draw("1>>hist", weight, "goff")

        mc_arr[row][sel] = sf * hist.Integral()
        mc_unc_arr[row][sel] = sf * np.sqrt(tree.GetEntries(cut))

        hist.Delete()

    mc[suff] = mc_arr[row]
    mc_unc[suff] = mc_unc_arr[row]
    row = row + 1
    inFile.Close()



##
##  ADD CHANNELS
##

for sample in [data, mc_arr, sig]:
    sample['2m2e']  = sample['2m2e'] + sample['2e2m']
    sample['4l']    = sample['4m'] + sample['2m2e'] + sample['4e']
    sample['2e2m']  = 0

# Handle uncertainty
for sample_unc in [mc_unc_arr, sig_unc]:
    sample_unc['2m2e']  = np.sqrt(sample_unc['2m2e'] ** 2 + sample_unc['2e2m'] ** 2)
    sample_unc['4l']    = np.sqrt(sample_unc['4m']**2 + sample_unc['2m2e']**2 + sample_unc['4e']**2)
    sample_unc['2e2m']  = 0



##
##  ADD SAMPLES
##

# Take average of ttbar (inclusive) and tt_2l2nu
if YEAR_STR != "2012":
    for sel in ["mumu", "ee"]:
        mc['ttbar'][sel] = (mc['ttbar'][sel] + mc['tt_2l2nu'][sel])
        mc_unc['ttbar'][sel] = np.sqrt(mc_unc['ttbar'][sel] ** 2 + mc_unc['tt_2l2nu'][sel] ** 2)
        mc['tt_2l2nu'][sel] = 0
        mc_unc['tt_2l2nu'][sel] = 0

# Get total expected and background events
if lumiUp:
    npt = npt_lumiUp
elif lumiDown:
    npt = npt_lumiDn

exp, exp_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bg, bg_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in selection:
    exp[sel]        = np.sum(mc_arr[sel]) + sig[sel] + npt[sel]
    exp_unc[sel]    = np.sqrt(np.sum(mc_unc_arr[sel] ** 2) + sig_unc[sel] ** 2 + npt_unc[sel] ** 2)
    bg[sel]         = np.sum(mc_arr[sel]) + npt[sel]
    bg_unc[sel]     = np.sqrt(np.sum(mc_unc_arr[sel] ** 2) + npt_unc[sel] ** 2)



##
##  WRITE TEX FILES
##

for sel in selection:
    if sel in ["mumu", "ee"]:
        fmt = '%.0f'
    elif sel in ["4l", "4m", "2m2e", "4e"]:
        fmt = '%.2f'
    else:
        continue

    fileName = "Yield" + selDef[sel] + YEAR_STR + ".tex"
    if lumiUp:
        fileName = "Yield" + selDef[sel] + YEAR_STR + "_lumiUp.tex"
    elif lumiDown:
        fileName = "Yield" + selDef[sel] + YEAR_STR + "_lumiDown.tex"
    f = open(fileName, "w")


    f.write(r"\begin{tabular}{l l l r c r}" + "\n")
    f.write(r"\toprule" + "\n")
    f.write("\t" + r"\multicolumn{3}{l}{Type} & \multicolumn{3}{l}{$N_{" 
            + selTeX[sel] + r"}$ (events)} \\" + "\n")
    f.write(r"\midrule" + "\n")
    f.write("\t" + r"\multicolumn{3}{l}{Observed} & \num{" + '%.0f' % np.squeeze(data[sel])
            + r"} \\" + "\n")
    f.write(r"\addlinespace\addlinespace" + "\n")
    f.write("\t" + r"\multicolumn{3}{l}{Expected} & \num{" + fmt % np.squeeze(exp[sel])
            + r"} & $\pm$ & \num{" + fmt % np.squeeze(exp_unc[sel]) + r"} \\" + "\n")
    f.write(r"\addlinespace" + "\n")
    f.write("\t&\t" + r"\multicolumn{2}{l}{Signal} & \num{" + fmt % np.squeeze(sig[sel])
            + r"} & $\pm$ & \num{" + fmt % np.squeeze(sig_unc[sel]) + r"} \\" + "\n")
    f.write(r"\addlinespace" + "\n")
    f.write("\t&\t" + r"\multicolumn{2}{l}{Background} & \num{" + fmt % np.squeeze(bg[sel])
            + r"} & $\pm$ & \num{" + fmt % np.squeeze(bg_unc[sel]) + r"} \\" + "\n")
    f.write(r"\addlinespace" + "\n")

    if sel in ["4l", "4m", "2m2e", "4e"]:
       f.write("\t&\t&\t" + r"Nonprompt & \num{" + fmt % npt[sel]
               + r"} & $\pm$ & \num{" + fmt % npt_unc[sel] + r"} \\" + "\n")

    for suff in MC_SUFF:
#       if suff == "zjets_m-50" and sel in ["mumu", "ee"]:
#           continue
       if suff in ["zjets_m-50", "ttbar"] and sel in ["4l", "4m", "2m2e", "4e"]:
           continue
       else:
           f.write("\t&\t&\t" + MC_TEX[suff] + r" & \num{" + fmt % np.squeeze(mc[suff][sel])
                   + r"} & $\pm$ & \num{" + fmt % np.squeeze(mc_unc[suff][sel]) + r"} \\" + "\n")

    f.write(r"\bottomrule" + "\n")
    f.write(r"\end{tabular}" + "\n")

    f.close()
    print("Wrote table to", fileName)
