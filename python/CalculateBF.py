from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

from Cuts2017 import *



##
##  SAMPLE INFO
##

selection   = ["mumu", "ee", "4l", "4m", "2m2e", "2e2m", "4e"]
selTeX      = { "mumu":r"\MM",  "ee":r"\EE",
                "4l":r"4\Pell", "4m":r"4\PGm",  "4e":r"4\Pe",   "2m2e":r"2\PGm 2\Pe"
                }
selDef      = { "mumu":"MM",    "ee":"EE",  "4l":"4L",  "4m":"4M",  "4e":"4E",  "2m2e":"2M2E"   }
channel     = { "mumu":"3",     "ee":"4",   "4m":"6",   "2m2e":"7", "2e2m":"8", "4e":"9"    }
T = np.dtype([(sel, 'f4') for sel in selection])



##
##  DATA
##

inPath = EOS_PATH + "/Selected/" + YEAR_STR + "_old/"
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
data, data_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    if sel == "4l":
        continue
    elif sel in ["mumu", "4m", "2m2e"]:
        tree = muFile.Get(sel + "_" + MU_SUFF)
    elif sel in ["ee", "4e", "2e2m"]:
        tree = elFile.Get(sel + "_" + EL_SUFF)

    data[sel] = tree.GetEntries()
    data_unc[sel] = np.sqrt(data[sel])

muFile.Close()
elFile.Close()

print("Data")
print(data)
print("")



##
##  MONTE CARLO
##

mc_arr, mc_unc_arr = np.zeros(N_MC, dtype=T), np.zeros(N_MC, dtype=T)
mc, mc_unc = {}, {}
j = 0

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
        
        sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

        hist = TH1D("hist", "", 1, 0, 2)
        tree = inFile.Get(sel + "_" + suff)
#       tree.Draw("1>>hist", "weight/trigWeight", "goff")
        tree.Draw("1>>hist", "weight", "goff")

        mc_arr[j][sel] = sf * hist.Integral()
        mc_unc_arr[j][sel] = sf * np.sqrt(tree.GetEntries())

        hist.Delete()

    mc[suff] = mc_arr[j]
    mc_unc[suff] = mc_unc_arr[j]
    j = j + 1
    inFile.Close()



##
##  PHASE SPACE
##

inPath = EOS_PATH + "/BLT/" + YEAR_STR + "_old/"
prefix = "genHardProc"

# ZZTo4L file
zzName = "gen_zz_4l/" + prefix + "_" + "zz_4l" + ".root"
zzFile = TFile.Open(inPath + zzName)
print("Opened", inPath + zzName)

# DYJetsToLL file
dyName = "gen_zjets_m-50/" + prefix + "_" + "zjets_m-50" + ".root"
dyFile = TFile.Open(inPath + dyName)
print("Opened", inPath + dyName)

# Get yields
ps, ps_unc = np.zeros(2, dtype=T), np.zeros(2, dtype=T)

for sel in selection:
    if sel == "4l":
        continue
    elif sel in ["mumu", "4m", "2m2e"]:
        lumi = MUON_TRIG_LUMI
    elif sel in ["ee", "4e", "2e2m"]:
        lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF
    
    elif sel in ["mumu", "ee"]:
        suff = "zjets_m-50"
        tree = dyFile.Get("tree_" + suff)
    elif sel in ["4m", "4e", "2e2m", "2m2e"]:
        suff = "zz_4l"
        tree = zzFile.Get("tree_" + suff)

    sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

    # JUST USE SELECTED EVENTS HIST
    hist = TH1D("hist", "", 1, 0, 2)
    tree.Draw("1>>hist", "genWeight * (decayChannel == " + channel[sel] + ")", "goff")

    ps[sel] = sf * hist.Integral()
    ps_unc[sel] = sf * np.sqrt(tree.GetEntries())

    hist.Delete()

dyFile.Close()
zzFile.Close()




##
##  ADD CHANNELS
##

for sample in [data, mc_arr, ps]:
    sample['2m2e']  = sample['2m2e'] + sample['2e2m']
    sample['4l']    = sample['4m'] + sample['2m2e'] + sample['4e']
    sample['2e2m']  = 0

# Handle uncertainty
for sample_unc in [data_unc, mc_unc_arr, ps_unc]:
    sample_unc['2m2e']  = np.sqrt(sample_unc['2m2e'] ** 2 + sample_unc['2e2m'] ** 2)
    sample_unc['4l']    = np.sqrt(sample_unc['4m']**2 + sample_unc['2m2e']**2 + sample_unc['4e']**2)
    sample_unc['2e2m']  = 0



##
##  ADD SAMPLES
##

# Get events from signal process
sig, sig_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in ["4l", "4m", "2m2e", "4e"]:
    sig[sel]        = mc['zz_4l'][sel]
    sig_unc[sel]    = mc_unc['zz_4l'][sel]
for sel in ["mumu", "ee"]:
    sig[sel]        = mc['zjets_m-50'][sel]
    sig_unc[sel]    = mc_unc['zjets_m-50'][sel]

# Calculated acceptance * efficiency
axe, axe_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in selection:
    if sel == "2e2m":
        continue
    axe[sel]    = sig['sel'] / ps['sel']
    # FIXME add uncertainty?
print("A x e")
print(axe)
print("")

# Get total expected and background events
exp, exp_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bg, bg_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in selection:
    exp[sel]        = np.sum(mc_arr[sel])
    exp_unc[sel]    = np.sqrt(np.sum(mc_unc_arr[sel] ** 2))
    bg[sel]         = exp[sel] - sig[sel]
    bg_unc[sel]     = np.sqrt(np.sum(mc_unc_arr[sel] ** 2) - sig_unc[sel] ** 2)


'''
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
    f = open(fileName, "w")


    f.write(r"\begin{tabular}{l l l r c r}" + "\n")
    f.write(r"\toprule" + "\n")
    f.write("\t" + r"\multicolumn{3}{l}{Type} & \multicolumn{3}{l}{$N_{" 
            + selTeX[sel] + r"}$ (events)} \\" + "\n")
    f.write(r"\midrule" + "\n")
    f.write("\t" + r"\multicolumn{3}{l}{Observed} & " + '%.0f' % np.squeeze(data[sel]) + " \\\\\n")
    f.write(r"\addlinespace\addlinespace" + "\n")
    f.write("\t" + r"\multicolumn{3}{l}{Expected} & " + fmt % np.squeeze(exp[sel])
            + r" & $\pm$ & " + fmt % np.squeeze(exp_unc[sel]) + " \\\\\n")
    f.write(r"\addlinespace" + "\n")
    f.write("\t&\t" + r"\multicolumn{2}{l}{Signal} & " + fmt % np.squeeze(sig[sel])
            + r" & $\pm$ & " + fmt % np.squeeze(sig_unc[sel]) + " \\\\\n")
    f.write(r"\addlinespace" + "\n")
    f.write("\t&\t" + r"\multicolumn{2}{l}{Background} & " + fmt % np.squeeze(bg[sel])
            + r" & $\pm$ & " + fmt % np.squeeze(bg_unc[sel]) + " \\\\\n")
    f.write(r"\addlinespace" + "\n")

    for suff in MC_SUFF:
        if suff == "zjets_m-50" and sel in ["mumu", "ee"]:
            continue
        elif suff == "zz_4l" and sel in ["4l", "4m", "2m2e", "4e"]:
            continue
        else:
            f.write("\t&\t&\t" + MC_TEX[suff] + " & " + fmt % np.squeeze(mc[suff][sel])
                    + r" & $\pm$ & " + fmt % np.squeeze(mc_unc[suff][sel]) + " \\\\\n")

    f.write(r"\bottomrule" + "\n")
    f.write(r"\end{tabular}" + "\n")

    f.close()
    print("Wrote table to", fileName)
'''