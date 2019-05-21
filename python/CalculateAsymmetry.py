from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TTree, TH1D
from secret_number import *

#from Cuts2018 import *
#from Cuts2017 import *
#from Cuts2016 import *
from Cuts2012 import *



##
##  SAMPLE INFO
##

selection   = [ "4l", "4m", "2m2e", "2e2m", "4e"]
selTeX      = { "4l":r"4\Pell", "4m":r"4\PGm",  "4e":r"4\Pe",   "2m2e":r"2\PGm 2\Pe"    }
T = np.dtype([(sel, 'f4') for sel in selection])

year = sys.argv[1]
if year != YEAR_STR:
    print("Wrong year in header file")



##
##  DATA
##

inPath = HOME_PATH + "/Boosted/" + YEAR_STR + "/"
prefix = "boosted"

# Muon file
muName = prefix + "_" + MU_SUFF + ".root"
muFile = TFile.Open(inPath + muName)
print("Opened", inPath + muName)

# Electron file
elName = prefix + "_" + EL_SUFF + ".root"
elFile = TFile.Open(inPath + elName)
print("Opened", inPath + elName)

# Get yields
pos, neg = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    if sel in ["4l"]:
        continue
    elif sel in ["4m", "2m2e"]:
        tree = muFile.Get(sel + "_" + MU_SUFF)
    elif sel in ["4e", "2e2m"]:
        tree = elFile.Get(sel + "_" + EL_SUFF)

    hist = TH1D("hist", "", 2, -1, 1)
    tree.Draw("sin_phi>>hist", "", "goff")

    pos[sel] = hist.GetBinContent(2)
    neg[sel] = hist.GetBinContent(1)

    hist.Delete()

muFile.Close()
elFile.Close()



##
##  MONTE CARLO
##

# Loop over all samples
for suff in MC_SUFF:
    if suff in ["zjets_m-50", "ttbar", "tt_2l2nu"]:
        continue

    inName = prefix + "_" + suff + ".root"
    inFile = TFile.Open(inPath + inName)
    print("Opened", inPath + inName)

    # Get histograms
    for sel in selection:
        if sel in [ "4l"]:
            continue
        elif sel in ["4m", "2m2e"]:
            lumi = MUON_TRIG_LUMI
        elif sel in ["ee", "4e", "2e2m"]:
            lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

        if suff in ["zz_4l"]:
            cut = "* (hasTauDecay)"
        else:
            cut = ""

        sf = lumi * 1000 * XSEC[suff] / NGEN[suff]
        weight = "weight/trigWeight/qtWeight"

        tree = inFile.Get(sel + "_" + suff)

        hist = TH1D("hist", "", 2, -1, 1)
        tree.Draw("sin_phi>>hist", weight + cut, "goff")
        
        pos[sel] = pos[sel] - sf * hist.GetBinContent(2)
        neg[sel] = neg[sel] - sf * hist.GetBinContent(1)

        hist.Delete()

    inFile.Close()



##
##  BACKGROUND
##

prefix = "boosted_bkg"

# Muon file
muName = prefix + "_" + MU_SUFF + ".root"
muFile = TFile.Open(inPath + muName)
print("Opened", inPath + muName)

# Electron file
elName = prefix + "_" + EL_SUFF + ".root"
elFile = TFile.Open(inPath + elName)
print("Opened", inPath + elName)

for sel in selection:
    if sel in ["4l"]:
        continue
    elif sel in ["4m", "2m2e"]:
        tree = muFile.Get(sel + "_" + MU_SUFF)
    elif sel in ["4e", "2e2m"]:
        tree = elFile.Get(sel + "_" + EL_SUFF)

    hist = TH1D("hist", "", 2, -1, 1)
    tree.Draw("sin_phi>>hist", "", "goff")

    sf = npt[sel] / hist.Integral()

    pos[sel] = pos[sel] - sf * hist.GetBinContent(2)
    neg[sel] = neg[sel] - sf * hist.GetBinContent(1)

    hist.Delete()

muFile.Close()
elFile.Close()




##
##  ADD CHANNELS
##

for sample in [pos, neg]:
    sample['2m2e']  = sample['2m2e'] + sample['2e2m']
    sample['4l']    = sample['4m'] + sample['2m2e'] + sample['4e']



##
##  CALCULATE
##

assy, unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
sig = np.zeros(1, dtype=T)

for sel in ["4l", "4m", "2m2e", "4e"]:
    assy[sel] = (pos[sel] - neg[sel]) / (pos[sel] + neg[sel])
    unc[sel] = 1 / np.sqrt(pos[sel] + neg[sel])
    sig[sel] = abs(assy[sel] / unc[sel])

print(assy)
print(unc)
print(sig)



##
##  WRITE TEX FILES
##

fileName = "Asymmetry" + YEAR_STR + ".tex"
f = open(fileName, "w")
fmt = '%.2f'

f.write(r"\begin{tabular}{l r r r c r c}" + "\n")
f.write(r"\toprule" + "\n")
f.write(r"Channel & \multicolumn{1}{l}{$N_{+}$} & \multicolumn{1}{l}{$N_{-}$}"
        + r"& \multicolumn{3}{l}{$A$ (\%)} & Significance ($\sigma$) \\" + "\n")
f.write(r"\midrule" + "\n")

for sel in ["4l", "4m", "2m2e", "4e"]:
    f.write("$" + selTeX[sel] + r"$ & " + fmt % np.squeeze(pos[sel])
            + r" & " + fmt % np.squeeze(neg[sel])
            + r" & $" + fmt % np.squeeze(100 * assy[sel])
            + r"$ & $\pm$ & $" + fmt % np.squeeze(100 * unc[sel])
            + r"$ & " + '%.2f' % np.squeeze(sig[sel]) + r" \\" + "\n")
    if sel == "4l":
        f.write(r"\addlinespace" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
