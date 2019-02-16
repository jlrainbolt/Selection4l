from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D
from secret_number import *

from Cuts2017 import *
#from Cuts2016 import *



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
    elif sel in ["mumu", "4m", "2m2e"]:
        tree = muFile.Get(sel + "_" + MU_SUFF)
    elif sel in ["ee", "4e", "2e2m"]:
        tree = elFile.Get(sel + "_" + EL_SUFF)

    hist = TH1D("hist", "", 2, -1, 1)
    tree.Draw("sin_phi>>hist", "", "goff")

    pos[sel] = hist.GetBinContent(2)
    neg[sel] = hist.GetBinContent(1)

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


'''
##
##  WRITE TEX FILES
##

fileName = "BranchingFrac" + YEAR_STR + ".tex"
f = open(fileName, "w")
fmt = '%.2f'

f.write(r"\begin{tabular}{l p{.5in} p{.75in} p{.75in}>{\leavevmode\color{purple}}p{.75in} c}"
        + "\n")
f.write(r"\toprule" + "\n")
f.write("\t" + r"& \multicolumn{4}{c}{$\BF$ ($\times 10^{-6}$)}"
        + r" & \multirow{2}{1in}[-1ex]{Statistical Precision (\%)} \\" + "\n")
f.write(r"\cmidrule{2-5}" + "\n")
f.write(r"Channel & & (stat.) & (syst.) & (MC stat.) \\" + "\n")
f.write(r"\midrule" + "\n")

for sel in ["4l", "4m", "2m2e", "4e"]:
    f.write("$" + selTeX[sel] + r"$ & " + fmt % np.squeeze(bf[sel])
            + r" & $\pm$ " + fmt % np.squeeze(bf_stat[sel])
            + r" & $\pm$ " + fmt % np.squeeze(bf_syst[sel])
            + r" & $\pm$ " + fmt % np.squeeze(mc_stat[sel])
            + r" & " + '%.1f' % np.squeeze(prec[sel]) + r" \\" + "\n")
    if sel == "4l":
        f.write(r"\addlinespace" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
'''
