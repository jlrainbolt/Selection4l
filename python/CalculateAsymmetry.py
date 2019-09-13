from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TTree, TH1D
from secret_number import *



##
##  SAMPLE INFO
##

selection   = [ "4l", "4m", "2m2e", "4e"]
selTeX      = { "4l":r"4\Pell", "4m":r"4\PGm",  "4e":r"4\Pe",   "2m2e":r"2\PGm 2\Pe"    }
T = np.dtype([(sel, 'f4') for sel in selection])

period      = ["2012", "2016", "2017", "2018"]


# Load header info
mc_suff, lumi, xsec, ngen = {}, {}, {}, {}

from Cuts2012 import *
mc_suff[YEAR_STR], lumi[YEAR_STR] = MC_SUFF, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN
from Cuts2016 import *
mc_suff[YEAR_STR], lumi[YEAR_STR] = MC_SUFF, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN
from Cuts2017 import *
mc_suff[YEAR_STR], lumi[YEAR_STR] = MC_SUFF, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN
from Cuts2018 import *
mc_suff[YEAR_STR], lumi[YEAR_STR] = MC_SUFF, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN

npt, npt_unc = {}, {}
for year in period:
    infile = "yields" + year + ".npz"
    npzfile = np.load(infile)
    npt[year], npt_unc[year] = npzfile["npt"], npzfile["npt_unc"]
print(npt)



##
##  DATA
##

prefix = "boosted"
pos, neg = {}, {}

for year in period:
    inPath = EOS_PATH + "/Boosted/" + year + "_new/"

    # Muon file
    muName = prefix + "_muon_" + year + ".root"
    muFile = TFile.Open(inPath + muName)
    print("Opened", inPath + muName)

    # Electron file
    elName = prefix + "_electron_" + year + ".root"
    elFile = TFile.Open(inPath + elName)
    print("Opened", inPath + elName)

    # Get yields
    pos_, neg_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

    for sel in selection:
        muTree = muFile.Get(sel + "_muon_" + year)
        elTree = elFile.Get(sel + "_electron_" + year)

        hist = TH1D("hist", "", 2, -1, 1)
        muTree.Draw("sin_phi>>hist", "", "goff")
        elTree.Draw("sin_phi>>+hist", "", "goff")

        pos_[sel] = hist.GetBinContent(2)
        neg_[sel] = hist.GetBinContent(1)

        hist.Delete()

    muFile.Close()
    elFile.Close()

    pos[year], neg[year] = pos_, neg_



##
##  MONTE CARLO
##

pos_bkg, neg_bkg, unc = {}, {}, {}

weight = "weight"
for year in period:
    pos_bkg_, neg_bkg_, unc_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T), np.zeros(1, dtype=T)
    inPath = EOS_PATH + "/Boosted/" + year + "_new/"

    for suff in mc_suff[year]:
        if suff in ["zjets_m-50", "ttbar", "tt_2l2nu"]:
            continue

        inName = prefix + "_" + suff + ".root"
        inFile = TFile.Open(inPath + inName)
        print("Opened", inPath + inName)

        sf = lumi[year] * 1000 * xsec[year][suff] / ngen[year][suff]

        if suff == "zz_4l":
            cut = "* (hasTauDecay)"
        else:
            cut = ""

        for sel in selection:
            tree = inFile.Get(sel + "_" + suff)
            hist = TH1D("hist", "", 2, -1, 1)

            tree.Draw("sin_phi>>hist", weight + cut, "goff")
            pos_bkg_[sel] += sf * hist.GetBinContent(2)
            neg_bkg_[sel] += sf * hist.GetBinContent(1)

            tree.Draw("0>>hist", weight + " * " + weight + cut, "goff")
            unc_[sel] += sf * hist.Integral()

            hist.Delete()

        inFile.Close()

    pos_bkg[year], neg_bkg[year], unc[year] = pos_bkg_, neg_bkg_, unc_



##
##  SAME SIGN
##

prefix = "boosted_bkg"

for year in period:
    inPath = EOS_PATH + "/Boosted/" + year + "_new/"

    # Muon file
    muName = prefix + "_muon_" + year + ".root"
    muFile = TFile.Open(inPath + muName)
    print("Opened", inPath + muName)

    # Electron file
    elName = prefix + "_electron_" + year + ".root"
    elFile = TFile.Open(inPath + elName)
    print("Opened", inPath + elName)

    # Get yields
    for sel in selection:
        muTree = muFile.Get(sel + "_muon_" + year)
        elTree = elFile.Get(sel + "_electron_" + year)

        hist = TH1D("hist", "", 2, -1, 1)
        muTree.Draw("sin_phi>>hist", "", "goff")
        elTree.Draw("sin_phi>>+hist", "", "goff")

        sf = npt[year][sel] / hist.Integral()

        pos_bkg[year][sel] += sf * hist.GetBinContent(2)
        neg_bkg[year][sel] += sf * hist.GetBinContent(1)

        unc[year][sel] += npt_unc[year][sel] ** 2
        unc[year][sel] += (DELTA_LAMBDA * npt[year][sel]) ** 2

        hist.Delete()

    muFile.Close()
    elFile.Close()



##
##  CALCULATE
##

assy, syst_unc, stat_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T), np.zeros(1, dtype=T)
pos_tot, neg_tot = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    for year in period:
        stat_unc[sel] += pos[year][sel] + neg[year][sel]

        pos[year][sel] -= pos_bkg[year][sel]
        neg[year][sel] -= neg_bkg[year][sel]
        pos_tot[sel] += pos[year][sel]
        neg_tot[sel] += neg[year][sel]

        syst_unc[sel] += unc[year][sel]
        unc[year][sel] = np.sqrt(unc[year][sel])

    total = pos_tot[sel] + neg_tot[sel]
    assy[sel] = (pos_tot[sel] - neg_tot[sel]) / total

    stat_unc[sel] = 1 / np.sqrt(total)
    syst_unc[sel] = np.sqrt(1 / (total - unc[year][sel]) - 1 / total)
#   syst_unc[sel] = max(np.abs(1 / np.sqrt(total) - 1 / np.sqrt(total + syst_unc[sel])),
#                           np.abs(1 / np.sqrt(total) - 1 / np.sqrt(total - syst_unc[sel])))

    stat_unc[sel] = 1 / np.sqrt(total)
    syst_unc[sel] = np.sqrt(1 / (total - unc[year][sel]) - 1 / total)

# FIXME
for year in period:
    pos[year]["4l"] = pos[year]["4m"] + pos[year]["2m2e"] + pos[year]["4e"]
    neg[year]["4l"] = neg[year]["4m"] + neg[year]["2m2e"] + neg[year]["4e"]

pos_tot["4l"] = pos_tot["4m"] + pos_tot["2m2e"] + pos_tot["4e"]
neg_tot["4l"] = neg_tot["4m"] + neg_tot["2m2e"] + neg_tot["4e"]
total = pos_tot["4l"] + neg_tot["4l"]
assy["4l"] = (pos_tot["4l"] - neg_tot["4l"]) / total

pos["Total"], neg["Total"] = pos_tot, neg_tot


print(assy)
print(stat_unc)
print(syst_unc)



##
##  WRITE TEX FILES
##

period      = ["2012", "2016", "2017", "2018", "Total"]

fileName = "Asymmetry.tex"
f = open(fileName, "w")
fmt = '%.2f'

f.write(r"\begin{tabular}{l c r r c r r c r r c r r c c r r c r@{\ $\pm$ }r@{\ $\pm$ }r}" + "\n")
f.write(r"\toprule" + "\n")
f.write("\t")

for year in period:
    if year == "Total":
        f.write("&")
    f.write(r"&& \multicolumn{2}{c}{" + year + "} ")

f.write(r"&& \multicolumn{3}{c}{$A_{\sin\phi}$ (\%)} \\" + "\n")
f.write(r"\rulepad \cline{3-4} \cline{6-7} \cline{9-10} \cline{12-13} \cline{16-17}"
        + r" \cline{19-21} \rulepad" + "\n")

f.write("Channel ")
for year in period:
    if year == "Total":
        f.write("&")
    f.write(r"&& \multicolumn{1}{l}{$N_{+}$} & \multicolumn{1}{l}{$N_{-}$} ")
f.write(r" && \multicolumn{3}{r}{\stat \ \syst} \\" + "\n")

f.write(r"\midrule" + "\n")
for sel in selection:
    f.write("$" + selTeX[sel] + "$ ")
    for year in period:
        if year == "Total":
            f.write("&")
        f.write("&& " + fmt % np.squeeze(pos[year][sel]) + " & " + fmt % np.squeeze(neg[year][sel])
                + " ")
    f.write("&& $" + fmt % np.squeeze(100 * assy[sel]) + "$ & " + fmt % np.squeeze(100 * stat_unc[sel])
            + " & " + fmt % np.squeeze(100 * syst_unc[sel]) + r" \\" + "\n")
    if sel == "4l":
        f.write(r"\addlinespace" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
