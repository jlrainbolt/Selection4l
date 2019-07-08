from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

#from Cuts2017 import *
from Cuts2016 import *



##
##  SAMPLE INFO
##

selection   = ["mumu", "4m"]
contingency = ["both", "single", "double", "neither"]
varexp      = { "both":"singleMuTrig && doubleMuTrig",
                "single":"singleMuTrig && !doubleMuTrig",
                "double":"!singleMuTrig && doubleMuTrig",
                "neither":"!singleMuTrig && !doubleMuTrig"}

T = np.dtype([(cont, 'f4') for cont in contingency])



##
##  MONTE CARLO
##

inPath = EOS_PATH + "/Systematics/" + YEAR_STR + "_new/"
prefix = "triggerEff"

# Drell-Yan file
dyName = prefix + "_zjets_m-50.root"
dyFile = TFile.Open(inPath + dyName)
print("Opened", inPath + dyName)

# ZZTo4L file
zzName = prefix + "_zz_4l.root"
zzFile = TFile.Open(inPath + zzName)
print("Opened", inPath + zzName)


# Get yields (without trigger scale factors)
mc_unw, mc_unw_frac = {}, {}

weight = "weight"

for sel in selection:
    mc_unw_, mc_unw_frac_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

    if sel in ["mumu", "ee"]:
        suff = "zjets_m-50"
        tree = dyFile.Get(sel + "_" + suff)
    elif sel in ["4m", "2e2m", "4e"]:
        suff = "zz_4l"
        tree = zzFile.Get(sel + "_" + suff)

    sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]

    total = 0
    for cont in contingency:
        cut = "(!hasTauDecay && " + varexp[cont] + ")"
        hist = TH1D("hist", "", 1, 0, 2)

        tree.Draw("1>>hist", cut + " * " + weight, "goff")
        mc_unw_[cont] = sf * hist.Integral()

#       tree.Draw("1>>hist", cut + " * " + weight + " * " + weight, "goff")
#       mc_unw_unc_[cont] = sf * np.sqrt(hist.Integral())

        if cont != "neither":
            total += mc_unw_[cont]

    for cont in contingency:
        mc_unw_frac_[cont] = mc_unw_[cont] / total
    mc_unw_frac_["neither"] = -1

    mc_unw[sel], mc_unw_frac[sel] = mc_unw_, mc_unw_frac_

dyFile.Close()
zzFile.Close()



##
##  Data
##

inPath = EOS_PATH + "/Selected/" + YEAR_STR + "_new/"
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
data, data_frac = {}, {}

for sel in selection:
    data_, data_frac_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

    muTree = muFile.Get(sel + "_" + MU_SUFF)
    elTree = elFile.Get(sel + "_" + EL_SUFF)

    total = 0
    for cont in contingency:
        data_[cont] = muTree.GetEntries(varexp[cont]) + elTree.GetEntries(varexp[cont])
#       data_unc_[cont] = np.sqrt(data_[cont])

        total += data_[cont]


    for cont in contingency:
        data_frac_[cont] = data_[cont] / total
#       data_unc_[cont] = data_unc_[cont] / total

    data[sel], data_frac[sel] = data_, data_frac_

muFile.Close()
elFile.Close()



##
##  PRINT
##

print("\n")
print("MC (no trigger SFs)")
print(mc_unw)
print("\nFraction")
print(mc_unw_frac)

print("\n")
print("Data")
print(data)
print("\nFraction")
print(data_frac)
