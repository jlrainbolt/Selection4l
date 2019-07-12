from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

#from Cuts2018 import *
#from Cuts2017 import *
#from Cuts2016 import *
from Cuts2012 import *



##
##  SAMPLE INFO
##

selection   = ["mumu", "ee", "4m", "2m2e", "2e2m", "4e"]
contingency = ["both", "single", "double", "neither"]

varexp_muon = { "both":"singleMuTrig && doubleMuTrig",
                "single":"singleMuTrig && !doubleMuTrig",
                "double":"!singleMuTrig && doubleMuTrig",
                "neither":"!singleMuTrig && !doubleMuTrig"}

varexp_elec = { "both":"singleElTrig && doubleElTrig",
                "single":"singleElTrig && !doubleElTrig",
                "double":"!singleElTrig && doubleElTrig",
                "neither":"!singleElTrig && !doubleElTrig"}

s = {"both":1, "single":1, "double":0, "neither":0}
d = {"both":1, "single":0, "double":1, "neither":0}

T = np.dtype([(sel, 'f4') for sel in selection])



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


# Get yields and efficiencies

mc_unw, mc_unw_frac = np.zeros((2,2), dtype=T), np.zeros((2,2), dtype=T)
mc_w, mc_w_frac = np.zeros((2,2), dtype=T), np.zeros((2,2), dtype=T)

weight, unweight = "weight", "weight/trigWeight"

for sel in selection:
    if sel in ["mumu", "ee"]:
        suff = "zjets_m-50"
        tree = dyFile.Get(sel + "_" + suff)
    elif sel in ["4m", "2m2e", "4e"]:
        suff = "zz_4l"
        tree = zzFile.Get(sel + "_" + suff)
    elif sel == "2e2m":
        tree = zzFile.Get("2m2e_zz_4l")


    if sel in ["mumu", "4m", "2m2e"]:
        varexp = varexp_muon
    elif sel in ["ee", "4e", "2e2m"]:
        varexp = varexp_elec


    sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]

    tot_w, tot_unw = 0, 0
    for cont in contingency:
        cut = "(!hasTauDecay && " + varexp[cont] + ")"
        hist = TH1D("hist", "", 1, 0, 2)

        tree.Draw("1>>hist", cut + " * " + unweight, "goff")
        mc_unw[s[cont], d[cont]][sel] = sf * hist.Integral()

        tree.Draw("1>>hist", cut + " * " + weight, "goff")
        mc_w[s[cont], d[cont]][sel] = sf * hist.Integral()

        if cont != "neither":
            tot_unw += mc_unw[s[cont], d[cont]][sel]
            tot_w += mc_w[s[cont], d[cont]][sel]

    for cont in contingency:
        mc_unw_frac[s[cont], d[cont]][sel] = mc_unw[s[cont], d[cont]][sel] / tot_unw
        mc_w_frac[s[cont], d[cont]][sel] = mc_w[s[cont], d[cont]][sel] / tot_w

    mc_unw_frac[0,0] = -1
    mc_w_frac[0,0] = -1

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
data, data_frac = np.zeros((2,2), dtype=T), np.zeros((2,2), dtype=T)

for sel in selection:
    if sel == "2e2m":
        muTree = muFile.Get("2m2e_" + MU_SUFF)
        elTree = elFile.Get("2m2e_" + EL_SUFF)
    else:
        muTree = muFile.Get(sel + "_" + MU_SUFF)
        elTree = elFile.Get(sel + "_" + EL_SUFF)


    if sel in ["mumu", "4m", "2m2e"]:
        varexp = varexp_muon
    elif sel in ["ee", "4e", "2e2m"]:
        varexp = varexp_elec


    total = 0
    for cont in contingency:
        data[s[cont], d[cont]][sel]  = muTree.GetEntries(varexp[cont]) + elTree.GetEntries(varexp[cont])
#       data_unc_[cont] = np.sqrt(data_[cont])

        total += data[s[cont], d[cont]][sel]


    for cont in contingency:
        data_frac[s[cont], d[cont]][sel] = data[s[cont], d[cont]][sel] / total
#       data_unc_[cont] = data_unc_[cont] / total

muFile.Close()
elFile.Close()



##
##  PRINT
##

print("\n")
print(YEAR_STR, end='\n\n')
print("Passed both", "Double only", sep='\t')
print("Single only", "Failed both", sep='\t')
print("\n\n")

for sel in ["mumu", "4m", "2m2e", "ee", "2e2m", "4e"]:
    if sel == "mumu":
        print("MUON TRIGGERS", end='\n\n')
    elif sel == "ee":
        print("ELECTRON TRIGGERS", end='\n\n')

    if sel in ["mumu", "ee"]:
        fmt = '%8.0f'
        sep='\t'
    else:
        fmt = '%6.2f'
        sep='\t\t'

    print('\t\t\t', sel, sep='', end='\n\n')
    print("\tMC", "Data", sep='\t\t\t\t')
    print(fmt % np.squeeze(mc_unw[1,1][sel]), fmt % np.squeeze(mc_unw[0,1][sel]), sep=sep, end=sep)
    print(fmt % np.squeeze(data[1,1][sel]), fmt % np.squeeze(data[0,1][sel]), sep=sep)
    print(fmt % np.squeeze(mc_unw[1,0][sel]), fmt % np.squeeze(mc_unw[0,0][sel]), sep=sep, end=sep)
    print(fmt % np.squeeze(data[1,0][sel]), fmt % np.squeeze(data[0,0][sel]), sep=sep)
    print("\n\n")



##
##  SAVE
##

print("\n")
outfile = "trigger_eff_" + YEAR_STR + ".py"
#np.savez(outfile, mc_noTW=mc_unw, mc_noTW_frac=mc_unw_frac, mc_tw=mc_w, mc_tw_frac=mc_w_frac,
#        data=data, data_frac=data_frac)
np.savez(outfile, mc=mc_unw, data=data)
