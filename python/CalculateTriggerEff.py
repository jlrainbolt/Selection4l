from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D

from Cuts2018 import *
#from Cuts2017 import *
#from Cuts2016 import *
#from Cuts2012 import *



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

varexp_pair = [ "diTrig12", "diTrig13", "diTrig14", "diTrig23", "diTrig24", "diTrig34"  ]
varexp_2m2e = [ "doubleMuTrig", "doubleElTrig"  ]

s = {"both":1, "single":1, "double":0, "neither":0}
d = {"both":1, "single":0, "double":1, "neither":0}

T = np.dtype([(sel, 'f4') for sel in selection])
V = np.dtype([(var, 'f4') for var in varexp_pair])
W = np.dtype([(var, 'f4') for var in varexp_2m2e])



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

mc, mc_var = np.zeros((2,2), dtype=T), np.zeros((2,2), dtype=T)

mc_pair_eff = {}

for sel in selection:
    if sel in ["mumu", "ee"]:
        suff = "zjets_m-50"
        tree = dyFile.Get(sel + "_" + suff)
        weight = "weight/trigWeight"
    elif sel in ["4m", "2m2e", "4e"]:
        suff = "zz_4l"
        tree = zzFile.Get(sel + "_" + suff)
        weight = "weight"
    elif sel == "2e2m":
        tree = zzFile.Get("2m2e_zz_4l")
        weight = "weight"


    if sel in ["mumu", "4m", "2m2e"]:
        varexp = varexp_muon
    elif sel in ["ee", "4e", "2e2m"]:
        varexp = varexp_elec


    sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]

    for cont in contingency:
        cut = "(!hasTauDecay && " + varexp[cont] + ")"
        hist = TH1D("hist", "", 1, 0, 2)

        tree.Draw("1>>hist", cut + " * " + weight, "goff")
        mc[s[cont], d[cont]][sel] = sf * hist.Integral()

        hist = TH1D("hist", "", 1, 0, 2)
        tree.Draw("1>>hist", cut + " * " + weight + " * " + weight, "goff")
        mc_var[s[cont], d[cont]][sel] = sf * hist.Integral()


    if sel in ["4m", "2m2e", "4e"]:
        if sel == "2m2e":
            varexp_list = varexp_2m2e
            X = W
        else:
            varexp_list = varexp_pair
            X = V

        mc_pair_eff_ = np.zeros(1, dtype=X)
        for varexp in varexp_list:
            hist = TH1D("hist", "", 2, -0.5, 1.5)

            cut = "(!hasTauDecay)"
            tree.Draw(varexp + ">>hist", cut + " * " + weight, "goff")
            mc_pair_eff_[varexp] = hist.GetBinContent(2) / hist.Integral()

        mc_pair_eff[sel] = mc_pair_eff_


    mc_var[0,0] = -1

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
        total += data[s[cont], d[cont]][sel]

    for cont in contingency:
        data_frac[s[cont], d[cont]][sel] = data[s[cont], d[cont]][sel] / total

muFile.Close()
elFile.Close()



##
##  CALCULATE
##

# Efficiency for single-lepton triggered events to also fire dilepton trigger
mc_ll_eff, data_ll_eff, ll_sf = {}, {}, {}
mc_ll_var, data_ll_var, ll_sf_unc = {}, {}, {}
for sel in ["mumu", "ee"]:
    mc_ll_eff[sel] = mc[1,1][sel] / (mc[1,1][sel] + mc[1,0][sel])
    mc_ll_var[sel] = mc_var[1,1][sel] * mc[1,0][sel] ** 2 + mc_var[1,0][sel] * mc[1,1][sel] ** 2
    mc_ll_var[sel] /= (mc[1,1][sel] + mc[1,0][sel]) ** 4

    data_ll_eff[sel] = data[1,1][sel] / (data[1,1][sel] + data[1,0][sel])
    data_ll_var[sel] = (data[1,1][sel] * data[1,0][sel]) / (data[1,1][sel] + data[1,0][sel]) ** 3

    ll_sf[sel] = data_ll_eff[sel] / mc_ll_eff[sel]
    ll_sf_unc[sel] = data_ll_var[sel] / data_ll_eff[sel] ** 2 + mc_ll_var[sel] / mc_ll_eff[sel] ** 2
    ll_sf_unc[sel] = np.sqrt(ll_sf_unc[sel]) * ll_sf[sel]


# Total 4l efficiency
mc_4l_eff = {}
for sel in ["4m", "2m2e", "4e"]:
    mc_4l_eff[sel] = mc[1,1][sel] + mc[0,1][sel] 
    if sel == "2m2e":
        mc_4l_eff[sel] += mc[1,1]["2e2m"] + mc[0,1]["2e2m"]
    mc_4l_eff[sel] /= mc[1,1][sel] + mc[1,0][sel] + mc[0,1][sel] + mc[0,0][sel]

# 4l efficiency from pairs
mc_4l_check = {}
for sel in ["4m", "2m2e", "4e"]:
    miss = 1
    if sel == "2m2e":
        varexp_list = varexp_2m2e
    else:
        varexp_list = varexp_pair

    for var in varexp_list:
        miss *= (1 - mc_pair_eff[sel][var])
    mc_4l_check[sel] = 1 - miss

# 4l efficiency in data
data_4l_check, data_4l_up, data_4l_down = {}, {}, {}

for data_4l in [data_4l_check, data_4l_up, data_4l_down]:
    for sel in ["4m", "2m2e", "4e"]:
        miss = 1

        if sel == "2m2e":
            for ll_sel in ["mumu", "ee"]:
                if ll_sel == "mumu":
                    var = "doubleMuTrig"
                elif ll_sel == "ee":
                    var = "doubleElTrig"

                sf = ll_sf[ll_sel]
                sf_unc = ll_sf_unc[ll_sel]

                if data_4l == data_4l_up:
                    sf += sf_unc
                elif data_4l == data_4l_down:
                    sf -= sf_unc

                miss *= (1 - sf * mc_pair_eff[sel][var])
            

        else:
            if sel == "4m":
                ll_sel = "mumu"
            elif sel == "4e":
                ll_sel = "ee"

            sf = ll_sf[ll_sel]
            sf_unc = ll_sf_unc[ll_sel]

            if data_4l == data_4l_up:
                sf += sf_unc
            elif data_4l == data_4l_down:
                sf -= sf_unc

            for var in varexp_pair:
                miss *= (1 - sf * mc_pair_eff[sel][var])

        data_4l[sel] = 1 - miss




##
##  PRINT
##

'''
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
    print(fmt % np.squeeze(mc[1,1][sel]), fmt % np.squeeze(mc[0,1][sel]), sep=sep, end=sep)
    print(fmt % np.squeeze(data[1,1][sel]), fmt % np.squeeze(data[0,1][sel]), sep=sep)
    print(fmt % np.squeeze(mc[1,0][sel]), fmt % np.squeeze(mc[0,0][sel]), sep=sep, end=sep)
    print(fmt % np.squeeze(data[1,0][sel]), fmt % np.squeeze(data[0,0][sel]), sep=sep)
    print("\n\n")

#print(mc_pair_eff)
'''


print("")
print("MC single-to-double efficiency")
for sel in ["mumu", "ee"]:
    print(sel, ": ", np.squeeze(mc_ll_eff[sel]), " +- ", np.squeeze(np.sqrt(mc_ll_var[sel])), sep='')

print("")
print("Data single-to-double efficiency")
for sel in ["mumu", "ee"]:
    print(sel, ": ", np.squeeze(data_ll_eff[sel]), " +- ", np.squeeze(np.sqrt(data_ll_var[sel])), sep='')

print("")
print("Scale factor")
for sel in ["mumu", "ee"]:
    print(sel, ": ", np.squeeze(ll_sf[sel]), " +- " , np.squeeze(ll_sf_unc[sel]), sep='')

print("")
print("Actual MC efficiency")
for sel in ["4m", "2m2e", "4e"]:
    print(sel, ": ", np.squeeze(mc_4l_eff[sel]), sep='')

print("")
print("Calculated MC efficiency")
for sel in ["4m", "2m2e", "4e"]:
    print(sel, ": ", np.squeeze(mc_4l_check[sel]), sep='')

print("")
print("Calculated data efficiency")
for sel in ["4m", "2m2e", "4e"]:
    print(sel, ": ", np.squeeze(data_4l_check[sel]), sep='')

print("")
print("Calculated data efficiency (up)")
for sel in ["4m", "2m2e", "4e"]:
    print(sel, ": ", np.squeeze(data_4l_up[sel]), sep='')

print("")
print("Calculated data efficiency (down)")
for sel in ["4m", "2m2e", "4e"]:
    print(sel, ": ", np.squeeze(data_4l_down[sel]), sep='')


'''
##
##  SAVE
##

print("\n")
outfile = "trigger_eff_" + YEAR_STR + ".py"
np.savez(outfile, mc=mc, data=data)
'''
