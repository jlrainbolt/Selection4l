from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D
from secret_number import *

#from Cuts2017 import *
#from Cuts2016 import *
from Cuts2012 import *



##
##  SAMPLE INFO
##

selection   = [ "ll", "mumu", "ee", "4l", "4m", "2m2e", "2e2m", "4e"]
selTeX      = { "mumu":r"\MM",  "ee":r"\EE",
                "4l":r"4\Pell", "4m":r"4\PGm",  "4e":r"4\Pe",   "2m2e":r"2\PGm 2\Pe"
                }
selDef      = { "mumu":"MM",    "ee":"EE",  "4l":"4L",  "4m":"4M",  "4e":"4E",  "2m2e":"2M2E"   }
channel     = { "mumu":3,       "ee":4,     "4m":6,     "2m2e":7,   "2e2m":8,   "4e":9  }
T = np.dtype([(sel, 'f4') for sel in selection])


if (YEAR_STR == "2012"):
    XSEC["zz_4l"] = 0.3305 / 0.31


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
data, data_stat = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    if sel in ["ll", "4l"]:
        continue
    elif sel in ["mumu", "4m", "2m2e"]:
        tree = muFile.Get(sel + "_" + MU_SUFF)
    elif sel in ["ee", "4e", "2e2m"]:
        tree = elFile.Get(sel + "_" + EL_SUFF)

    data[sel] = tree.GetEntries()
    data_stat[sel] = np.sqrt(data[sel])

muFile.Close()
elFile.Close()



##
##  MONTE CARLO
##

mc_arr, mc_stat_arr = np.zeros(N_MC, dtype=T), np.zeros(N_MC, dtype=T)
sig, sig_stat = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
mc, mc_stat = {}, {}
j = 0

# Loop over all samples
for suff in MC_SUFF:
    inName = prefix + "_" + suff + ".root"
    inFile = TFile.Open(inPath + inName)
    print("Opened", inPath + inName)

    # Get histograms
    for sel in selection:
        if sel in ["ll", "4l"]:
            continue
        elif sel in ["mumu", "4m", "2m2e"]:
            lumi = MUON_TRIG_LUMI
        elif sel in ["ee", "4e", "2e2m"]:
            lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

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
            sig_stat[sel] = sf * np.sqrt(tree.GetEntries("!hasTauDecay"))
            cut = "hasTauDecay"
            weight = cut + " * " + weight


        tree.Draw("1>>hist", weight, "goff")

        mc_arr[j][sel] = sf * hist.Integral()
        mc_stat_arr[j][sel] = sf * np.sqrt(tree.GetEntries(cut))

        hist.Delete()

    mc[suff] = mc_arr[j]
    mc_stat[suff] = mc_stat_arr[j]
    j = j + 1
    inFile.Close()



##
##  PHASE SPACE
##
if (YEAR_STR == "2012"):
    XSEC["zz_4l"] = 0.3305 / 0.15

inPath = EOS_PATH + "/BLT/" + YEAR_STR + "/"
prefix, hname = "gen", "PhaseSpaceEvents"

# ZZTo4L file
zzName = prefix + "_zz_4l_" + "0" + ".root"
zzFile = TFile.Open(inPath + zzName)
print("Opened", inPath + zzName)
zzHist = zzFile.Get(hname + "_zz_4l")
zzHist.SetDirectory(0)
zzFile.Close()

# DYJetsToLL file
for i in range(N_DY):
    dyName = prefix + "_zjets_m-50_" + str(i) + ".root"
    dyFile = TFile.Open(inPath + dyName)
    print("Opened", inPath + dyName)

    if i == 0:
        dyHist = dyFile.Get(hname + "_zjets_m-50")
        dyHist.SetDirectory(0)
    else:
        dyHist_ = dyFile.Get(hname + "_zjets_m-50")
        dyHist.Add(dyHist_)
        dyHist_.Delete()

    dyFile.Close()


# Get yields
ps, ps_stat = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    if sel in ["ll", "4l"]:
        continue
    elif sel in ["mumu", "ee"]:
        suff = "zjets_m-50"
        hist = dyHist
    elif sel in ["4m", "4e", "2e2m", "2m2e"]:
        suff = "zz_4l"
        hist = zzHist

    if sel in ["mumu", "4m", "2m2e"]:
        lumi = MUON_TRIG_LUMI
    elif sel in ["ee", "4e", "2e2m"]:
        lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

    sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

    ps[sel] = sf * hist.GetBinContent(channel[sel])
    ps_stat[sel] = sf * hist.GetBinError(channel[sel])

dyHist.Delete()
zzHist.Delete()



##
##  ADD CHANNELS
##

for sample in [data, mc_arr, ps, sig]:
    sample['ll']    = sample['mumu'] + sample['ee']
    sample['2m2e']  = sample['2m2e'] + sample['2e2m']
    sample['4l']    = sample['4m'] + sample['2m2e'] + sample['4e']
    sample['2e2m']  = 0

print("")
print("")
print("Total observed")
print(data)
print("")

# Handle uncertainty
for sample_stat in [data_stat, mc_stat_arr, ps_stat, sig_stat]:
    sample_stat['ll']    = np.sqrt(sample_stat['mumu'] ** 2 + sample_stat['ee'] ** 2)
    sample_stat['2m2e']  = np.sqrt(sample_stat['2m2e'] ** 2 + sample_stat['2e2m'] ** 2)
    sample_stat['4l']    = np.sqrt(sample_stat['4m']**2 + sample_stat['2m2e']**2 + sample_stat['4e']**2)
    sample_stat['2e2m']  = 0



##
##  ADD SAMPLES
##

# Get total expected and background events
exp, exp_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bg, bg_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
diff, diff_stat = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in selection:
    exp[sel]        = np.sum(mc_arr[sel]) + sig[sel] + npt[sel]
    exp_unc[sel]    = np.sqrt(np.sum(mc_stat_arr[sel] ** 2) + sig_stat[sel] ** 2 + npt_unc[sel] ** 2)
    bg[sel]         = np.sum(mc_arr[sel]) + npt[sel]
    bg_unc[sel]     = np.sqrt(np.sum(mc_stat_arr[sel] ** 2) + npt_unc[sel] ** 2)
    diff[sel]       = data[sel] - bg[sel]
    diff_stat[sel]  = np.sqrt(diff[sel])

print("Total background")
print(bg)
print("")

print("Observed - background")
print(diff)
print("")

print("Total expected")
print(exp)
print("")

print("Expected signal")
print(sig)
print("")
print("")


print("Statistical uncertainty in obs - bkg")
print(diff_stat)
print("")


# Calculate acceptance * efficiency
axe, axe_stat = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
for sel in selection:
    if sel == "2e2m":
        continue
    axe[sel] = sig[sel] / ps[sel]
    # FIXME add uncertainty?

print("Acc * eff")
print(axe)
print("")



##
##  CALCULATE
##

ll_factor = (1 - F_NR) * BF_LL * axe['ll'] / (data['ll'] - bg['ll'])
print("Constant for diff dists")
print(ll_factor)
print("")

bf, bf_stat = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
bf_syst, bg_stat = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
prec = np.zeros(1, dtype=T)

for sel in ["4l", "4m", "2m2e", "4e"]:
    factor = 1000000 * (1 + read_secret_number()) * ll_factor / axe[sel]
    bf[sel] = diff[sel] * factor
    bf_stat[sel] = diff_stat[sel] * factor
    bg_stat[sel] = bg_unc[sel] * factor

    for src in [mu_id, el_id, el_reco, mu_pt, el_pt]:
        bf_syst[sel] = bf_syst[sel] + bf[sel] * src[sel]

    prec[sel] = 1 / diff_stat[sel] * 100

print("Branching fractions")
print(bf)
print("")

print("Statistical uncertainty")
print(bf_stat)
print("")
print("")



##
##  WRITE TEX FILES
##

fileName = "BranchingFrac" + YEAR_STR + ".tex"
f = open(fileName, "w")
fmt = '%.3f'

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
            + r" & $\pm$ " + fmt % np.squeeze(bg_stat[sel])
            + r" & " + '%.1f' % np.squeeze(prec[sel]) + r" \\" + "\n")
    if sel == "4l":
        f.write(r"\addlinespace" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)

