from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D



##
##  SAMPLE INFO
##

selection   = ["4l", "4m", "2m2e", "4e"]
#selection   = ["mumu", "ee", "4l", "4m", "2m2e", "4e"]
selTeX      = {"mumu":r"\MM", "ee":r"\EE", "4l":r"\fL", "4m":r"\fM", "2m2e":r"\tMtE", "4e":r"\fE", "ll":r"\LL"}
selDef      = {"mumu":"MM", "ee":"EE", "4l":"4L", "4m":"4M", "4e":"4E", "2m2e":"2M2E", "ll":"LL"}
T = np.dtype([(sel, 'f4') for sel in selection])

period      = ["2012", "2016", "2017", "2018"]

# Get list of MC samples
mc_suff = {}
from Cuts2012 import MC_SUFF
mc_suff["2012"] = MC_SUFF
from Cuts2016 import MC_SUFF
mc_suff["2016"] = MC_SUFF
from Cuts2017 import MC_SUFF
mc_suff["2017"] = MC_SUFF
from Cuts2018 import MC_SUFF, MC_TEX
mc_suff["2018"] = MC_SUFF

from Cuts2018 import UNC_DIBOSON, UNC_TTBAR, UNC_TAUTAU, UNC_OTHER

sel_4l = ["4l", "4m", "2m2e", "4e"]
#sel_2l = ["ll", "mumu", "ee"]
sel_2l = []

cutstr = "Pt7"



##
##  DATA
##

data, exp, sig, bg, mc, npt = {}, {}, {}, {}, {}, {}
exp_unc, sig_unc, bg_unc, mc_unc, npt_unc = {}, {}, {}, {}, {}
pur, pur_unc, sf, sf_unc = {}, {}, {}, {}

for year in period:
    infile = "yields" + year + "_" + cutstr + ".npz"
    npzfile = np.load(infile)
    data[year], exp[year], sig[year] = npzfile["data"], npzfile["exp"], npzfile["sig"]
    bg[year], npt[year] = npzfile["bg"], npzfile["npt"]
    exp_unc[year], sig_unc[year] = npzfile["exp_unc"], npzfile["sig_unc"]
    bg_unc[year], npt_unc[year] = npzfile["bg_unc"], npzfile["npt_unc"]

    mc_arr, mc_unc_arr = npzfile["mc"], npzfile["mc_unc"]
    mc_, mc_unc_ = {}, {}
    row = 0
    for suff in mc_suff[year]:
        mc_[suff], mc_unc_[suff] = mc_arr[row], mc_unc_arr[row]
        row += 1
    mc[year], mc_unc[year] = mc_, mc_unc_

    pur_, pur_unc_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
    sf_, sf_unc_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T)
    for sel in selection:
        pur_[sel]       = sig[year][sel] / exp[year][sel] * 100
        pur_unc_[sel]   = pur_[sel] * np.sqrt(sig_unc[year][sel] / sig[year][sel] ** 2 + exp_unc[year][sel] / exp[year][sel] ** 2)
        sf_[sel]        = data[year][sel] / exp[year][sel] * 100
        sf_unc_[sel]    = sf_[sel] * np.sqrt(1 / data[year][sel] + exp_unc[year][sel] / exp[year][sel] ** 2)
    pur[year], pur_unc[year] = pur_, pur_unc_
    sf[year], sf_unc[year] = sf_, sf_unc_




##
##  WRITE TEX FILES
##

# Individual channels

for sel in selection:
    if sel in sel_2l:
        fmt = '{:,.0f}'
    else:
        fmt = '{:,.2f}'

    fileName = "Yield" + selDef[sel] + cutstr + ".tex"
    f = open(fileName, "w")

    f.write(r"\begin{tabular}{lll r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r r@{ $\pm$ }r}" + "\n")
    f.write(r"\toprule" + "\n")

    f.write("\t" + r"& & & \multicolumn{8}{c}{$N_{" + selTeX[sel] + r"}$ (events)} \\" + "\n")
    f.write(r"\cmidrule{4-11}" + "\n")
    f.write("\t" + r"\multicolumn{3}{l}{Type}")
    for year in period:
        f.write(r" & \multicolumn{2}{l}{" + year + "}")
    f.write(r" \\" + "\n")

    f.write(r"\midrule" + "\n")

    f.write("\t" + r"\multicolumn{3}{l}{Observed}")
    for year in period:
        f.write(r" & \multicolumn{2}{l}{" + '{:,.0f}'.format(np.squeeze(data[year][sel])) + "}")
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace\addlinespace" + "\n")

    f.write("\t" + r"\multicolumn{3}{l}{Expected}")
    for year in period:
        f.write(" & " + fmt.format(np.squeeze(exp[year][sel])) + " & "
                + fmt.format(np.squeeze(exp_unc[year][sel])))
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    f.write("\t&\t" + r"\multicolumn{2}{l}{Signal}")
    for year in period:
        f.write(" & " + fmt.format(np.squeeze(sig[year][sel])) + " & "
                + fmt.format(np.squeeze(sig_unc[year][sel])))
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    f.write("\t&\t" + r"\multicolumn{2}{l}{Background}")
    for year in period:
        f.write(" & " + fmt.format(np.squeeze(bg[year][sel])) + " & "
                + fmt.format(np.squeeze(bg_unc[year][sel])))
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    if sel in sel_4l:
        f.write("\t&\t&\t" + "Nonprompt")
        for year in period:
            f.write(" & " + fmt.format(np.squeeze(npt[year][sel])) + " & "
                    + fmt.format(np.squeeze(npt_unc[year][sel])))
        f.write(r" \\" + "\n")

    for suff in MC_SUFF:
        if suff == "tt_2l2nu":
            continue
        elif suff in ["zjets_m-50", "ttbar"] and sel in sel_4l:
            continue
        else:
            f.write("\t&\t&\t" + MC_TEX[suff])
            for year in period:
                if suff in mc_suff[year]:
                    f.write(" & " + fmt.format(np.squeeze(mc[year][suff][sel])) + " & "
                            + fmt.format(np.squeeze(mc_unc[year][suff][sel])))
                else:
                    f.write(r" & \multicolumn{2}{l}{}")
            f.write(r" \\" + "\n")
     
    f.write(r"\midrule" + "\n")

    if sel in sel_2l:
        fmt = '%.4f'
    else:
        fmt = '%.2f'

    f.write("\t" + r"\multicolumn{3}{l}{Signal purity (\%)}")
    for year in period:
        f.write(" & " + fmt % np.squeeze(pur[year][sel]) + " & " + fmt % np.squeeze(pur_unc[year][sel]))
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    f.write("\t" + r"\multicolumn{3}{l}{$N^\text{obs} / N^\text{exp}$ (\%)}")
    for year in period:
        f.write(" & " + fmt % np.squeeze(sf[year][sel]) + " & " + fmt % np.squeeze(sf_unc[year][sel]))
    f.write(r" \\" + "\n")

    f.write(r"\bottomrule" + "\n")
    f.write(r"\end{tabular}" + "\n")

    f.close()
    print("Wrote table to", fileName)




##
##  ADD SAMPLES
##

category = ["Non", "VV", "tt", "tau", "VVV", "H", "ttZ"]
name = {"Non":"Nonprompt", "VV":"Diboson", "VVV":"Triboson", "H":"Higgs", "ttZ":r"$\ttbar\PZ$",
        "tt":r"$\ttbar$", "tau":R"$\ZtoTT$"}

V = np.dtype([(sel, 'f4') for sel in sel_2l + sel_4l])

obs_, exp_ = np.zeros(1, dtype=V), np.zeros(1, dtype=V)
sig_, bg_ = np.zeros(1, dtype=V), np.zeros(1, dtype=V)
exp_unc_, sig_unc_, bg_unc_ = np.zeros(1, dtype=V), np.zeros(1, dtype=V), np.zeros(1, dtype=V)
pur_, pur_unc_ = np.zeros(1, dtype=V), np.zeros(1, dtype=V)
sf_, sf_unc_ = np.zeros(1, dtype=V), np.zeros(1, dtype=V)

cat_arr, cat_unc_arr = np.zeros(len(category), dtype=V), np.zeros(len(category), dtype=V)
cat_, cat_unc_ = {}, {}

row = 0
for cat in category:
    cat_[cat] = cat_arr[row]
    cat_unc_[cat] = cat_unc_arr[row]
    row += 1


for year in period:
    for sel in selection:
        if sel == "4l":
            continue

        obs_[sel] += data[year][sel]
        exp_[sel] += exp[year][sel]
        sig_[sel] += sig[year][sel]
        bg_[sel] += bg[year][sel]

        exp_unc_[sel] += exp_unc[year][sel] ** 2
        sig_unc_[sel] += sig_unc[year][sel] ** 2
        bg_unc_[sel] += bg_unc[year][sel] ** 2

        for suff in mc_suff[year]:
            if suff in ["zz_4l", "ww_2l2nu", "wz_2l2q", "wz_3lnu", "zz_2l2nu", "zz_2l2q"]:
                cat_["VV"][sel] += mc[year][suff][sel]
                cat_unc_["VV"][sel] += mc_unc[year][suff][sel] ** 2 
            elif suff in ["wwz_4l2nu", "wzz_4l2nu", "zzz_4l2nu", "zzg_4l2nu"]:
                cat_["VVV"][sel] += mc[year][suff][sel]
                cat_unc_["VVV"][sel] += mc_unc[year][suff][sel] ** 2
            elif suff in ["ggH_zz_4l", "vbfH_zz_4l"]:
                cat_["H"][sel] += mc[year][suff][sel]
                cat_unc_["H"][sel] += mc_unc[year][suff][sel] ** 2
            elif suff in ["ttbar", "tt_2l2nu"]:
                cat_["tt"][sel] += mc[year][suff][sel]
                cat_unc_["tt"][sel] += mc_unc[year][suff][sel] ** 2
            elif suff == "ttz_2l2nu":
                cat_["ttZ"][sel] += mc[year][suff][sel]
                cat_unc_["ttZ"][sel] += mc_unc[year][suff][sel] ** 2
            elif suff == "zjets_m-50":
                cat_["tau"][sel] += mc[year][suff][sel]
                cat_unc_["tau"][sel] += mc_unc[year][suff][sel] ** 2

    for sel in sel_4l:
        if sel == "4l":
            continue

        cat_["Non"][sel] += npt[year][sel]
        cat_unc_["Non"][sel] += npt_unc[year][sel] ** 2


for sel in sel_4l:
    if sel == "4l":
        continue

    obs_["4l"] += obs_[sel]
    exp_["4l"] += exp_[sel]
    sig_["4l"] += sig_[sel]
    bg_["4l"] += bg_[sel]
    
    exp_unc_["4l"] += exp_unc_[sel]
    sig_unc_["4l"] += sig_unc_[sel]
    bg_unc_["4l"] += bg_unc_[sel]

    for cat in category:
        cat_[cat]["4l"] += cat_[cat][sel]
        cat_unc_[cat]["4l"] += cat_unc_[cat][sel]



for sel in sel_2l:
    if sel == "ll":
        continue

    obs_["ll"] += obs_[sel]
    exp_["ll"] += exp_[sel]
    sig_["ll"] += sig_[sel]
    bg_["ll"] += bg_[sel]
    
    exp_unc_["ll"] += exp_unc_[sel]
    sig_unc_["ll"] += sig_unc_[sel]
    bg_unc_["ll"] += bg_unc_[sel]

    for cat in category:
        cat_[cat]["ll"] += cat_[cat][sel]
        cat_unc_[cat]["ll"] += cat_unc_[cat][sel]


for sel in (sel_2l + sel_4l):
    exp_unc_[sel] = np.sqrt(exp_unc_[sel])
    sig_unc_[sel] = np.sqrt(sig_unc_[sel])
    bg_unc_[sel] = np.sqrt(bg_unc_[sel])

    pur_[sel]       = sig_[sel] / exp_[sel] * 100
    pur_unc_[sel]   = pur_[sel] * np.sqrt(sig_unc_[sel] / sig_[sel] ** 2 + exp_unc_[sel] / exp_[sel] ** 2)
    sf_[sel]        = obs_[sel] / exp_[sel] * 100
    sf_unc_[sel]    = sf_[sel] * np.sqrt(1 / obs_[sel] + exp_unc_[sel] / exp_[sel] ** 2)

    for cat in category:
        cat_unc_[cat][sel] = np.sqrt(cat_unc_[cat][sel])



##
##  WRITE TEX FILES
##

for sel_ in ["ll", "4l"]:
    if sel_ in sel_2l:
        fmt = '{:,.0f}'
        selection = sel_2l
        categories = ["VV", "tt", "tau", "VVV", "H", "ttZ"]
    else:
        fmt = '{:,.1f}'
        selection = sel_4l
        categories = ["Non", "VV", "VVV", "H", "ttZ"]

    fileName = "Yield" + selDef[sel_] + "Sum" + cutstr + ".tex"
    f = open(fileName, "w")

    if sel_ == "ll":
        f.write(r"\begin{tabular}{l l l r@{ $\pm$ }r c r@{ $\pm$ }r r@{ $\pm$ }r}" + "\n")
    else:
        f.write(r"\begin{tabular}{l l l r@{ $\pm$ }l c r@{ $\pm$ }l r@{ $\pm$ }l r@{ $\pm$ }l}" + "\n")
    f.write(r"\toprule" + "\n")

    f.write("\t" + r"\multicolumn{3}{l}{Type}")
    for sel in selection:
        f.write(r" & \multicolumn{2}{l}{$N_{" + selTeX[sel] + r"}$")
        if sel == sel_:
            f.write(" (total)} &")
        else:
            f.write("}")
    f.write(r" \\" + "\n")

    f.write(r"\midrule" + "\n")

    f.write("\t" + r"\multicolumn{3}{l}{Observed}")
    for sel in selection:
        f.write(r" & \multicolumn{2}{l}{" + '{:,.0f}'.format(np.squeeze(obs_[sel])) + "}")
        if sel == sel_:
            f.write(" &")
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace\addlinespace" + "\n")

    f.write("\t" + r"\multicolumn{3}{l}{Expected}")
    for sel in selection:
        f.write(" & " + fmt.format(np.squeeze(exp_[sel])) + " & "
                + fmt.format(np.squeeze(exp_unc_[sel])))
        if sel == sel_:
            f.write(" &")
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    f.write("\t&\t" + r"\multicolumn{2}{l}{Signal}")
    for sel in selection:
        f.write(" & " + fmt.format(np.squeeze(sig_[sel])) + " & "
                + fmt.format(np.squeeze(sig_unc_[sel])))
        if sel == sel_:
            f.write(" &")
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    f.write("\t&\t" + r"\multicolumn{2}{l}{Background}")
    for sel in selection:
        f.write(" & " + fmt.format(np.squeeze(bg_[sel])) + " & "
                + fmt.format(np.squeeze(bg_unc_[sel])))
        if sel == sel_:
            f.write(" &")
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    for cat in categories:
        f.write("\t&\t&\t" + name[cat])
        for sel in selection:
            f.write(" & " + fmt.format(np.squeeze(cat_[cat][sel])) + " & "
                    + fmt.format(np.squeeze(cat_unc_[cat][sel])))
            if sel == sel_:
                f.write(" &")
        f.write(r" \\" + "\n")
     
    f.write(r"\midrule" + "\n")

    if sel in sel_2l:
        fmt = '%.4f'
    else:
        fmt = '%.1f'

    f.write("\t" + r"\multicolumn{3}{l}{Signal purity (\%)}")
    for sel in selection:
        f.write(" & " + fmt % np.squeeze(pur_[sel]) + " & " + fmt % np.squeeze(pur_unc_[sel]))
        if sel == sel_:
            f.write(" &")
    f.write(r" \\" + "\n")

    f.write(r"\addlinespace" + "\n")

    f.write("\t" + r"\multicolumn{3}{l}{$N^\text{obs} / N^\text{exp}$ (\%)}")
    for sel in selection:
        f.write(" & " + fmt % np.squeeze(sf_[sel]) + " & " + fmt % np.squeeze(sf_unc_[sel]))
        if sel == sel_:
            f.write(" &")
    f.write(r" \\" + "\n")

    f.write(r"\bottomrule" + "\n")
    f.write(r"\end{tabular}" + "\n")

    f.close()
    print("Wrote table to", fileName)
