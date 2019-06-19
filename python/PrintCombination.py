from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D
from secret_number import *



##
##  SAMPLE INFO
##

selTeX      = { "4m":r"\fM", "4e":r"\fE", "2m2e":r"\tMtE", "4l":r"\fL" }
selection   = [ "4m",   "2m2e", "4e"    ]
I           = { "4m":0, "2m2e":1, "4e":2, "2012":0, "2016":1, "2017":2, "2018":3}
period      = [ "2012", "2016", "2017", "2018"  ]

P, M = len(period), len(selection)
n_params = [1, M, P, P*M]
name = ["SM", "BSM", "Yearly", "Indiv."]
name = dict(zip(n_params, name))



##
##  DATA
##

bf, bf_stat, bf_syst = {}, {}, {}

for year in period:
    infile = "bf_measurements.npz"
    npzfile = np.load(infile)

    bf["2012"], bf_stat["2012"] = npzfile["bf_2012"], npzfile["bf_stat_2012"]
    bf_syst["2012"] = npzfile["bf_syst_2012"]

    bf["2016"], bf_stat["2016"] = npzfile["bf_2016"], npzfile["bf_stat_2016"]
    bf_syst["2016"] = npzfile["bf_syst_2016"]

    bf["2017"], bf_stat["2017"] = npzfile["bf_2017"], npzfile["bf_stat_2017"]
    bf_syst["2017"] = npzfile["bf_syst_2017"]

    bf["2018"], bf_stat["2018"] = npzfile["bf_2018"], npzfile["bf_stat_2018"]
    bf_syst["2018"] = npzfile["bf_syst_2018"]


alpha_stat, alpha_total = {}, {}
delta_stat, delta_syst, delta_total = {}, {}, {}
chi_sq_stat, chi_sq_total = {}, {}
bf_pred = {}

for N in n_params:
    infile = "combination_" + str(N) + ".npz"
    npzfile = np.load(infile)

    alpha_stat[N], alpha_total[N] = npzfile["alpha_stat"], npzfile["alpha_total"]
    delta_stat[N], delta_total[N] = npzfile["delta_stat"], npzfile["delta_total"]
    delta_syst[N] = npzfile["delta_syst"]
    chi_sq_stat[N], chi_sq_total[N] = npzfile["chi_sq_stat"], npzfile["chi_sq_total"]

    if N == 1:
        bf_pred = npzfile["bf_pred"]

bf_comb = bf_pred * alpha_total[1]
bf_comb_stat = bf_comb * delta_stat[1]
bf_comb_syst = bf_comb * delta_syst[1]


##
##  WRITE TEX FILES
##

fileName = "BFMeasurements.tex"
f = open(fileName, "w")
fmt, fmt2 = '%.2f', '%.1f'

f.write(r"\begin{tabular}{l l@{$\quad\quad$} r @{$\quad\pm$ }l @{$\quad\pm$ }l l r r}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t" + r"& & \multicolumn{3}{c}{$\BF$ ($\ten{$-6$}$)}"
        + r" && \multicolumn{2}{l}{Precision (\%)} \\" + "\n")
f.write(r"\rulepad \cline{3-5} \cline{7-8} \rulepad" + "\n")
f.write(r"Year & Channel & \multicolumn{3}{c}{\quad \quad \stat \quad \syst}"
        + r" && $\stat$ & $\syst$ \\" + "\n")
f.write(r"\midrule" + "\n")

for year in period:
    for sel in selection:
        if sel == "4m":
            f.write(year + " ")
        else:
            f.write("\t")
        f.write(" & $" + selTeX[sel] + "$ & " + fmt % np.squeeze(bf[year][sel]) + " & "
                + fmt % np.squeeze(bf_stat[year][sel]) + " & "
                + fmt % np.squeeze(bf_syst[year][sel]) + " && "
                + fmt2 % np.squeeze(bf_stat[year][sel] / bf[year][sel] * 100) + " & "
                + fmt2 % np.squeeze(bf_syst[year][sel] / bf[year][sel] * 100) + r" \\" + "\n")
    if year == "2018":
        f.write(r"\midrule" + "\n")
    else:
        f.write(r"\addlinespace" + "\n")

for sel in selection:
    i = I[sel]
    if sel == "4m":
        f.write("Combined ")
    else:
        f.write("\t")
    f.write("& $" + selTeX[sel] + "$ & "
            + fmt % np.squeeze(bf_comb[i]) + " & " + fmt % np.squeeze(bf_comb_stat[i]) + " & "
            + fmt % np.squeeze(bf_comb_syst[i]) + " && "
            + fmt2 % np.squeeze(delta_stat[1] * 100) + " & "
            + fmt2 % np.squeeze(delta_syst[1]* 100) + r" \\" + "\n")

f.write(r"\addlinespace" + "\n")
f.write("\t" + "& $" + selTeX["4l"] + "$ & "
            + fmt % np.sum(bf_comb) + " & " + fmt % np.squeeze(np.sum(bf_comb_stat)) + " & "
            + fmt % np.sum(bf_comb_syst) + " && " + fmt2 % np.squeeze(delta_stat[1] * 100) + " & "
            + fmt2 % np.squeeze(delta_syst[1]* 100) + r" \\" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)



fileName = "BFCombination.tex"
f = open(fileName, "w")
fmt = '%.2f'


f.write(r"\begin{tabular}{l ll c rr c rrr}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t" + r"& & && \multicolumn{2}{c}{Signal strength} && \multicolumn{3}{c}{Precision (\%)} \\"
        + "\n")
f.write(r"\rulepad \cline{5-6} \cline{8-10} \rulepad" + "\n")
f.write("\t" + r"Case & \multicolumn{2}{l}{Category} && \stat & \total && \stat & \syst & \total \\"
        + "\n")
f.write(r"\midrule" + "\n")

for N in n_params:
    if N == 1:
        f.write("\t" + name[N] + " & & && " + fmt % np.squeeze(alpha_stat[N]) + " & "
                + fmt % np.squeeze(alpha_total[N]) + " && " + fmt % np.squeeze(100 * delta_stat[N])
                + " & " + fmt % np.squeeze(100 * delta_syst[N]) + " & "
                + fmt % np.squeeze(100 * delta_total[N]) + r" \\" + "\n")

    elif N == M:
        for sel in selection:
            i = I[sel]
            if sel == "4m":
                f.write("\t" + name[N])
            else:
                f.write("\t\t")
            f.write(r" & \multicolumn{2}{l}{$" + selTeX[sel] + "$} && "
                    + fmt % np.squeeze(alpha_stat[N][i]) + " & "
                    + fmt % np.squeeze(alpha_total[N][i]) + " && "
                    + fmt % np.squeeze(100 * delta_stat[N][i]) + " & "
                    + fmt % np.squeeze(100 * delta_syst[N][i]) + " & "
                    + fmt % np.squeeze(100 * delta_total[N][i]) + r" \\" + "\n")

    elif N == P:
        for year in period:
            i = I[year]
            if year == "2012":
                f.write("\t" + name[N])
            else:
                f.write("\t\t")
            f.write(r" & \multicolumn{2}{l}{" + year + "} && "
                    + fmt % np.squeeze(alpha_stat[N][i]) + " & "
                    + fmt % np.squeeze(alpha_total[N][i]) + " && "
                    + fmt % np.squeeze(100 * delta_stat[N][i]) + " & "
                    + fmt % np.squeeze(100 * delta_syst[N][i]) + " & "
                    + fmt % np.squeeze(100 * delta_total[N][i]) + r" \\" + "\n")

    elif N == P*M:
        for year in period:
            if year == "2012":
                f.write("\t" + name[N])
            else:
                f.write("\t\t")
            for sel in selection:
                i = I[sel] + M * I[year]
                if sel == "4m":
                    f.write(" & " + year)
                else:
                    f.write("\t\t\t" + "&")
                f.write(" & $" + selTeX[sel] + "$ && "
                    + fmt % np.squeeze(alpha_stat[N][i]) + " & "
                    + fmt % np.squeeze(alpha_total[N][i]) + " && "
                    + fmt % np.squeeze(100 * delta_stat[N][i]) + " & "
                    + fmt % np.squeeze(100 * delta_syst[N][i]) + " & "
                    + fmt % np.squeeze(100 * delta_total[N][i]) + r" \\" + "\n")

            if year != "2018":
                f.write(r"\addlinespace" + "\n")

    if N == P*M:
        f.write(r"\bottomrule" + "\n")
    else:
        f.write(r"\addlinespace\addlinespace" + "\n")

f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
