from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D



##
##  SAMPLE INFO
##

period = ["2012", "2016", "2017", "2018"]



##
##  DATA
##

cond_num, n_iter, chisq_smr, chisq_unf, bins = {}, {}, {}, {}, {}

for year in period:
    infile = "unfolding" + year + ".npz"
    npzfile = np.load(infile)

    cond_num[year], n_iter[year], bins[year] = npzfile["cond_num"], npzfile["n_iter"], npzfile['bins']
    chisq_smr[year], chisq_unf[year] = npzfile["chisq_smr"], npzfile["chisq_unf"]

names = npzfile['names']

nameTeX = {"b_z1m":r"$m_{\mathrm{Z}_{1}}$", "b_z2m":r"$m_{\mathrm{Z}_{2}}$",
            "b_l1p":r"$p_{\ell_{1}}$", "b_ttm":r"$m_{\ell_{2,3,4}}$",
            "angle_z1l2_z2":r"$\beta$", "angle_z1leps":r"$\alpha_{\mathrm{Z}_{1}}$",
            "angle_z2leps":r"$\alpha_{\mathrm{Z}_{2}}$",
            "cos_theta_z1":r"$\cos\theta_{\mathrm{Z}_{1}}$",
            "cos_theta_z2":r"$\cos\theta_{\mathrm{Z}_{2}}$",
            "sin_phi":r"$\sin\phi$", "sin_phi_10":r"$\sin\phi$"}




##
##  WRITE TEX FILES
##

fileName = "Unfolding.tex"
f = open(fileName, "w")

f.write(r"\begin{tabular}{l l r l r r r}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t\t" + r"$x$ & Year & \multicolumn{1}{l}{$\nu$} & $\cond(R)$ & "
        + r"\multicolumn{1}{l}{$n$} & \multicolumn{1}{l}{$\chi^2_\smr$} & "
        + r"\multicolumn{1}{l}{$\chi^2_\unf$} \\" + "\n")

f.write(r"\midrule" + "\n")

for h in range(len(names)):
    if names[h] == "sin_phi":
        continue

    for year in period:
        if year == "2012":
            f.write("\t" + nameTeX[names[h]] + " ")
        else:
            f.write("\t\t")

        f.write("& " + year + " & {:,.0f}".format(np.squeeze(bins[year][h])) 
                + " & {:,.2f}".format(np.squeeze(cond_num[year][h]))
                + " & {:,.0f}".format(np.squeeze(n_iter[year][h]))
                + " & {:,.2f}".format(np.squeeze(chisq_smr[year][h])))
        if chisq_unf[year][h] < 50:
            f.write(" & {:,.2f}".format(np.squeeze(chisq_unf[year][h])))
        else:
            f.write(" & \multicolumn{1}{l}{$> 50$}")

        f.write(r" \\" + "\n")
        if year == "2018":
            f.write(r"\addlinespace" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
