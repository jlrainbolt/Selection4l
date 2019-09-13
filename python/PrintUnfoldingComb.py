from __future__ import print_function
from __future__ import division

import numpy as np



##
##  DATA
##

infile = "unfoldingcomb.npz"
npzfile = np.load(infile)

cond_num, n_iter, bins = npzfile["cond_num"], npzfile["n_iter"], npzfile['bins']
chisq_smr, chisq_unf = npzfile["chisq_smr"], npzfile["chisq_unf"]
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

f.write(r"\begin{tabular}{l r l r r r}" + "\n")
f.write(r"\toprule" + "\n")

f.write("\t\t" + r"$x$ & \multicolumn{1}{l}{$\nu$} & $\cond(R)$ & "
        + r"\multicolumn{1}{l}{$n$} & \multicolumn{1}{l}{$\chi^2_\smr$} & "
        + r"\multicolumn{1}{l}{$\chi^2_\unf$} \\" + "\n")

f.write(r"\midrule" + "\n")

for h in range(len(names)):
#   if names[h] in ["sin_phi", "sin_phi_10"]:
#       continue

    f.write("\t" + nameTeX[names[h]] + " & {:,.0f}".format(np.squeeze(bins[h])) 
            + " & {:,.2f}".format(np.squeeze(cond_num[h]))
            + " & {:,.0f}".format(np.squeeze(n_iter[h]))
            + " & {:,.2f}".format(np.squeeze(chisq_smr[h])))
    if chisq_unf[h] < 50:
        f.write(" & {:,.2f}".format(np.squeeze(chisq_unf[h])))
    else:
        f.write(" & \multicolumn{1}{l}{$> 50$}")

    f.write(r" \\" + "\n")

f.write(r"\bottomrule" + "\n")
f.write(r"\end{tabular}" + "\n")

f.close()
print("Wrote table to", fileName)
