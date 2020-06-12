from __future__ import print_function
from __future__ import division

import numpy as np
from numpy.linalg import multi_dot, pinv

from ROOT import TFile, TH1, TH2, TH2D

from PlotUtils import *

from Cuts2018 import *


np.set_printoptions(precision=6, suppress=True, linewidth=200)

##
##  SAMPLE INFO
##

period  = ["2018", "2017", "2016", "2012"]
hnames  = ["b_z1m", "b_z2m", "b_l1p", "b_ttm", "cos_theta_z1", "cos_theta_z2",
        "angle_z1leps", "angle_z2leps", "angle_z1l2_z2", "sin_phi"]
uncertainty = ["MuonID", "ElecID", "ElecReco"]
delta = {"MuonID":0.01213, "ElecID":0.00931, "ElecReco":0.00791, "Others":0.01397}
total_syst = 0.02217
# FIXME load from npz

H = len(hnames)
T = np.dtype([(sel, object) for sel in ["4l"]])
V = np.dtype([(unc, object) for unc in uncertainty])

sel = "4l"
suff = "zz_4l"




# Get lists of lumi, xsec, ngen for signal
lumi, xsec, ngen, mc_suff = {}, {}, {}, {}
from Cuts2012 import INT_LUMI, XSEC, NGEN, MC_SUFF_4L
lumi["2012"], xsec["2012"], ngen["2012"], mc_suff["2012"] = INT_LUMI, XSEC, NGEN, MC_SUFF_4L
from Cuts2016 import INT_LUMI, XSEC, NGEN, MC_SUFF_4L
lumi["2016"], xsec["2016"], ngen["2016"], mc_suff["2016"] = INT_LUMI, XSEC, NGEN, MC_SUFF_4L
from Cuts2017 import INT_LUMI, XSEC, NGEN, MC_SUFF_4L
lumi["2017"], xsec["2017"], ngen["2017"], mc_suff["2017"] = INT_LUMI, XSEC, NGEN, MC_SUFF_4L
from Cuts2018 import INT_LUMI, XSEC, NGEN, MC_SUFF_4L
lumi["2018"], xsec["2018"], ngen["2018"], mc_suff["2018"] = INT_LUMI, XSEC, NGEN, MC_SUFF_4L



##
##  INPUT FILES
##

inPath = ""

# Migration and MC matrices
nom = np.empty(H, dtype=object)
up, dn = np.empty(H, dtype=V), np.empty(H, dtype=V)

for year in period:
    inName = inPath + "ddr_" + year + "_" + suff + ".root"
    inFile = TFile(inName, "READ")
    print("Opened", inName)
    h = 0

    sf = lumi[year] * 1000 * xsec[year][suff] / ngen[year][suff]

    for hname in hnames:
        nom_ = inFile.Get("weight/" + sel + "/" + hname + "_" + suff)
        nom_.SetDirectory(0)
        nom_.Scale(sf)

        if year == "2018":
            nom[h] = nom_
        else:
            nom[h].Add(nom_)

        for unc in uncertainty:
            up_ = inFile.Get("wt" + unc + "Up/" + sel + "/" + hname + "_" + suff)
            up_.SetDirectory(0)
            up_.Scale(sf)

            dn_ = inFile.Get("wt" + unc + "Dn/" + sel + "/" + hname + "_" + suff)
            dn_.SetDirectory(0)
            dn_.Scale(sf)

            if year == "2018":
                up[h][unc] = up_
                dn[h][unc] = dn_
            else:
                up[h][unc].Add(up_)
                dn[h][unc].Add(dn_)

        h = h + 1
    h = 0
    inFile.Close()

print("Got nominal and up/down histograms", "\n")


# Unfolding covariance matrices and scaling

infile = "unfoldingcomb.npz"
npzfile = np.load(infile)

cov_unf = []
for hname in hnames:
    cov_unf.append(npzfile["cov_mat_" + hname])

infile = "scaling.npz"
npzfile = np.load(infile)

scl_vec = []
for hname in hnames:
    scl_vec.append(npzfile["scl_vec_" + hname])





####
####
####    LOOP OVER DISTS
####
####


# Arrays
cov, cov_tot, scl_mat = [], [], []

# Histograms
hcov_new, hcov_unf, hcov_tot = np.empty(H, dtype=T), np.empty(H, dtype=T), np.empty(H, dtype=T)


for h in range(H):
    print(hnames[h])


    ##
    ##  GET BIN CONTENT
    ##
 
    nbins = nom[h].GetNbinsX() + 2
 
    f = np.zeros(nbins)
    for i in range(nbins):
        f[i] = nom[h].GetBinContent(i)

    # Slicing
    if hnames[h] == "b_z1m":
        s = slice(2, -1)
    elif hnames[h] == "b_z2m":
        s = slice(1, None)
    elif hnames[h] == "b_l1p":
        s = slice(0, -1)
    elif hnames[h] == "angle_z1leps":
        s = slice(2, -1)
    elif hnames[h] in ["angle_z1l2_z2", "cos_theta_z1", "cos_theta_z2", "angle_z2leps", "sin_phi", "sin_phi_10"]:
        s = slice(1, -1)
    else:
        s = slice(0, None)

    # Slice nominal histogram to visible bins
    f = f[s]
    nvbins = len(f)

    # Get systematics from "other" sources (just fraction of bin count)
    cov.append(delta["Others"] ** 2 * np.tensordot(f, f.T, axes=0))


    # Difference of upper and lower uncertainty variations
    for unc in uncertainty:
        delta_f = np.zeros(nbins)

        for i in range(nbins):
            delta_f[i] = (up[h][unc].GetBinContent(i) - dn[h][unc].GetBinContent(i)) / 2

        # Slice to visible bins
        delta_f = delta_f[s]
     
        cov_h = np.tensordot(delta_f, delta_f.T, axes=0)
     
        # Find total uncertainty such that change in chi-squared = 1
        chi_sq_1 = multi_dot([np.zeros(nvbins), pinv(cov_h), np.zeros(nvbins).T])
        chi_sq_unc = multi_dot([f, pinv(cov_h), f.T])
     
        sigma_sq = (chi_sq_unc - chi_sq_1) * delta[unc] ** 2
        cov_h = sigma_sq * cov_h

        cov[-1] = cov[-1] + cov_h



    # Scaling
    scl_vec[h] = scl_vec[h][s]
    scl_mat.append(np.tensordot(scl_vec[h], scl_vec[h].T, axes=0))

    cov[h] = np.multiply(scl_mat[h], cov[h])
    cov_unf[h] = np.multiply(scl_mat[h], cov_unf[h])
    cov_tot.append(cov[h] + cov_unf[h])



    ##
    ##  DRAW RESULTS
    ##

    o = s.start
    bmin = o - 0.5
    bmax = bmin + nvbins

    # New covariance matrix
    hcov_new[h][sel] = TH2D(hnames[h] + "_cov_syst", "", nvbins, bmin, bmax, nvbins, bmin, bmax)
    hcov_new[h][sel].SetXTitle("i (bin)");
    hcov_new[h][sel].SetYTitle("j (bin)");
    for i in range(nvbins):
        for j in range(nvbins):
            hcov_new[h][sel].SetBinContent(i + 1, j + 1, cov[-1][i][j])


    # Unfolding covariance matrix
    hcov_unf[h][sel] = TH2D(hnames[h] + "_cov_unf", "", nvbins, bmin, bmax, nvbins, bmin, bmax)
    hcov_unf[h][sel].SetXTitle("i (bin)");
    hcov_unf[h][sel].SetYTitle("j (bin)");
    for i in range(nvbins):
        for j in range(nvbins):
            hcov_unf[h][sel].SetBinContent(i + 1, j + 1, cov_unf[h][i][j])


    # Total covariance matrix
    hcov_tot[h][sel] = TH2D(hnames[h] + "_cov_tot", "", nvbins, bmin, bmax, nvbins, bmin, bmax)
    hcov_tot[h][sel].SetXTitle("i (bin)");
    hcov_tot[h][sel].SetYTitle("j (bin)");
    for i in range(nvbins):
        for j in range(nvbins):
            hcov_tot[h][sel].SetBinContent(i + 1, j + 1, cov_tot[-1][i][j])



##
##  OUTPUT FILE
##

# Create file
outName = "ddr_uncertainty.root"
outFile = TFile(outName, "RECREATE")


for h in range(H):
    outFile.mkdir(hnames[h])
    outFile.cd(hnames[h])

    # Write histograms
    hcov_new[h]['4l'].Write()
    hcov_unf[h]['4l'].Write()
    hcov_tot[h]['4l'].Write()


outFile.Close()
print("Wrote output to", outName)



##
##  PRINT TEXT
##

filePref = "total_cov"

for h in range(H):
    fileName = filePref + "_" + hnames[h] + ".tex"
    np.savetxt(fileName, cov_tot[h], fmt='% .2e', comments='', header=r'\begin{verbatim}',
            footer=r'\end{verbatim}')

print("Printed arrays to", filePref + "_*.tex")



##
##  SAVE MATRICES
##

names = np.array(hnames, dtype=object)   # strings of (ugly) names

print("")
outfile = "diffuncertainties.npz"
np.savez(outfile, names=names,
        cov_tot_b_z1m=cov_tot[0], cov_tot_b_z2m=cov_tot[1], cov_tot_b_l1p=cov_tot[2],
        cov_tot_b_ttm=cov_tot[3], cov_tot_cos_theta_z1=cov_tot[4],
        cov_tot_cos_theta_z2=cov_tot[5], cov_tot_angle_z1leps=cov_tot[6],
        cov_tot_angle_z2leps=cov_tot[7], cov_tot_angle_z1l2_z2=cov_tot[8],
        cov_tot_sin_phi=cov_tot[9])


print("Wrote arrays to", outfile)

