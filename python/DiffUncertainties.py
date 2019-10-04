from __future__ import print_function
from __future__ import division

import numpy as np
from numpy.linalg import multi_dot, pinv

from ROOT import TFile, TH1

from PlotUtils import *

from Cuts2018 import *


np.set_printoptions(precision=5, suppress=True, linewidth=110)

##
##  SAMPLE INFO
##

period  = ["2018", "2017", "2016", "2012"]
hnames  = ["b_z1m", "b_z2m", "b_l1p", "b_ttm", "cos_theta_z1", "cos_theta_z2",
        "angle_z1leps", "angle_z2leps", "angle_z1l2_z2", "sin_phi"]
uncertainty = ["MuonID", "ElecID", "ElecReco"]
delta = {"MuonID":0.01213, "ElecID":0.00931, "ElecReco":0.00791}  # FIXME load from npz

H = len(hnames)
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




####
####
####    LOOP OVER DISTS
####
####

for unc in uncertainty:
    for h in range(H):

        ##
        ##  GET BIN CONTENT
        ##
     
        # Diff
        nbins = nom[h].GetNbinsX()
     
        f, delta_f = np.zeros(nbins), np.zeros(nbins)
        for i in range(nbins):
            f[i] = nom[h].GetBinContent(i+1)
            delta_f[i] = (up[h][unc].GetBinContent(i+1) - dn[h][unc].GetBinContent(i+1)) / 2
     
        # Slicing
        if hnames[h] == "b_z1m":
            s = slice(2, -1)
        elif hnames[h] == "b_z2m":
            s = slice(1, None)
        elif hnames[h] == "b_l1p":
            s = slice(0, -1)
        elif hnames[h] == "angle_z1leps":
            s = slice(2, -1)
        elif hnames[h] == "angle_z1l2_z2":
            s = slice(1, -1)
        elif hnames[h] in ["cos_theta_z1", "cos_theta_z2", "angle_z2leps", "sin_phi", "sin_phi_10"]:
            s = slice(1, -1)
        else:
            s = slice(0, None)
     
        f = f[s]
        delta_f = delta_f[s]
        nbins = len(f)
     
        cov = np.tensordot(delta_f, delta_f.T, axes=0)
     
        print(hnames[h])
        print(unc)
#       print("f:", f, sep="\n", end="\n\n")
#       print("delta f:", delta_f, sep="\n", end="\n\n")
#       print("w:", cov, sep="\n", end="\n\n")
     
        chi_sq_1 = multi_dot([np.zeros(nbins), pinv(cov), np.zeros(nbins).T])
        chi_sq_unc = multi_dot([f, pinv(cov), f.T])
     
        sigma_sq = (chi_sq_unc - chi_sq_1) * delta[unc] ** 2
#       print("chi^2(1):", chi_sq_1, sep="\n", end="\n\n")
#       print("chi^2(1+unc):", chi_sq_unc, sep="\n", end="\n\n")
#       print("sigma^2:", sigma_sq)
        print("sigma:", np.sqrt(sigma_sq))
        print("")
    print("\n\n\n")
