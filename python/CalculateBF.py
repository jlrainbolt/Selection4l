from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D
from secret_number import *



##
##  SAMPLE INFO
##

selection   = [ "mumu", "ee", "4m", "2m2e", "4e"]
selTeX      = { "4m":r"\fM", "4e":r"\fE", "2m2e":r"\tMtE" }
selDef      = { "mumu":"MM",    "ee":"EE",  "4l":"4L",  "4m":"4M",  "4e":"4E",  "2m2e":"2M2E"   }
channel     = { "mumu":3,       "ee":4,     "4m":6,     "2m2e":7,   "2e2m":8,   "4e":9  }
T = np.dtype([(sel, 'f4') for sel in selection])

period      = ["2012", "2016", "2017", "2018"]

sel_4l = ["4m", "2m2e", "4e"]
sel_2l = ["mumu", "ee"]



##
##  DATA
##

data, sig, bg, npt, bg_unc, diff = {}, {}, {}, {}, {}, {}

for year in period:
    infile = "yields" + year + ".npz"
    npzfile = np.load(infile)
    data[year], sig[year], bg[year] = npzfile["data"], npzfile["sig"], npzfile["bg"]
    bg_unc[year], npt[year] = npzfile["bg_unc"], npzfile["npt"]

    diff_ = np.zeros(1, dtype=T)
    for sel in selection:
        diff_[sel] = data[year][sel] - bg[year][sel]
    diff[year] = diff_

n_dy, lumi, xsec, ngen = {}, {}, {}, {}
mu_id_unc, el_id_unc, el_reco_unc, ecal_unc = {}, {}, {}, {}
qcd_unc, pdf_unc, pu_unc = {}, {}, {}

from Cuts2012 import *
n_dy[YEAR_STR], lumi[YEAR_STR] = N_DY, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN 
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup
from Cuts2016 import *
n_dy[YEAR_STR], lumi[YEAR_STR] = N_DY, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN 
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup
from Cuts2017 import *
n_dy[YEAR_STR], lumi[YEAR_STR] = N_DY, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN 
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup
from Cuts2018 import *
n_dy[YEAR_STR], lumi[YEAR_STR] = N_DY, INT_LUMI
xsec[YEAR_STR], ngen[YEAR_STR] = XSEC, NGEN 
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup

print("Total background")
print(bg)
print("")

print("Observed - background")
print(diff)
print("")

print("Expected signal")
print(sig)
print("")
print("")



##
##  ACCEPTANCE AND EFFICIENCY
##

acc, eff, axe = {}, {}, {}

for year in period:
    fid, ps = np.empty(1, dtype=T), np.empty(1, dtype=T)
    acc_, eff_, axe_ = np.empty(1, dtype=T), np.empty(1, dtype=T), np.empty(1, dtype=T)

    # Get number of gen signal events
    inPath, pref = EOS_PATH + "/BLT/" + year + "_new/", "gen"
    fidHistName, psHistName = "FiducialEvents", "PhaseSpaceEvents"

    for suff in ["zz_4l", "zjets_m-50"]:
        inName = pref + "_" + suff + "_0.root"
        inFile = TFile.Open(inPath + inName)
        print("Opened", inPath + inName)
#       fidHist = inFile.Get(fidHistName + "_" + suff)
#       fidHist.SetDirectory(0)
        psHist = inFile.Get(psHistName + "_" + suff)
        psHist.SetDirectory(0)
        inFile.Close()

        if suff == "zjets_m-50":
            for i in range(1, n_dy[year]):
                inName = pref + "_" + suff + "_" + str(i) + ".root"
                inFile = TFile.Open(inPath + inName)
                print("Opened", inPath + inName)
#               fidHist.Add(inFile.Get(fidHistName + "_" + suff))
                psHist.Add(inFile.Get(psHistName + "_" + suff))
                inFile.Close()

        sf = lumi[year] * 1000 * xsec[year][suff] / ngen[year][suff]
#       fidHist.Scale(sf)
        psHist.Scale(sf)

        for sel in selection:
            if (suff == "zz_4l" and sel in sel_4l) or (suff == "zjets_m-50" and sel in sel_2l):
#               fid[sel] = fidHist.GetBinContent(channel[sel])
                ps[sel] = psHist.GetBinContent(channel[sel])
            if sel == "2m2e":
                ps[sel] += psHist.GetBinContent(channel["2e2m"])

    for sel in selection:
#       acc_[sel] = fid[sel] / ps[sel]
#       eff_[sel] = sig[year][sel] / fid[sel]
        axe_[sel] = sig[year][sel] / ps[sel]

    acc[year], eff[year], axe[year] = acc_, eff_, axe_

print("(A * e)")
print(axe)
print("")



##
##  CALCULATE
##

adj = {}

for year in period:
    adj_ = np.zeros(1, dtype=T)
    for sel in selection:
        adj_[sel] = diff[year][sel] / axe[year][sel]
    adj[year] = adj_

bf, bf_stat, bf_syst, pr_stat, pr_syst = {}, {}, {}, {}, {}
for year in period:
    bf_, bf_stat_, bf_syst_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T), np.zeros(1, dtype=T)
    pr_stat_, pr_syst_ = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

    for sel in sel_4l:
        if sel == "4m":
            denom = adj[year]["mumu"]
        elif sel == "2m2e":
            denom = (adj[year]["mumu"] + adj[year]["ee"]) / 2
        elif sel == "4e":
            denom = adj[year]["ee"]

        f = 1000000 * (1 + read_secret_number()) * (1 - F_NR) * BF_LL / denom

        bf_[sel] = f * adj[year][sel]

        pr_stat_[sel] = np.sqrt(1 / diff[year][sel])
        bf_stat_[sel] = bf_[sel] * pr_stat_[sel]

        pr_syst_[sel] = (bg_unc[year][sel] / diff[year][sel]) ** 2
        pr_syst_[sel] += (DELTA_LAMBDA * npt[year][sel] / diff[year][sel]) ** 2

        for src in [mu_id_unc, el_id_unc, el_reco_unc, ecal_unc]:
            pr_syst_[sel] += src[year][sel] ** 2
        for src in [qcd_unc, pdf_unc, pu_unc]:
            pr_syst_[sel] += src[year] ** 2

        pr_syst_[sel] = np.sqrt(pr_syst_[sel])
        bf_syst_[sel] = bf_[sel] * pr_syst_[sel]

    bf[year], bf_stat[year], bf_syst[year] = bf_, bf_stat_, bf_syst_
    pr_stat[year], pr_syst[year] = pr_stat_, pr_syst_

print("Branching fractions")
print(bf)
print("")

print("Statistical uncertainty")
print(bf_stat)
print(pr_stat)
print("")

print("Systematic uncertainty")
print(bf_syst)
print(pr_syst)
print("")
print("")



##
##  SAVE
##

print("")
outfile = "bf_measurements.npz"
np.savez(outfile, bf_2012=bf["2012"], bf_2016=bf["2016"], bf_2017=bf["2017"], bf_2018=bf["2018"],
        bf_stat_2012=bf_stat["2012"], bf_stat_2016=bf_stat["2016"], bf_stat_2017=bf_stat["2017"],
        bf_stat_2018=bf_stat["2018"], bf_syst_2012=bf_syst["2012"], bf_syst_2016=bf_syst["2016"],
        bf_syst_2017=bf_syst["2017"], bf_syst_2018=bf_syst["2018"], diff_2012=diff["2012"],
        diff_2016=diff["2016"], diff_2017=diff["2017"], diff_2018=diff["2018"],
        bg_unc_2012=bg_unc["2012"], bg_unc_2016=bg_unc["2016"], bg_unc_2017=bg_unc["2017"],
        bg_unc_2018=bg_unc["2018"],
        sig_2012=sig["2012"], sig_2016=sig["2016"], sig_2017=sig["2017"], sig_2018=sig["2018"],
        npt_2012=npt["2012"], npt_2016=npt["2016"], npt_2017=npt["2017"], npt_2018=npt["2018"])


print("Wrote arrays to", outfile)
