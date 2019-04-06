from __future__ import print_function
from __future__ import division
import sys

import numpy as np
from pyunfold import iterative_unfold
from pyunfold.callbacks import Logger
from scipy.stats import chi2

from ROOT import TFile, TH1, TH2, TH2D, TCanvas, TLegend

from PlotUtils import *
from Cuts2017 import *
#from Cuts2016 import *


np.set_printoptions(precision=1, suppress=True)

##
##  SAMPLE INFO
##

selection = ["4l", "4m", "2m2e", "2e2m", "4e"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f8'), ("y", 'f8'), ("ex", 'f8'), ("ey", 'f8'), ("b", 'f8')])

hnames = ["b_l1p", "b_ttm", "cos_theta_z1", "cos_theta_z2", "angle_z1leps", "angle_z2leps",
            "angle_z1l2_z2"]
H = len(hnames)

year = sys.argv[1]
if year != YEAR_STR:
    print("Wrong year in header file")



##
##  INPUT FILES
##

# Migration and MC matrices
inName = "migration_" + YEAR_STR + "_zz_4l.root"
inFile = TFile(inName, "READ")
print("Opened", inName)

gen, reco = np.empty(H, dtype=T), np.empty(H, dtype=T)
mig = np.empty(H, dtype=T)
h = 0

for sel in selection:
    if sel == "4l":
        continue
    elif sel in ["4m", "2m2e"]:
        lumi = MUON_TRIG_LUMI
    elif sel in ["4e", "2e2m"]:
        Lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

    sf = lumi * 1000 * XSEC['zz_4l'] / NGEN['zz_4l']

    for hname in hnames:
        gen[h][sel] = inFile.Get(sel + "/" + hname + "_gen")
        gen[h][sel].SetDirectory(0)
        gen[h][sel].Scale(sf)
    
        reco[h][sel] = inFile.Get(sel + "/" + hname + "_reco")
        reco[h][sel].SetDirectory(0)
        reco[h][sel].Scale(sf)
    
        mig[h][sel] = inFile.Get(sel + "/" + hname + "_2d")
        mig[h][sel].SetDirectory(0)
        mig[h][sel].Scale(sf)
    
        h = h + 1
    h = 0

inFile.Close()
print("Got migration and MC histograms", "\n")



##
##  DATA
##

prefix = "4l"

# Muon file
muName = prefix + "_" + YEAR_STR + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + YEAR_STR + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)

# Get histograms
data = np.empty(H, dtype=T)
h = 0

for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        if sel in ["4m", "2m2e"]:
            data[h][sel] = muFile.Get(sel + "/" + hname + "_" + MU_SUFF)
        elif sel in ["4e", "2e2m"]:
            data[h][sel] = elFile.Get(sel + "/" + hname + "_" + EL_SUFF)

        data[h][sel].SetName(hnames[h] + "_data")
        data[h][sel].SetDirectory(0)
        h = h + 1
    h = 0

muFile.Close()
elFile.Close()
print("Got data histograms")
print("")



##
##  BACKGROUND
##

prefix = "bkg_all"

# Muon file
muName = prefix + "_" + YEAR_STR + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + YEAR_STR + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)

# Get histograms
bkg = np.empty(H, dtype=T)
h = 0

for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        if sel in ["4m", "2m2e"]:
            bkg[h][sel] = muFile.Get(sel + "/" + hname + "_" + MU_SUFF)
        elif sel in ["4e", "2e2m"]:
            bkg[h][sel] = elFile.Get(sel + "/" + hname + "_" + EL_SUFF)

        bkg[h][sel].SetName(hnames[h] + "_bkg")
        bkg[h][sel].SetDirectory(0)
        h = h + 1
    h = 0

muFile.Close()
elFile.Close()
print("Got nonprompt background histograms")
print("")


# Monte Carlo
prefix = "4l"
h = 0

# Loop over all samples
for suff in MC_SUFF_4L:
    if suff in ["zz_4l", "zjets_m-50", "tt_2l2nu", "ttbar"]:
        continue

    inName = prefix + "_" +year + "_" + suff + ".root"
    inFile = TFile.Open(inName)
    print("Opened", inName)

    # Get histograms
    for sel in selection:
        if sel == "4l":
            continue
        elif sel in ["4m", "2m2e"]:
            lumi = MUON_TRIG_LUMI
        elif sel in ["4e", "2e2m"]:
            lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

        sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

        for hname in hnames:
            hist = inFile.Get(sel + "/" + hname + "_" + suff)
            hist.SetDirectory(0)
            hist.Scale(sf)
            bkg[h][sel].Add(hist)

            h = h + 1
        h = 0

    inFile.Close()

print("Got background MC histograms")
print("")



##
##  ADD CHANNELS
##

# Get 4l and 2m2e
for h in range(H):
    for sample in [data, gen, reco, mig, bkg]:
        sample[h]['2m2e'].Add(sample[h]['2e2m'])
        sample[h]['4l'] = sample[h]['2m2e'].Clone()
        sample[h]['4l'].Add(sample[h]['4m'])
        sample[h]['4l'].Add(sample[h]['4e'])

# Background subtraction
for sel in ["4l"]:
    for h in range(H):
        sf = npt[sel] / bkg[h][sel].Integral()
        bkg[h][sel].Scale(sf)
        data[h][sel].Add(bkg[h][sel], -1)

        err_sf = npt_unc[sel] / np.sqrt(bkg[h][sel].GetSumw2().GetSum())
        for i in range(1, bkg[h][sel].GetNbinsX()):
            bkg[h][sel].SetBinError(i, err_sf * bkg[h][sel].GetBinError(i))






####
####
####    LOOP OVER DISTS
####
####


# Store results in histograms (for now)
result = np.empty(H, dtype=T)
resp, unf, cov = np.empty(H, dtype=T), np.empty(H, dtype=T), np.empty(H, dtype=T)
var = np.empty(H, dtype=T)

for sel in ["4l"]:
    for h in range(H):

        ##
        ##  GET BIN CONTENT
        ##

        # Data
        v_data = np.zeros(data[h][sel].GetNbinsX() + 2, dtype=V)
        for i in range(len(v_data)):
            v_data[i]['x']  = data[h][sel].GetBinCenter(i)
            v_data[i]['y']  = np.maximum(data[h][sel].GetBinContent(i), 0)
            v_data[i]['ey'] = data[h][sel].GetBinError(i)

        # Gen MC
        v_gen = np.zeros(gen[h][sel].GetNbinsX() + 2, dtype=V)
        for i in range(len(v_gen)):
            v_gen[i]['x']   = gen[h][sel].GetBinLowEdge(i)
            v_gen[i]['y']   = gen[h][sel].GetBinContent(i)
            v_gen[i]['ey']  = gen[h][sel].GetBinError(i)

        # Reco MC (only used for tests)
        v_reco = np.zeros(reco[h][sel].GetNbinsX() + 2, dtype=V)
        for i in range(len(v_gen)):
            v_reco[i]['x']  = reco[h][sel].GetBinLowEdge(i)
            v_reco[i]['y']  = reco[h][sel].GetBinContent(i)
            v_reco[i]['ey'] = reco[h][sel].GetBinError(i)

        # Migrations
        v_mig = np.zeros([mig[h][sel].GetNbinsX() + 2, mig[h][sel].GetNbinsY() + 2], dtype=V)
        for i in range(np.size(v_mig, 0)):
            for j in range(np.size(v_mig, 1)):
                v_mig[i][j]['y']    = mig[h][sel].GetBinContent(i, j)
                v_mig[i][j]['ey']   = mig[h][sel].GetBinError(i, j)



        ##
        ##  CREATE MATRICES
        ##

        # Slicing
        if hnames[h] == "b_l1p":
            s = slice(0, -1)
        elif hnames[h] == "b_ttm":
            s = slice(1, None)
        elif hnames[h] in ["cos_theta_z1", "cos_theta_z2", "angle_z2leps"]:
            s = slice(1, -1)
        elif hnames[h] == "angle_z1leps":
            s = slice(3, -1)
        elif hnames[h] == "angle_z1l2_z2":
            s = slice(1, -2)

        v_data = v_data[s]
        v_gen = v_gen[s]
        v_reco = v_reco[s]
        v_mig = v_mig[s, s]


        # Efficiencies (assume 1)
        v_eff       = np.zeros_like(v_gen, dtype=V)
        v_eff['y']  = np.ones_like(v_gen['y'])

        # Response matrix
        col_sums    = v_mig['y'].sum(axis=0)
        norm_factor = v_eff['y'] / col_sums

        v_resp = np.zeros_like(v_mig, dtype=V)
        v_resp['y']     = v_mig['y'] * norm_factor
        v_resp['ey']    = v_mig['ey'] * norm_factor

        # Test response matrix
#       print("Reco:", v_reco['y'])
#       print("Response * gen:", np.dot(v_resp['y'], v_gen['y'].T))

#       v_eff['y'] = np.dot(v_reco['y'] / v_resp['y'], v_gen['y'].T)
#       v_resp['y'] = v_resp['y'] * v_eff['y']
#       print("Efficiencies:", v_eff['y'])




        ##
        ##  CONDITION NUMBER
        ##

        sigma = np.linalg.svd(v_resp['y'], compute_uv = False)
        cond_num = np.max(sigma) / np.min(sigma)
        print("Condition number for", hnames[h], "is", cond_num)



        ##
        ##  PERFORM UNFOLDING
        ##

        chisq_nu = chi2.isf(0.997, len(v_resp['y'])) / len(v_resp['y'])
        results = iterative_unfold( data = v_data['y'],         data_err = v_data['ey'],
                                    response = v_resp['y'],     response_err = v_resp['ey'],
                                    efficiencies = v_eff['y'],  efficiencies_err = v_eff['ey'],
                                    ts = 'chi2',                ts_stopping = chisq_nu,
                                    callbacks = [Logger()]
                                    )



        ##
        ##  GET RESULTS
        ##

        # Unfolded result
        v_result        = np.zeros_like(v_gen, dtype=V)
        v_result['y']   = results['unfolded']
        v_result['ey']  = np.sqrt(results['stat_err']**2, results['sys_err']**2)

        o = s.start

        result[h][sel] = data[h][sel].Clone(hnames[h] + "_result");
        result[h][sel].SetYTitle("");
        for i in range(len(v_result)):
            result[h][sel].SetBinContent(i + o, v_result[i]['y'])
            result[h][sel].SetBinError(i + o, v_result[i]['ey'])


        # Bin histograms in terms of indices
        xbins, ybins = len(v_result), len(v_result)
        xmin, ymin = o - 0.5, o - 0.5
        xmax, ymax = xmin + xbins, ymin + ybins


        # Response matrix
        resp[h][sel] = TH2D(hnames[h] + "_response", "", xbins, xmin, xmax, ybins, ymin, ymax);
        resp[h][sel].SetXTitle("i (bin)");
        resp[h][sel].SetYTitle("j (bin)");
        for i in range(xbins):
            for j in range(ybins):
                resp[h][sel].SetBinContent(i + 1, j + 1, v_resp[i][j]['y'])
                resp[h][sel].SetBinError(i + 1, j + 1, v_resp[i][j]['ey'])

        # Unfolding matrix
        v_unf   = results['unfolding_matrix']

        unf[h][sel] = TH2D(hnames[h] + "_unfolding", "", xbins, xmin, xmax, ybins, ymin, ymax);
        unf[h][sel].SetXTitle("i (bin)");
        unf[h][sel].SetYTitle("j (bin)");
        for i in range(xbins):
            for j in range(ybins):
                unf[h][sel].SetBinContent(i + 1, j + 1, v_unf[i][j])

        # Covariance matrix
        v_cov   = results['covariance_matrix']

        cov[h][sel] = TH2D(hnames[h] + "_covariance", "", xbins, xmin, xmax, ybins, ymin, ymax);
        cov[h][sel].SetXTitle("i (bin)");
        cov[h][sel].SetYTitle("j (bin)");
        for i in range(xbins):
            for j in range(ybins):
                cov[h][sel].SetBinContent(i + 1, j + 1, v_cov[i][j])

        # Data covariance matrix
        var[h][sel] = TH2D(hnames[h] + "_data_cov", "", xbins, xmin, xmax, ybins, ymin, ymax);
        var[h][sel].SetXTitle("i (bin)");
        var[h][sel].SetYTitle("j (bin)");
        for i in range(xbins):
            var[h][sel].SetBinContent(i + 1, i + 1, v_data[i]['y'])


        print("")



##
##  OUTPUT FILE
##

# Create file
outName = "unfolding_" + YEAR_STR + ".root"
outFile = TFile(outName, "RECREATE")


for h in range(H):
    outFile.mkdir(hnames[h])
    outFile.cd(hnames[h])

    # Write histograms
    mig[h]['4l'].Write()
    resp[h]['4l'].Write()
    unf[h]['4l'].Write()
    cov[h]['4l'].Write()
    var[h]['4l'].Write()

    gen[h]['4l'].Write()
    reco[h]['4l'].Write()
    data[h]['4l'].Write()
    bkg[h]['4l'].Write()
    result[h]['4l'].Write()



#   ##
#   ##  DRAW
#   ##

#   # Data/MC comparison
#   c_comp = TCanvas(hnames[h] + "_comparison", "", 1000, 1000)

#   data[h]['4l'].SetLineColor(1)
#   data[h]['4l'].SetLineWidth(2)
#   data[h]['4l'].SetMarkerColor(1)
#   data[h]['4l'].SetMarkerStyle(20)
#   data[h]['4l'].SetMarkerSize(2)

#   reco[h]['4l'].SetLineColor(8)
#   reco[h]['4l'].SetLineWidth(2)

#   gen[h]['4l'].SetLineColor(4)
#   gen[h]['4l'].SetLineWidth(2)

#   result[h]['4l'].SetLineColor(2)
#   result[h]['4l'].SetLineWidth(2)
#   result[h]['4l'].SetMarkerColor(2)
#   result[h]['4l'].SetMarkerStyle(22)
#   result[h]['4l'].SetMarkerSize(2)

#   if hnames[h] in ["b_l1p", "angle_z1leps"]:
#       l = TLegend(0.18, 0.68, 0.48, 0.98)
#   elif hnames[h] in ["cos_theta_z1"]:
#       l = TLegend(0.35, 0.68, 0.65, 0.98)
#   else:
#       l = TLegend(0.68, 0.68, 0.98, 0.98)
#   l.AddEntry(reco[h]['4l'], "Reco MC", "L")
#   l.AddEntry(gen[h]['4l'], "Gen MC", "L")
#   l.AddEntry(data[h]['4l'], "Smeared data", "LP")
#   l.AddEntry(result[h]['4l'], "Unfolded data", "LP")

#   c_comp.cd()
#   result[h]['4l'].SetMinimum(0);
#   result[h]['4l'].Draw("E1")
#   reco[h]['4l'].Draw("SAME")
#   gen[h]['4l'].Draw("SAME")
#   data[h]['4l'].Draw("E1 SAME")
#   result[h]['4l'].Draw("E1 SAME")
#   l.Draw()

#   c_comp.Write()


#   # Unfolded covariance matrix
#   c_cov = TCanvas(hnames[h] + "_comparison", "", 1000, 800)
#   c_cov.cd()
#   cov.Draw("COLZ")


outFile.Close()
print("Wrote output to", outName)
