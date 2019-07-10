from __future__ import print_function
from __future__ import division
import sys

import numpy as np
from pyunfold import iterative_unfold
from pyunfold.callbacks import Logger
from scipy.stats import chi2

from ROOT import TFile, TH1, TH2, TH2D, TCanvas, TLegend

from PlotUtils import *

from Cuts2012 import *


np.set_printoptions(precision=4, suppress=True)

##
##  SAMPLE INFO
##

selection = ["4l"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f8'), ("y", 'f8'), ("ex", 'f8'), ("ey", 'f8'), ("b", 'f8')])

hnames = ["b_l1p"]
H = len(hnames)

year = sys.argv[1]
if year != YEAR_STR:
    print("Wrong year in header file")



##
##  INPUT FILES
##

#inPath = EOS_PATH + "/Unfolding/new/"
inPath = ""

# Migration and MC matrices
inName = inPath + "migration_" + YEAR_STR + "_zz_4l.root"
inFile = TFile(inName, "READ")
print("Opened", inName)

gen, reco = np.empty(H, dtype=T), np.empty(H, dtype=T)
mig = np.empty(H, dtype=T)
h = 0

for sel in selection:
    sf = INT_LUMI * 1000 * XSEC['zz_4l'] / NGEN['zz_4l']

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
#inPath = EOS_PATH + "/Histograms/new/"

# Muon file
muName = inPath + prefix + "_" + YEAR_STR + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = inPath + prefix + "_" + YEAR_STR + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)

# Get histograms
data = np.empty(H, dtype=T)
h = 0

for sel in selection:
    for hname in hnames:
        data[h][sel] = muFile.Get(sel + "/" + hname + "_" + MU_SUFF)
        data[h][sel].Add(elFile.Get(sel + "/" + hname + "_" + EL_SUFF))

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
muName = inPath + prefix + "_" + YEAR_STR + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = inPath + prefix + "_" + YEAR_STR + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)

# Get histograms
bkg = np.empty(H, dtype=T)
h = 0

for sel in selection:
    for hname in hnames:
        bkg[h][sel] = muFile.Get(sel + "/" + hname + "_" + MU_SUFF)
        bkg[h][sel].Add(elFile.Get(sel + "/" + hname + "_" + EL_SUFF))

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

    inName = inPath + prefix + "_" + year + "_" + suff + ".root"
    inFile = TFile.Open(inName)
    print("Opened", inName)

    # Get histograms
    for sel in selection:
        sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]

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

# Get nonprompt background
infile = "nonprompt" + YEAR_STR + ".npz"
npzfile = np.load(infile)
npt, npt_unc = npzfile['npt'], npzfile['npt_unc']

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
result, stat = np.empty(H, dtype=T), np.empty(H, dtype=T)
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
                                    ts = 'chi2',                ts_stopping = 0,
                                    max_iter = 200,             return_iterations = True,
#                                   max_iter = 20,              return_iterations = True,
                                    callbacks = [Logger()]
                                    )

#       iters = [5, 20, 50, 99, 100]
#       iters = [3, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        iters = range(1, 200)
        for it in iters:
            print("Iteration", it, end='\n\n')
#           print("Result", results.loc[it - 1, 'unfolded'], sep='\n', end='\n\n')
#           print("Covariance matrix", results.loc[it - 1, 'covariance_matrix'], sep='\n', end='\n\n\n')
            print(np.allclose(results.loc[it - 1, 'covariance_matrix'], results.loc[it, 'covariance_matrix'], rtol=1e-04, atol=1e-07))


'''
        ##
        ##  GET RESULTS
        ##

        # Unfolded result (total unc. and stat. only)
        v_result, v_stat            = np.zeros_like(v_gen, dtype=V), np.zeros_like(v_gen, dtype=V)
        v_result['y'], v_stat['y']  = results['unfolded'], results['unfolded']
        v_result['ey']  = np.sqrt(results['stat_err']**2, results['sys_err']**2)
        v_stat['ey']    = results['stat_err']

        o = s.start

        result[h][sel] = data[h][sel].Clone(hnames[h] + "_result")
        result[h][sel].SetYTitle("")
        for i in range(len(v_result)):
            result[h][sel].SetBinContent(i + o, v_result[i]['y'])
            result[h][sel].SetBinError(i + o, v_result[i]['ey'])

        stat[h][sel] = data[h][sel].Clone(hnames[h] + "_stat")
        stat[h][sel].SetYTitle("")
        for i in range(len(v_stat)):
            stat[h][sel].SetBinContent(i + o, v_stat[i]['y'])
            stat[h][sel].SetBinError(i + o, v_stat[i]['ey'])

        # Probability of final chi-squared
        chisq_prob = chi2.sf(results['ts_iter'] * len(v_resp['y']), len(v_resp['y']))
        print("Probability for", hnames[h], "is", chisq_prob)



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
    stat[h]['4l'].Write()


outFile.Close()
print("Wrote output to", outName)
'''
