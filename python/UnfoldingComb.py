from __future__ import print_function
from __future__ import division
import sys

import numpy as np
from numpy.linalg import multi_dot, pinv
from pyunfold import iterative_unfold
#from pyunfold.callbacks import Logger
#from scipy.stats import chi2

from ROOT import TFile, TH1, TH1D, TH2, TH2D, TCanvas, TLegend

from Cuts2018 import *
from PlotUtils import *


np.set_printoptions(precision=2, suppress=True, linewidth=99)

##
##  SAMPLE INFO
##

selection   = ["4l"]
period      = ["2018", "2017", "2016", "2012"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f8'), ("y", 'f8'), ("ex", 'f8'), ("ey", 'f8'), ("b", 'f8')])

hnames = ["b_z1m", "b_z2m", "b_l1p", "b_ttm", "cos_theta_z1", "cos_theta_z2",
            "angle_z1leps", "angle_z2leps", "angle_z1l2_z2", "sin_phi", "sin_phi_10"]
H = len(hnames)


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

# Get nonprompt background
npt, npt_unc = {}, {}
for year in period:
    infile = "nonprompt" + year + ".npz"
    npzfile = np.load(infile)
    npt[year], npt_unc[year] = npzfile['npt'], npzfile['npt_unc']



##
##  INPUT FILES
##

inPath = ""

# Migration and MC matrices
gen, reco = np.empty(H, dtype=T), np.empty(H, dtype=T)
mig = np.empty(H, dtype=T)

for year in period:
    inName = inPath + "migration_" + year + "_zz_4l.root"
    inFile = TFile(inName, "READ")
    print("Opened", inName)
    h = 0

    for sel in selection:
        sf = lumi[year] * 1000 * xsec[year]["zz_4l"] / ngen[year]["zz_4l"]

        for hname in hnames:
            gen_ = inFile.Get(sel + "/" + hname + "_gen")
            gen_.SetDirectory(0)
            gen_.Scale(sf)
    
            reco_ = inFile.Get(sel + "/" + hname + "_reco")
            reco_.SetDirectory(0)
            reco_.Scale(sf)
    
            mig_ = inFile.Get(sel + "/" + hname + "_2d")
            mig_.SetDirectory(0)
            mig_.Scale(sf)

            if year == "2018":
                gen[h][sel] = gen_
                reco[h][sel] = reco_
                mig[h][sel] = mig_
            else:
                gen[h][sel].Add(gen_)
                reco[h][sel].Add(reco_)
                mig[h][sel].Add(mig_)
    
            h = h + 1
        h = 0
    inFile.Close()

print("Got migration and MC histograms", "\n")



##
##  DATA
##

prefix = "4l"
data = np.empty(H, dtype=T)

for year in period:

    # Muon file
    muName = inPath + prefix + "_" + year + "_muon_" + year + ".root"
    muFile = TFile(muName, "READ")
    print("Opened", muName)

    # Electron file
    elName = inPath + prefix + "_" + year + "_electron_" + year + ".root"
    elFile = TFile(elName, "READ")
    print("Opened", elName)

    h = 0
    for sel in selection:
        for hname in hnames:
            data_ = muFile.Get(sel + "/" + hname + "_muon_" + year)
            data_.Add(elFile.Get(sel + "/" + hname + "_electron_" + year))

            data_.SetName(hnames[h] + "_data")
            data_.SetDirectory(0)

            if year == "2018":
                data[h][sel] = data_
            else:
                data[h][sel].Add(data_)

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
bkg = np.empty(H, dtype=T)

for year in period:

    # Muon file
    muName = inPath + prefix + "_" + year + "_muon_" + year + ".root"
    muFile = TFile(muName, "READ")
    print("Opened", muName)

    # Electron file
    elName = inPath + prefix + "_" + year + "_electron_" + year + ".root"
    elFile = TFile(elName, "READ")
    print("Opened", elName)

    h = 0
    for sel in selection:
        for hname in hnames:
            bkg_ = muFile.Get(sel + "/" + hname + "_muon_" + year)
            bkg_.Add(elFile.Get(sel + "/" + hname + "_electron_" + year))
            bkg_.SetName(hnames[h] + "_bkg")
            bkg_.SetDirectory(0)
        
            sf = npt[year][sel] / bkg_.Integral()
            bkg_.Scale(sf)
            err_sf = npt_unc[year][sel] / np.sqrt(bkg_.GetSumw2().GetSum())
            for i in range(1, bkg_.GetNbinsX()):
                bkg_.SetBinError(i, err_sf * bkg_.GetBinError(i))

            if year == "2018":
                bkg[h][sel] = bkg_
            else:
                bkg[h][sel].Add(bkg_)

            h = h + 1
        h = 0
    muFile.Close()
    elFile.Close()
print("Got nonprompt background histograms")
print("")


# Monte Carlo
prefix = "4l"
h = 0

for year in period:

    # Loop over all samples
    for suff in mc_suff[year]:
        if suff in ["zz_4l", "zjets_m-50", "tt_2l2nu", "ttbar"]:
            continue

        inName = inPath + prefix + "_" + year + "_" + suff + ".root"
        inFile = TFile.Open(inName)
        print("Opened", inName)

        # Get histograms
        for sel in selection:
            sf = lumi[year] * 1000 * xsec[year][suff] / ngen[year][suff]

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






####
####
####    LOOP OVER DISTS
####
####


# Store results in histograms
result, stat = np.empty(H, dtype=T), np.empty(H, dtype=T)
resp, unf, cov = np.empty(H, dtype=T), np.empty(H, dtype=T), np.empty(H, dtype=T)
var = np.empty(H, dtype=T)


# Store unfolding info in arrays
bins, n_iter, cond_num = np.empty(H), np.empty(H), np.empty(H)
chisq_smr, chisq_unf = np.empty(H), np.empty(H)

cov_mat = []


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
            v_data[i]['ex'] = data[h][sel].GetBinContent(i) # dummy for unfolding
            if v_data[i]['ex'] < 1:
                v_data[i]['ex'] = reco[h][sel].GetBinContent(i)
            v_data[i]['ey'] = data[h][sel].GetBinErrorUp(i)

        # Gen MC
        v_gen = np.zeros(gen[h][sel].GetNbinsX() + 2, dtype=V)
        for i in range(len(v_gen)):
            v_gen[i]['x']   = gen[h][sel].GetBinLowEdge(i)
            v_gen[i]['y']   = gen[h][sel].GetBinContent(i)
            v_gen[i]['ey']  = gen[h][sel].GetBinErrorUp(i)

        # Reco MC (only used for tests)
        v_reco = np.zeros(reco[h][sel].GetNbinsX() + 2, dtype=V)
        for i in range(len(v_gen)):
            v_reco[i]['x']  = reco[h][sel].GetBinLowEdge(i)
            v_reco[i]['y']  = reco[h][sel].GetBinContent(i)
            v_reco[i]['ey'] = reco[h][sel].GetBinErrorUp(i)

        # Migrations
        v_mig = np.zeros([mig[h][sel].GetNbinsX() + 2, mig[h][sel].GetNbinsY() + 2], dtype=V)
        for i in range(np.size(v_mig, 0)):
            for j in range(np.size(v_mig, 1)):
                v_mig[i][j]['y']    = mig[h][sel].GetBinContent(i, j)
                v_mig[i][j]['ey']   = mig[h][sel].GetBinErrorUp(i, j)



        ##
        ##  CREATE MATRICES
        ##

        # Efficiencies (assume 1)
        v_eff       = np.zeros_like(v_gen, dtype=V)
        v_eff['y']  = np.ones_like(v_gen['y'])
        col_sums    = v_mig['y'].sum(axis=0)

#       print(hnames[h])
#       print(v_data['y'])
#       print(v_reco['y'])
#       print(v_eff['y'])
#       print(v_gen['y'])

        # Slicing
        if hnames[h] == "b_z1m":
            s = slice(2, -1)
            v_eff['y'][3] = col_sums[2] / (col_sums[2] + col_sums[1] + col_sums[0])
            v_eff['y'][-2] = col_sums[-2] / (col_sums[-2] + col_sums[-1])
        elif hnames[h] == "b_z2m":
            s = slice(1, None)
            v_eff['y'][1] = col_sums[1] / (col_sums[1] + col_sums[0])
        elif hnames[h] == "b_l1p":
            s = slice(0, -1)
            v_eff['y'][-2] = col_sums[-2] / (col_sums[-2] + col_sums[-1])
        elif hnames[h] == "angle_z1leps":
            s = slice(2, -1)
            v_eff['y'][2] = col_sums[2] / (col_sums[2] + col_sums[1] + col_sums[0])
            v_eff['y'][-2] = col_sums[-2] / (col_sums[-2] + col_sums[-1])
        elif hnames[h] == "angle_z1l2_z2":
            s = slice(1, -1)
            v_eff['y'][1] = col_sums[1] / (col_sums[1] + col_sums[0])
            v_eff['y'][-2] = col_sums[-2] / (col_sums[-2] + col_sums[-1])
        elif hnames[h] in ["cos_theta_z1", "cos_theta_z2", "angle_z2leps", "sin_phi", "sin_phi_10"]:
            s = slice(1, -1)
            v_eff['y'][1] = col_sums[1] / (col_sums[1] + col_sums[0])
            v_eff['y'][-2] = col_sums[-2] / (col_sums[-2] + col_sums[-1])
        else:
            s = slice(0, None)

        v_data = v_data[s]
        v_gen = v_gen[s]
        v_reco = v_reco[s]
        v_mig = v_mig[s, s]
        v_eff = v_eff[s]

        bins[h] = len(v_data)
#       print(v_eff['y'])
#       print(v_data['y'])
#       print(v_reco['y'])
#       print(v_gen['y'])
#       print("")


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
        cond_num[h] = np.max(sigma) / np.min(sigma)
        print("Condition number for", hnames[h], "is", cond_num[h])



        ##
        ##  PERFORM UNFOLDING
        ##

        results = iterative_unfold( data = v_data['ex'],        data_err = v_data['ey'],
#       results = iterative_unfold( data = v_data['y'],         data_err = v_data['ey'],
                                    response = v_resp['y'],     response_err = v_resp['ey'],
                                    efficiencies = v_eff['y'],  efficiencies_err = v_eff['ey'],
                                    ts = 'conv',                ts_stopping = 1e-2,
                                    max_iter=1000,              #callbacks = [Logger()]
                                    )
        n_iter[h] = results['num_iterations']
        print("Iterations:", n_iter[h])

        cov_mat.append(results['covariance_matrix'])



        ##
        ##  GET RESULTS
        ##

        # Unfolded result (total unc. and stat. only)
        v_result, v_stat= np.zeros_like(v_gen, dtype=V), np.zeros_like(v_gen, dtype=V)
        v_result['y']   = np.dot(v_data['y'], results['unfolding_matrix'])
        v_stat['y']     = results['unfolded']
        v_result['ey']  = np.sqrt(results['stat_err']**2, results['sys_err']**2)
        v_stat['ey']    = results['stat_err']

        print(v_stat['y'])
        print(v_result['y'])

        o = s.start

        if hnames[h] in ["b_z1m", "angle_z1leps"]:
            xmin = data[h][sel].GetBinLowEdge(o)
            xmax = data[h][sel].GetBinLowEdge(data[h][sel].GetNbinsX() + 1)
            print(xmin, xmax)

            rebin = TH1D(hnames[h] + "_rebin", "", bins[h], xmin, xmax)

            for i in range(len(v_data)):
                rebin[h][sel].SetBinContent(i + 1, data[h][se].GetBinContent(i + o))
                rebin[h][sel].DetBinError(i + 1, data[h][se].GetBinError(i + o))

            data[h][sel].Delete()
            data[h][sel] = rebin[h][sel]
            rebin[h][sel].SetName(hnames[h] + "_data")

            o = 1
            

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


        # Bottom line test 
        diff_unf = v_result['y'] - v_gen['y']
        chisq_unf[h] = multi_dot([diff_unf.T, pinv(v_cov), diff_unf])
        print("Unfolded chi^2:", chisq_unf[h])

        diff_smr = v_data['y'] - v_reco['y']
        cov_smr = np.diag(v_data['y'])
        chisq_smr[h] = multi_dot([diff_smr.T, pinv(cov_smr), diff_smr])
        print("Smeared chi^2:", chisq_smr[h])
        print("Pass bottom line test:", chisq_smr[h] > chisq_unf[h])


        print("")



##
##  OUTPUT FILE
##

# Create file
outName = "unfolding_comb.root"
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


# Save arrays
names = np.array(hnames, dtype=object)   # strings of (ugly) names

print("")
outfile = "unfoldingcomb.npz"
np.savez(outfile, names=names, bins=bins, cond_num=cond_num, n_iter=n_iter,
        chisq_smr=chisq_smr, chisq_unf=chisq_unf,
        cov_mat_b_z1m=cov_mat[0], cov_mat_b_z2m=cov_mat[1], cov_mat_b_l1p=cov_mat[2],
        cov_mat_b_ttm=cov_mat[3], cov_mat_cos_theta_z1=cov_mat[4],
        cov_mat_cos_theta_z2=cov_mat[5], cov_mat_angle_z1leps=cov_mat[6],
        cov_mat_angle_z2leps=cov_mat[7], cov_mat_angle_z1l2_z2=cov_mat[8],
        cov_mat_sin_phi=cov_mat[9], cov_mat_sin_phi_10=cov_mat[10])


print("Wrote arrays to", outfile)

