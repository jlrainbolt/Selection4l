from __future__ import print_function
from __future__ import division

import numpy as np
from pyunfold import iterative_unfold
from pyunfold.callbacks import Logger
from scipy.stats import chi2

from ROOT import TFile, TH1, TH2, TH2D, TCanvas, TLegend

from PlotUtils import *

from Cuts2018 import *


np.set_printoptions(precision=5, suppress=True, linewidth=110)

##
##  SAMPLE INFO
##

selection = ["4l"]
period = ["2018", "2017", "2016", "2012"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f8'), ("y", 'f8'), ("ex", 'f8'), ("ey", 'f8'), ("b", 'f8')])

hnames = ["b_z1m", "b_z2m", "b_l1p", "b_ttm", "cos_theta_z1", "cos_theta_z2",
        "angle_z1leps", "angle_z2leps", "angle_z1l2_z2", "sin_phi"]
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


draw_opt = ['unc', 'min', 'max']



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

prefix = "bkg"
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
            v_data[i]['ex'] = data[h][sel].GetBinContent(i) # dummy for unfolding
            if v_data[i]['ex'] < 1:
                v_data[i]['ex'] = reco[h][sel].GetBinContent(i)
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
        if hnames[h] == "b_z2m":
            s = slice(1, None)
            v_eff['y'][1] = col_sums[1] / (col_sums[1] + col_sums[0])
        elif hnames[h] == "b_l1p":
            s = slice(0, -1)
            v_eff['y'][-2] = col_sums[-2] / (col_sums[-2] + col_sums[-1])
        elif hnames[h] == "angle_z1leps":
            s = slice(2, -1)
            v_eff['y'][2] = col_sums[2] / (col_sums[2] + col_sums[1] + col_sums[0])
            v_eff['y'][-2] = col_sums[-2] / (col_sums[-2] + col_sums[-1])
        elif hnames[h] in ["b_z1m", "angle_z1l2_z2", "cos_theta_z1", "cos_theta_z2", "angle_z2leps", "sin_phi", "sin_phi_10"]:
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



        ##
        ##  CONDITION NUMBER
        ##

        sigma = np.linalg.svd(v_resp['y'], compute_uv = False)
        cond_num = np.max(sigma) / np.min(sigma)
        print("")
        print("Condition number for", hnames[h], "is", cond_num)



        ##
        ##  PERFORM UNFOLDING
        ##

        max_iter = 100

        results = iterative_unfold( data = v_data['ex'],        data_err = v_data['ey'],
                                    response = v_resp['y'],     response_err = v_resp['ey'],
                                    efficiencies = v_eff['y'],  efficiencies_err = v_eff['ey'],
                                    ts = 'conv',                ts_stopping = 0,
                                    max_iter = max_iter,        return_iterations = True,
#                                   callbacks = [Logger()]
                                    )

        y = {}
        true_x, true_y = [], {opt:[] for opt in draw_opt}
        false_x, false_y = [], {opt:[] for opt in draw_opt}

        for it in range(max_iter):
#           print("Iteration", it + 1, end='\n\n')
#           print("Result", results.loc[it, 'unfolded'], sep='\n', end='\n\n')
#           print("Covariance matrix", results.loc[it, 'covariance_matrix'], sep='\n', end='\n\n\n')
#           if it > 0:
#               print(np.isclose(results.loc[it - 1, 'covariance_matrix'], results.loc[it, 'covariance_matrix'], rtol=1e-2, atol=1e-2))

            result = np.dot(v_data['y'], results.loc[it, 'unfolding_matrix'])

            unc = np.sqrt(results.loc[it, 'stat_err']**2, results.loc[it, 'sys_err']**2)
#           unc = np.mean(unc / result)
            y['unc'] = np.mean(unc)

            y['max'] = np.max(result)
            y['min'] = np.min(result)

            if it == 0:
                converged = np.allclose(np.zeros_like(results.loc[it, 'covariance_matrix']),
                        results.loc[it, 'covariance_matrix'], rtol=1e-2, atol=1e-2)
            else:
                converged = np.allclose(results.loc[it - 1, 'covariance_matrix'],
                        results.loc[it, 'covariance_matrix'], rtol=1e-2, atol=1e-2)

            if converged:
                true_x.append(it+1)
                for opt in draw_opt:
                    true_y[opt].append(y[opt])
            else:
                false_x.append(it+1)
                for opt in draw_opt:
                    false_y[opt].append(y[opt])


        ##
        ##  PLOT
        ##

        if hnames[h] == "b_l1p":
            xtitle = r"$p_{\ell_{1}}$"
        elif hnames[h] == "b_z1m":
            xtitle = r"$m_{\mathrm{Z}_{1}}$"
        elif hnames[h] == "b_z2m":
            xtitle = r"$m_{\mathrm{Z}_{2}}$"
        elif hnames[h] == "b_ttm":
            xtitle = r"$m_{\ell_{2,3,4}}$"
        elif hnames[h] == "angle_z1l2_z2":
            xtitle = r"$\beta$"
        elif hnames[h] == "angle_z1leps":
            xtitle = r"$\alpha_{\mathrm{Z}_{1}}$"
        elif hnames[h] == "angle_z2leps":
            xtitle = r"$\alpha_{\mathrm{Z}_{2}}$"
        elif hnames[h] == "cos_theta_z1":
            xtitle = r"$\cos\theta_{\mathrm{Z}_{1}}$"
        elif hnames[h] == "cos_theta_z2":
            xtitle = r"$\cos\theta_{\mathrm{Z}_{2}}$"
        elif hnames[h] == "sin_phi":
            xtitle = r"$\sin\phi$"



        for opt in draw_opt:
            fig = plt.figure(figsize = (8,3))
            ax = plt.axes()

            p_false = ax.scatter(false_x, false_y[opt], marker='o', c = 'r')
            p_true = ax.scatter(true_x, true_y[opt], marker='o', c = 'b')

            ax.set_title(xtitle)
            ax.set_title(r"cond$(K) = {:,.2f}$".format(cond_num), loc='right')

            y_min, y_max = ax.get_ylim()
            if opt == 'min':
                ax.set_ylim(y_min - (y_max - y_min)/4, y_max)
            else:
                ax.set_ylim(y_min, y_max + (y_max - y_min)/4)
                
            y_offset = (y_max - y_min)/4
            if opt != 'min':
                y_offset *= -1

            ax.annotate(r"$n = {:,.0f}$".format(true_x[0]), xy=(true_x[0], true_y[opt][0]),
                    xytext=(true_x[0], true_y[opt][0] + y_offset), fontsize='large',
                    horizontalalignment='center', arrowprops=dict(arrowstyle='->'))

            ax.set_xlabel(r"Iteration $n$")
            if opt == 'unc':
                ax.set_ylabel('Mean unc. in result')
            elif opt == 'max':
                ax.set_ylabel('Max. bin of result')
            elif opt == 'min':
                ax.set_ylabel('Min. bin of result')

            if opt == 'min':
                leg_loc = 'upper right'
            else:
                leg_loc = 'lower right'

            ax.legend([r"Stopping criterion not satisfied for $n$ and $n-1$",
                r"Stopping criterion satisfied for $n$ and $n-1$"],
                loc = leg_loc, numpoints = 1, frameon = False)

            fig.tight_layout()
            figname = opt + "_" + hnames[h] + "_comb.pdf"
            fig.savefig(figname)

            print("Saved plot as", figname)
