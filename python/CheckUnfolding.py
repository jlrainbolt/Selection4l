from __future__ import print_function
from __future__ import division

import numpy as np
from pyunfold import iterative_unfold
from pyunfold.callbacks import Logger
from scipy.stats import chi2

from ROOT import TFile, TH1, TH2, TH2D, TCanvas, TLegend

from PlotUtils import *

from Cuts2012 import *


np.set_printoptions(precision=5, suppress=True, linewidth=110)

##
##  SAMPLE INFO
##

selection = ["4l"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f8'), ("y", 'f8'), ("ex", 'f8'), ("ey", 'f8'), ("b", 'f8')])

hnames = ["b_z1m", "b_z2m", "b_l1p", "b_ttm", "cos_theta_z1", "cos_theta_z2",
        "angle_z1leps", "angle_z2leps", "angle_z1l2_z2", "sin_phi"]
H = len(hnames)

draw_opt = ['unc', 'min', 'max']



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

    inName = inPath + prefix + "_" + YEAR_STR + "_" + suff + ".root"
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

        if hnames[h] == "b_z1m":
            if YEAR_STR in ["2012", "2016"]:
                s = slice(2, -2)
            else:
                s = slice(2, -1)
        elif hnames[h] == "b_z2m":
            if YEAR_STR == "2012":
                s = slice(1, -1)
            else:
                s = slice(1, None)
        elif hnames[h] in ["b_l1p"]:
            if YEAR_STR in ["2012", "2017"]:
                s = slice(1, -1)
            else:
                s = slice(0, -1)
        elif hnames[h] == "b_ttm":
            if YEAR_STR in ["2017", "2018"]:
                s = slice(1, -1)
            else:
                s = slice(1, None)
        elif hnames[h] == "angle_z1l2_z2":
            s = slice(1, -2)
        elif hnames[h] == "angle_z1leps":
            if YEAR_STR == "2012":
                s = slice(4, -1)
            else:
                s = slice(3, -1)
        elif hnames[h] in ["cos_theta_z1", "cos_theta_z2", "angle_z2leps", "sin_phi"]:
            s = slice(1, -1)

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
        print("")
        print("Condition number for", hnames[h], "is", cond_num)



        ##
        ##  PERFORM UNFOLDING
        ##

        max_iter = 100

        chisq_nu = chi2.isf(0.997, len(v_resp['y'])) / len(v_resp['y'])
        results = iterative_unfold( data = v_data['y'],         data_err = v_data['ey'],
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

            unc = np.sqrt(results.loc[it, 'stat_err']**2, results.loc[it, 'sys_err']**2)
#           unc = np.mean(unc / results.loc[it, 'unfolded'])
            y['unc'] = np.mean(unc)

            y['max'] = np.max(results.loc[it, 'unfolded'])
            y['min'] = np.min(results.loc[it, 'unfolded'])

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
            figname = opt + "_" + hnames[h] + "_" + YEAR_STR + ".pdf"
            fig.savefig(figname)

            print("Saved plot as", figname)
