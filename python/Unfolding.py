from __future__ import print_function
from __future__ import division
import sys

import numpy as np
from pyunfold import iterative_unfold
from pyunfold.callbacks import Logger

from ROOT import TFile, TH1, TH2, TCanvas, TLegend

from PlotUtils import *
from Cuts2017 import *


np.set_printoptions(precision=1, suppress=True)

##
##  SAMPLE INFO
##

selection = ["4l", "4m", "2m2e", "2e2m", "4e"]

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
        lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

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

        data[h][sel].SetDirectory(0)
        h = h + 1
    h = 0

muFile.Close()
elFile.Close()
print("Got data histograms")
print("")


'''
##
##  ACC * EFF
##

# Unscaled signal events
zzName = prefix + "_" + year + "_zz_4l.root"
zzFile = TFile(zzName, "READ")
print("Opened", zzName)

axe = np.empty(H, dtype=T)
h = 0

for sel in selection:
    if sel == "4l":
        continue

        for hname in hnames:
            axe[h][sel] = zzFile.Get(sel + "/" + hname + "_zz_4l")
            axe[h][sel].SetDirectory(0)
            axe[h][sel].SetName(hname + "_acc_x_eff")

            h = h + 1
        h = 0

zzFile.Close()


# Phase space events
ps = np.empty(H, dtype=T)
h = 0

psName = prefix + "_" + year + "_phase_space.root"
psFile = TFile(psName, "READ")
print("Opened", psName)

for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        ps[h][sel] = psFile.Get(sel + "/" + hname + "_phase_space")
        ps[h][sel].SetDirectory(0)

        h = h + 1
    h = 0

psFile.Close()

print("Got acc * eff histograms", "\n")
'''


##
##  ADD CHANNELS
##

# Get 4l and 2m2e, rebin 4e
for h in range(H):
#   for sample in [data, gen, reco, mig, ps, axe]:
    for sample in [data, gen, reco, mig]:
        sample[h]['2m2e'].Add(sample[h]['2e2m'])
        sample[h]['4l'] = sample[h]['2m2e'].Clone()
        sample[h]['4l'].Add(sample[h]['4m'])
        sample[h]['4l'].Add(sample[h]['4e'])


'''
##
##  SCALING
##

for sel in selection:
    if sel == "2e2m":
        continue

    for h in range(H):
        axe[h][sel].Divide(ps[h][sel])

        for sample in [data, gen, reco, mig]:
'''






####
####
####    LOOP OVER DISTS
####
####


# Store results in histograms (for now)
result = np.empty(H, dtype=T)
resp, unf, cov = np.empty(H, dtype=T), np.empty(H, dtype=T), np.empty(H, dtype=T)


for sel in ["4l"]:
    for h in range(H):

        ##
        ##  GET BIN CONTENT
        ##

        # Data
        v_data = np.zeros(data[h][sel].GetNbinsX() + 2, dtype=V)
        for i in range(len(v_data)):
            v_data[i]['x']  = data[h][sel].GetBinCenter(i)
            v_data[i]['y']  = data[h][sel].GetBinContent(i)
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
        print("Condition number of", hname[h], "is", cond_num)



        ##
        ##  PERFORM UNFOLDING
        ##

        results = iterative_unfold( data = v_data['y'],         data_err = v_data['ey'],
                                    response = v_resp['y'],     response_err = v_resp['ey'],
                                    efficiencies = v_eff['y'],  efficiencies_err = v_eff['ey'],
#                                   ts = 'chi2',                ts_stopping = 0.05,
                                    ts = 'chi2',                ts_stopping = 0.1,
                                    callbacks = [Logger()]
                                    )



        ##
        ##  GET RESULTS
        ##

        # Unfolded result
        v_result        = np.zeros_like(v_data, dtype=V)
        v_result['y']   = results['unfolded']
        v_result['ey']  = np.sqrt(results['stat_err']**2, results['sys_err']**2)

        result[h][sel] = data[h][sel].Clone(hnames[h] + "_result");
        for i in range(len(v_result)):
            result[h][sel].SetBinContent(i, v_result[i]['y'])
            result[h][sel].SetBinError(i, v_result[i]['ey'])


        # Response matrix
        resp[h][sel] = mig[h][sel].Clone(hnames[h] + "_response");
        resp[h][sel].Reset()
        resp[h][sel].SetTitle("")
        for i in range(np.size(v_resp, 0)):
            for j in range(np.size(v_resp, 1)):
                resp[h][sel].SetBinContent(i, j, v_resp[i][j]['y'])
                resp[h][sel].SetBinError(i, j, v_resp[i][j]['ey'])

        # Unfolding matrix
        v_unf   = results['unfolding_matrix']

        unf[h][sel] = mig[h][sel].Clone(hnames[h] + "_unfolding");
        unf[h][sel].Reset()
        unf[h][sel].SetTitle("")
        for i in range(np.size(v_unf, 0)):
            for j in range(np.size(v_unf, 1)):
                unf[h][sel].SetBinContent(i, j, v_unf[i][j])

        # Covariance matrix
        v_cov   = results['covariance_matrix']

        cov[h][sel] = mig[h][sel].Clone(hnames[h] + "_covariance");
        cov[h][sel].Reset()
        cov[h][sel].SetTitle("")
        for i in range(np.size(v_cov, 0)):
            for j in range(np.size(v_cov, 1)):
                cov[h][sel].SetBinContent(i, j, v_cov[i][j])



        ##
        ##  BOTTOM-LINE TEST
        ##

        dof = len(v_data[1:-1])
        chi2_smr = np.sum(((v_data[1:-1]['y'] - v_reco[1:-1]['y']) / v_data[1:-1]['ey']) ** 2) / dof
        print("Smeared chi square is", chi2_smr)

        num = v_result[1:-1]['y'] - v_gen[1:-1]['y']
        chi2_unf = np.dot(num, np.dot(np.linalg.pinv(v_cov[1:-1,1:-1]), num.T)) / dof
        print("Unfolded chi square is", chi2_unf)
        




##
##  DRAW
##

c = TCanvas("canvas", "", 800, 600)

data[0]['4l'].SetLineColor(1)
data[0]['4l'].SetLineWidth(2)
data[0]['4l'].SetMarkerColor(1)
data[0]['4l'].SetMarkerStyle(20)
data[0]['4l'].SetMarkerSize(2)

reco[0]['4l'].SetLineColor(8)
reco[0]['4l'].SetLineWidth(2)

gen[0]['4l'].SetLineColor(4)
gen[0]['4l'].SetLineWidth(2)

result[0]['4l'].SetLineColor(2)
result[0]['4l'].SetLineWidth(2)
result[0]['4l'].SetMarkerColor(2)
result[0]['4l'].SetMarkerStyle(22)
result[0]['4l'].SetMarkerSize(2)

l = TLegend(0.78, 0.68, 0.98, 0.98)
l.AddEntry(reco[0]['4l'], "Reco", "L")
l.AddEntry(gen[0]['4l'], "Gen", "L")
l.AddEntry(data[0]['4l'], "Data", "LP")
l.AddEntry(result[0]['4l'], "Result", "LP")

c.cd()
result[0]['4l'].SetMinimum(0);
result[0]['4l'].Draw("E1")
reco[0]['4l'].Draw("SAME")
gen[0]['4l'].Draw("SAME")
data[0]['4l'].Draw("E1 SAME")
result[0]['4l'].Draw("E1 SAME")
l.Draw()



##
##  OUTPUT FILE
##

# Create file
outName = "unfolding_py_" + YEAR_STR + ".root"
outFile = TFile(outName, "RECREATE")
outFile.cd()

mig[0]['4l'].Write()
resp[0]['4l'].Write()
unf[0]['4l'].Write()
cov[0]['4l'].Write()

gen[0]['4l'].Write()
reco[0]['4l'].Write()
data[0]['4l'].Write()
result[0]['4l'].Write()
c.Write()

outFile.Close()
print("Wrote output to", outName)
