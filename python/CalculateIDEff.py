from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TH1D
from Cuts2016 import *



##
##  SAMPLE INFO
##

selection   = [ "mumu", "ee", "4l", "4m", "2m2e", "4e"]
selTeX      = { "4m":r"\fM", "4e":r"\fE", "2m2e":r"\tMtE" }
channel     = { "mumu":3,       "ee":4,     "4l":5,     "4m":6,     "2m2e":7,   "4e":9  }
T = np.dtype([(sel, 'f4') for sel in selection])

year = "2016"



##
##  GET HISTOGRAMS
##

inPath, pref = EOS_PATH + "/BLT/" + year + "_update/", "eff2"
selHistName, matHistName = "SelectedEvents", "MatchedEvents"
selHist, matHist = {}, {}

for suff in ["zz_4l", "zjets_m-50"]:
    inName = pref + "_" + suff + ".root"
    inFile = TFile.Open(inPath + inName)
    print("Opened", inPath + inName)

    selHist[suff] = inFile.Get(selHistName + "_" + suff)
    selHist[suff].SetDirectory(0)
    matHist[suff] = inFile.Get(matHistName + "_" + suff)
    matHist[suff].SetDirectory(0)
    inFile.Close()
    print("Closed", inPath + inName)

    sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]
    selHist[suff].Scale(sf)
    matHist[suff].Scale(sf)

    


##
##  CALCULATE
##

matH = matHist["zz_4l"].Clone()
matH.Add(matHist["zjets_m-50"])
selH = selHist["zz_4l"].Clone()
selH.Add(selHist["zjets_m-50"])
effH = matH.Clone()
effH.Divide(selH)

eff = np.zeros(1, dtype=T)
for sel in selection:
    eff[sel] = effH.GetBinContent(channel[sel])

print("Channel", "\t", "Efficiency")
for sel in selection:
    print(sel, "\t\t", np.squeeze(eff[sel]))
