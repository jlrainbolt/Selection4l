from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TTree, TH1D
from secret_number import *

#from Cuts2017 import *
from Cuts2016 import *



##
##  SAMPLE INFO
##

selection   = [ "mumu", "ee", "4m", "2m2e", "2e2m", "4e"]
selTeX      = { "mumu":r"\MM",  "ee":r"\EE",
                "4l":r"4\Pell", "4m":r"4\PGm",  "4e":r"4\Pe",   "2m2e":r"2\PGm 2\Pe"
                }
selDef      = { "mumu":"MM",    "ee":"EE",  "4l":"4L",  "4m":"4M",  "4e":"4E",  "2m2e":"2M2E"   }
channel     = { "mumu":3,       "ee":4,     "4m":6,     "2m2e":7,   "2e2m":8,   "4e":9  }
T = np.dtype([(sel, 'f4') for sel in selection])



##
##  NUMERATOR
##

inPath = EOS_PATH + "/Selected/" + YEAR_STR + "/"
prefix = "selected"

# Drell-Yan file
dyName = prefix + "_zjets_m-50.root"
dyFile = TFile.Open(inPath + dyName)
print("Opened", inPath + dyName)

# ZZTo4L file
zzName = prefix + "_zz_4l.root"
zzFile = TFile.Open(inPath + zzName)
print("Opened", inPath + zzName)

# Get yields
num, num_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    if sel in ["mumu", "ee"]:
        tree = dyFile.Get(sel + "_zjets_m-50")
    elif sel in ["4e", "4m", "2m2e", "2e2m"]:
        tree = zzFile.Get(sel + "_zz_4l")

    hist = TH1D("hist", "", 1, 0, 2)
    weight = "weight/trigWeight/qtWeight"

    tree.Draw("1>>hist", "!hasTauDecay * " + weight, "goff")

    num[sel] = hist.Integral()
    num_unc[sel] = np.sqrt(tree.GetEntries())

dyFile.Close()
zzFile.Close()

num["2m2e"] = num["2m2e"] + num["2e2m"]
num["2e2m"] = 0
num_unc["2m2e"] = np.sqrt(num_unc["2m2e"] ** 2 + num_unc["2e2m"] ** 2)
num_unc["2e2m"] = 0



##
##  DENOMINATOR
##

inPath = EOS_PATH + "/Systematics/" + YEAR_STR + "/"
prefix = "triggerEff"

# Drell-Yan file
dyName = prefix + "_zjets_m-50.root"
dyFile = TFile.Open(inPath + dyName)
print("Opened", inPath + dyName)

# ZZTo4L file
zzName = prefix + "_zz_4l.root"
zzFile = TFile.Open(inPath + zzName)
print("Opened", inPath + zzName)

# Get yields
denom, denom_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

for sel in selection:
    if sel in ["2e2m"]:
        continue
    elif sel in ["mumu", "ee"]:
        tree = dyFile.Get(sel + "_zjets_m-50")
    elif sel in ["4e", "4m", "2m2e"]:
        tree = zzFile.Get(sel + "_zz_4l")

    hist = TH1D("hist", "", 1, 0, 2)
    weight = "weight/trigWeight/qtWeight"

    tree.Draw("1>>hist", "!hasTauDecay * " + weight, "goff")

    denom[sel] = hist.Integral()
    denom_unc[sel] = np.sqrt(tree.GetEntries())

dyFile.Close()
zzFile.Close()



##
##  CALCULATE
##

ratio, ratio_unc = np.zeros(1, dtype=T), np.zeros(1, dtype=T)

print(YEAR_STR)

for sel in selection:
    if sel == "2e2m":
        continue

    ratio[sel] = num[sel] / denom[sel]
    ratio_unc[sel] = ratio[sel] * np.sqrt((num_unc[sel] / num[sel]) ** 2 + (denom_unc[sel] / denom[sel]) ** 2)

    print(sel, "\t", ratio[sel], "+-", ratio_unc[sel], "pass trigger")
