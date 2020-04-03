from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TH1D, TCanvas



##
##  LOAD DATA
##

# Single-parameter result
infile = "combination_1.npz"
npzfile = np.load(infile)

alpha_total_1 = npzfile['alpha_total']
sigma_total_1 = npzfile['sigma_total']


# Three-parameter result
infile = "combination_3.npz"
npzfile = np.load(infile)

alpha_total_3 = npzfile['alpha_total']
sigma_total_3 = npzfile['sigma_total']


# Four-parameter result
infile = "combination_4.npz"
npzfile = np.load(infile)

alpha_total_4 = npzfile['alpha_total']
sigma_total_4 = npzfile['sigma_total']


# Twelve-parameter result ("data")
infile = "combination_12.npz"
npzfile = np.load(infile)

data = npzfile['alpha_total']
unc = npzfile['sigma_total']



##
##  PULLS
##

print("2012, 2016, 2017, 2018")
print("4m, 2m2e, 4e")
print("")

print("Data")
print(data)
print("")


# Single parameter

fit_1 = np.squeeze(alpha_total_1)
res_1 = data - fit_1
pulls_1 = res_1 / unc

print("One parameter")
print(pulls_1)
print("")


# Three parameters
fit_3 = np.tile(alpha_total_3, 4)
res_3 = data - fit_3
pulls_3 = res_3 / unc

print("Three parameters")
print(pulls_3)
print("")


# Four parameters
fit_4 = np.repeat(alpha_total_4, 3)
res_4 = data - fit_4
pulls_4 = res_4 / unc

print("Four parameters")
print(pulls_4)
print("")



##
##  HISTOGRAMS
##

outname = "pulls.root"
outfile = TFile(outname, "RECREATE")
outfile.cd()

h1 = TH1D("hPulls_SM", "SM parameterization", 12, -3, 3)
h3 = TH1D("hPulls_BSM", "BSM parameterization", 12, -3, 3)
h4 = TH1D("hPulls_yearly", "Yearly parameterization", 12, -3, 3)

for i in range(len(data)):
    h1.Fill(np.squeeze(pulls_1[i]))
    h3.Fill(np.squeeze(pulls_3[i]))
    h4.Fill(np.squeeze(pulls_4[i]))

for h in [h1, h3, h4]:
    h.SetFillColor(9)
    h.SetLineColor(9)
    h.SetStats(0)
    h.SetXTitle("Pull")

c1 = TCanvas("pulls_SM", "", 500, 500)
c3 = TCanvas("pulls_BSM", "", 500, 500)
c4 = TCanvas("pulls_yearly", "", 500, 500)

c1.cd()
h1.Draw("HIST")
c3.cd()
h3.Draw("HIST")
c4.cd()
h4.Draw("HIST")

for h in [h1, h3, h4]:
    h.Write()
for c in [c1, c3, c4]:
    c.SaveAs(".pdf")
    c.Write()

outfile.Close()

print("Wrote histograms to", outname)
