from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TH1D, TCanvas, gStyle



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

alpha_total_12 = npzfile['alpha_total']
sigma_total_12 = npzfile['sigma_total']


# "Observed" branching fractions
infile = "bf_measurements.npz"
npzfile = np.load(infile)

bf, bf_stat, bf_syst = {}, {}, {}

bf["2012"], bf_stat["2012"], bf_syst["2012"] = npzfile["bf_2012"], npzfile["bf_stat_2012"], npzfile["bf_syst_2012"]
bf["2016"], bf_stat["2016"], bf_syst["2016"] = npzfile["bf_2016"], npzfile["bf_stat_2016"], npzfile["bf_syst_2016"]
bf["2017"], bf_stat["2017"], bf_syst["2017"] = npzfile["bf_2017"], npzfile["bf_stat_2017"], npzfile["bf_syst_2017"]
bf["2018"], bf_stat["2018"], bf_syst["2018"] = npzfile["bf_2018"], npzfile["bf_stat_2018"], npzfile["bf_syst_2018"]

bf_meas = np.array([
                    [   bf["2012"]["4m"],       bf["2012"]["2m2e"],     bf["2012"]["4e"],   ],
                    [   bf["2016"]["4m"],       bf["2016"]["2m2e"],     bf["2016"]["4e"],   ],
                    [   bf["2017"]["4m"],       bf["2017"]["2m2e"],     bf["2017"]["4e"],   ],
                    [   bf["2018"]["4m"],       bf["2018"]["2m2e"],     bf["2018"]["4e"],   ],
                    ])
bf_meas = np.squeeze(bf_meas)
bf_meas = bf_meas.flatten()

unc_stat = np.array([
                    [   bf_stat["2012"]["4m"], bf_stat["2012"]["2m2e"], bf_stat["2012"]["4e"],  ],
                    [   bf_stat["2016"]["4m"], bf_stat["2016"]["2m2e"], bf_stat["2016"]["4e"],  ],
                    [   bf_stat["2017"]["4m"], bf_stat["2017"]["2m2e"], bf_stat["2017"]["4e"],  ],
                    [   bf_stat["2018"]["4m"], bf_stat["2018"]["2m2e"], bf_stat["2018"]["4e"],  ],
                    ])
unc_stat = np.squeeze(unc_stat)
unc_stat = unc_stat.flatten()

unc_syst = np.array([
                    [   bf_syst["2012"]["4m"], bf_syst["2012"]["2m2e"], bf_syst["2012"]["4e"],  ],
                    [   bf_syst["2016"]["4m"], bf_syst["2016"]["2m2e"], bf_syst["2016"]["4e"],  ],
                    [   bf_syst["2017"]["4m"], bf_syst["2017"]["2m2e"], bf_syst["2017"]["4e"],  ],
                    [   bf_syst["2018"]["4m"], bf_syst["2018"]["2m2e"], bf_syst["2018"]["4e"],  ],
                    ])
unc_syst = np.squeeze(unc_syst)
unc_syst = unc_syst.flatten()

bf_unc = np.sqrt(unc_stat ** 2 + unc_syst ** 2)


# Predictions
bf_pred = np.array([    1.195,  2.310,  1.195   ])
pred = np.tile(bf_pred, 4)

# Scaling
data = bf_meas / pred
unc = bf_unc / pred


print("Observed signal strengths")
print(data)
print("")

print("Individual signal strengths")
print(alpha_total_12)
print("")



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


# Twelve parameters
fit_12 = alpha_total_12
res_12 = data - fit_12
pulls_12 = res_12 / unc

print("Twelve parameters")
print(pulls_12)
print("")



#
#   Comparison of fit parameterizations
#

# Single parameter

fit_1 = np.squeeze(alpha_total_1)
res_1 = alpha_total_12 - fit_1
pulls_fit_1 = res_1 / sigma_total_12

print("FIT: One parameter")
print(pulls_fit_1)
print("")


# Three parameters
fit_3 = np.tile(alpha_total_3, 4)
res_3 = alpha_total_12 - fit_3
pulls_fit_3 = res_3 / sigma_total_12

print("FIT: Three parameters")
print(pulls_fit_3)
print("")


# Four parameters
fit_4 = np.repeat(alpha_total_4, 3)
res_4 = alpha_total_12 - fit_4
pulls_fit_4 = res_4 / sigma_total_12

print("FIT: Four parameters")
print(pulls_fit_4)
print("")



##
##  HISTOGRAMS
##

gStyle.SetOptStat(1100)
gStyle.SetStatW(0.4)
gStyle.SetStatH(0.2)
gStyle.SetStatY(0.9)

outname = "pulls.root"
outfile = TFile(outname, "RECREATE")
outfile.cd()

h1 = TH1D("hPulls_SM", "SM parameterization", 12, -3, 3)
h3 = TH1D("hPulls_BSM", "BSM parameterization", 12, -3, 3)
h4 = TH1D("hPulls_yearly", "Yearly parameterization", 12, -3, 3)
h12 = TH1D("hPulls_indiv", "Individual parameterization", 12, -3, 3)

for i in range(len(data)):
    h1.Fill(np.squeeze(pulls_1[i]))
    h3.Fill(np.squeeze(pulls_3[i]))
    h4.Fill(np.squeeze(pulls_4[i]))
    h12.Fill(np.squeeze(pulls_12[i]))

for h in [h1, h3, h4, h12]:
    h.SetFillColor(9)
    h.SetLineColor(9)
#   h.SetStats(1100)
    h.SetXTitle("Pull")

c1 = TCanvas("pulls_SM", "", 500, 500)
c3 = TCanvas("pulls_BSM", "", 500, 500)
c4 = TCanvas("pulls_yearly", "", 500, 500)
c12 = TCanvas("pulls_indiv", "", 500, 500)

c1.cd()
h1.Draw("HIST")
c3.cd()
h3.Draw("HIST")
c4.cd()
h4.Draw("HIST")
c12.cd()
h12.Draw("HIST")


hf1 = TH1D("hPulls_Fit_SM", "SM parameterization", 12, -3, 3)
hf3 = TH1D("hPulls_Fit_BSM", "BSM parameterization", 12, -3, 3)
hf4 = TH1D("hPulls_Fit_yearly", "Yearly parameterization", 12, -3, 3)

for i in range(len(data)):
    hf1.Fill(np.squeeze(pulls_fit_1[i]))
    hf3.Fill(np.squeeze(pulls_fit_3[i]))
    hf4.Fill(np.squeeze(pulls_fit_4[i]))

for h in [hf1, hf3, hf4]:
    h.SetFillColor(9)
    h.SetLineColor(9)
#   h.SetStats(1100)
    h.SetXTitle("Pull")

cf1 = TCanvas("pulls_fit_SM", "", 500, 500)
cf3 = TCanvas("pulls_fit_BSM", "", 500, 500)
cf4 = TCanvas("pulls_fit_yearly", "", 500, 500)

cf1.cd()
hf1.Draw("HIST")
cf3.cd()
hf3.Draw("HIST")
cf4.cd()
hf4.Draw("HIST")

for h in [h1, h3, h4, h12, hf1, hf3, hf4]:
    h.Write()
for c in [c1, c3, c4, c12, cf1, cf3, cf4]:
    c.SaveAs(".pdf")
    c.Write()

outfile.Close()

print("Wrote histograms to", outname)
