from __future__ import print_function
from __future__ import division

import numpy as np
from numpy.linalg import multi_dot, pinv

from iminuit import Minuit



##
##  OPTIONS
##

N = 12   # number of parameters


##
##  SAMPLE INFO
##

selection   = [ "4m",   "2m2e", "4e"    ]
period      = [ "2012", "2016", "2017", "2018"  ]

M, P = len(selection), len(period)

np.set_printoptions(precision=5, suppress=True, linewidth=100)
#np.set_printoptions(precision=2, linewidth=100)



##
##  DATA
##

bf, bf_stat, diff, bg_unc, npt, sig = {}, {}, {}, {}, {}, {}

infile = "bf_measurements.npz"
npzfile = np.load(infile)

bf["2012"], bf_stat["2012"] = npzfile["bf_2012"], npzfile["bf_stat_2012"]
diff["2012"], bg_unc["2012"] = npzfile["diff_2012"], npzfile["bg_unc_2012"]
npt["2012"], sig["2012"] = npzfile["npt_2012"], npzfile["sig_2012"]

bf["2016"], bf_stat["2016"] = npzfile["bf_2016"], npzfile["bf_stat_2016"]
diff["2016"], bg_unc["2016"] = npzfile["diff_2016"], npzfile["bg_unc_2016"]
npt["2016"], sig["2016"] = npzfile["npt_2016"], npzfile["sig_2016"]

bf["2017"], bf_stat["2017"] = npzfile["bf_2017"], npzfile["bf_stat_2017"]
diff["2017"], bg_unc["2017"] = npzfile["diff_2017"], npzfile["bg_unc_2017"]
npt["2017"], sig["2017"] = npzfile["npt_2017"], npzfile["sig_2017"]

bf["2018"], bf_stat["2018"] = npzfile["bf_2018"], npzfile["bf_stat_2018"]
diff["2018"], bg_unc["2018"] = npzfile["diff_2018"], npzfile["bg_unc_2018"]
npt["2018"], sig["2018"] = npzfile["npt_2018"], npzfile["sig_2018"]


mu_id_unc, el_id_unc, el_reco_unc, ecal_unc = {}, {}, {}, {}
qcd_unc, pdf_unc, pu_unc = {}, {}, {}

from Cuts2012 import *
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup
from Cuts2016 import *
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup
from Cuts2017 import *
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup
from Cuts2018 import *
mu_id_unc[YEAR_STR], el_id_unc[YEAR_STR], el_reco_unc[YEAR_STR] = mu_id, el_id, el_reco
ecal_unc[YEAR_STR], qcd_unc[YEAR_STR], pdf_unc[YEAR_STR], pu_unc[YEAR_STR] = ecal, qcd, pdf, pileup


# Trigger efficiency uncertainties
mumu_trig_unc, ee_trig_unc = {}, {}
for year in period:
    infile = "trigger_eff_" + year + ".npz"
    npzfile = np.load(infile)

    ll_sf = npzfile["sf"]
    mumu_trig_unc[year] = abs(1 - np.squeeze(ll_sf["mumu"])) / 2
    ee_trig_unc[year] = abs(1 - np.squeeze(ll_sf["ee"])) / 2



##
##  FUNCTIONS
##

# Repeat yearly uncertainty
def get_unc_year(unc):
    return np.tile(unc, (1, M))

# Repeat yearly rho
def get_rho_year(rho):
    rho = np.repeat(rho, M, axis=1)
    rho = np.repeat(rho, M, axis=0)
    return rho


# Uncorrelated (diagonal) covariance matrix
def get_cov_uncorr(unc):
    unc = unc.flatten()                     # puts rows next to each other
    return np.diag(unc ** 2)


# Fully correlated covariance matrix
def get_cov_corr(unc):
    unc = unc.flatten()                     # puts rows next to each other
    return np.tensordot(unc, unc.T, axes=0)


# Target function to minimize
def least_squares(vec, cov):
    return multi_dot([vec, pinv(cov), vec.T])
#   return np.matmul(vec, np.matmul(pinv(cov), vec.T))



def get_vec(meas, pred, params):
    pred = np.tile(pred, P)       # flat predictions for every row (year)
    meas = meas.flatten()

    if len(params) == 1:                # one alpha (SM result)
        alpha = np.squeeze(params)
    elif len(params) == M:              # three alphas (BSM result)
        alpha = np.tile(params, P)
    elif len(params) == P:              # four alphas (CMS result)
        alpha = np.repeat(params, M)
    else:
        alpha = params

#   print("Parameters", "\n", alpha, "\n")

    return meas - alpha * pred



##
##  CENTRAL VALUES
##

print("\n")

# SM predictions
                    #   4m      2m2e    4e
bf_pred = np.array([    1.195,  2.310,  1.195   ])

print("Predictions", "\n",          bf_pred.flatten(), "\n")


# Measured central values
bf_meas = np.array([
                    [   bf["2012"]["4m"],       bf["2012"]["2m2e"],     bf["2012"]["4e"],   ],
                    [   bf["2016"]["4m"],       bf["2016"]["2m2e"],     bf["2016"]["4e"],   ],
                    [   bf["2017"]["4m"],       bf["2017"]["2m2e"],     bf["2017"]["4e"],   ],
                    [   bf["2018"]["4m"],       bf["2018"]["2m2e"],     bf["2018"]["4e"],   ],
                    ])
bf_meas = np.squeeze(bf_meas)

print("Measurements", "\n",         bf_meas.flatten(), "\n")



##
##  UNCERTAINTIES
##


# Statistical
unc_stat = np.array([
                    [   bf_stat["2012"]["4m"], bf_stat["2012"]["2m2e"], bf_stat["2012"]["4e"],  ],
                    [   bf_stat["2016"]["4m"], bf_stat["2016"]["2m2e"], bf_stat["2016"]["4e"],  ],
                    [   bf_stat["2017"]["4m"], bf_stat["2017"]["2m2e"], bf_stat["2017"]["4e"],  ],
                    [   bf_stat["2018"]["4m"], bf_stat["2018"]["2m2e"], bf_stat["2018"]["4e"],  ],
                    ])
bf_stat = np.squeeze(bf_stat)
print("Statistical uncertainties", "\n",  unc_stat.flatten(), "\n")

cov_stat = get_cov_uncorr(unc_stat)
print("Statistical covariance", "\n",     cov_stat, "\n")


# Pileup
unc_pu = bf_meas * get_unc_year(np.array([
                    [   pu_unc["2012"]  ],
                    [   pu_unc["2016"]  ],
                    [   pu_unc["2017"]  ],
                    [   pu_unc["2018"]  ],
                    ]))
print("Pileup uncertainties", "\n", unc_pu.flatten(), "\n")

cov_pu = get_cov_corr(unc_pu)
print("Pileup covariance", "\n",    cov_pu, "\n")


# PDFs
unc_pdf = bf_meas * get_unc_year(np.array([
                    [   pdf_unc["2012"] ],
                    [   pdf_unc["2016"] ],
                    [   pdf_unc["2017"] ],
                    [   pdf_unc["2018"] ],
                    ]))
print("PDF uncertainties", "\n",    unc_pdf.flatten(), "\n")

cov_pdf = get_cov_corr(unc_pdf)
print("PDF covariance", "\n",       cov_pdf * 100)
print("(* 1e-2)", "\n")


# QCD scales
unc_qcd = bf_meas * get_unc_year(np.array([
                    [   qcd_unc["2012"] ],
                    [   qcd_unc["2016"] ],
                    [   qcd_unc["2017"] ],
                    [   qcd_unc["2018"] ],
                    ]))
print("QCD uncertainties", "\n",    unc_qcd.flatten(), "\n")

cov_qcd = get_cov_corr(unc_qcd)
print("QCD covariance", "\n",       cov_qcd, "\n")


# Prefiring weight
#ecal_unc *= 1.0016
unc_ecal = bf_meas * np.array([
                    [ ecal_unc["2012"]["4m"], ecal_unc["2012"]["2m2e"], ecal_unc["2012"]["4e"], ],
                    [ ecal_unc["2016"]["4m"], ecal_unc["2016"]["2m2e"], ecal_unc["2016"]["4e"], ],
                    [ ecal_unc["2017"]["4m"], ecal_unc["2017"]["2m2e"], ecal_unc["2017"]["4e"], ],
                    [ ecal_unc["2018"]["4m"], ecal_unc["2018"]["2m2e"], ecal_unc["2018"]["4e"], ],
                    ])
print("Prefiring uncertainties", "\n",  unc_ecal.flatten(), "\n")

cov_ecal = get_cov_corr(unc_ecal)
print("Prefiring covariance", "\n",     cov_ecal * 100)
print("(* 1e-2)", "\n")


# Dimuon trigger efficiency
#unc_mutr = np.zeros_like(bf_meas)
unc_mutr = bf_meas * np.array([
                                [   mumu_trig_unc["2012"],  0.5 * mumu_trig_unc["2012"],    0,  ],
                                [   mumu_trig_unc["2016"],  0.5 * mumu_trig_unc["2016"],    0,  ],
                                [   mumu_trig_unc["2017"],  0.5 * mumu_trig_unc["2017"],    0,  ],
                                [   mumu_trig_unc["2018"],  0.5 * mumu_trig_unc["2018"],    0,  ],
                                ])
print("Dimuon trigger uncertainties", "\n",  unc_mutr.flatten(), "\n")

cov_mutr = get_cov_uncorr(unc_mutr)
print("Dimuon trigger covariance", "\n",     cov_mutr, "\n")


# Dielectron trigger efficiency
#unc_eltr = np.zeros_like(bf_meas)
unc_eltr = bf_meas * np.array([
                                [   0,  0.5 * ee_trig_unc["2012"],      ee_trig_unc["2012"],    ],
                                [   0,  0.5 * ee_trig_unc["2016"],      ee_trig_unc["2016"],    ],
                                [   0,  0.5 * ee_trig_unc["2017"],      ee_trig_unc["2017"],    ],
                                [   0,  0.5 * ee_trig_unc["2018"],      ee_trig_unc["2018"],    ],
                                ])
print("Dielectron trigger uncertainties", "\n",  unc_eltr.flatten(), "\n")

cov_eltr = get_cov_uncorr(unc_eltr)
print("Dielectron trigger covariance", "\n",     cov_eltr, "\n")


# Muon ID
unc_muid = bf_meas * np.array([
                [ mu_id_unc["2012"]["4m"], mu_id_unc["2012"]["2m2e"], mu_id_unc["2012"]["4e"], ],
                [ mu_id_unc["2016"]["4m"], mu_id_unc["2016"]["2m2e"], mu_id_unc["2016"]["4e"], ],
                [ mu_id_unc["2017"]["4m"], mu_id_unc["2017"]["2m2e"], mu_id_unc["2017"]["4e"], ],
                [ mu_id_unc["2018"]["4m"], mu_id_unc["2018"]["2m2e"], mu_id_unc["2018"]["4e"], ],
                ])
print("Muon ID uncertainties", "\n",     unc_muid.flatten(), "\n")

mu_13, mu_8_13 = 1, 0

rho_muid = np.array([#  2012        2016        2017        2018
                    [   1,          mu_8_13,    mu_8_13,    mu_8_13,    ],  # 2012
                    [   mu_8_13,    1,          mu_13,      mu_13,      ],  # 2016
                    [   mu_8_13,    mu_13,      1,          mu_13,      ],  # 2017
                    [   mu_8_13,    mu_13,      mu_13,      1,          ],  # 2017
                    ])
rho_muid = get_rho_year(rho_muid)
print("Muon ID correlations", "\n",    rho_muid, "\n")

cov_muid = rho_muid * get_cov_corr(unc_muid) 
print("Muon ID cov.", "\n",     cov_muid, "\n")


# Electron ID
unc_elid = bf_meas * np.array([
                [ el_id_unc["2012"]["4m"], el_id_unc["2012"]["2m2e"], el_id_unc["2012"]["4e"], ],
                [ el_id_unc["2016"]["4m"], el_id_unc["2016"]["2m2e"], el_id_unc["2016"]["4e"], ],
                [ el_id_unc["2017"]["4m"], el_id_unc["2017"]["2m2e"], el_id_unc["2017"]["4e"], ],
                [ el_id_unc["2018"]["4m"], el_id_unc["2018"]["2m2e"], el_id_unc["2018"]["4e"], ],
                ])
print("Electron ID uncertainties", "\n", unc_elid.flatten(), "\n")

el_13, el_8_13 = 1, 0.5

rho_elid = np.array([#  2012        2016        2017        2018
                    [   1,          el_8_13,    el_8_13,    el_8_13,    ],  # 2012
                    [   el_8_13,    1,          el_13,      el_13,      ],  # 2016
                    [   el_8_13,    el_13,      1,          el_13,      ],  # 2017
                    [   el_8_13,    el_13,      el_13,      1,          ],  # 2017
                    ])
rho_elid = get_rho_year(rho_elid)
print("Electron ID correlations", "\n",    rho_elid, "\n")

cov_elid = rho_elid * get_cov_corr(unc_elid) 
print("Electron ID cov.", "\n",     cov_elid, "\n")


# Electron reco
unc_reco = bf_meas * np.array([
           [ el_reco_unc["2012"]["4m"], el_reco_unc["2012"]["2m2e"], el_reco_unc["2012"]["4e"], ],
           [ el_reco_unc["2016"]["4m"], el_reco_unc["2016"]["2m2e"], el_reco_unc["2016"]["4e"], ],
           [ el_reco_unc["2017"]["4m"], el_reco_unc["2017"]["2m2e"], el_reco_unc["2017"]["4e"], ],
           [ el_reco_unc["2018"]["4m"], el_reco_unc["2018"]["2m2e"], el_reco_unc["2018"]["4e"], ],
           ])
print("Electron reco uncertainties", "\n", unc_reco.flatten(), "\n")

el_13, el_8_13 = 1, 0

rho_reco = np.array([#  2012        2016        2017        2018
                    [   1,          el_8_13,    el_8_13,    el_8_13,    ],  # 2012
                    [   el_8_13,    1,          el_13,      el_13,      ],  # 2016
                    [   el_8_13,    el_13,      1,          el_13,      ],  # 2017
                    [   el_8_13,    el_13,      el_13,      1,          ],  # 2017
                    ])
rho_reco = get_rho_year(rho_reco)
print("Electron reco correlations", "\n",  rho_reco, "\n")

cov_reco = rho_reco * get_cov_corr(unc_reco) 
print("Electron reco covariance", "\n",  cov_reco, "\n")


# Background statistical
diff_evts = np.array([
                    [   diff["2012"]["4m"],     diff["2012"]["2m2e"],   diff["2012"]["4e"]  ],
                    [   diff["2016"]["4m"],     diff["2016"]["2m2e"],   diff["2016"]["4e"]  ],
                    [   diff["2017"]["4m"],     diff["2017"]["2m2e"],   diff["2017"]["4e"]  ],
                    [   diff["2018"]["4m"],     diff["2018"]["2m2e"],   diff["2018"]["4e"]  ],
                    ])
diff_evts = np.squeeze(diff_evts)

unc_bsta = np.array([
                    [   bg_unc["2012"]["4m"],   bg_unc["2012"]["2m2e"], bg_unc["2012"]["4e"]    ],
                    [   bg_unc["2016"]["4m"],   bg_unc["2016"]["2m2e"], bg_unc["2016"]["4e"]    ],
                    [   bg_unc["2017"]["4m"],   bg_unc["2017"]["2m2e"], bg_unc["2017"]["4e"]    ],
                    [   bg_unc["2018"]["4m"],   bg_unc["2018"]["2m2e"], bg_unc["2018"]["4e"]    ],
                    ])
unc_bsta = np.squeeze(unc_bsta)
unc_bsta = bf_meas * unc_bsta / diff_evts
print("Background statistical uncertainties", "\n",  unc_bsta.flatten(), "\n")

cov_bsta = get_cov_uncorr(unc_bsta)
print("Background statistical covariance", "\n", cov_bsta, "\n")


# Background systematic
unc_bsys = np.array([
                    [   npt["2012"]["4m"],      npt["2012"]["2m2e"],    npt["2012"]["4e"]   ],
                    [   npt["2016"]["4m"],      npt["2016"]["2m2e"],    npt["2016"]["4e"]   ],
                    [   npt["2017"]["4m"],      npt["2017"]["2m2e"],    npt["2017"]["4e"]   ],
                    [   npt["2018"]["4m"],      npt["2018"]["2m2e"],    npt["2018"]["4e"]   ],
                    ])
unc_bsys = np.squeeze(unc_bsys)
unc_bsys = bf_meas * DELTA_LAMBDA * unc_bsys / diff_evts
print("Background systematic uncertainties", "\n",  unc_bsys.flatten(), "\n")

cov_bsys = get_cov_corr(unc_bsys)
print("Background systematic covariance", "\n", cov_bsys, "\n")


# Total covariance
cov_syst = cov_pu + cov_pdf + cov_qcd + cov_ecal + cov_mutr + cov_eltr
cov_syst += cov_muid + cov_elid + cov_reco + cov_bsta + cov_bsys
cov_total = cov_stat + cov_syst

print("Total covariance", "\n", cov_total, "\n")



##
##  PARAMETERS
##

alpha_0 = np.ones(N)

print("Initial parameters", "\n", alpha_0, "\n")
print("\n")



##
##  STATISTICAL ONLY RESULT
##

def target_func_stat(alpha):
    vec = get_vec(bf_meas, bf_pred, alpha)
    cov = cov_stat
    return least_squares(vec, cov)

print("STATISTICAL COVARIANCE", "\n")
minuit = Minuit.from_array_func(target_func_stat, alpha_0, error=0, errordef=1)
minuit.migrad()

alpha_stat = minuit.np_values()
sigma_stat = minuit.np_errors()
delta_stat = sigma_stat / alpha_stat
chi_sq_stat = minuit.fval

print("\n\n")



##
##  TOTAL RESULT
##

def target_func_total(alpha):
    vec = get_vec(bf_meas, bf_pred, alpha)
    cov = cov_total
    return least_squares(vec, cov)

print("TOTAL COVARIANCE", "\n")
minuit = Minuit.from_array_func(target_func_total, alpha_0, error=0, errordef=1)
minuit.migrad()

alpha_total = minuit.np_values()
sigma_total = minuit.np_errors()
delta_total = sigma_total / alpha_total
chi_sq_total = minuit.fval

delta_syst = np.sqrt(delta_total ** 2 - delta_stat ** 2)

if len(alpha_total) == M: 
    alpha_f = np.tile(alpha_total, P)
elif len(alpha_total) == P: 
    alpha_f = np.repeat(alpha_total, M)
else:
    alpha_f = alpha_total

bf_comb = alpha_f * np.tile(bf_pred, P)

print("\n\n")



##
##  PRINT
##

print("Statistical-only parameters", "\n", alpha_stat, "\n")
print("Total covariance parameters", "\n", alpha_total, "\n")

print("Fractional statistical errors", "\n", delta_stat, "\n")
print("Fractional systematic errors", "\n", delta_syst, "\n")
print("Fractional total errors", "\n", delta_total, "\n")

print("Statistical-only chi-squared", "\n", chi_sq_stat, "\n")
print("Total chi-squared", "\n", chi_sq_total, "\n")

print("Results from total parameters", "\n", bf_comb, "\n")

print("\n")



##
##  LATEX
##





##
##  SAVE
##

sigma_stat = delta_stat * alpha_total
sigma_syst = delta_syst * alpha_total

outfile = "combination_" + str(len(alpha_0)) + ".npz"
np.savez(outfile, bf_pred=bf_pred, alpha_total=alpha_total, alpha_stat=alpha_stat,
        sigma_stat=sigma_stat, sigma_syst=sigma_syst, sigma_total=sigma_total,
        delta_stat=delta_stat, delta_syst=delta_syst, delta_total=delta_total,
        chi_sq_total=chi_sq_total, chi_sq_stat=chi_sq_stat)

print("Wrote arrays to", outfile, "\n")
