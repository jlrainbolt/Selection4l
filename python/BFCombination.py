from __future__ import print_function
from __future__ import division

import numpy as np
from numpy.linalg import multi_dot, pinv

from iminuit import Minuit

from CutsComb import *



##
##  SAMPLE INFO
##

selection   = [ "4m",   "2m2e", "4e"    ]
period      = [ "2012", "2016", "2017", "2018"  ]

M, P = len(selection), len(period)

np.set_printoptions(precision=5, suppress=True, linewidth=100)
#np.set_printoptions(precision=2, linewidth=100)


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

    print("Parameters", "\n", alpha, "\n")

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
bf_meas = np.array([#   4m      2m2e    4e
                    [   1.348,  1.797,  1.000   ],  # 2012
                    [   1.372,  2.169,  1.277   ],  # 2016
                    [   1.368,  2.401,  1.109   ],  # 2017
                    [   1.384,  2.458,  1.240   ]   # 2018
                    ])

print("Measurements", "\n",         bf_meas.flatten(), "\n")



##
##  UNCERTAINTIES
##


# Statistical
unc_stat = np.array([#  4m      2m2e    4e
                    [   0.166,  0.274,  0.246   ],  # 2012
                    [   0.094,  0.196,  0.187   ],  # 2016
                    [   0.086,  0.185,  0.157   ],  # 2017
                    [   0.072,  0.161,  0.144   ],  # 2018
                    ])
print("Statistical uncertainties", "\n",  unc_stat.flatten(), "\n")

cov_stat = get_cov_uncorr(unc_stat)
print("Statistical covariance", "\n",     cov_stat, "\n")


# Pileup
unc_pu = np.array([
                    [   pileup[2012]],  # 2012
                    [   pileup[2016]],  # 2016
                    [   pileup[2017]],  # 2017
                    [   pileup[2018]],  # 2018
                    ])
unc_pu = bf_meas * get_unc_year(unc_pu)
print("Pileup uncertainties", "\n", unc_pu.flatten(), "\n")

cov_pu = get_cov_corr(unc_pu)
print("Pileup covariance", "\n",    cov_pu, "\n")


# PDFs
unc_pdf = np.array([
                    [   pdf[2012]   ],  # 2012
                    [   pdf[2016]   ],  # 2016
                    [   pdf[2017]   ],  # 2017
                    [   pdf[2018]   ],  # 2018
                    ])
unc_pdf = bf_meas * get_unc_year(unc_pdf)
print("PDF uncertainties", "\n",    unc_pdf.flatten(), "\n")

cov_pdf = get_cov_corr(unc_pdf)
print("PDF covariance", "\n",       cov_pdf * 100)
print("(* 1e-2)", "\n")


# QCD scales
unc_qcd = np.array([
                    [   qcd[2012]   ],  # 2012
                    [   qcd[2016]   ],  # 2016
                    [   qcd[2017]   ],  # 2017
                    [   qcd[2018]   ],  # 2018
                    ])
unc_qcd = bf_meas * get_unc_year(unc_qcd)
print("QCD uncertainties", "\n",    unc_qcd.flatten(), "\n")

cov_qcd = get_cov_corr(unc_qcd)
print("QCD covariance", "\n",       cov_qcd, "\n")


# Prefiring weight
unc_ecal = np.array([#  4m                  2m2e                4e
                    [   0,                  0,                  0,              ],  # 2012
                    [   ecal_2016["4m"],    ecal_2016["2m2e"],  ecal_2016["4e"],],  # 2016
                    [   ecal_2017["4m"],    ecal_2017["2m2e"],  ecal_2017["4e"],],  # 2017
                    [   0,                  0,                  0,              ],  # 2018
                    ])
print("Prefiring uncertainties", "\n",  unc_ecal.flatten(), "\n")

cov_ecal = get_cov_corr(unc_ecal)
#cov_ecal = get_cov_uncorr(unc_ecal)
print("Prefiring covariance", "\n",     cov_ecal * 100)
print("(* 1e-2)", "\n")


# Trigger efficiency (FIXME)

unc_trig = np.zeros_like(bf_meas)
print("Trigger uncertainties", "\n",  unc_trig.flatten(), "\n")

cov_trig = get_cov_corr(unc_trig)
print("Trigger covariance", "\n",     cov_trig, "\n")


# Muon ID (FIXME)
                    #   4m                  2m2e                4e
unc_muid = np.array([   mu_id["4m"],        mu_id["2m2e"],      mu_id["4e"],    ])
unc_muid = bf_meas * np.tile(unc_muid, (P, 1))
print("Muon ID uncertainties", "\n",     unc_muid.flatten(), "\n")

mu_13, mu_8_13 = 0, 0

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


# Electron ID (FIXME)
                    #   4m                  2m2e                4e
unc_elid = np.array([   el_id["4m"],        el_id["2m2e"],      el_id["4e"],    ])
unc_elid = bf_meas * np.tile(unc_elid, (P, 1))
print("Electron ID uncertainties", "\n", unc_elid.flatten(), "\n")

el_13, el_8_13 = 0, 0

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


# Electron reco (FIXME)
                    #   4m                  2m2e                4e
unc_reco = np.array([   el_reco["4m"],      el_reco["2m2e"],    el_reco["4e"],    ])
unc_reco = bf_meas * np.tile(unc_reco, (P, 1))
print("Electron reco uncertainties", "\n", unc_reco.flatten(), "\n")

#el_13, el_8_13 = 1, 1

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
unc_bsta = np.array([#  4m                      2m2e                    4e
                    [   npt_unc_2012["4m"],     npt_unc_2012["2m2e"],   npt_unc_2012["4e"]], # 2012
                    [   npt_unc_2016["4m"],     npt_unc_2016["2m2e"],   npt_unc_2016["4e"]], # 2016
                    [   npt_unc_2017["4m"],     npt_unc_2017["2m2e"],   npt_unc_2017["4e"]], # 2017
                    [   npt_unc_2018["4m"],     npt_unc_2018["2m2e"],   npt_unc_2018["4e"]], # 2018
                    ])
unc_bsta = bf_meas * unc_bsta
print("Background statistical uncertainties", "\n",  unc_bsta.flatten(), "\n")

cov_bsta = get_cov_uncorr(unc_bsta)
print("Background statistical covariance", "\n", cov_bsta, "\n")


# Background systematic (FIXME)
unc_bsys = np.full_like(bf_meas, 0.005)
unc_bsys = bf_meas * unc_bsys
print("Background systematic uncertainties", "\n",  unc_bsys.flatten(), "\n")

cov_bsys = get_cov_corr(unc_bsys)
print("Background systematic covariance", "\n", cov_bsys, "\n")


# Total covariance
cov_syst = cov_pu + cov_pdf + cov_qcd + cov_ecal + cov_trig
cov_syst += cov_muid + cov_elid + cov_reco + cov_bsta + cov_bsys
cov_total = cov_stat + cov_syst

print("Total covariance", "\n", cov_total, "\n")


##
##  PARAMETERS
##

alpha_0 = np.ones(1)

print("Initial parameters", "\n", alpha_0, "\n")
print("\n")



##
##  SOLVE
##

def target_func(alpha):
    vec = get_vec(bf_meas, bf_pred, alpha)
#   cov = cov_stat
    cov = cov_total
    return least_squares(vec, cov)

minuit = Minuit.from_array_func(target_func, alpha_0, error=0, errordef=1)

minuit.migrad()



##
##  RESULTS
##

print("\n")

alpha_f = minuit.np_values()
print("Parameters", "\n", alpha_f, "\n")

sigma_f = minuit.np_errors()
print("Hesse errors", "\n", sigma_f, "\n")

delta_f = sigma_f / alpha_f
print("Fractional Hesse errors", "\n", delta_f, "\n")

chi_sq = target_func(alpha_f) / (P * M)
print("Reduced chi-squared", "\n", chi_sq, "\n")

if len(alpha_f) == P:
    alpha_f = np.tile(alpha_f.T, (M, 1)).T
    bf_pred = np.tile(bf_pred, (P, 1))

bf_comb = alpha_f * bf_pred
print("Combined result", "\n", bf_comb, "\n")


# Find systematic uncertainty
print("\n")

alpha_stat = np.full_like(alpha_f, 1.0654246397258216)
print("Stat. only parameters", "\n", alpha_stat, "\n")

sigma_stat = np.full_like(sigma_f, 0.026306929883156572)
print("Stat. only Hesse errors", "\n", sigma_stat, "\n")

delta_stat = sigma_stat / alpha_stat
print("Stat. only fractional errors", "\n", delta_stat, "\n")

delta_syst = np.sqrt(delta_f ** 2 - delta_stat ** 2)
print("Fractional systematic uncertainty", "\n", delta_syst, "\n")
