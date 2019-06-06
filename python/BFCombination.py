from __future__ import print_function
from __future__ import division

import numpy as np
from numpy.linalg import *
from numpy.matlib import *

from iminuit import Minuit



##
##  SAMPLE INFO
##

selection   = [ "4m",   "2m2e", "4e"]
period      = [ "2012", "2016", "2017", "2018"]

M, P = len(selection), len(period)


##
##  CENTRAL VALUES
##

# SM predictions
                    #   4m      2m2e    4e
bf_pred = np.array([    1.195,  2.310,  1.195   ])


# Measured central values
                    #   4m      2m2e    4e
bf_meas = np.array([
                    [   1.348,  1.797,  1.000   ],  # 2012
                    [   1.372,  2.169,  1.277   ],  # 2016
                    [   1.368,  2.401,  1.109   ],  # 2017
                    [   1.384,  2.458,  1.240   ]   # 2018
                    ])



##
##  UNCERTAINTIES
##


# Statistical
                    #   4m      2m2e    4e
bf_stat = np.array([
                    [   0.166,  0.274,  0.246   ],  # 2012
                    [   0.094,  0.196,  0.187   ],  # 2016
                    [   0.086,  0.185,  0.157   ],  # 2017
                    [   0.072,  0.161,  0.144   ],  # 2018
                    ])

# Systematic
                    #   4m      2m2e    4e
bf_syst = np.array([
                    [   0.034,  0.084,  0.069   ],  # 2012
                    [   0.037,  0.071,  0.058   ],  # 2016
                    [   0.037,  0.065,  0.048   ],  # 2017
                    [   0.037,  0.067,  0.044   ],  # 2018
                    ])


##
##  CORRELATIONS
##
'''
rho_stat = np.zeros((P, P))

                    #   2012    2016    2017    2018
rho_syst = np.array([
                    [   .5,     .5,     .5,    .5   ],  # 2012
                    [   .5,     .5,     .5,    .5   ],  # 2016
                    [   .5,     .5,     .5,    .5   ],  # 2017
                    [   .5,     .5,     .5,    .5   ],  # 2018
                    ])
'''
rho_stat = 0
rho_syst = 0.5



##
##  PARAMETERS
##

alpha_0 = np.ones_like(bf_meas)



##
##  FUNCTIONS
##

def least_squares(vec, cov):
    return multi_dot([vec, pinv(cov), vec.T])
#   return np.matmul(vec, np.matmul(pinv(cov), vec.T))


def get_cov_mat(unc, rho):
    unc = unc.flatten()                     # puts rows next to each other
    cov = np.tensordot(unc, unc.T, axes=0)

    for i, j in np.ndindex(cov.shape):
       if i != j:
           cov[i,j] *= rho

    return cov


def get_vec(meas, pred, alpha):
    pred = repmat(pred, 1, meas.shape[0])   # flat predictions for every row (year)
    meas = meas.flatten()
    alpha = alpha.flatten()
    return meas - np.multiply(alpha, pred)



##
##  TEST
##

# Only look at 4mu

bf_meas = bf_meas[:,0]
bf_stat = bf_stat[:,0]
bf_syst = bf_syst[:,0]
bf_pred = bf_pred[0]
alpha_0 = alpha_0[:,0]
error_0 = np.full_like(alpha_0, 0.1)

cov_stat = get_cov_mat(bf_stat, rho_stat)
cov_syst = get_cov_mat(bf_syst, rho_syst)
total_cov = cov_stat + cov_syst

print("\n")

print("Predictions", "\n",          bf_pred, "\n")
print("Measurements", "\n",         bf_meas, "\n")
print("Stat. uncertainties", "\n",  bf_stat, "\n")
print("Stat. covariance", "\n",     cov_stat, "\n")
print("Syst. uncertainties", "\n",  bf_syst, "\n")
print("Syst. covariance", "\n",     cov_syst, "\n")
print("Total covariance", "\n",     total_cov, "\n")
print("Initial parameters", "\n",   alpha_0, "\n")

print("\n")



##
##  SOLVE
##

def target_func(alpha):
    vec = get_vec(bf_meas, bf_pred, alpha)
    cov = total_cov
    return least_squares(vec, cov)

minuit = Minuit.from_array_func(target_func, alpha_0, error=error_0, errordef=1)

minuit.migrad()



##
##  INTERPRET
##

alpha_f = minuit.np_values()
#alpha_f = alpha_f * (len(alpha_f) / np.sum(alpha_f))

chi_sq = np.squeeze(target_func(alpha_f)) / len(alpha_0)

bf_comb = np.multiply(alpha_f, bf_pred)
bf_comb = np.sum(bf_comb) / len(alpha_0)

print("\n")

print("Normalized parameters", "\n",    alpha_f, "\n")
print("Reduced chi-squared", "\n",      chi_sq, "\n")
print("Combined result", "\n",          bf_comb, "\n")

