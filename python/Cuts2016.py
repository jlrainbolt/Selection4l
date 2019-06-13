from __future__ import division

from PlotUtils import *

##
##  BASICS
##

YEAR_STR    = "2016"
EOS_PATH    = "root://cmseos.fnal.gov//store/user/jrainbol"
HOME_PATH   = "/uscms/home/jrainbol/nobackup"



##
##  SYSTEMATICS
##

F_NR = 0.04
BF_LL = 0.033658
CAP_K = 5.085509e-10
GAMMA_Z = 2.4952e6

mu_id = {   "4l":0.0171,    "4m":0.0213,    "2m2e":0.0095,      "4e":0          }
el_id = {   "4l":0.0062,    "4m":0,         "2m2e":0.0138,      "4e":0.0265     }
el_reco = { "4l":0.0027,    "4m":0,         "2m2e":0.0093,      "4e":0.0154     }
mu_pt = {   "4l":0.000042,  "4m":0.000071,  "2m2e":0.000042,    "4e":0          }
el_pt = {   "4l":0.001017,  "4m":0,         "2m2e":0.002930,    "4e":0.000169   }
ecal =  {   "4l":0.000799,  "4m":0,         "2m2e":0.000536,    "4e":0.001595   }
qcd = 0.008592
pdf = 0.001185
pileup = 0.008885



##
##  SAMPLE INFO
##

ELEC_TRIG_SF = 1
INT_LUMI, LUMI_UNC = 36.42, .025
SQRT_S  = 13
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

NGEN_ZZ_4L,         XSEC_ZZ_4L          = 6669988,      1.256     #1.212
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      = 80924255,     5765.4    #6225.4
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      = 999800,       0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     = 499262,       0.001034
NGEN_TTBAR,         XSEC_TTBAR          = 15173839,     831.76 / 2
NGEN_TT_2L2NU,      XSEC_TT_2L2NU       = 79140880,     87.31 / 2
NGEN_TTZ_2L2NU,     XSEC_TTZ_2L2NU      = 2654179,      0.2529
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       = 1999000,      12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        = 15879472,     5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        = 7387013,      4.42965
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        = 496436,       3.22
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       = 48655100,     0.564

N_MC = 12
N_DY = 10

MC_SUFF = [ "zz_4l",            "zjets_m-50",
            "ttbar",            "tt_2l2nu",         "ttz_2l2nu",        "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
            "ggH_zz_4l",        "vbfH_zz_4l",
            ]

NGEN_   = [ NGEN_ZZ_4L,         NGEN_ZJETS_M50,
            NGEN_TTBAR,         NGEN_TT_2L2NU,      NGEN_TTZ_2L2NU,     NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
            NGEN_GGH_ZZ_4L,     NGEN_VBFH_ZZ_4L,
            ]

XSEC_   = [ XSEC_ZZ_4L,         XSEC_ZJETS_M50,
            XSEC_TTBAR,         XSEC_TT_2L2NU,      XSEC_TTZ_2L2NU,     XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
            XSEC_GGH_ZZ_4L,     XSEC_VBFH_ZZ_4L,
            ]

COLOR_  = [ lLightBlue,         lYellow,
            lGreen,             lGreen,             lGreen,             lOrange,
            lOrange,            lOrange,            lOrange,            lOrange,
            lPurple,            lPurple,
            ]

MC_TEX_ = [ r"$\ZZtofL$",       r"$\ZtoLL$",
            r"$\ttbar$",        r"$\TTtoLLNuNu$",   r"$\TTZtoLLNuNu$",  r"$\WWtoLLNuNu$",
            r"$\WZtoLLQQ$",     r"$\WZtoLLLNu$",    r"$\ZZtoLLNuNu$",   r"$\ZZtoLLQQ$",
            r"$\ggF\HtoZZtofL$",r"$\VBF\HtoZZtofL$",
            ]

NGEN    = dict(zip(MC_SUFF, NGEN_))
XSEC    = dict(zip(MC_SUFF, XSEC_))
COLOR   = dict(zip(MC_SUFF, COLOR_))
MC_TEX  = dict(zip(MC_SUFF, MC_TEX_))

MC_SUFF_4L = MC_SUFF
MC_SUFF_2L = list(MC_SUFF)
MC_SUFF_2L[0], MC_SUFF_2L[1] = MC_SUFF[1], MC_SUFF[0]
