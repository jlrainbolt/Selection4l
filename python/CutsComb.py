from __future__ import division

from PlotUtils import *
import numpy as np

##
##  BASICS
##

EOS_PATH    = "root://cmseos.fnal.gov//store/user/jrainbol"
HOME_PATH   = "/uscms/home/jrainbol/nobackup"



##
##  SYSTEMATICS
##

F_NR = 0.04
BF_LL = 0.033658 * 2
GAMMA_Z = 2.4952e6
BF_4L = 5.14e-6

# Same for all channels (fractional)
pileup = {  2012:0.009704,  2016:0.008885,  2017:0.005765,  2018:0.005730   }
qcd =   {   2012:0.008592,  2016:0.008592,  2017:0.008592,  2018:0.008592   }
pdf =   {   2012:0.001185,  2016:0.001185,  2017:0.001185,  2018:0.001185   }

# Same for all years (FIXME, should not really be the same)
mu_id = {   "4l":0.0171,    "4m":0.0213,    "2m2e":0.0095,      "4e":0          }
el_id = {   "4l":0.0062,    "4m":0,         "2m2e":0.0138,      "4e":0.0265     }
el_reco = { "4l":0.0027,    "4m":0,         "2m2e":0.0093,      "4e":0.0154     }



##
##  FIT
##

alpha_stat_sm = np.array([1.06543])
sigma_stat_sm = np.array([0.02631])

alpha_stat_bsm = np.array([1.1497,  0.99137, 0.98628])
sigma_stat_bsm = np.array([0.03829, 0.04182, 0.07231])

alpha_stat_cms = np.array([0.91126, 1.0535,  1.07354, 1.10787])
sigma_stat_cms = np.array([0.08262, 0.05412, 0.04957, 0.04263])

alpha_stat_all = np.array([1.12803, 0.77792, 0.83682, 1.14812, 0.93896, 1.06862, 1.14477, 1.03939, 0.92803, 1.15816, 1.06407, 1.03766])
sigma_stat_all = np.array([0.13891, 0.11861, 0.20586, 0.07866, 0.08485, 0.15649, 0.07197, 0.08009, 0.13138, 0.06025, 0.0697,  0.1205 ])


##
##  2012
##

ecal_2012 = {   "4m":0,         "2m2e":0,           "4e":0      }
sig_2012 = {"4m":65.85,    "2m2e":43.14,  "4e":16.56}
npt_2012 = { "4l":3.14, "4m":0.05, "2m2e":2.86/2, "4e":0.33, "mumu":0, "ee":0, "2e2m":2.86/2, "ll":0}   
npt_unc_2012 = { "4l":2.00,      "4m":0.06,      "2m2e":1.73,    "4e":1.00}

MUON_TRIG_LUMI_2012, ELEC_TRIG_LUMI_2012, ELEC_TRIG_SF_2012 = 19.712, 19.712, 1

NGEN_ZZ_4L  = 1499064 + 1499093 + 1497445 + 823922 + 823911 + 824466
XSEC_ZZ_4L  = 3 * 0.1767 + 3 * 0.07691
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      = 30459503,       3531.9
NGEN_ZJETS_M10,     XSEC_ZJETS_M10      = 33648307,       11050
NGEN_TTBAR,         XSEC_TTBAR          = 12011428,       25.81
NGEN_TTZ_2L2NU,     XSEC_TTZ_2L2NU      = 210160,         0.2057
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       = 10000431,       57.25
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        = 3215990,        5.09
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        = 2017979,        1.086
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       = 1936727,        2.47
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        = 954911,         0.71

MC_SUFF_2012 = [ "zz_4l",            "zjets_m-50",
            "ttbar",            "ttz_2l2nu",        "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
            ]

NGEN_   = [ NGEN_ZZ_4L,         NGEN_ZJETS_M50,
            NGEN_TTBAR,         NGEN_TTZ_2L2NU,     NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
            ]

XSEC_   = [ XSEC_ZZ_4L,         XSEC_ZJETS_M50,
            XSEC_TTBAR,         XSEC_TTZ_2L2NU,     XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
            ]

NGEN_2012    = dict(zip(MC_SUFF_2012, NGEN_))
XSEC_2012    = dict(zip(MC_SUFF_2012, XSEC_))



##
##  2016
##

CAP_K_2016 = 5.15258169e-10 / 2

ecal_2016 = {   "4m":0,         "2m2e":0.000536,    "4e":0.001595   }
sig_2016 = {"4m":215.08,    "2m2e":122.79,  "4e":46.43}
npt_2016 = { "4l":14.19,     "4m":5.68,      "2m2e":7.45/2,    "4e":1.06, "mumu":0, "ee":0, "2e2m":7.45/2,  "ll":0}   
npt_unc_2016 = { "4l":4.00,      "4m":2.45,      "2m2e":2.83,    "4e":1.42}

MUON_TRIG_LUMI_2016, ELEC_TRIG_LUMI_2016, ELEC_TRIG_SF_2016 = 36.42, 36.42, 1

NGEN_ZZ_4L,         XSEC_ZZ_4L          = 6762740,        1.256     #1.212
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      = 80924255,       5765.4    #6225.4
NGEN_ZJETS_M10,     XSEC_ZJETS_M10      = 29374008,       18610
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      = 999800,         0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     = 499262,         0.001034
NGEN_TTBAR,         XSEC_TTBAR          = 15173839,       831.76 / 2
NGEN_TT_2L2NU,      XSEC_TT_2L2NU       = 65899840,       87.31 / 2
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       = 1832358,        12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        = 15879472,       5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        = 7387013,        4.42965
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        = 496436,         3.22
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       = 48623080,       0.564

MC_SUFF_2016 = [ "zz_4l",            "zjets_m-50",     # "zjets_m-10",
            "ttbar",            "tt_2l2nu",       # "ttz_2l2nu",        
            "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
#           "wwz_4l2nu",        "wzz_4l2nu",        "zzz_4l2nu",        "zzg_4l2nu",
            "ggH_zz_4l",        "vbfH_zz_4l",
            ]

NGEN_   = [ NGEN_ZZ_4L,         NGEN_ZJETS_M50,   # NGEN_ZJETS_M10,
            NGEN_TTBAR,         NGEN_TT_2L2NU,    # NGEN_TTZ_2L2NU,
            NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
#           NGEN_WWZ_4L2NU,     NGEN_WZZ_4L2NU,     NGEN_ZZZ_4L2NU,     NGEN_ZZG_4L2NU,
            NGEN_GGH_ZZ_4L,     NGEN_VBFH_ZZ_4L,
            ]

XSEC_   = [ XSEC_ZZ_4L,         XSEC_ZJETS_M50,   # XSEC_ZJETS_M10,
            XSEC_TTBAR,         XSEC_TT_2L2NU,    # XSEC_TTZ_2L2NU,
            XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
#           XSEC_WWZ_4L2NU,     XSEC_WZZ_4L2NU,     XSEC_ZZZ_4L2NU,     XSEC_ZZG_4L2NU,
            XSEC_GGH_ZZ_4L,     XSEC_VBFH_ZZ_4L,
            ]

NGEN_2016    = dict(zip(MC_SUFF_2016, NGEN_))
XSEC_2016    = dict(zip(MC_SUFF_2016, XSEC_))




##
##  2017
##

CAP_K_2017 = 4.64972283e-10 / 2

ecal_2017 = {   "4m":0,         "2m2e":0.001136,    "4e":0.002990   }
sig_2017 = {"4m":253.09,    "2m2e":167.67,  "4e":49.65}
npt_2017 = { "4l":17.58,     "4m":9.59,      "2m2e":6.54/2,    "4e":1.45,  "mumu":0, "ee":0, "2e2m":6.54/2, "ll":0}
npt_unc_2017 = { "4l":4.36,  "4m":3.16,  "2m2e":2.65,    "4e":1.42}

MUON_TRIG_LUMI_2017, ELEC_TRIG_LUMI_2017, ELEC_TRIG_SF_2017 = 41.37, 41.37, 0.991

NGEN_ZZ_4L,         XSEC_ZZ_4L          =   6893887,        1.256      #1.212
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      =   123485957,      5765.4    #6225.43
NGEN_ZJETS_M10,     XSEC_ZJETS_M10      =   39505301,       71310
NGEN_TTBAR,         XSEC_TTBAR          =   57584555,       831.76 / 2
NGEN_TT_2L2NU,      XSEC_TT_2L2NU       =   8926992,        87.31 / 2
NGEN_TTZ_2L2NU,     XSEC_TTZ_2L2NU      =   3570720,        0.2529
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       =   1992522,        12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        =   16664610,       5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        =   6887413,        4.42965
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       =   8733658,        0.564
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        =   17768294,       3.22
NGEN_WWZ_4L2NU,     XSEC_WWZ_4L2NU      =   1707572,        0.1651
NGEN_WZZ_4L2NU,     XSEC_WZZ_4L2NU      =   1690058,        0.05565
NGEN_ZZZ_4L2NU,     XSEC_ZZZ_4L2NU      =   1673322,        0.01398
NGEN_ZZG_4L2NU,     XSEC_ZZG_4L2NU      =   1455362,        0.01398
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      =   955384,         0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     =   984662,         0.001034

N_MC = 14

MC_SUFF = [ "zz_4l",            "zjets_m-50",       "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
            "wwz_4l2nu",        "wzz_4l2nu",        "zzz_4l2nu",        "zzg_4l2nu",
            "ggH_zz_4l",        "vbfH_zz_4l",       "ttz_2l2nu",
            ]

NGEN_   = [ NGEN_ZZ_4L,         NGEN_ZJETS_M50,     NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
            NGEN_WWZ_4L2NU,     NGEN_WZZ_4L2NU,     NGEN_ZZZ_4L2NU,     NGEN_ZZG_4L2NU,
            NGEN_GGH_ZZ_4L,     NGEN_VBFH_ZZ_4L,    NGEN_TTZ_2L2NU,
            ]

XSEC_   = [ XSEC_ZZ_4L,         XSEC_ZJETS_M50,     XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
            XSEC_WWZ_4L2NU,     XSEC_WZZ_4L2NU,     XSEC_ZZZ_4L2NU,     XSEC_ZZG_4L2NU,
            XSEC_GGH_ZZ_4L,     XSEC_VBFH_ZZ_4L,    XSEC_TTZ_2L2NU,
            ]

COLOR_  = [ lLightBlue,         lYellow,            lOrange,
            lOrange,            lOrange,            lOrange,            lOrange,
            lRed,               lRed,               lRed,               lRed,
            lPurple,            lPurple,            lGreen,
            ]

NGEN_2017   = dict(zip(MC_SUFF, NGEN_))
XSEC_2017   = dict(zip(MC_SUFF, XSEC_))



##
##  2018
##

YEAR_STR = "2018"

CAP_K_2018 = 0

ecal_2018 = {   "4m":0,         "2m2e":0,           "4e":0      }
sig_2018 = {"4m":368.52,    "2m2e":233.40,  "4e":74.09}
npt_2018 = {"4l":34.51,     "4m":22.10,     "2m2e":13.33,   "4e":0.91,  "mumu":0, "ee":0, "2e2m":0, "ll":0}
npt_unc_2018 = { "4l":6.09,      "4m":4.80,      "2m2e":3.74,    "4e":1.02}

MUON_TRIG_LUMI_2018, ELEC_TRIG_LUMI_2018, ELEC_TRIG_SF_2018 = 58.83, 58.83, 1

NGEN_ZZ_4L,         XSEC_ZZ_4L          =   6622242,        1.256      #1.212
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      =   100114403,      5765.4    #6225.43
NGEN_TTBAR,         XSEC_TTBAR          =   53887126,       831.76 / 2
NGEN_TT_2L2NU,      XSEC_TT_2L2NU       =   63791484,       87.31 / 2
NGEN_TTZ_2L2NU,     XSEC_TTZ_2L2NU      =   6274046,        0.2529
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       =   7729266,        12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        =   17048434,       5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        =   6739437,        4.42965
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       =   8371752,        0.564
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        =   17815223,       3.22
NGEN_WWZ_4L2NU,     XSEC_WWZ_4L2NU      =   1692397,        0.1651
NGEN_WZZ_4L2NU,     XSEC_WZZ_4L2NU      =   1722646,        0.05565
NGEN_ZZZ_4L2NU,     XSEC_ZZZ_4L2NU      =   1639163,        0.01398
NGEN_ZZG_4L2NU,     XSEC_ZZG_4L2NU      =   1527024,        0.01398
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      =   948128,         0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     =   499348,         0.001034

NGEN_   = [ NGEN_ZZ_4L,         NGEN_ZJETS_M50,     NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
            NGEN_WWZ_4L2NU,     NGEN_WZZ_4L2NU,     NGEN_ZZZ_4L2NU,     NGEN_ZZG_4L2NU,
            NGEN_GGH_ZZ_4L,     NGEN_VBFH_ZZ_4L,    NGEN_TTZ_2L2NU,
            ]

XSEC_   = [ XSEC_ZZ_4L,         XSEC_ZJETS_M50,     XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
            XSEC_WWZ_4L2NU,     XSEC_WZZ_4L2NU,     XSEC_ZZZ_4L2NU,     XSEC_ZZG_4L2NU,
            XSEC_GGH_ZZ_4L,     XSEC_VBFH_ZZ_4L,    XSEC_TTZ_2L2NU,
            ]

NGEN    = dict(zip(MC_SUFF, NGEN_))
XSEC    = dict(zip(MC_SUFF, XSEC_))
COLOR   = dict(zip(MC_SUFF, COLOR_))
