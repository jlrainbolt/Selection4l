from __future__ import division

from PlotUtils import *

##
##  BASICS
##

YEAR_STR    = "2012"
EOS_PATH    = "root://cmseos.fnal.gov//store/user/jrainbol"
HOME_PATH   = "/uscms/home/jrainbol/nobackup"



##
##  SYSTEMATICS
##

F_NR = 0.04
BF_LL = 0.033658 * 2
CAP_K = 5.10487652e-10
GAMMA_Z = 2.4952e6

mu_id = {   "4l":0,         "4m":0,         "2m2e":0,       "4e":0 }
el_id = {   "4l":0,         "4m":0,         "2m2e":0,       "4e":0 }
el_reco = { "4l":0,         "4m":0,         "2m2e":0,       "4e":0 }
mu_pt = {   "4l":0,         "4m":0,         "2m2e":0,       "4e":0 }
el_pt = {   "4l":0,         "4m":0,         "2m2e":0,       "4e":0 }

npt     = { "4l":0,     "4m":0,      "2m2e":0,    "4e":0, "mumu":0, "ee":0, "2e2m":0   }   
npt_unc = { "4l":0,      "4m":0,      "2m2e":0,    "4e":0, "mumu":0, "ee":0, "2e2m":0   }



##
##  SAMPLE INFO
##

MUON_TRIG_LUMI, ELEC_TRIG_LUMI, ELEC_TRIG_SF = 19.712, 19.712, 1
SQRT_S  = 8
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

NGEN_ZZ,            XSEC_ZZ             = 9799908,        17.0
NGEN_ZZ_4L,         XSEC_ZZ_4L          = 4807893,        0.3305
NGEN_ZZ_4M,         XSEC_ZZ_4M          = 1499064,        0.07691 * 5
NGEN_ZZ_2M2E,       XSEC_ZZ_2M2E        = 1488445,        0.1767 * 5
NGEN_ZZ_4E,         XSEC_ZZ_4E          = 1499093,        0.07691 * 5
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      = 30459503,       3531.9
NGEN_ZJETS_M10,     XSEC_ZJETS_M10      = 33648307,       11050
NGEN_TTBAR,         XSEC_TTBAR          = 12011428,       25.81
NGEN_TTZ_2L2NU,     XSEC_TTZ_2L2NU      = 210160,         0.2057
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       = 10000431,       57.25
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        = 3215990,        5.09
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        = 2017979,        1.086
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       = 1936727,        2.47
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        = 954911,         0.71

N_MC = 11
N_DY = 1

MC_SUFF = [ "zz_4m",            "zz_2m2e",          "zz_4e",            "zjets_m-50",
            "ttbar",            "ttz_2l2nu",        "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
            ]

NGEN_   = [ NGEN_ZZ_4M,         NGEN_ZZ_2M2E,       NGEN_ZZ_4E,         NGEN_ZJETS_M50,
            NGEN_TTBAR,         NGEN_TTZ_2L2NU,     NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
            ]

XSEC_   = [ XSEC_ZZ_4M,         XSEC_ZZ_2M2E,       XSEC_ZZ_4E,         XSEC_ZJETS_M50,
            XSEC_TTBAR,         XSEC_TTZ_2L2NU,     XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
            ]

COLOR_  = [ lLightBlue,         lLightBlue,         lLightBlue,         lYellow,
            lGreen,             lGreen,             lOrange,
            lOrange,            lOrange,            lOrange,            lOrange,
            ]

MC_TEX_ = [ r"$\ZZtofM$",       r"$\ZZtotMtE$",     r"$\ZZtofE$",       r"$\ZtoLL$",
            r"$\ttbar$",        r"$\TTZtoLLNuNu$",  r"$\WWtoLLNuNu$",
            r"$\WZtoLLQQ$",     r"$\WZtoLLLNu$",    r"$\ZZtoLLNuNu$",   r"$\ZZtoLLQQ$",
            ]

NGEN    = dict(zip(MC_SUFF, NGEN_))
XSEC    = dict(zip(MC_SUFF, XSEC_))
COLOR   = dict(zip(MC_SUFF, COLOR_))
MC_TEX  = dict(zip(MC_SUFF, MC_TEX_))

MC_SUFF_4L = MC_SUFF
MC_SUFF_2L = list(MC_SUFF)
MC_SUFF_2L[0], MC_SUFF_2L[3] = MC_SUFF[3], MC_SUFF[0]
