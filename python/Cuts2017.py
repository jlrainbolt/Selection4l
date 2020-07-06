from __future__ import division

from PlotUtils import *

##
##  BASICS
##

YEAR_STR    = "2017"
EOS_PATH    = "root://cmseos.fnal.gov//store/user/jrainbol"
HOME_PATH   = "/uscms/home/jrainbol/nobackup"



##
##  SYSTEMATICS
##

F_NR = 0.04
BF_LL = 0.033658
CAP_K = 4.6494342e-10
GAMMA_Z = 2.4952

mu_id = {   "4l":0.0154,    "4m":0.0208,    "2m2e":0.0093,      "4e":0          }
el_id = {   "4l":0.0077,    "4m":0,         "2m2e":0.0138,      "4e":0.0274     }
el_reco = { "4l":0.0176,    "4m":0,         "2m2e":0.0247,      "4e":0.0588     }
ecal =  {   "4l":0.0000,    "4m":0.000400,  "2m2e":0.001126,    "4e":0.002811 + 0.0016 }
qcd = 0.008592
pdf = 0.001185
pileup = 0.005765

# Systematic uncertainty for nonprompt
DELTA_LAMBDA = 0.3

# Systematic uncertainties for MC
UNC_DIBOSON,        UNC_TTBAR,      UNC_TAUTAU,     UNC_OTHER    =  0.1,    0.1,    0.05,   0.2




##
##  SAMPLE INFO
##

# Muon trigger lumi doesn't include 2017B
INT_LUMI, LUMI_UNC = 41.53, .023
SQRT_S  = 13
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

# Event numbers from TotalEvents histogram in "selected" ntuples (bin1 - 2*bin10)
#
# Cross sections from
#      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
NGEN_ZZ_4L,         XSEC_ZZ_4L          =   6893887,        1.256
NGEN_ZZ_4L_AMC,     XSEC_ZZ_4L_AMC      =   18530335,       1.212
NGEN_ZZ_4L_M1,      XSEC_ZZ_4L_M1       =   68016559,       7.957       # est
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      =   123485957,      5765.4
NGEN_TTBAR,         XSEC_TTBAR          =   57584555,       831.76 / 2
NGEN_TT_2L2NU,      XSEC_TT_2L2NU       =   8926992,        87.31 / 2
NGEN_TTZ_2L2NU,     XSEC_TTZ_2L2NU      =   3570720,        0.2529
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       =   1992522,        12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        =   16664610,       5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        =   6887413,        4.42965
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       =   8733658,        0.564
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        =   17768294,       3.22
NGEN_WWZ_4L2NU,     XSEC_WWZ_4L2NU      =   1707572,        0.0006024   #+- 5.770e-07 pb
NGEN_WZZ_4L2NU,     XSEC_WZZ_4L2NU      =   1690058,        0.0002692   #+- 2.656e-07 pb
NGEN_ZZZ_4L2NU,     XSEC_ZZZ_4L2NU      =   1673322,        0.0001907   #+- 1.778e-07 pb
NGEN_ZZG_4L2NU,     XSEC_ZZG_4L2NU      =   1455362,        0.001042    #+- 1.080e-06 pb
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      =   955384,         0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     =   984662,         0.001034

N_MC = 16
N_DY = 15

MC_SUFF = [ "zz_4l",            "zz_4l_aMC",        "zz_4l_m-1",        "zjets_m-50",
            "ttbar",            "tt_2l2nu",                             "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
            "wwz_4l2nu",        "wzz_4l2nu",        "zzz_4l2nu",        "zzg_4l2nu",
            "ggH_zz_4l",        "vbfH_zz_4l",                           "ttz_2l2nu",
            ]

NGEN_   = [ NGEN_ZZ_4L,         NGEN_ZZ_4L_AMC,     NGEN_ZZ_4L_M1,      NGEN_ZJETS_M50,
            NGEN_TTBAR,         NGEN_TT_2L2NU,                          NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
            NGEN_WWZ_4L2NU,     NGEN_WZZ_4L2NU,     NGEN_ZZZ_4L2NU,     NGEN_ZZG_4L2NU,
            NGEN_GGH_ZZ_4L,     NGEN_VBFH_ZZ_4L,                        NGEN_TTZ_2L2NU,
            ]

XSEC_   = [ XSEC_ZZ_4L,         XSEC_ZZ_4L_AMC,     XSEC_ZZ_4L_M1,      XSEC_ZJETS_M50,
            XSEC_TTBAR,         XSEC_TT_2L2NU,                          XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
            XSEC_WWZ_4L2NU,     XSEC_WZZ_4L2NU,     XSEC_ZZZ_4L2NU,     XSEC_ZZG_4L2NU,
            XSEC_GGH_ZZ_4L,     XSEC_VBFH_ZZ_4L,                        XSEC_TTZ_2L2NU,
            ]

COLOR_  = [ lLightBlue,         lBlue,              lLightBlue,         lYellow,
            lGreen,             lGreen,                                 lOrange,
            lOrange,            lOrange,            lOrange,            lOrange,
            lRed,               lRed,               lRed,               lRed,
            lPurple,            lPurple,                                lGreen,
            ]

MC_TEX_ = [ r"$\ZZtofL$",       r"$\ZZtofL$ (a\textsc{mc@nlo})",        r"$\ZZtofL$",
            r"$\ZtoLL$",        r"$\ttbar$",        r"$\TTtoLLNuNu$",   r"$\WWtoLLNuNu$",
            r"$\WZtoLLQQ$",     r"$\WZtoLLLNu$",    r"$\ZZtoLLNuNu$",   r"$\ZZtoLLQQ$",
            r"$\WWZtofLtNu$",   r"$\WZZtofLtNu$",   r"$\ZZZtofLtNu$",   r"$\ZZGtofLtNu$",
            r"$\ggF\HtoZZtofL$",r"$\VBF\HtoZZtofL$",                    r"$\TTZtoLLNuNu$",
            ]

MC_UNC_ = [ UNC_DIBOSON,        UNC_DIBOSON,        UNC_DIBOSON,        UNC_TAUTAU,
            UNC_TTBAR,          UNC_TTBAR,                              UNC_DIBOSON,
            UNC_DIBOSON,        UNC_DIBOSON,        UNC_DIBOSON,        UNC_DIBOSON,
            UNC_OTHER,          UNC_OTHER,          UNC_OTHER,          UNC_OTHER,
            UNC_OTHER,          UNC_OTHER,                              UNC_OTHER,
            ]

NGEN    = dict(zip(MC_SUFF, NGEN_))
XSEC    = dict(zip(MC_SUFF, XSEC_))
COLOR   = dict(zip(MC_SUFF, COLOR_))
MC_TEX  = dict(zip(MC_SUFF, MC_TEX_))
MC_UNC  = dict(zip(MC_SUFF, MC_UNC_))

MC_SUFF_4L = list(MC_SUFF)
MC_SUFF_4L.remove("zz_4l_aMC")
MC_SUFF_4L.remove("zz_4l_m-1")

MC_SUFF_2L = list(MC_SUFF)
MC_SUFF_2L.remove("zz_4l_aMC")
MC_SUFF_2L.remove("zz_4l_m-1")
MC_SUFF_2L[0], MC_SUFF_2L[1] = MC_SUFF_2L[1], MC_SUFF_2L[0]

MC_SUFF_AMC = list(MC_SUFF)
MC_SUFF_AMC.remove("zz_4l")
MC_SUFF_AMC.remove("zz_4l_m-1")
MC_SUFF_AMC.append("zz_4l")

MC_SUFF_MLL = list(MC_SUFF)
MC_SUFF_MLL.remove("zz_4l")
MC_SUFF_MLL.remove("zz_4l_aMC")

MC_SUFF.remove("zz_4l_aMC")
MC_SUFF.remove("zz_4l_m-1")

MC_SUFF_HZZ = [ "ggH_zz_4l",        "vbfH_zz_4l",       "zz_4l",            "zjets_m-50",
                "ttbar",            "tt_2l2nu",                             "ww_2l2nu",
                "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
                "wwz_4l2nu",        "wzz_4l2nu",        "zzz_4l2nu",        "zzg_4l2nu",
                "ttz_2l2nu",
                ]
