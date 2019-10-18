from __future__ import division

from PlotUtils import *

##
##  BASICS
##

YEAR_STR    = "2018"
EOS_PATH    = "root://cmseos.fnal.gov//store/user/jrainbol"
HOME_PATH   = "/uscms/home/jrainbol/nobackup"



##
##  SYSTEMATICS
##

F_NR = 0.04
BF_LL = 0.033658
GAMMA_Z = 2.4952

mu_id = {   "4l":0.0157,    "4m":0.0208,    "2m2e":0.0092,      "4e":0          }
el_id = {   "4l":0.0109,    "4m":0,         "2m2e":0.0188,      "4e":0.0406     }
el_reco = { "4l":0.0037,    "4m":0,         "2m2e":0.0093,      "4e":0.0152     }
ecal =  {   "4l":0,         "4m":0,         "2m2e":0,           "4e":0          }
qcd = 0.008592
pdf = 0.001185
pileup = 0.005730

DELTA_LAMBDA = 0.3




##
##  SAMPLE INFO
##

ELEC_TRIG_SF = 1
INT_LUMI, LUMI_UNC = 59.74, .025
SQRT_S  = 13
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

# Event numbers from TotalEvents histogram in "selected" ntuples (bin1 - 2*bin10)
NGEN_ZZ_4L,         XSEC_ZZ_4L          =   6622242,        1.256
NGEN_ZZ_4L_AMC,     XSEC_ZZ_4L_AMC      =   18530083,       1.212
NGEN_ZZ_4L_M1,      XSEC_ZZ_4L_M1       =   68654388,       7.957       #est
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      =   100114403,      5765.4
NGEN_TTBAR,         XSEC_TTBAR          =   53887126,       831.76 / 2
NGEN_TT_2L2NU,      XSEC_TT_2L2NU       =   63791484,       87.31 / 2
NGEN_TTZ_2L2NU,     XSEC_TTZ_2L2NU      =   6274046,        0.2529
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       =   7729266,        12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        =   17048434,       5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        =   6739437,        4.42965
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       =   8371752,        0.564
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        =   17815223,       3.22
NGEN_WWZ_4L2NU,     XSEC_WWZ_4L2NU      =   1692397,        0.0006024   #+- 5.770e-07 pb
NGEN_WZZ_4L2NU,     XSEC_WZZ_4L2NU      =   1722646,        0.0002692   #+- 2.656e-07 pb
NGEN_ZZZ_4L2NU,     XSEC_ZZZ_4L2NU      =   1639163,        0.0001907   #+- 1.778e-07 pb
NGEN_ZZG_4L2NU,     XSEC_ZZG_4L2NU      =   1527024,        0.001042    #+- 1.080e-06 pb
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      =   948128,         0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     =   499348,         0.001034

N_MC = 16
N_DY = 10   # number of gen zjets_m-50 files

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

NGEN    = dict(zip(MC_SUFF, NGEN_))
XSEC    = dict(zip(MC_SUFF, XSEC_))
COLOR   = dict(zip(MC_SUFF, COLOR_))
MC_TEX  = dict(zip(MC_SUFF, MC_TEX_))

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
