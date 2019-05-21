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
CAP_K = 3.3016956e-10
GAMMA_Z = 2.4952e6

mu_id = {   "4l":0.0484,    "4m":0.0676,    "2m2e":0.0350,      "4e":0          }
el_id = {   "4l":0.0219,    "4m":0,         "2m2e":0.0371,      "4e":0.0804     }
el_reco = { "4l":0.0093,    "4m":0,         "2m2e":0.0156,      "4e":0.0361     }
mu_pt = {   "4l":0.000042,  "4m":0.000071,  "2m2e":0.000042,    "4e":0          }
el_pt = {   "4l":0.001017,  "4m":0,         "2m2e":0.002930,    "4e":0.000169   }
ecal =  {   "4l":0,         "4m":0,         "2m2e":0,           "4e":0          }
qcd = 0.008592
pdf = 0.001185
pileup = 0.005730

npt     = { "4l":34.51,     "4m":22.10,     "2m2e":13.33,   "4e":0.91,  "mumu":0, "ee":0, "2e2m":0, "ll":0}
npt_unc = { "4l":6.09,      "4m":4.80,      "2m2e":3.74,    "4e":1.02,  "mumu":0, "ee":0, "2e2m":0, "ll":0}




##
##  SAMPLE INFO
##

# Muon trigger lumi doesn't include 2018B
MUON_TRIG_LUMI, ELEC_TRIG_LUMI, ELEC_TRIG_SF = 58.83, 58.83, 1
LUMI_UNC = .025;
SQRT_S  = 13
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

# Event numbers from TotalEvents histogram in "selected" ntuples (bin1 - 2*bin10)
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

N_MC = 16
N_DY = 10

MC_SUFF = [ "zz_4l",            "zjets_m-50",
            "ttbar",            "tt_2l2nu",         "ttz_2l2nu",        "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2nu",         "zz_2l2q",
            "wwz_4l2nu",        "wzz_4l2nu",        "zzz_4l2nu",        "zzg_4l2nu",
            "ggH_zz_4l",        "vbfH_zz_4l",
            ]

NGEN_   = [ NGEN_ZZ_4L,         NGEN_ZJETS_M50,
            NGEN_TTBAR,         NGEN_TT_2L2NU,      NGEN_TTZ_2L2NU,     NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2NU,      NGEN_ZZ_2L2Q,
            NGEN_WWZ_4L2NU,     NGEN_WZZ_4L2NU,     NGEN_ZZZ_4L2NU,     NGEN_ZZG_4L2NU,
            NGEN_GGH_ZZ_4L,     NGEN_VBFH_ZZ_4L,
            ]

XSEC_   = [ XSEC_ZZ_4L,         XSEC_ZJETS_M50,
            XSEC_TTBAR,         XSEC_TT_2L2NU,      XSEC_TTZ_2L2NU,     XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2NU,      XSEC_ZZ_2L2Q,
            XSEC_WWZ_4L2NU,     XSEC_WZZ_4L2NU,     XSEC_ZZZ_4L2NU,     XSEC_ZZG_4L2NU,
            XSEC_GGH_ZZ_4L,     XSEC_VBFH_ZZ_4L,
            ]

COLOR_  = [ lLightBlue,         lYellow,
            lGreen,             lGreen,             lGreen,             lOrange,
            lOrange,            lOrange,            lOrange,            lOrange,
            lRed,               lRed,               lRed,               lRed,
            lPurple,            lPurple,
            ]

MC_TEX_ = [ r"$\ZZtofL$",       r"$\ZtoLL$",
            r"$\ttbar$",        r"$\TTtoLLNuNu$",   r"$\TTZtoLLNuNu$",  r"$\WWtoLLNuNu$",
            r"$\WZtoLLQQ$",     r"$\WZtoLLLNu$",    r"$\ZZtoLLNuNu$",   r"$\ZZtoLLQQ$",
            r"$\WWZtofLtNu$",   r"$\WZZtofLtNu$",   r"$\ZZZtofLtNu$",   r"$\ZZGtofLtNu$",
            r"$\ggF\HtoZZtofL$",r"$\VBF\HtoZZtofL$",
            ]

NGEN    = dict(zip(MC_SUFF, NGEN_))
XSEC    = dict(zip(MC_SUFF, XSEC_))
COLOR   = dict(zip(MC_SUFF, COLOR_))
MC_TEX  = dict(zip(MC_SUFF, MC_TEX_))

MC_SUFF_4L = MC_SUFF
MC_SUFF_2L = list(MC_SUFF)
MC_SUFF_2L[0], MC_SUFF_2L[1] = MC_SUFF[1], MC_SUFF[0]
