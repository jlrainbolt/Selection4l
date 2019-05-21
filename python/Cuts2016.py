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

mu_id = {   "4l":0.0484,    "4m":0.0676,    "2m2e":0.0350,      "4e":0          }
el_id = {   "4l":0.0219,    "4m":0,         "2m2e":0.0371,      "4e":0.0804     }
el_reco = { "4l":0.0093,    "4m":0,         "2m2e":0.0156,      "4e":0.0361     }
mu_pt = {   "4l":0.000042,  "4m":0.000071,  "2m2e":0.000042,    "4e":0          }
el_pt = {   "4l":0.001017,  "4m":0,         "2m2e":0.002930,    "4e":0.000169   }
ecal =  {   "4l":0.000799,  "4m":0,         "2m2e":0.000536,    "4e":0.001595   }
qcd = 0.008592
pdf = 0.001185
pileup = 0.008885

npt_unc = { "4l":4.00,      "4m":2.45,      "2m2e":2.83,    "4e":1.42, "mumu":0, "ee":0, "2e2m":0,  "ll":0}
npt     = { "4l":14.19,     "4m":5.68,      "2m2e":7.45,    "4e":1.06, "mumu":0, "ee":0, "2e2m":0,  "ll":0}   



##
##  SAMPLE INFO
##

# Muon trigger lumi doesn't include 2017B
MUON_TRIG_LUMI, ELEC_TRIG_LUMI, ELEC_TRIG_SF = 36.42, 36.42, 1
LUMI_UNC = .025;
SQRT_S  = 13
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

# Event numbers from DAS, negative fractions from XSDB (FIXME?)
#
# Cross sections from
#      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
#      https://twiki.cern.ch/twiki/bin/view/ns
NGEN_ZZ_4L,         XSEC_ZZ_4L          = 5932142,        1.256     #1.212
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      = 74513923,       5765.4    #6225.4
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

N_MC = 11
N_DY = 10

MC_SUFF = [ "zz_4l",            "zjets_m-50",     # "zjets_m-10",
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

COLOR_  = [ lLightBlue,         lYellow,          # lYellow,
            lGreen,             lGreen,           # lGreen,             
            lOrange,
            lOrange,            lOrange,            lOrange,            lOrange,
#           lRed,               lRed,               lRed,               lRed,
            lPurple,            lPurple,
            ]

MC_TEX_ = [ r"$\ZZtofL$",       r"$\ZtoLL$", # r"$\ZtoLL$ $(10 < \mll < 50\GeV)$",
            r"$\ttbar$",        r"$\TTtoLLNuNu$", # r"$\TTZtoLLNuNu$",
            r"$\WWtoLLNuNu$",
            r"$\WZtoLLQQ$",     r"$\WZtoLLLNu$",    r"$\ZZtoLLNuNu$",   r"$\ZZtoLLQQ$",
#           r"$\WWZtofLtNu$",   r"$\WZZtofLtNu$",   r"$\ZZZtofLtNu$",   r"$\ZZGtofLtNu$",
            r"$\ggF\HtoZZtofL$",r"$\VBF\HtoZZtofL$",
            ]

NGEN    = dict(zip(MC_SUFF, NGEN_))
XSEC    = dict(zip(MC_SUFF, XSEC_))
COLOR   = dict(zip(MC_SUFF, COLOR_))
MC_TEX  = dict(zip(MC_SUFF, MC_TEX_))

MC_SUFF_4L = MC_SUFF
MC_SUFF_2L = list(MC_SUFF)
MC_SUFF_2L[0], MC_SUFF_2L[1] = MC_SUFF[1], MC_SUFF[0]
