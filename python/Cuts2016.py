from __future__ import division

from PlotUtils import *

##
##  BASICS
##

YEAR_STR    = "2016"
EOS_PATH    = "root://cmseos.fnal.gov//store/user/jrainbol"
HOME_PATH   = "/uscms/home/jrainbol/nobackup"



##
##  SAMPLE INFO
##

# Muon trigger lumi doesn't include 2017B
MUON_TRIG_LUMI, ELEC_TRIG_LUMI, ELEC_TRIG_SF = 35.9, 35.9, 1
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

# Event numbers from DAS, negative fractions from XSDB (FIXME?)
#
# Cross sections from
#      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
#      https://twiki.cern.ch/twiki/bin/view/ns
NGEN_ZZ_4L,         XSEC_ZZ_4L          = 6762740,        1.212
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      = 80924255,       5765.4
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      = 999800,         0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     = 499262,         0.001034
NGEN_TTBAR,         XSEC_TTBAR          = 12284545,       831.76
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       = 1832358,        12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        = 15879472,       5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        = 7387013,        4.42965
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        = 496436,         3.22
NGEN_ZZ_2L2NU,      XSEC_ZZ_2L2NU       = 48623080,       0.564

N_MC = 10
MC_SUFF = [ "zjets_m-50",               "zz_4l",                        "ttbar",
            "ww_2l2nu",                 "wz_2l2q",                      "wz_3lnu",
            "zz_2l2nu",                 "zz_2l2q",
            "vbfH_zz_4l",               "ggH_zz_4l"
            ]
MC_SUFF_4L = [ "zz_4l",                 "zjets_m-50",                   "ttbar",
            "ww_2l2nu",                 "wz_2l2q",                      "wz_3lnu",
            "zz_2l2nu",                 "zz_2l2q",
            "vbfH_zz_4l",               "ggH_zz_4l"
            ]
XSEC    = { "zz_4l":XSEC_ZZ_4L,         "zjets_m-50":XSEC_ZJETS_M50,    "ggH_zz_4l":XSEC_GGH_ZZ_4L,
            "vbfH_zz_4l":XSEC_VBFH_ZZ_4L,   "ttbar":XSEC_TTBAR,         "ww_2l2nu":XSEC_WW_2L2NU,
            "wz_2l2q":XSEC_WZ_2L2Q,     "wz_3lnu":XSEC_WZ_3LNU,         "zz_2l2q":XSEC_ZZ_2L2Q,
            "zz_2l2nu":XSEC_ZZ_2L2NU
            }
NGEN    = { "zz_4l":NGEN_ZZ_4L,         "zjets_m-50":NGEN_ZJETS_M50,    "ggH_zz_4l":NGEN_GGH_ZZ_4L,
            "vbfH_zz_4l":NGEN_VBFH_ZZ_4L,   "ttbar":NGEN_TTBAR,         "ww_2l2nu":NGEN_WW_2L2NU,
            "wz_2l2q":NGEN_WZ_2L2Q,     "wz_3lnu":NGEN_WZ_3LNU,         "zz_2l2q":NGEN_ZZ_2L2Q,
            "zz_2l2nu":NGEN_ZZ_2L2NU
            }
COLOR   = { "zz_4l":lLightBlue,         "zjets_m-50":lYellow,           "ggH_zz_4l":lPurple,
            "vbfH_zz_4l":lPurple,       "ttbar":lGreen,                 "ww_2l2nu":lOrange,
            "wz_2l2q":lOrange,          "wz_3lnu":lOrange,              "zz_2l2q":lOrange,
            "zz_2l2nu":lOrange
            }
MC_TEX  = { "zz_4l":r"$\ZZtofl$",                   "zjets_m-50":r"$\Ztoll$",
            "vbfH_zz_4l":r"VBF $\PH\to\ZZtofl$",    "ggH_zz_4l":r"$\Pg\Pg$F $\PH \to \ZZtofl$",
            "ttbar":r"$\ttbar$",                    "ww_2l2nu":r"$\WW \to 2\Pell2\PGn$",
            "wz_2l2q":r"$\WZ \to 2\Pell2\PQq$",     "wz_3lnu":r"$\WZ \to 3\Pell\PGn$",
            "zz_2l2q":r"$\ZZ \to 2\Pell2\PQq$",     "zz_2l2nu":r"$\ZZ \to 2\Pell2\PGn$"
            }
