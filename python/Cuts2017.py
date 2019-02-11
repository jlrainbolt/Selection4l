from PlotUtils import *

##
##  BASICS
##

YEAR_STR    = "2017"
EOS_PATH    = "root://cmseos.fnal.gov//store/user/jrainbol"
HOME_PATH   = "/uscms/home/jrainbol/nobackup"



##
##  SAMPLE INFO
##

# Muon trigger lumi doesn't include 2017B
MUON_TRIG_LUMI, ELEC_TRIG_LUMI, ELEC_TRIG_SF = 36.735, 45.529, 0.991
MU_SUFF, EL_SUFF = "muon_" + YEAR_STR, "electron_" + YEAR_STR

# Event numbers from DAS, negative fractions from XSDB (FIXME?)
#
# Cross sections from
#      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
#      https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections
NGEN_ZZ_4L,         XSEC_ZZ_4L          =   6967853   * (1 - 2 * 0.005244),           1.212
NGEN_ZJETS_M50,     XSEC_ZJETS_M50      =   27413121  * (1 - 2 * 0.1624),             5765.4
NGEN_GGH_ZZ_4L,     XSEC_GGH_ZZ_4L      =   1000000   * (1 - 2 * 0.004958),           0.01212
NGEN_VBFH_ZZ_4L,    XSEC_VBFH_ZZ_4L     =   234800    * (1 - 2 * 0.00073),            0.001034
NGEN_TTBAR,         XSEC_TTBAR          =   (153531390-626036) * (1 - 2 * 0.3163),    831.76
NGEN_WW_2L2NU,      XSEC_WW_2L2NU       =   2000000   * (1 - 2 * 0.001928),           12.178
NGEN_WZ_2L2Q,       XSEC_WZ_2L2Q        =   27582164  * (1 - 2 * 0.2006),             5.595
NGEN_WZ_3LNU,       XSEC_WZ_3LNU        =   10881896  * (1 - 2 * 0.1879),             4.42965
NGEN_ZZ_2L2Q,       XSEC_ZZ_2L2Q        =   27840918  * (1 - 2 * 0.1804),             3.22

N_MC, ZZ, DY = 9, 0, 1
MC_SUFF = [ "zz_4l",            "zjets_m-50",       "ggH_zz_4l",
            "vbfH_zz_4l",       "ttbar",            "ww_2l2nu",
            "wz_2l2q",          "wz_3lnu",          "zz_2l2q"       ]
NGEN    = [ NGEN_ZZ_4L,         NGEN_ZJETS_M50,     NGEN_GGH_ZZ_4L,
            NGEN_VBFH_ZZ_4L,    NGEN_TTBAR,         NGEN_WW_2L2NU,
            NGEN_WZ_2L2Q,       NGEN_WZ_3LNU,       NGEN_ZZ_2L2Q    ]
XSEC    = [ XSEC_ZZ_4L,         XSEC_ZJETS_M50,     XSEC_GGH_ZZ_4L,
            XSEC_VBFH_ZZ_4L,    XSEC_TTBAR,         XSEC_WW_2L2NU,
            XSEC_WZ_2L2Q,       XSEC_WZ_3LNU,       XSEC_ZZ_2L2Q    ]
COLOR   = [ lLightBlue,         lYellow,            lPurple,
            lPurple,            lGreen,             lOrange,
            lOrange,            lOrange,            lOrange         ]
