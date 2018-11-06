#ifndef CUTS2017_HH
#define CUTS2017_HH


// ROOT
#include "TString.h"



//
//  BASICS
//

const TString   YEAR_STR = "2017";
const TString   EOS_PATH = "root://cmseos.fnal.gov//store/user/jrainbol";
const TString   HOME_PATH = "/uscms/home/jrainbol/nobackup";



//
//  PHYSICAL CONSTANTS
//

const float     MUON_MASS = 0.105658369,                        ELEC_MASS = 0.000511;



//
//  PHASE SPACE
//

const float     M_MIN = 80,             M_MAX = 100,            MLL_MIN = 4;



//
//  FIDUCIAL REGION
//

const float     FID_PT1_MIN = 20,       FID_PT2_MIN = 10,       FID_PT_MIN = 5;
const float     FID_ETA_MAX = 2.5;



//
//  RECO SELECTION
//

const float     Z1_M_MIN = 12,          Z_M_MAX = 120;
const float     MUON_PT1_MIN = 20,      MUON_PT2_MIN = 10,      MUON_PT_MIN = 5;
const float     ELEC_PT1_MIN = 25,      ELEC_PT2_MIN = 15,      ELEC_PT_MIN = 7;
const float     MUON_ETA_MAX = 2.4,     ELEC_ETA_MAX = 2.5;
const float     MUON_ISO_MAX = 0.35,    ELEC_ISO_MAX = 0.35;
const float     SF_DR_MIN = 0.02,       OF_DR_MIN = 0.05;



//
//  OBJECT MATCHING
//

const float     MATCH_DR_MAX = 0.02;



//
//  SAMPLES
//

// Muon trigger lumi doesn't include 2017B
const float     MUON_TRIG_LUMI = 36.735,                        ELEC_TRIG_LUMI = 41.529;
const float                                                     ELEC_TRIG_SF = 0.991;
const TString   MU_SUFF = "muon_" + YEAR_STR,                   EL_SUFF = "electron_" + YEAR_STR;

// Event numbers from DAS, negative fractions from XSDB (FIXME?)
//
// Cross sections from
//      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
//      https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections
const float     NGEN_ZZ_4L      = 6967853   * (1 - 2 * 0.005244),       XSEC_ZZ_4L = 1.212;
const float     NGEN_ZJETS_M50  = 27413121  * (1 - 2 * 0.1624),         XSEC_ZJETS_M50 = 5765.4;
const float     NGEN_GGH_ZZ_4L  = 1000000   * (1 - 2 * 0.004958),       XSEC_GGH_ZZ_4L = 0.01212;
const float     NGEN_VBFH_ZZ_4L = 234800    * (1 - 2 * 0.00073),        XSEC_VBFH_ZZ_4L = 0.001034;
const float     NGEN_TTBAR      = 153531390 * (1 - 2 * 0.3163),         XSEC_TTBAR = 831.76;
const float     NGEN_WW_2L2NU   = 2000000   * (1 - 2 * 0.001928),       XSEC_WW_2L2NU = 12.178;
const float     NGEN_WZ_2L2Q    = 27582164  * (1 - 2 * 0.2006),         XSEC_WZ_2L2Q = 5.595;
const float     NGEN_WZ_3LNU    = 10881896  * (1 - 2 * 0.1879),         XSEC_WZ_3LNU = 4.42965;
const float     NGEN_ZZ_2L2Q    = 27840918  * (1 - 2 * 0.1804),         XSEC_ZZ_2L2Q = 3.22;

const unsigned  N_MC = 9,       ZZ = 0,     DY = 1;
const TString   MC_SUFF[N_MC] = {   "zz_4l",        "zjets_m-50",   "ggH_zz_4l",    "vbfH_zz_4l",
                                    "ttbar",        "ww_2l2nu",     "wz_2l2q",      "wz_3lnu",
                                    "zz_2l2q"};
const float     NGEN[N_MC] = {  NGEN_ZZ_4L,     NGEN_ZJETS_M50, NGEN_GGH_ZZ_4L, NGEN_VBFH_ZZ_4L,
                                NGEN_TTBAR,     NGEN_WW_2L2NU,  NGEN_WZ_2L2Q,   NGEN_WZ_3LNU,
                                NGEN_ZZ_2L2Q};
const float     XSEC[N_MC] = {  XSEC_ZZ_4L,     XSEC_ZJETS_M50, XSEC_GGH_ZZ_4L, XSEC_VBFH_ZZ_4L,
                                XSEC_TTBAR,     XSEC_WW_2L2NU,  XSEC_WZ_2L2Q,   XSEC_WZ_3LNU,
                                XSEC_ZZ_2L2Q};
const int       COLOR[N_MC] = { lLightBlue,     lYellow,        lPurple,        lPurple,
                                lGreen,         lOrange,        lOrange,        lOrange,
                                lOrange};
const TString   MC_TEX[N_MC] = {"\\Z\\Zto",     "\\DY",         "\\text{VBF }\\Higgs\\to\\Z\\Zto",
                                "\\text{ggF }\\Higgs\\to\\Z\\Zto",      "\\tq\\tbar",
                                "\\W\\W\\to 2\\ell 2\\nu",      "\\W\\Z\\to 2\\ell 2\\qq",
                                "\\W\\Z\\to 3\\ell \\nu",       "\\Z\\Z\\to 2\\ell 2\\qq"};




#endif // CUTS2017_HH
