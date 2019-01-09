#ifndef CUTS2017_HH
#define CUTS2017_HH


// ROOT
#include "TString.h"



//
//  BASICS
//

const TString   YEAR_STR = "2016";
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
const float     MUON_PT1_MIN = 25,      MUON_PT2_MIN = 10,      MUON_PT_MIN = 5;
const float     ELEC_PT1_MIN = 30,      ELEC_PT2_MIN = 10,      ELEC_PT_MIN = 7;
const float     MUON_ETA_MAX = 2.4,     ELEC_ETA_MAX = 2.5;
const float     MUON_ISO_MAX = 0.35,    ELEC_ISO_MAX = 0.35;
const float     SF_DR_MIN = 0.02,       OF_DR_MIN = 0.05;



//
//  TIGHT ID
//

const float     MUON_TIGHT_ISO_MAX = 0.15;
const float     EB_TIGHT_ISO_MAX = 0.0588,                      EE_TIGHT_ISO_MAX = 0.0571;
const float     EB_ETA_MAX = 1.479;



//
//  OBJECT MATCHING
//

const float     MATCH_DR_MAX = 0.02;



//
//  SYSTEMATICS
//

const float     MUON_PT_SHIFT = 0.002,  ELEC_PT_SHIFT = 0.005;



//
//  SAMPLES
//

// Muon trigger lumi doesn't include 2017B
const float     MUON_TRIG_LUMI = 35.9,                          ELEC_TRIG_LUMI = 35.9;
const float                                                     ELEC_TRIG_SF = 1;
const TString   MU_SUFF = "muon_" + YEAR_STR,                   EL_SUFF = "electron_" + YEAR_STR;

// Event numbers from hTotalEvents with negative events subtracted (notebook p 64)
//
// Cross sections from
//      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
//      https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections
const float     NGEN_ZZ_4L      = 6838244,                      XSEC_ZZ_4L      = 1.212;
const float     NGEN_ZJETS_M50  = 75087060,                     XSEC_ZJETS_M50  = 5765.4;
//const float     NGEN_ZJETS_M50  = 81616008,                     XSEC_ZJETS_M50  = 5765.4 * 1.18;
const float     NGEN_ZJETS_M10  = 99482541,                     XSEC_ZJETS_M10  = 18610;
const float     NGEN_GGH_ZZ_4L  = 999778,                       XSEC_GGH_ZZ_4L  = 0.01212;
const float     NGEN_VBFH_ZZ_4L = 499293,                       XSEC_VBFH_ZZ_4L = 0.001034;
const float     NGEN_TTBAR      = 155211361,                    XSEC_TTBAR      = 831.76;
const float     NGEN_TTZ_2L2NU  = 3717852,                      XSEC_TTZ_2L2NU  = 0.2529;
const float     NGEN_WW_2L2NU   = 1998956,                      XSEC_WW_2L2NU   = 12.178;
const float     NGEN_WZ_3LNU    = 1993154,                      XSEC_WZ_3LNU    = 4.42965;

const unsigned  N_MC = 9,       ZZ = 0,     DY = 1;
const TString   MC_SUFF[N_MC]={ "zz_4l",        "dy_m-50",      "dy_m-10to50",  "ggH_zz_4l",
                                "H_zz_4l",      "ttbar",        "ttz_2l2nu",    "ww_2l2nu",
                                "wz_3lnu"};
const float     NGEN[N_MC] = {  NGEN_ZZ_4L,     NGEN_ZJETS_M50, NGEN_ZJETS_M10, NGEN_GGH_ZZ_4L,
                                NGEN_VBFH_ZZ_4L,NGEN_TTBAR,     NGEN_TTZ_2L2NU, NGEN_WW_2L2NU,
                                NGEN_WZ_3LNU};
const float     XSEC[N_MC] = {  XSEC_ZZ_4L,     XSEC_ZJETS_M50, XSEC_ZJETS_M10, XSEC_GGH_ZZ_4L,
                                XSEC_VBFH_ZZ_4L,XSEC_TTBAR,     XSEC_TTZ_2L2NU, XSEC_WW_2L2NU,
                                XSEC_WZ_3LNU};
const int       COLOR[N_MC] = { lLightBlue,     lYellow,        lYellow,        lPurple,
                                lPurple,        lGreen,         lGreen,         lOrange,
                                lOrange};
const TString   MC_TEX[N_MC] = {"\\Z\\Zto",     "\\DY",         "\\DY",
                                "\\text{VBF }\\Higgs\\to\\Z\\Zto",
                                "\\text{ggF }\\Higgs\\to\\Z\\Zto",
                                "\\tq\\tbar",                   "\\tq\\tq\\to 2\\ell 2\\nu",
                                "\\W\\W\\to 2\\ell 2\\nu",      "\\W\\Z\\to 3\\ell \\nu"};




#endif // CUTS2017_HH
