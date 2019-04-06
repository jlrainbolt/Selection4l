#ifndef CUTS2012_HH
#define CUTS2012_HH


// ROOT
#include "TString.h"



//
//  BASICS
//

const TString   YEAR_STR = "2012";
const TString   EOS_PATH = "root://cmseos.fnal.gov//store/user/jrainbol";
const TString   HOME_PATH = "/uscms/home/jrainbol/nobackup";

const unsigned  RNG_SEED = 2012;


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
const float     ELEC_PT1_MIN = 20,      ELEC_PT2_MIN = 10,      ELEC_PT_MIN = 7;
const float     MUON_ETA_MAX = 2.4,     ELEC_ETA_MAX = 2.5;
const float     MUON_ISO_MAX = 0.4,     ELEC_ISO_MAX = 0.4;
const float     SF_DR_MIN = 0.02,       OF_DR_MIN = 0.05;



//
//  LEPTON ID
//

const float     MUON_D0_MAX = 0.5,      MUON_DZ_MAX = 1,        MUON_SIP_MAX = 4;
const int       MUON_BAD_TRACK = 2;



//
//  OBJECT MATCHING
//

const float     MATCH_DR_MAX = 0.02;



//
//  SAMPLES
//

const float     MUON_TRIG_LUMI = 19.7,              ELEC_TRIG_LUMI = 19.7;
const float                                         ELEC_TRIG_SF = 1;
const TString   MU_SUFF = "muon_" + YEAR_STR,       EL_SUFF = "electron_" + YEAR_STR;

// Event numbers from TotalEvents histogram in "selected" ntuples (bin1 - 2*bin10)
//
// Cross sections from
//      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
const float     NGEN_ZZ         = 9799908,          XSEC_ZZ         = 17.0;
const float     NGEN_ZJETS_M50  = 30459503,         XSEC_ZJETS_M50  = 3531.9;
const float     NGEN_ZJETS_M10  = 33648307,         XSEC_ZJETS_M10  = 11050;
const float     NGEN_TTBAR      = 12011428,         XSEC_TTBAR      = 25.81;
const float     NGEN_TTZ_2L2NU  = 210160,           XSEC_TTZ_2L2NU  = 0.2057;
const float     NGEN_WW_2L2NU   = 10000431,         XSEC_WW_2L2NU   = 57.25;
const float     NGEN_WZ_2L2Q    = 3215990,          XSEC_WZ_2L2Q    = 5.09;
const float     NGEN_WZ_3LNU    = 2017979,          XSEC_WZ_3LNU    = 1.086;
const float     NGEN_ZZ_2L2NU   = 954911,           XSEC_ZZ_2L2NU   = 0.71;
const float     NGEN_ZZ_2L2Q    = 1936727,          XSEC_ZZ_2L2Q    = 2.47;

const unsigned  N_MC = 9,      ZZ = 0,     DY = 1;
const TString   MC_SUFF[N_MC] = {   "zz",           "zjets_m-50",   // "zjets_m-10",
                                    "ttbar",        "ttz_2l2nu",    "ww_2l2nu",
                                    "wz_2l2q",      "wz_3lnu",      "zz_2l2nu",     "zz_2l2q"
                                };
const float     NGEN[N_MC] = {  NGEN_ZZ,        NGEN_ZJETS_M50, // NGEN_ZJETS_M10,
                                NGEN_TTBAR,     NGEN_TTZ_2L2NU, NGEN_WW_2L2NU,
                                NGEN_WZ_2L2Q,   NGEN_WZ_3LNU,   NGEN_ZZ_2L2NU,  NGEN_ZZ_2L2Q
                                };
const float     XSEC[N_MC] = {  XSEC_ZZ,        XSEC_ZJETS_M50, // XSEC_ZJETS_M10,
                                XSEC_TTBAR,     XSEC_TTZ_2L2NU, XSEC_WW_2L2NU,
                                XSEC_WZ_2L2Q,   XSEC_WZ_3LNU,   XSEC_ZZ_2L2NU,  XSEC_ZZ_2L2Q
                                };
const int       COLOR[N_MC] = { lLightBlue,     lYellow,        // lYellow,
                                lGreen,         lGreen,         lOrange,
                                lOrange,        lOrange,        lOrange,        lOrange
                                };

#endif // CUTS2012_HH
