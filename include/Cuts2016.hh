#ifndef CUTS2016_HH
#define CUTS2016_HH


// ROOT
#include "TString.h"



//
//  BASICS
//

const TString   YEAR_STR = "2016";
const TString   EOS_PATH = "root://cmseos.fnal.gov//store/user/jrainbol";
const TString   HOME_PATH = "/uscms/home/jrainbol/nobackup";
const TString   BLT_PATH = HOME_PATH + "/lacey_legacy/CMSSW_10_2_13/src/BLT";

const unsigned  RNG_SEED = 6102;


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
//  TRIGGER REQUIREMENTS
//

const float     MUON_SINGLE_PT = 24,    ELEC_SINGLE_PT = 27;
const float     MUON_LEG1_PT = 17,      ELEC_LEG1_PT = 23;
const float     MUON_LEG2_PT = 8,       ELEC_LEG2_PT = 12;



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
//  SYSTEMATICS
//

const float     MUON_PT_SHIFT = 0.002,  ELEC_PT_SHIFT = 0.005;

const unsigned  N_HISTS = 100;
//  ll          mumu        ee
const float     AXE_2L[3] = {   0.320909,   0.376059,   0.265793};
//  4l          4m          2m2e        4e
const float     AXE_4L[4] = {   0.043777,   0.0868125,  0.0335567,  0.0204522};



//
//  SAMPLES
//

const float     MUON_TRIG_LUMI = 36.42,                         ELEC_TRIG_LUMI = 36.42;
const float                                                     ELEC_TRIG_SF = 1;
const TString   MU_SUFF = "muon_" + YEAR_STR,                   EL_SUFF = "electron_" + YEAR_STR;

// FIXME
const float     NGEN_ZZ_4L      = 6669988,          XSEC_ZZ_4L      = 1.256;    //1.212;
const float     NGEN_ZJETS_M50  = 80924255,         XSEC_ZJETS_M50  = 6225.42;  //5765.4;
const float     NGEN_ZJETS_M10  = 29374008,         XSEC_ZJETS_M10  = 18610;
const float     NGEN_GGH_ZZ_4L  = 999800,           XSEC_GGH_ZZ_4L  = 0.01212;
const float     NGEN_VBFH_ZZ_4L = 499262,           XSEC_VBFH_ZZ_4L = 0.001034;
const float     NGEN_TTBAR      = 15173839,         XSEC_TTBAR      = 831.76 * 0.5;
const float     NGEN_TT_2L2NU   = 79140880,         XSEC_TT_2L2NU   = 87.31 * 0.5;
const float     NGEN_TTZ_2L2NU  = 2654179,          XSEC_TTZ_2L2NU  = 0.2529;
const float     NGEN_WW_2L2NU   = 1999000,          XSEC_WW_2L2NU   = 12.178;
const float     NGEN_WZ_2L2Q    = 15879472,         XSEC_WZ_2L2Q    = 5.595;
const float     NGEN_WZ_3LNU    = 7387013,          XSEC_WZ_3LNU    = 4.42965;
const float     NGEN_ZZ_2L2Q    = 496436,           XSEC_ZZ_2L2Q    = 3.22;
const float     NGEN_ZZ_2L2NU   = 48655100,         XSEC_ZZ_2L2NU   = 0.564;

const unsigned  N_MC = 11,      ZZ = 0,     DY = 1,     TT = 4;
const TString   MC_SUFF[N_MC] = {   "zz_4l",        "zjets_m-50",   "zjets_m-10",
                                    "ttbar",        "ww_2l2nu",     "wz_2l2q",      "wz_3lnu",
                                    "zz_2l2q",      "zz_2l2nu",     "ggH_zz_4l",    "vbfH_zz_4l"
                                };
const float     NGEN[N_MC] = {  NGEN_ZZ_4L,     NGEN_ZJETS_M50, NGEN_ZJETS_M10,
                                NGEN_TTBAR,     NGEN_WW_2L2NU,  NGEN_WZ_2L2Q,   NGEN_WZ_3LNU,
                                NGEN_ZZ_2L2Q,   NGEN_ZZ_2L2NU,  NGEN_GGH_ZZ_4L, NGEN_VBFH_ZZ_4L
                                };
const float     XSEC[N_MC] = {  XSEC_ZZ_4L,     XSEC_ZJETS_M50, XSEC_ZJETS_M10,
                                XSEC_TTBAR,     XSEC_WW_2L2NU,  XSEC_WZ_2L2Q,   XSEC_WZ_3LNU,
                                XSEC_ZZ_2L2Q,   XSEC_ZZ_2L2NU,  XSEC_GGH_ZZ_4L, XSEC_VBFH_ZZ_4L
                                };
const int       COLOR[N_MC] = { lLightBlue,     lYellow,        kRed,
                                lGreen,         lOrange,        lOrange,        lOrange,
                                lOrange,        lOrange,        lPurple,        lPurple
                                };




#endif // CUTS2016_HH
