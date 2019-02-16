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

const unsigned  RNG_SEED = 16;


//
//  PHYSICAL CONSTANTS
//

const float     MUON_MASS = 0.105658369,                        ELEC_MASS = 0.000511;



//
////  PHASE SPACE
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

const float     MUON_TRIG_LUMI = 35.9,                          ELEC_TRIG_LUMI = 35.9;
const float                                                     ELEC_TRIG_SF = 1;
const TString   MU_SUFF = "muon_" + YEAR_STR,                   EL_SUFF = "electron_" + YEAR_STR;

// Event numbers from DAS, negative fractions from XSDB (FIXME?)
//
// Cross sections from
//      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
//      https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections
const float     NGEN_ZZ_4L      = 6762740,          XSEC_ZZ_4L      = 1.212;
const float     NGEN_ZJETS_M50  = 80924255,         XSEC_ZJETS_M50  = 5765.4;
const float     NGEN_GGH_ZZ_4L  = 999800,           XSEC_GGH_ZZ_4L  = 0.01212;
const float     NGEN_VBFH_ZZ_4L = 499262,           XSEC_VBFH_ZZ_4L = 0.001034;
const float     NGEN_TTBAR      = 12284545,         XSEC_TTBAR      = 831.76;
const float     NGEN_WW_2L2NU   = 1832358,          XSEC_WW_2L2NU   = 12.178;
const float     NGEN_WZ_2L2Q    = 15879472,         XSEC_WZ_2L2Q    = 5.595;
const float     NGEN_WZ_3LNU    = 7387013,          XSEC_WZ_3LNU    = 4.42965;
const float     NGEN_ZZ_2L2Q    = 496436,           XSEC_ZZ_2L2Q    = 3.22;
const float     NGEN_ZZ_2L2NU   = 48623080,         XSEC_ZZ_2L2NU   = 0.564;

const unsigned  N_MC = 10,      ZZ = 0,     DY = 1;
const TString   MC_SUFF[N_MC] = {   "zz_4l",        "zjets_m-50",   "ggH_zz_4l",    "vbfH_zz_4l",
                                    "ttbar",        "ww_2l2nu",     "wz_2l2q",      "wz_3lnu",
                                    "zz_2l2q",      "zz_2l2nu"};
const float     NGEN[N_MC] = {  NGEN_ZZ_4L,     NGEN_ZJETS_M50, NGEN_GGH_ZZ_4L, NGEN_VBFH_ZZ_4L,
                                NGEN_TTBAR,     NGEN_WW_2L2NU,  NGEN_WZ_2L2Q,   NGEN_WZ_3LNU,
                                NGEN_ZZ_2L2Q,   NGEN_ZZ_2L2NU   };
const float     XSEC[N_MC] = {  XSEC_ZZ_4L,     XSEC_ZJETS_M50, XSEC_GGH_ZZ_4L, XSEC_VBFH_ZZ_4L,
                                XSEC_TTBAR,     XSEC_WW_2L2NU,  XSEC_WZ_2L2Q,   XSEC_WZ_3LNU,
                                XSEC_ZZ_2L2Q,   XSEC_ZZ_2L2NU   };
const int       COLOR[N_MC] = { lLightBlue,     lYellow,        lPurple,        lPurple,
                                lGreen,         lOrange,        lOrange,        lOrange,
                                lOrange,        lOrange };




#endif // CUTS2016_HH
