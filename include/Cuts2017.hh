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
const TString   BLT_PATH = HOME_PATH + "/lacey_legacy/CMSSW_10_2_13/src/BLT";

const unsigned  RNG_SEED = 7102;


//
//  PHYSICAL CONSTANTS
//

const float     MUON_MASS = 0.105658369,                        ELEC_MASS = 0.000511;



//
//  PHASE SPACE
//

//const float     M_MIN = 80,             M_MAX = 100,            MLL_MIN = 4;
const float     M_MIN = 80,             M_MAX = 100,            MLL_MIN = 1;
//const float     M_MIN = 50,             M_MAX = TMath::Infinity(),      MLL_MIN = 4;



//
//  FIDUCIAL REGION
//

const float     FID_PT1_MIN = 20,       FID_PT2_MIN = 10,       FID_PT_MIN = 5;
const float     FID_ETA_MAX = 2.5;



//
//  RECO SELECTION
//

const float     Z1_M_MIN = 1,          Z_M_MAX = 120;
//const float     Z1_M_MIN = 12,          Z_M_MAX = 120;
const float     MUON_PT1_MIN = 20,      MUON_PT2_MIN = 10,      MUON_PT_MIN = 5;
const float     ELEC_PT1_MIN = 25,      ELEC_PT2_MIN = 15,      ELEC_PT_MIN = 7;
const float     MUON_ETA_MAX = 2.4,     ELEC_ETA_MAX = 2.5;
const float     MUON_ISO_MAX = 0.35,    ELEC_ISO_MAX = 0.35;
const float     SF_DR_MIN = 0.02,       OF_DR_MIN = 0.05;



//
//  TRIGGERS
//

const float     MUON_SINGLE_PT = 27,    ELEC_SINGLE_PT = 35;
const float     MUON_LEG1_PT = 17,      ELEC_LEG1_PT = 23;
const float     MUON_LEG2_PT = 8,       ELEC_LEG2_PT = 12;
const float     ELEC_ETA_PREF = 2.1;

const unsigned  N_TRIG = 1;             // number of trigger sf files
const float     TRIG_FRAC = 1;          // fraction of events for first file
const TString   TRIG_NAME = "IsoMu27";
const TString   TRIG_PD[2] = {"BCDEF", ""};



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

const float     INT_LUMI = 41.37;
const float                                         ELEC_TRIG_SF = 0.991;
const TString   MU_SUFF = "muon_" + YEAR_STR,       EL_SUFF = "electron_" + YEAR_STR;

// Event numbers from TotalEvents histogram in "selected" ntuples (bin1 - 2*bin10)
//
// Cross sections from
//      https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
const float     NGEN_ZZ_4L      = 6893887,          XSEC_ZZ_4L      = 1.256;    //1.212
const float     NGEN_ZJETS_M50  = 123485957,        XSEC_ZJETS_M50  = 6225.43;  //5765.4
const float     NGEN_GGH_ZZ_4L  = 955384,           XSEC_GGH_ZZ_4L  = 0.01212;
const float     NGEN_VBFH_ZZ_4L = 984662,           XSEC_VBFH_ZZ_4L = 0.001034;
const float     NGEN_TTBAR      = 57584555,         XSEC_TTBAR      = 831.76 * 0.5;
const float     NGEN_TT_2L2NU   = 8926992,          XSEC_TT_2L2NU   = 87.31 * 0.5;
const float     NGEN_TTZ_2L2NU  = 3570720,          XSEC_TTZ_2L2NU  = 0.2529;
const float     NGEN_WW_2L2NU   = 1992522,          XSEC_WW_2L2NU   = 12.178;
const float     NGEN_WZ_2L2Q    = 16664610,         XSEC_WZ_2L2Q    = 5.595;
const float     NGEN_WZ_3LNU    = 6887413,          XSEC_WZ_3LNU    = 4.42965;
const float     NGEN_ZZ_2L2NU   = 8733658,          XSEC_ZZ_2L2NU   = 0.564;
const float     NGEN_ZZ_2L2Q    = 17768294,         XSEC_ZZ_2L2Q    = 3.22;
const float     NGEN_WWZ_4L2NU  = 1707572,          XSEC_WWZ_4L2NU  = 0.1651;
const float     NGEN_WZZ_4L2NU  = 1690058,          XSEC_WZZ_4L2NU  = 0.05565;
const float     NGEN_ZZZ_4L2NU  = 1673322,          XSEC_ZZZ_4L2NU  = 0.01398;
const float     NGEN_ZZG_4L2NU  = 1455362,          XSEC_ZZG_4L2NU  = 0.01398;

const unsigned  N_MC = 16,      ZZ = 0,     DY = 1;
const unsigned  N_DY = 15;

const TString   MC_SUFF[N_MC] = {   "zz_4l",        "zjets_m-50",   "ggH_zz_4l",    "vbfH_zz_4l",
                                    "ttbar",        "tt_2l2nu",     "ttz_2l2nu",    "ww_2l2nu",
                                    "wz_2l2q",      "wz_3lnu",      "zz_2l2nu",     "zz_2l2q",
                                    "wwz_4l2nu",    "wzz_4l2nu",    "zzz_4l2nu",    "zzg_4l2nu"
                                };
const float     NGEN[N_MC] = {  NGEN_ZZ_4L,     NGEN_ZJETS_M50, NGEN_GGH_ZZ_4L, NGEN_VBFH_ZZ_4L,
                                NGEN_TTBAR,     NGEN_TT_2L2NU,  NGEN_TTZ_2L2NU, NGEN_WW_2L2NU,
                                NGEN_WZ_2L2Q,   NGEN_WZ_3LNU,   NGEN_ZZ_2L2NU,  NGEN_ZZ_2L2Q,
                                NGEN_WWZ_4L2NU, NGEN_WZZ_4L2NU, NGEN_ZZZ_4L2NU, NGEN_ZZG_4L2NU
                                };
const float     XSEC[N_MC] = {  XSEC_ZZ_4L,     XSEC_ZJETS_M50, XSEC_GGH_ZZ_4L, XSEC_VBFH_ZZ_4L,
                                XSEC_TTBAR,     XSEC_TT_2L2NU,  XSEC_TTZ_2L2NU, XSEC_WW_2L2NU,
                                XSEC_WZ_2L2Q,   XSEC_WZ_3LNU,   XSEC_ZZ_2L2NU,  XSEC_ZZ_2L2Q,
                                XSEC_WWZ_4L2NU, XSEC_WZZ_4L2NU, XSEC_ZZZ_4L2NU, XSEC_ZZG_4L2NU
                                };
const int       COLOR[N_MC] = { lLightBlue,     lYellow,        lPurple,        lPurple,
                                lGreen,         lGreen,         lGreen,         lOrange,
                                lOrange,        lOrange,        lOrange,        lOrange,
                                lBlue,          lBlue,          lBlue,          lBlue
                                };

#endif // CUTS2017_HH
