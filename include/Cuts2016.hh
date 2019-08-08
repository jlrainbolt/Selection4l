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

const unsigned  RNG_SEED = 2016;


//
//  PHYSICAL CONSTANTS
//

const float     MUON_MASS = 0.105658369,                        ELEC_MASS = 0.000511;



//
//  SAMPLES
//

const float     INT_LUMI = 35.921,                  SQRT_S = 13;
const float     ELEC_TRIG_SF = 1;
const TString   MU_SUFF = "muon_" + YEAR_STR,       EL_SUFF = "electron_" + YEAR_STR;


const float     NGEN_ZZ_4L = 6669988,               XSEC_ZZ_4L = 1.256;



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
//  TRIGGERS
//

const float     MUON_SINGLE_PT = 24,    ELEC_SINGLE_PT = 27;
const float     MUON_LEG1_PT = 17,      ELEC_LEG1_PT = 23;
const float     MUON_LEG2_PT = 8,       ELEC_LEG2_PT = 12;
const float     ELEC_ETA_PREF = 2.1;

const unsigned  N_TRIG = 2;                             // number of trigger sf files
const float     TRIG_FRAC = 19.695 / INT_LUMI;          // fraction of events for first file
const TString   TRIG_NAME = "IsoMu24_OR_IsoTkMu24";
const TString   TRIG_PD[2] = {"BCDEF", "GH"};



//
//  LEPTON ID
//

const float     MUON_D0_MAX = 0.5,      MUON_DZ_MAX = 1,        MUON_SIP_MAX = 4;
const int       MUON_BAD_TRACK = 2;



//
//  OBJECT MATCHING
//

const float     MATCH_DR_MAX = 0.02;


#endif // CUTS2016_HH
