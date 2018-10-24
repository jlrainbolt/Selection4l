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

const float     MATCH_DR_MAX = 0.15;



//
//  SAMPLES
//

const float     MUON_TRIG_LUMI = 36.735,                        ELEC_TRIG_LUMI = 41.529;
const float     NGEN_ZZ_4L = 6967853 * (1 - 2 * 0.005244),      XSEC_ZZ_4L = 1.212;



#endif // CUTS2017_HH
