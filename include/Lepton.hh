#ifndef LEPTON_HH
#define LEPTON_HH


// STL
#include <vector>

// ROOT
#include "TLorentzVector.h"


using namespace std;




//////////////////
//    STRUCT    //
//////////////////


struct Lepton
{

    ////  DATA MEMBERS

    TLorentzVector  p4;         // Lab frame 4-momentum
    TLorentzVector  b_p4;       // Z CM frame 4-momentum
    int             q;          // Charge
    int             pdg;        // PDG ID (signed)
    unsigned        mother;     // Mother Z index
};




/////////////////////
//    UTILITIES    //
/////////////////////


////  SORTING

// Sort conditions for STD container of leptons
bool    DecreasingLabPt(        const Lepton&,  const Lepton&);
bool    DecreasingBoostedP(     const Lepton&,  const Lepton&);



////  VECTOR SUMS

// Get TLorentzVector sum of leptons in lab frame
TLorentzVector  LabP4Sum(       const vector<Lepton>&);

// Get TLorentzVector sum of leptons in boosted frame
TLorentzVector  BoostedP4Sum(   const vector<Lepton>&);




#endif  // LEPTON_HH
