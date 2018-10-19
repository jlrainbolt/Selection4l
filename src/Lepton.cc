// STL
#include <iostream>
#include <vector>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"


using namespace std;


/////////////////////
//    UTILITIES    //
/////////////////////


////  SORTING

// Sort std container by Pt in lab frame
bool    DecreasingLabPt(    const Lepton &i,    const Lepton &j)
{   
    return (i.p4.Pt() > j.p4.Pt());
}


// Sort std container by P in Z rest ("boosted") frame
bool    DecreasingBoostedP( const Lepton &i,    const Lepton &j)
{   
    return (i.b_p4.P() > j.b_p4.P());
}



////  VECTOR SUMS

// Get TLorentzVector sum of leptons in lab frame
TLorentzVector  LabP4Sum(       const vector<Lepton> &leps)
{
    TLorentzVector p4sum;
    for (unsigned i = 0; i < leps.size(); i++)
        p4sum += leps[i].p4;
    return p4sum;
}


// Get TLorentzVector sum of leptons in boosted frame
TLorentzVector  BoostedP4Sum(   const vector<Lepton> &leps)
{
    TLorentzVector b_p4sum;
    for (unsigned i = 0; i < leps.size(); i++)
        b_p4sum += leps[i].b_p4;
    return b_p4sum;
}
