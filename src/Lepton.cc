// STL
#include <iostream>
#include <vector>

// ROOT
#include "TLorentzVector.h"

// Custom
#include "Lepton.hh"


using namespace std;




//
//  MEMBERS
//
/*
// SetBoostedP4
bool Lepton :: SetBoostedP4(const TVector3 &beta)
{
    // b_p4 is in rest frame of "system" where beta = system.BoostVector()
    // (which is why a minus sign is needed)
    b_p4 = p4;
    b_p4.Boost(-beta);
    b_p3 = b_p4.Vect();

    return kTRUE;
} 
*/


//
//    "FRIENDS"
//

// DecreasingPt
bool DecreasingPt(const Lepton &i, const Lepton &j)
{   
    return (i.p4.Pt() > j.p4.Pt());
}
/*
// DecreasingBoostedPt
bool DecreasingBoostedP(const Lepton &i, const Lepton &j)
{   
    return (i.b_p4.P() > j.b_p4.P());
}
*/

// TotalP4
TLorentzVector TotalP4(const vector<Lepton> &leps)
{
    TLorentzVector p4;
    for (unsigned i = 0; i < leps.size(); i++)
        p4 += leps[i].p4;
    return p4;
}

/*
// TotalBoostedP4
TLorentzVector  TotalBoostedP4( const vector<Lepton> &leps)
{
    TLorentzVector b_p4sum;
    for (unsigned i = 0; i < leps.size(); i++)
        b_p4sum += leps[i].b_p4;
    return b_p4sum;
}
*/
