// STL
#include <iostream>
#include <vector>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"


using namespace std;




//
//  MEMBERS
//

// SetBoostedP4
void Lepton :: SetBoostedP4(const TVector3 &beta)
{
    b_p4 = p4;          // in rest frame of "system" where beta = system.BoostVector()
    b_p4.Boost(-beta);  // so a minus sign is needed here
    b_v3 = b_p4.Vect();
}

// SetMatch
void Lepton :: SetMatch(const Lepton &match)
{
    m_p4 = match.p4;
    dr = p4.DeltaR(m_p4);

    m_b_p4 = match.b_p4;
    m_b_v3 = m_b_p4.Vect();
}



//
//    "FRIENDS"
//

// DecreasingPt
bool DecreasingPt(const Lepton &i, const Lepton &j)
{   
    return (i.p4.Pt() > j.p4.Pt());
}

// DecreasingBoostedPt
bool DecreasingBoostedP(const Lepton &i, const Lepton &j)
{   
    return (i.b_p4.P() > j.b_p4.P());
}


// TotalP4
TLorentzVector TotalP4(const vector<Lepton> &leps)
{
    TLorentzVector p4;
    for (unsigned i = 0; i < leps.size(); i++)
        p4 += leps[i].p4;
    return p4;
}
