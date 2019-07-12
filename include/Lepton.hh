#ifndef LEPTON_HH
#define LEPTON_HH


// STL
#include <vector>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;



//
//  STRUCT
//

struct Lepton
{
    // Basics
    TLorentzVector      p4;                 // Lab frame 4-momentum
    int                 q;                  // Charge
    int                 pdg;                // PDG ID (signed)
    unsigned            mother;             // Mother Z index
    bool                tight;              // Lepton passes tight HZZ ID

    // Reco quantities
    TLorentzVector      u_p4;               // 4-momentum before scale factor ("uncorrected")
    float               iso;                // Relative isolation
    pair<float, float>  id_sf;              // Scale factors for ID/reco
    pair<bool, bool>    di_hlt;             // Fired trigger leg1/leg2
    bool                si_hlt;             // Fired single-lepton trigger
    float               trig_sf;            // Scale factor for single-lepton trigger

    // Boosted quantities
    TLorentzVector      b_p4;               // boosted 4-momentum
    TVector3            b_v3;               // 3-momentum corresponding to b_p4

    // Matched quantities
    TLorentzVector      m_p4;               // 4-momentum of matched gen lepton
    float               dr;                 // DeltaR between p4 and m_p4
    TLorentzVector      m_b_p4;             // Boosted m_p4
    TVector3            m_b_v3;             // 3-momentum corresponding to m_b_p4


    // Calculate boosted p4 (and p3) using p4 and given boost vector
    void SetBoostedP4(const TVector3&);

    // Calculate dr, m_b_p4, and m_b_v3 given matched lepton
    void SetMatch(const Lepton&);
};



//
//  "FRIENDS"
//

// Sort std container by Pt (in lab frame)
bool DecreasingPt(const Lepton&, const Lepton&);

// Sort std container by P in Z rest ("boosted") frame
bool DecreasingBoostedP(const Lepton&, const Lepton&);



// Get sum of p4 of all leptons in vector
TLorentzVector  TotalP4(const vector<Lepton>&);
/*
// Get sum of b_p4 of all leptons in vector
TLorentzVector  TotalBoostedP4( const vector<Lepton>&);
*/


#endif  // LEPTON_HH
