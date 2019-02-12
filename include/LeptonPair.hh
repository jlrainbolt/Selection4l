#ifndef LEPTONPAIR_HH
#define LEPTONPAIR_HH

// STL
#include <utility>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"

using namespace std;



//
//    STRUCT
//

struct LeptonPair
{
    // Basics
    TLorentzVector  p4;                                 // Sum of lepton p4
    int             pdg;                                // PDG ID (absolute value)
    unsigned        mother;                             // Mother Z index of member leptons

    // Boosted quantities
    TLorentzVector  b_p4;                               // Sum of lepton b_p4
    TVector3        b_v3;                               // 3-momentum corresponding to b_p4

    // Matches quantities
    TLorentzVector  m_p4;                               // Sum of lepton m_p4
    TLorentzVector  m_b_p4;                             // Sum of lepton m_b_p4
    TVector3        m_b_v3;                             // 3-momentum corresponding to m_b_p4


    // Constructors
    LeptonPair(const Lepton&, const Lepton&);           // Constructor calling SetMembers
    LeptonPair(){};                                     // Default constructor


    // "Getters"
    Lepton          First() const,  Second() const;     // Return first, second Pt lepton (from p4)
    Lepton          Plus() const,   Minus() const;      // Return positive, negative lepton
    Lepton          BFirst() const, BSecond() const;    // Return first, second P lepton (b_p4)
    vector<Lepton>  GetMembers() const;                 // Return vector of leptons (in pair order)


    // "Setters"
    void SetMembers(const Lepton&, const Lepton&);      // Copy leptons to self; call Initialize()
    void SetMothers(unsigned);                          // Set mother of self and member leptons
    void SetBoostedP4(const TVector3&);                 // Set b_p4, b_v3 of self and members
    void SetBoostedP4(const TVector3&, const TVector3&);// Also sets matched boost


    // Other
    bool BlindCharges(const float);                     // Swap member charges (for blinding)


    // Protected attributes
    private:        pair<Lepton, Lepton> leptons;       // STD pair of member leptons
};



#endif  // LEPTONPAIR_HH
