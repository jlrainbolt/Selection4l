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
    // Public data members
    TLorentzVector  p4;                                 // Sum of lepton p4
    TLorentzVector  b_p4;                               // Sum of lepton b_p4
    TVector3        b_v3;                               // 3-momentum corresponding to b_p4
    int             pdg;                                // PDG ID (absolute value)
    unsigned        mother;                             // Mother Z index of member leptons


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
    void SetBoostedP4(const TVector3&);                 // Bet b_p4, b_v3 of self and members


    // Protected attributes
    private:        pair<Lepton, Lepton> leptons;       // STD pair of member leptons
};



#endif  // LEPTONPAIR_HH
