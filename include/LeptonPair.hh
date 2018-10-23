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
    TLorentzVector  p4;                                 // Sum of lab frame 4-momenta
    int             pdg;                                // PDG ID (absolute value)
    unsigned        mother;                             // Mother Z index


    // Constructors
    LeptonPair(const Lepton&, const Lepton&);           // Constructor calling SetMembers
    LeptonPair(){};                                     // Default constructor


    // "Getters"
    Lepton          First() const,      Second() const; // Return first, second Pt lepton
    Lepton          Plus() const,       Minus() const;  // Return positive, negative lepton
    vector<Lepton>  GetMembers() const;                 // Return vector of leptons (in pair order)


    // "Setters"
    void    SetMembers(const Lepton&, const Lepton&);   // COPY leptons to self; call Initialize()
    void    SetMothers(unsigned);                       // Set mother of self and member leptons


    // Protected attributes
    private:        pair<Lepton, Lepton> leptons;       // STD pair of member leptons
};



//
//    "FRIENDS"
//



#endif  // LEPTONPAIR_HH
