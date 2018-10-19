#ifndef LEPTONPAIR_HH
#define LEPTONPAIR_HH


// STL
#include <utility>

// ROOT
#include "TLorentzVector.h"

// Custom
#include "Lepton.hh"


using namespace std;




//////////////////
//    STRUCT    //
//////////////////


struct LeptonPair
{
    ////  FUNCTIONS

    // Constructors
    LeptonPair( Lepton&,  Lepton&);                 // Constructor with leptons
    LeptonPair(){};

    // Members
    bool            SetMembers( Lepton&,  Lepton&); // Set leptons



    ////  DATA MEMBERS

    TLorentzVector  p4;                         // Sum of lab frame 4-momenta
    TLorentzVector  b_p4;                       // Sum of Z CM frame 4-momenta

    int             pdg;                        // PDG ID (absolute value)
    unsigned        mother;                     // Mother Z index

    Lepton          *firstPt,   *secondPt;      // Pointers to leptons ordered by Pt
    Lepton          *firstP,    *secondP;       // Pointers to leptons ordered by P
    Lepton          *plus,      *minus;         // Pointers to leptons ordered by q


    private:
        pair<Lepton*, Lepton*>  members;        // STD pair of pointers to member leptons
};




/////////////////////
//    UTILITIES    //
/////////////////////


////  PAIRING

// Make pairs using mother indices
bool    MakePairsFromMother(vector<Lepton>&,    LeptonPair*,    LeptonPair*);




#endif  // LEPTONPAIR_HH
