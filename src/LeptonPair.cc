// STL
#include <iostream>
#include <utility>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "LeptonPair.hh"
#include "Lepton.hh"

using namespace std;



//
//    CONSTRUCTOR
//

LeptonPair :: LeptonPair(const Lepton& lep1, const Lepton& lep2)
{
    SetMembers(lep1, lep2);
}



//
//    "GETTERS"
//


// First
Lepton LeptonPair :: First() const
{
    if (leptons.first.p4.Pt() > leptons.second.p4.Pt())
        return leptons.first;
    else
        return leptons.second;
}

// Second
Lepton LeptonPair :: Second() const
{
    if (leptons.first.p4.Pt() < leptons.second.p4.Pt())
        return leptons.first;
    else
        return leptons.second;
}

// Plus
Lepton LeptonPair :: Plus() const
{
    if      (leptons.first.q > 0 && leptons.second.q < 0)
        return leptons.first;

    else if (leptons.second.q > 0 && leptons.first.q < 0)
        return leptons.second;

    else    // return an empty lepton :(
    {
        Lepton lep;
        return lep;
    }
}

// Minus
Lepton LeptonPair :: Minus() const
{
    if      (leptons.first.q < 0 && leptons.second.q > 0)
        return leptons.first;

    else if (leptons.second.q < 0 && leptons.first.q > 0)
        return leptons.second;

    else    // return an empty lepton :(
    {
        Lepton lep;
        return lep;
    }
}


// GetMembers
vector<Lepton> LeptonPair :: GetMembers() const
{
    vector<Lepton> vect = {leptons.first, leptons.second};
    return vect;   
}



//
//  "SETTERS"
//

//  SetMembers
void LeptonPair :: SetMembers(const Lepton& lep1, const Lepton& lep2)
{
    leptons = make_pair(lep1, lep2);


    // Momentum
    p4 = leptons.first.p4 + leptons.second.p4;


    // PDG ID (if it's a match)
    if (abs(leptons.first.pdg) == abs(leptons.second.pdg))
        pdg = abs(leptons.first.pdg);
    else
        pdg = 0;


    // Mother (if it's a match)
    if (leptons.first.mother == leptons.second.mother)
        mother = leptons.first.mother;
    else
        mother = 0;
}



//  SetMothers
void LeptonPair :: SetMothers(unsigned mom)
{
    mother                  = mom;
    leptons.first.mother    = mom;
    leptons.second.mother   = mom;

    return;
}



//
//  "FRIENDS"
//
