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
        lep.q = 0;
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
        lep.q = 0;
        return lep;
    }
}


// BFirst
Lepton LeptonPair :: BFirst() const
{
    if (leptons.first.b_p4.P() > leptons.second.b_p4.P())
        return leptons.first;
    else
        return leptons.second;
}

// Second
Lepton LeptonPair :: BSecond() const
{
    if (leptons.first.b_p4.P() < leptons.second.b_p4.P())
        return leptons.first;
    else
        return leptons.second;
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


    // Momentum & boost
    p4 = leptons.first.p4 + leptons.second.p4;
    b_p4 = leptons.first.b_p4 + leptons.second.b_p4;
    b_v3 = b_p4.Vect();

    m_p4 = leptons.first.m_p4 + leptons.second.m_p4;
    m_b_p4 = leptons.first.m_b_p4 + leptons.second.m_b_p4;
    m_b_v3 = m_b_p4.Vect();


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


// SetBoostedP4
void LeptonPair :: SetBoostedP4(const TVector3& beta)
{
    leptons.first.SetBoostedP4(beta);
    leptons.second.SetBoostedP4(beta);
    b_p4 = leptons.first.b_p4 + leptons.second.b_p4;
    b_v3 = b_p4.Vect();
}
void LeptonPair :: SetBoostedP4(const TVector3& beta, const TVector3& m_beta)
{
    SetBoostedP4(beta);

    // Make intermediate leptons to propagate matched boost
    Lepton first_, second_;
    first_.p4   = leptons.first.m_p4;       first_.SetBoostedP4(m_beta);
    second_.p4  = leptons.second.m_p4;      second_.SetBoostedP4(m_beta);
    leptons.first.SetMatch(first_);         leptons.second.SetMatch(second_);

    // Recalculate matched boost for pair
    m_b_p4 = leptons.first.m_b_p4 + leptons.second.m_b_p4;
    m_b_v3 = m_b_p4.Vect();
}
