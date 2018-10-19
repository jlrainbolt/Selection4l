// STL
#include <iostream>
#include <utility>

// ROOT
#include "TLorentzVector.h"

// Custom
#include "LeptonPair.hh"
#include "Lepton.hh"


using namespace std;




/////////////////////
//    FUNCTIONS    //
/////////////////////


////  CONSTRUCTORS

// Constructor with lepton member initialization
LeptonPair::LeptonPair(Lepton& lep1, Lepton& lep2)
{
    this->SetMembers(lep1, lep2);
}



////  MEMBERS

// Initialize member leptons
bool LeptonPair::SetMembers(Lepton& lep1, Lepton& lep2)
{
    bool isGoodMatch = kTRUE;

    members = make_pair(&lep1, &lep2);


    // Add momenta
    p4      = lep1.p4   + lep2.p4;
    b_p4    = lep1.b_p4 + lep2.b_p4;


    // Set PDG ID (if it's a match)
    if (abs(lep1.pdg) == abs(lep2.pdg))
        pdg = abs(lep1.pdg);
    else
    {
        cout << "Flavor mismatch" << endl;
        isGoodMatch = kFALSE;
    }


    // Set mother (if it's a match)
    if (lep1.mother == lep2.mother)
        mother = lep1.mother;
    else
    {
        cout << "Mother mismatch" << endl;
        isGoodMatch = kFALSE;
    }


    // Set order pointers
    if (lep1.p4.Pt() > lep2.p4.Pt())
    {
        firstPt     = members.first;
        secondPt    = members.second;
    }
    else
    {
        firstPt     = members.second;
        secondPt    = members.first;
    }

    if (lep1.b_p4.P() > lep2.b_p4.P())
    {
        firstP      = members.first;
        secondP     = members.second;
    }
    else
    {
        firstP      = members.second;
        secondP     = members.first;
    }

    if      (lep1.q > lep2.q)
    {
        plus        = members.first;
        minus       = members.second;
    }
    else if (lep1.q < lep2.q)
    {
        plus        = members.second;
        minus       = members.first;
    }
    else
    {
        cout << "Charge mismatch" << endl;
        isGoodMatch = kFALSE;
    }


    return isGoodMatch;
}



/////////////////////
//    UTILITIES    //
/////////////////////


////  PAIRING

// Sort input vector of four leptons into two std pairs, returned as pointers
// Selection is determined by matching mother Z index
bool    MakePairsFromMother(vector<Lepton> &leps, LeptonPair *z1, LeptonPair*z2)
{

    // Get unique mother Z indices
    vector<unsigned> mothers;

    for (unsigned i = 0; i < leps.size(); i++)
    {
        if (find(mothers.begin(), mothers.end(), leps[i].mother) == mothers.end())
            mothers.push_back(leps[i].mother);
    }

    if (mothers.size() != 2)
    {
        cout << "Wrong number of mothers" << endl;
        return kFALSE;
    }


    // Sort leptons by mother
    vector<Lepton> pair1leps, pair2leps;

    for (unsigned i = 0; i < leps.size(); i++)
    {
        if      (leps[i].mother == mothers[0])
            pair1leps.push_back(leps[i]);

        else if (leps[i].mother == mothers[1])
            pair2leps.push_back(leps[i]);
    }

    if (pair1leps.size() != 2 || pair2leps.size() != 2)
    {
        cout << "Wrong number of daughters" << endl;
        return kFALSE;
    }


    // Sort pairs by mass
    LeptonPair  pair1(pair1leps[0], pair1leps[1]);
    LeptonPair  pair2(pair2leps[0], pair2leps[1]);

    if (pair1.p4.M() > pair2.p4.M())
    {
        *z1 = pair1;
        *z2 = pair2;
    }
    else
    {
        *z1 = pair2;
        *z2 = pair1;
    }


    return kTRUE;
}

