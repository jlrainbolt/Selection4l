// STL
#include <vector>
#include <iostream>

// ROOT
#include "TLorentzVector.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"
#include "SelectionTools.hh"

using namespace std;






////
////
////    PAIRING
////
////



//
//  FROM MOTHER
//

bool MakePairsFromMother(const vector<Lepton> &leps, LeptonPair *z1, LeptonPair *z2)
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
    if ((pair1leps.size() != 2) || (pair2leps.size() != 2))
    {
        cout << "Wrong number of daughters" << endl;
        return kFALSE;
    }


    // Sort pairs by mass
    TLorentzVector pair1p4 = TotalP4(pair1leps),    pair2p4 = TotalP4(pair2leps);

    if (pair1p4.M() > pair2p4.M())
    {
        z1->SetMembers(pair1leps[0], pair1leps[1]);
        z2->SetMembers(pair2leps[0], pair2leps[1]);
    }
    else
    {
        z1->SetMembers(pair2leps[0], pair2leps[1]);
        z2->SetMembers(pair1leps[0], pair1leps[1]);
    }


    // Set mothers (can always be changed outside function)
    z1->SetMothers(1);
    z2->SetMothers(2);


    return kTRUE;
}



//
//  MAXIMUM MASS DIFFERENCE
//

bool MakePairsMaxDiff(const vector<Lepton> &leps, LeptonPair *z1, LeptonPair *z2)
{
    if (leps.size() != 4)
        return kFALSE;


    // Configuration 1: opposite-sign leptons with largest momenta are paired

    // Assign first lepton (assumed highest Pt) as lepA1
    Lepton lepA1 = leps[0];

    // Pair with lepA2, subleading lepton of zA
    unsigned j;     // index of lepA2
    Lepton lepA2;
    for (unsigned i = 1; i < leps.size(); i++)
    {
        if (leps[i].q != lepA1.q)
        {
            j = i;
            lepA2 = leps[j];
            break;
        }
    }
    LeptonPair c1_zA(lepA1, lepA2);
    if ((c1_zA.Plus().q != 1) || (c1_zA.Minus().q != -1))
        return kFALSE;

    // Pair remaining leptons
    LeptonPair c1_zB;

    if      (j == 1)
        c1_zB.SetMembers(leps[2], leps[3]);
    else if (j == 2)
        c1_zB.SetMembers(leps[1], leps[3]);
    else if (j == 3)
        c1_zB.SetMembers(leps[1], leps[2]);

    if ((c1_zB.Plus().q != 1) || (c1_zB.Minus().q != -1))
        return kFALSE;


    // Configuration 2: swap negatively-charged leptons

    LeptonPair c2_zA(c1_zA.Plus(), c1_zB.Minus());
    if ((c2_zA.Plus().q != 1) || (c2_zA.Minus().q != -1))
        return kFALSE;

    LeptonPair c2_zB(c1_zB.Plus(), c1_zA.Minus());
    if ((c2_zB.Plus().q != 1) || (c2_zB.Minus().q != -1))
        return kFALSE;


    // Compare mass differences and assign selected configuration

    float diff_c1 = c1_zA.p4.M() - c1_zB.p4.M(),        diff_c2 = c2_zA.p4.M() - c2_zB.p4.M();

    if (fabs(diff_c1) > fabs(diff_c2))
    {
        if (diff_c1 > 0)
        {
            z1->SetMembers(c1_zA.First(), c1_zA.Second());
            z2->SetMembers(c1_zB.First(), c1_zB.Second());
        }
        else
        {
            z1->SetMembers(c1_zB.First(), c1_zB.Second());
            z2->SetMembers(c1_zA.First(), c1_zA.Second());
        }
    }
    else
    {
        if (diff_c2 > 0)
        {
            z1->SetMembers(c2_zA.First(), c2_zA.Second());
            z2->SetMembers(c2_zB.First(), c2_zB.Second());
        }
        else
        {
            z1->SetMembers(c2_zB.First(), c2_zB.Second());
            z2->SetMembers(c2_zA.First(), c2_zA.Second());
        }
    }


    // Set mothers (can always be changed outside function)
    z1->SetMothers(1);
    z2->SetMothers(2);


    return kTRUE;
}



////
////
////    TRIGGER SF
////
////


//
//  SINGLE LEPTON
//

float GetSingleTriggerSF(const Lepton& lep1, const Lepton& lep2)
{
    float sf = 1;

    // Only lepton 1 triggered
    if      (lep1.fired.first && !lep2.fired.first)
        sf = lep1.te_data.first / lep1.te_mc.first;

    // Only lepton 2 triggered
    else if (lep2.fired.first && !lep1.fired.first)
        sf = lep2.te_data.first / lep2.te_mc.first;

    // Both leptons triggered
    else if (lep1.fired.first && lep2.fired.first)
    {
        sf = 1. - (1. - lep1.te_data.first) * (1. - lep2.te_data.first);
        sf /= 1. - (1. - lep1.te_mc.first) * (1. - lep2.te_mc.first);
    }

    return sf;
}
