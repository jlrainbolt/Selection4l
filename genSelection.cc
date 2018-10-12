#include <vector>
#include <iostream>
#include <utility>
#include <tuple>
#include <algorithm>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

using namespace std;




//--- HELPERS ---//

bool            SortDecPt(  const pair<TLorentzVector, Short_t> &i_,
                            const pair<TLorentzVector, Short_t> &j_ );

bool            SortDecP(   const pair<TLorentzVector, Short_t> &i_,
                            const pair<TLorentzVector, Short_t> &j_ );

TLorentzVector  GetBoosted( const TLorentzVector    &p4_,
                            const TVector3          &beta   );

pair<TLorentzVector, Short_t>  GetBoosted( const pair<TLorentzVector, Short_t> &lep,
                            const TVector3          &beta   );

TLorentzVector GetP4Sum(const vector<pair<TLorentzVector, Short_t>> &leps);
TLorentzVector GetP4Sum(const pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> &leps);

int GetMassPairs(     vector<pair<TLorentzVector, Short_t>> leps_list,
    pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>>  *large_pair,
    pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>>  *small_pair );




void genSelection(const TString suffix)
{



    ///////////////////////
    //      OPTIONS      //
    ///////////////////////


    //--- OUTPUT ---//

    // File
    TString output  = "trees_gen_" + suffix + ".root";
    TFile *outFile  = new TFile(output, "RECREATE");


    // Selection
    const unsigned N = 6;
    unsigned                MM = 0, EE = 1, M4 = 2, ME = 3, EM = 4, E4 = 5; // Indices
    TString selection[N] = {"mumu", "ee",   "4m",   "2m2e", "2e2m", "4e"};


    // Trees
    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_gen_" + suffix, suffix);


    // Initialize branch variables
    UInt_t          evtNum;
    Float_t         weight;

    TLorentzVector  z1p4, z2p4, zzp4, l1p4, l2p4, l3p4, l4p4;
    Short_t         z1pdg, z2pdg, l1pdg, l2pdg, l3pdg, l4pdg;

    TLorentzVector  b_z1p4, b_z2p4, b_ttp4;
    TLorentzVector  b_l1p4, b_l2p4, b_l3p4, b_l4p4;
    Short_t         b_l1pdg, b_l2pdg, b_l3pdg, b_l4pdg;

    Float_t         b_theta, b_phi, b_z1alpha, b_z2alpha;
    Float_t         bb_z1theta, bb_z2theta;


    // Branches
    for (unsigned i = 0; i < N; i++)
    {
        // Event info
        tree[i]->Branch("evtNum",   &evtNum);       tree[i]->Branch("weight",   &weight);


        // Pair/group momenta
        tree[i]->Branch("zzp4",     &zzp4);
        tree[i]->Branch("z1p4",     &z1p4);         tree[i]->Branch("z2p4",     &z2p4);


        // Lepton momenta, id, iso
        tree[i]->Branch("l1p4",     &l1p4);         tree[i]->Branch("l1pdg",    &l1pdg);
        tree[i]->Branch("l2p4",     &l2p4);         tree[i]->Branch("l2pdg",    &l2pdg);
        tree[i]->Branch("l3p4",     &l3p4);         tree[i]->Branch("l3pdg",    &l3pdg);
        tree[i]->Branch("l4p4",     &l4p4);         tree[i]->Branch("l4pdg",    &l4pdg);


        // Boosted quantities for differential distributions
        tree[i]->Branch("b_z1p4",   &b_z1p4);       tree[i]->Branch("b_z2p4",   &b_z2p4);
        tree[i]->Branch("b_ttp4",   &b_ttp4);

        tree[i]->Branch("b_l1p4",   &b_l1p4);       tree[i]->Branch("b_l1pdg",  &b_l1pdg);
        tree[i]->Branch("b_l2p4",   &b_l2p4);       tree[i]->Branch("b_l2pdg",  &b_l2pdg);
        tree[i]->Branch("b_l3p4",   &b_l3p4);       tree[i]->Branch("b_l3pdg",  &b_l3pdg);
        tree[i]->Branch("b_l4p4",   &b_l4p4);       tree[i]->Branch("b_l4pdg",  &b_l4pdg);

        tree[i]->Branch("b_theta",  &b_theta);      tree[i]->Branch("b_phi",    &b_phi);
        tree[i]->Branch("b_z1alpha", &b_z1alpha);   tree[i]->Branch("b_z2alpha", &b_z2alpha);
        tree[i]->Branch("bb_z1theta", &bb_z1theta); tree[i]->Branch("bb_z2theta", &bb_z2theta);
    }




    //////////////////////
    //    INPUT FILE    //
    //////////////////////


    // Path to directory
    TString path = "";      // guess it's gotta be in this directory


    // Input ROOT file
    TString input = "trees_" + suffix + ".root";
    cout << input << endl;
    TFile *file = TFile::Open(input);



    /////////////////////
    //    SELECTION    //
    /////////////////////

    for (unsigned i = 0; i < N; i++)
    {
        TTreeReader reader(selection[i] + "_" + suffix, file);

        // Event
        TTreeReaderValue<UInt_t>            evtNum_(reader,             "evtNum");
        TTreeReaderValue<Float_t>           weight_(reader,             "weight");          //"genWeight");

        TTreeReaderValue<UShort_t>          nMuons_(reader,             "nGenMuons");
        TTreeReaderValue<UShort_t>          nElecs_(reader,             "nGenElectrons");


        // Leptons
        TTreeReaderArray<TLorentzVector>    muonP4_(reader,             "genMuonP4");
        TTreeReaderValue<vector<Short_t>>   muonQ_(reader,              "genMuonQ");

        TTreeReaderArray<TLorentzVector>    elecP4_(reader,             "genElectronP4");
        TTreeReaderValue<vector<Short_t>>   elecQ_(reader,              "genElectronQ");




        //////////////////////
        //    EVENT LOOP    //
        //////////////////////


        unsigned event = 0;
        while (reader.Next() && event < 100)
        {
//          event++;



            //--- INITIALIZE ---//

            weight  = *weight_;
            UShort_t    nMuons  = *nMuons_,         nElecs  = *nElecs_;

            vector<pair<TLorentzVector, Short_t>> muons, elecs;



            //--- FILL CONTAINERS ---//

            unsigned i_m = 0, i_e = 0;

            // Muons
            for (const TLorentzVector& muonP4: muonP4_)
            {
                muons.push_back(make_pair(muonP4, (*muonQ_)[i_m] * 13));
                i_m++;
            }


            // Electrons
            for (const TLorentzVector& elecP4: elecP4_)
            {
                elecs.push_back(make_pair(elecP4, (*elecQ_)[i_e] * 11));
                i_e++;
            }



            //--- CATEGORIZE ---//

            vector<pair<TLorentzVector, Short_t>>   allLeps;
            pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> z1pair, z2pair;
            Short_t lPDG, tPDG;

            if      (nElecs == 0 && nMuons == 4)
            {
                z1pdg   = 13;           z2pdg   = 13;
                allLeps = muons;
                GetMassPairs(allLeps, &z1pair, &z2pair);
            }
            else if (nElecs == 2 && nMuons == 2)
            {
                sort(muons.begin(),     muons.end(),    SortDecPt);
                sort(elecs.begin(),     elecs.end(),    SortDecPt);
                TLorentzVector  muon_pair = GetP4Sum(muons), elec_pair = GetP4Sum(elecs);

                allLeps = muons;            allLeps.insert(allLeps.end(), elecs.begin(), elecs.end());

                if (muon_pair.M() > elec_pair.M())
                {
                    z1pdg   = 13;           z2pdg   = 13;
                    z1pair  = make_pair(muons[0], muons[1]);
                    z2pair  = make_pair(elecs[0], elecs[1]);
                }
                else
                {
                    z1pdg   = 11;           z2pdg   = 11;
                    z1pair  = make_pair(elecs[0], elecs[1]);
                    z2pair  = make_pair(muons[0], muons[1]);
                }
            }
            else if (nElecs == 4 && nMuons == 0)
            {
                z1pdg   = 11;           z2pdg   = 11;
                allLeps = elecs;
                GetMassPairs(allLeps, &z1pair, &z2pair);
            }
            else
                continue;




            /////////////////////
            //    FILL TREE    //
            /////////////////////


            if (allLeps.size() != 4)
                continue;

            sort(allLeps.begin(), allLeps.end(), SortDecPt);



            //--- CONTAINERS ---//


            evtNum  = *evtNum_;                     weight  = *weight_;
            z1p4    = GetP4Sum(z1pair);             z2p4    = GetP4Sum(z2pair); 
            zzp4    = z1p4 + z2p4;
            tie(l1p4, l1pdg) = allLeps[0];          tie(l2p4, l2pdg) = allLeps[1];
            tie(l3p4, l3pdg) = allLeps[2];          tie(l4p4, l4pdg) = allLeps[3];



            //--- BOOST ---//

            // Sort leptons by magnitude of momentum (P) in Z COM frame

            TVector3 zzb3 = zzp4.BoostVector();

            vector<pair<TLorentzVector, Short_t>> b_allLeps;
            for (unsigned i = 0; i < allLeps.size(); i++)
                b_allLeps.push_back(GetBoosted(allLeps[i], zzb3));
            sort(b_allLeps.begin(), b_allLeps.end(), SortDecP);


            // Sort lepton pairs by P

            pair<TLorentzVector, Short_t> b_z1l1, b_z1l2, b_z2l1, b_z2l2;
            b_z1l1 = GetBoosted(z1pair.first, zzb3);
            b_z1l2 = GetBoosted(z1pair.second, zzb3);
            if (b_z1l1.first.P() < b_z1l2.first.P())
                swap(b_z1l1, b_z1l2);

            b_z2l1 = GetBoosted(z2pair.first, zzb3);
            b_z2l2 = GetBoosted(z2pair.second, zzb3);
            if (b_z2l1.first.P() < b_z2l2.first.P())
                swap(b_z2l1, b_z2l2);

            pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> b_z1pair, b_z2pair;
            b_z1pair = make_pair(b_z1l1, b_z1l2);
            b_z2pair = make_pair(b_z2l1, b_z2l2);


            // Fill boosted containers

            tie(b_l1p4, b_l1pdg) = b_allLeps[0];    tie(b_l2p4, b_l2pdg) = b_allLeps[1];
            tie(b_l3p4, b_l3pdg) = b_allLeps[2];    tie(b_l4p4, b_l4pdg) = b_allLeps[3];

            TLorentzVector b_z1l1p4 = b_z1l1.first, b_z1l2p4 = b_z1l2.first;
            TLorentzVector b_z2l1p4 = b_z2l1.first, b_z2l2p4 = b_z2l2.first;


            // Pairs
            b_z1p4  = GetBoosted(z1p4, zzb3);               b_z2p4  = GetBoosted(z2p4, zzb3);

            // "Trailing trio"
            b_ttp4  = b_l2p4 + b_l3p4 + b_l4p4;



            //--- ANGLES ---//

            // "phi":   angle between pairs
            TLorentzVector  b_z1lpp4 = b_z1l1.second > 0    ? b_z1l1.first  : b_z1l2.first;
            TLorentzVector  b_z1lmp4 = b_z1l1.second < 0    ? b_z1l1.first  : b_z1l2.first;
            TLorentzVector  b_z2lpp4 = b_z2l1.second > 0    ? b_z2l1.first  : b_z2l2.first;
            TLorentzVector  b_z2lmp4 = b_z2l1.second < 0    ? b_z2l1.first  : b_z2l2.first;

            // (normal to z1,z2 decay plane)
            TVector3        b_z1n   = b_z1lpp4.Vect().Cross(b_z1lmp4.Vect());
            TVector3        b_z2n   = b_z2lpp4.Vect().Cross(b_z2lmp4.Vect());
            b_phi   = b_z1n.Angle(b_z2n);

            // "theta": angle between trailing pair 1 lepton and low-mass pair
            b_theta = b_z1l2p4.Angle(b_z2p4.Vect());

            // "alpha": angle between paired leptons
            b_z1alpha = b_z1l1p4.Angle(b_z1l2p4.Vect());
            b_z2alpha = b_z2l1p4.Angle(b_z2l2p4.Vect());



            // OTHERS
            TVector3        b_z1b3    = b_z1p4.BoostVector(),       b_z2b3  = b_z2p4.BoostVector();

            TLorentzVector  bb_z1p4 = GetBoosted(b_z1p4, b_z2b3),   bb_z2p4 = GetBoosted(b_z2p4, b_z1b3);
            TLorentzVector  bb_z1lpp4 = GetBoosted(b_z1lpp4, b_z1b3);
            TLorentzVector  bb_z2lpp4 = GetBoosted(b_z2lpp4, b_z2b3);

            bb_z1theta  = bb_z1lpp4.Angle(bb_z2p4.Vect());
            bb_z2theta  = bb_z2lpp4.Angle(bb_z1p4.Vect());




            tree[i]->Fill();
//          event++;

        }   // END event loop
    }   // END selection loop



    //--- HISTOGRAMS ---//

    TH1D *hTotalEvents, *hPhaseSpaceEvents, *hFiducialEvents, *hSelectedEvents;

    file->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);
    hTotalEvents->SetName("TotalEvents_gen_" + suffix);

    file->GetObject("PhaseSpaceEvents_" + suffix, hPhaseSpaceEvents);
    hPhaseSpaceEvents->SetDirectory(outFile);
    hPhaseSpaceEvents->SetName("PhaseSpaceEvents_gen_" + suffix);

    file->GetObject("FiducialEvents_" + suffix, hFiducialEvents);
    hFiducialEvents->SetDirectory(outFile);
    hFiducialEvents->SetName("FiducialEvents_gen_" + suffix);

    file->GetObject("SelectedEvents_" + suffix, hSelectedEvents);
    hSelectedEvents->SetDirectory(outFile);
    hSelectedEvents->SetName("SelectedEvents_gen_" + suffix);



    //--- CLEANUP ---//

    // Write to file
    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();
    hTotalEvents->Write();
    hSelectedEvents->Write();
    hFiducialEvents->Write();
    hPhaseSpaceEvents->Write();

    outFile->Purge();
    outFile->Close();
    file->Close();



//  delete outFile;
//  delete tree;
//  delete hSelectedEvents;
}




////////////////////////////
//    HELPER FUNCTIONS    //
////////////////////////////


//--- SORT_DEC_? ---//

// Sort std container by TLorentzVector quantity

bool SortDecPt( const pair<TLorentzVector, Short_t> &i_,
                const pair<TLorentzVector, Short_t> &j_)
{
    TLorentzVector i = i_.first, j = j_.first;
    return (i.Pt() > j.Pt());
}


bool SortDecP ( const pair<TLorentzVector, Short_t> &i_,
                const pair<TLorentzVector, Short_t> &j_)
{
    TLorentzVector i = i_.first, j = j_.first;
    return (i.P() > j.P());
}



//--- GET_BOOSTED ---//

// Returns TLorentzVector boosted in frame of "system"
// where beta = system.BoostVector()

TLorentzVector GetBoosted(const TLorentzVector &p4_, const TVector3 &beta)
{
    TLorentzVector p4(p4_);
    p4.Boost(-beta);
    return p4;
}

pair<TLorentzVector,Short_t> GetBoosted(const pair<TLorentzVector, Short_t> &lep,
                                        const TVector3 &beta)
{
    TLorentzVector p4 = GetBoosted(lep.first, beta);
    return make_pair(p4, lep.second);
}



//--- GET_P4_SUM ---//

// Get sum of momentum from "leptons" stored in container

TLorentzVector GetP4Sum(const vector<pair<TLorentzVector, Short_t>> &leps)
{
    TLorentzVector p4sum;
    for (unsigned i = 0; i < leps.size(); i++)
        p4sum += leps[i].first;
    return p4sum;
}


TLorentzVector GetP4Sum(const pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> &leps)
{
    vector<pair<TLorentzVector, Short_t>> leps_vec = {leps.first, leps.second};
    return GetP4Sum(leps_vec);
}


// Sorts input vector of four "leptons" into two std pairs, returned as pointers
// Selection is determined by maximizing the invariant mass difference of the pairs

int GetMassPairs(vector<pair<TLorentzVector, Short_t>> leps_list,
    pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> *_large_pair,
    pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> *_small_pair)
{


    //--- INITIALIZE ---//

    // Declare pairs
    pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> pair1, pair2;
    pair<pair<TLorentzVector, Short_t>, pair<TLorentzVector, Short_t>> large_pair, small_pair;

    // Sort leptons by decreasing TRANSVERSE momentum
    sort(leps_list.begin(), leps_list.end(), SortDecPt);

    // Assign highest-Pt lepton as lep1b, which is always in pair 1
    pair<TLorentzVector, Short_t> lep1a = leps_list[0];
    leps_list.erase(leps_list.begin());
    
    unsigned config = 1;



    //--- CONFIG 1 ---//

    // OS leptons with largest momenta are Z1

    // Find lep1b, subleading lepton of pair 1
    pair<TLorentzVector, Short_t> lep1b;
    for (unsigned i = 0; i < leps_list.size(); i++)
    {
        if (leps_list[i].second != lep1a.second)
        {
            lep1b = leps_list[i];
            leps_list.erase(leps_list.begin() + i);
            break;
        }
    }
    pair1 = make_pair(lep1a, lep1b);


    // Assemble pair 2 from remaining leptons
    pair2 = make_pair(leps_list[0], leps_list[1]);


    // Find difference between masses of pairs and assign them
    TLorentzVector  pair1p4     = GetP4Sum(pair1),      pair2p4     = GetP4Sum(pair2); 
    Float_t         mass_diff_1 = pair1p4.M() - pair2p4.M();

    if (mass_diff_1 > 0)
    {
        large_pair = pair1;
        small_pair = pair2;
    }
    else
    {
        large_pair = pair2;
        small_pair = pair1;
    }



    //--- CONFIG 2 ---//

    // Swap leptons that are not SS as lep1
    // (i.e. look for pair 2 lepton with SS as trailing lep of pair 1)

    if (pair1.first.second == pair2.second.second)          // leading pair 2 lepton
        swap(pair1.first, pair2.second);
    else if (pair1.second.second == pair2.second.second)    // trailing pair 2 lepton
        swap(pair1.second, pair2.second);


    // Reassign pair momenta
    pair1p4     = GetP4Sum(pair1);      pair2p4     = GetP4Sum(pair2);
    Float_t mass_diff_2 = pair1p4.M() - pair2p4.M();


    // If this config has larger mass difference, reassign returned values
    if (fabs(mass_diff_2) > fabs(mass_diff_1))
    {
        config = 2;

        if (mass_diff_2 > 0)
        {
            large_pair = pair1;
            small_pair = pair2;
        }
        else
        {
            large_pair = pair2;
            small_pair = pair1;
        }
    }


    //--- DEBUG ---//
/*
    //   if (False):
    if (small_pair.M() < 4):
    {
        cout << "Configuration " << config << endl;
        cout << mass_diff_1 << " (1) " << mass_diff_2 << " (2)" << endl;
        cout << "\nPair 1: mass " << large_pair.M() << endl;
        cout << pair1_list[0], q_dict[pair1_list[0]])
        cout << pair1_list[1], q_dict[pair1_list[1]])
        cout << "\nPair 2: mass", small_pair.M())
        cout << pair2_list[0], q_dict[pair2_list[0]])
        cout << pair2_list[1], q_dict[pair2_list[1]])
        cout << "\n")
        sys.stdout.flush()
    }
*/


    //--- RETURN---//

    *_large_pair = large_pair;
    *_small_pair = small_pair;

    return config;
}
