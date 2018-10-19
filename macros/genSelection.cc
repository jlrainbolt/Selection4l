// STL
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"
#include "Utilities.hh"


using namespace std;




void genSelection(const TString suffix)
{
    // Writes out a gen tree that is "identical" to the reco downstream tree.
    // Reads from post-BLT analyzer (PhaseSpaceAnalyzer) ntuples.
    // So it only works for signal (4-lepton) events...



    ///////////////////
    //    OPTIONS    //
    ///////////////////


    bool debug = kFALSE;



    ////  DATASET

    int year = 2017;
    TString yearStr = TString::Format("%i", year);
    TString dir = yearStr;



    ////  SELECTION

    const unsigned N = 5;
    unsigned                L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N] = {"4l",   "4m",   "2m2e", "2e2m", "4e"};



    ////  FIDUCIAL CUTS

    const float PT1_MIN = 20,   PT2_MIN = 10,   PT_MIN = 5;
    const float ETA_MAX = 2.5;




    ///////////////////////
    //    OUTPUT FILE    //
    ///////////////////////


    TString output  = "trees_gen_" + suffix + ".root";
    TFile *outFile  = new TFile(output, "RECREATE");

    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);



    ////  BRANCHES

    UInt_t          evtNum;
    Float_t         weight;
    Bool_t          isFiducial;

    TLorentzVector  z1p4, z2p4, zzp4, l1p4, l2p4, l3p4, l4p4;
    Short_t         z1pdg, z2pdg, l1pdg, l2pdg, l3pdg, l4pdg;

    TLorentzVector  b_z1p4, b_z2p4, b_ttp4;
    TLorentzVector  b_l1p4, b_l2p4, b_l3p4, b_l4p4;
    Short_t         b_l1pdg, b_l2pdg, b_l3pdg, b_l4pdg;

    Float_t         b_theta, b_phi, b_z1alpha, b_z2alpha;
    Float_t         bb_z1theta, bb_z2theta;


    for (unsigned i = 0; i < N; i++)
    {
        // Event info
        tree[i]->Branch("evtNum",   &evtNum);           tree[i]->Branch("weight",   &weight);
        tree[i]->Branch("isFiducial", &isFiducial);


        // Pair/group momenta, id
        tree[i]->Branch("zzp4",     &zzp4);
        tree[i]->Branch("z1p4",     &z1p4);             tree[i]->Branch("z1pdg",    &z1pdg);
        tree[i]->Branch("z2p4",     &z2p4);             tree[i]->Branch("z2pdg",    &z2pdg);


        // Lepton momenta, id
        tree[i]->Branch("l1p4",     &l1p4);             tree[i]->Branch("l1pdg",    &l1pdg);
        tree[i]->Branch("l2p4",     &l2p4);             tree[i]->Branch("l2pdg",    &l2pdg);
        tree[i]->Branch("l3p4",     &l3p4);             tree[i]->Branch("l3pdg",    &l3pdg);
        tree[i]->Branch("l4p4",     &l4p4);             tree[i]->Branch("l4pdg",    &l4pdg);


        // Boosted quantities for differential distributions
        tree[i]->Branch("b_z1p4",   &b_z1p4);           tree[i]->Branch("b_z2p4",   &b_z2p4);
        tree[i]->Branch("b_ttp4",   &b_ttp4);

        tree[i]->Branch("b_l1p4",   &b_l1p4);           tree[i]->Branch("b_l1pdg",  &b_l1pdg);
        tree[i]->Branch("b_l2p4",   &b_l2p4);           tree[i]->Branch("b_l2pdg",  &b_l2pdg);
        tree[i]->Branch("b_l3p4",   &b_l3p4);           tree[i]->Branch("b_l3pdg",  &b_l3pdg);
        tree[i]->Branch("b_l4p4",   &b_l4p4);           tree[i]->Branch("b_l4pdg",  &b_l4pdg);

        // Observables
        tree[i]->Branch("b_theta",  &b_theta);          tree[i]->Branch("b_phi",    &b_phi);
        tree[i]->Branch("b_z1alpha", &b_z1alpha);       tree[i]->Branch("b_z2alpha", &b_z2alpha);
        tree[i]->Branch("bb_z1theta", &bb_z1theta);     tree[i]->Branch("bb_z2theta", &bb_z2theta);
    }




    //////////////////////
    //    INPUT FILE    //
    //////////////////////


    TString path    = "root://cmseos.fnal.gov//store/user/jrainbol/Trees/";
    TString input   = path + dir + "/gen_" + suffix + "/hard_" + suffix + ".root";

    cout << input << endl;

    TFile *inFile   = TFile::Open(input);



    ////  HISTOGRAMS

    TH1D *hTotalEvents, *hPhaseSpaceEvents, *hFiducialEvents;

    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);

    inFile->GetObject("PhaseSpaceEvents_" + suffix, hPhaseSpaceEvents);
    hPhaseSpaceEvents->SetDirectory(outFile);

    hFiducialEvents = new TH1D("FiducialEvents_" + suffix, "FiducialEvents", 10, 0.5, 10.5);
    hFiducialEvents->SetDirectory(outFile);
    hFiducialEvents->Sumw2();



    ////  BRANCHES

    TTreeReader reader("tree_" + suffix, inFile);

    // Event
    TTreeReaderValue<UInt_t>            evtNum_(reader,     "evtNumber.eventNumber");
    TTreeReaderValue<Float_t>           weight_(reader,     "genWeight");
    TTreeReaderValue<UInt_t>            channel_(reader,    "decayChannel");


    // Leptons
    TTreeReaderValue<UShort_t>          nMuons_(reader,     "nHardProcMuons");
    TTreeReaderValue<UShort_t>          nElecs_(reader,     "nHardProcElectrons");

    TTreeReaderArray<TLorentzVector>    muonP4_(reader,     "hardProcMuonP4");
    TTreeReaderValue<vector<Short_t>>   muonQ_(reader,      "hardProcMuonQ");
    TTreeReaderValue<vector<UShort_t>>  muonZ_(reader,      "hardProcMuonZIndex");

    TTreeReaderArray<TLorentzVector>    elecP4_(reader,     "hardProcElectronP4");
    TTreeReaderValue<vector<Short_t>>   elecQ_(reader,      "hardProcElectronQ");
    TTreeReaderValue<vector<UShort_t>>  elecZ_(reader,      "hardProcElectronZIndex");

    TTreeReaderValue<TLorentzVector>    lepsP4_(reader,     "hardProcLeptonsP4");




    //////////////////////
    //    EVENT LOOP    //
    //////////////////////


    unsigned event = 0;
    while (reader.Next() && event < 10000)
    {
//      event++;


        if (debug)
        {
            event++;
            cout << endl << event << endl;
        }




        /////////////////////
        //    FILL INFO    //
        /////////////////////


        ////  EVENT

        // Event info 
        evtNum  = *evtNum_;
        weight  = *weight_;

        // Counters
        UShort_t        channel = *channel_;
        UShort_t        nMuons  = *nMuons_,         nElecs  = *nElecs_;

        // Z frame boost
        TLorentzVector  lepsP4  = *lepsP4_;
        TVector3        z_boost = lepsP4.BoostVector();



        ////  LEPTONS

        vector<Lepton> muons, elecs;
        unsigned i_m = 0, i_e = 0;


        // Muons
        for (const TLorentzVector& muonP4: muonP4_)
        {
            Lepton      muon;

            muon.p4     = muonP4;
            muon.b_p4   = BoostP4(muonP4, z_boost);
            muon.q      = (*muonQ_)[i_m];
            muon.pdg    = -13 * muon.q;
            muon.mother = (*muonZ_)[i_m];

            muons.push_back(muon);
            i_m++;
        }
        sort(muons.begin(), muons.end(), DecreasingLabPt);


        // Electrons
        for (const TLorentzVector& elecP4: elecP4_)
        {
            Lepton      elec;

            elec.p4     = elecP4;
            elec.b_p4   = BoostP4(elecP4, z_boost);
            elec.q      = (*elecQ_)[i_e];
            elec.pdg    = -11 * elec.q;
            elec.mother = (*elecZ_)[i_e];

            elecs.push_back(elec);
            i_e++;
        }
        sort(elecs.begin(), elecs.end(), DecreasingLabPt);


        // All leptons
        vector<Lepton> leps = muons;
        leps.insert(leps.end(), elecs.begin(), elecs.end());

        if (leps.size() != 4)
        {
            cout << "Event does not have four leptons" << endl;
            continue;
        }


        // Sort by Pt in lab frame
        sort(leps.begin(), leps.end(), DecreasingLabPt);

        if (debug)
        {
            cout << "Lab Pt:\t";
            for (unsigned i = 0; i < leps.size(); i++)
                cout << leps[i].p4.Pt() << "\t";
            cout << endl;
        }


        // Sort by P in Z CM frame
        vector<Lepton> b_leps = leps;
        sort(b_leps.begin(), b_leps.end(), DecreasingBoostedP);




        //////////////////////
        //    CATEGORIZE    //
        //////////////////////


        ////  DECAY CHANNEL

        unsigned    C;                  // channel index
        LeptonPair  z1, z2;             // first pair member has higher Pt
        
        if (debug)
            cout << channel << endl;

        if      (channel == 6)  // 4m
        {
            C = M4;
            MakePairsFromMother(muons, &z1, &z2);
        }
        else if (channel == 7)  // 2m2e
        {
            C = ME;
            z1.SetMembers(muons[0], muons[1]);
            z2.SetMembers(elecs[0], elecs[1]);
        }
        else if (channel == 8)  // 2e2m
        {
            C = EM;
            z1.SetMembers(elecs[0], elecs[1]);
            z2.SetMembers(muons[0], muons[1]);
        }
        else if (channel == 9)  // 4e
        {
            C = E4;
            MakePairsFromMother(elecs, &z1, &z2);
        }
        else
        {
            cout << "Invalid decay channel" << endl;
            continue;
        }



        ////  FIDUCIAL CHECK

        isFiducial = kTRUE;     // innocent until proven guilty?

        // Lepton Pt
        if (leps[0].p4.Pt() < PT1_MIN)
            isFiducial = kFALSE;

        if (leps[1].p4.Pt() < PT2_MIN)
            isFiducial = kFALSE;

        if (leps[2].p4.Pt() < PT_MIN)
            isFiducial = kFALSE;

        if (leps[3].p4.Pt() < PT_MIN)
            isFiducial = kFALSE;


        // Lepton Eta
        for (unsigned i = 0; i < leps.size(); i++)
        {
            if (fabs(leps[i].p4.Eta()) > ETA_MAX)
                isFiducial = kFALSE;
        }


        hFiducialEvents->Fill(1, weight);

        if (isFiducial)
            hFiducialEvents->Fill(channel, weight);




        ///////////////////////
        //    OBSERVABLES    //
        ///////////////////////


        ////  Z FRAME

        // Get positive and negative lepton P3s
        TVector3    z1plus_p3   = z1.plus->b_p4.Vect();
        TVector3    z1minus_p3  = z1.minus->b_p4.Vect();

        TVector3    z2plus_p3   = z2.plus->b_p4.Vect();
        TVector3    z2minus_p3  = z2.minus->b_p4.Vect();


        // "alpha": angle between paired leptons
        b_z1alpha   = z1plus_p3.Angle(z1minus_p3);
        b_z2alpha   = z2plus_p3.Angle(z2minus_p3);


        // Find normals to z1, z2 decay planes
        TVector3    z1norm  = z1plus_p3.Cross(z1minus_p3);
        TVector3    z2norm  = z2plus_p3.Cross(z2minus_p3);


        // "phi": angle between decay planes
        b_phi   = z1norm.Angle(z2norm);


        // "theta": angle between trailing pair 1 lepton and low-mass pair
        TVector3    z1low_p3    = z1.secondP->b_p4.Vect();

        b_theta = z2.b_p4.Angle(z1low_p3);



        ////  OTHER FRAMES

        // Get boost vectors for each pair
        TVector3    z1_boost = z1.p4.BoostVector(),     z2_boost = z2.p4.BoostVector();


        // Boost positive lepton of each pair into its pair's CM frame
        TVector3    b1_z1plus_p3    = BoostP3(z1.plus->p4,  z1_boost);
        TVector3    b2_z2plus_p3    = BoostP3(z2.plus->p4,  z2_boost);


        // Boost each pair into the other pair's CM frame
        TVector3    b1_z2_p3    = BoostP3(z2.p4, z1_boost);
        TVector3    b2_z1_p3    = BoostP3(z1.p4, z2_boost);


        // "theta_Zx": angle between positive pair x lepton and pair y in pair x CM frame
        bb_z1theta  = b1_z1plus_p3.Angle(b1_z2_p3);
        bb_z2theta  = b2_z2plus_p3.Angle(b2_z1_p3);




        /////////////////////
        //    FILL TREE    //
        /////////////////////


        // Pairs, etc.
        zzp4    = lepsP4;

        z1p4    = z1.p4;                    z1pdg   = z1.pdg;
        z2p4    = z2.p4;                    z2pdg   = z2.pdg;
        
        b_z1p4  = z1.b_p4;                  b_z2p4  = z2.b_p4;

        b_ttp4  = b_leps[1].b_p4 + b_leps[2].b_p4 + b_leps[3].b_p4;


        // Leptons
        l1p4    = leps[0].p4;               l1pdg   = leps[0].pdg;
        l2p4    = leps[1].p4;               l2pdg   = leps[1].pdg;
        l3p4    = leps[2].p4;               l3pdg   = leps[2].pdg;
        l4p4    = leps[3].p4;               l4pdg   = leps[3].pdg;

        b_l1p4  = b_leps[0].b_p4;           b_l1pdg = b_leps[0].pdg;
        b_l2p4  = b_leps[1].b_p4;           b_l2pdg = b_leps[1].pdg;
        b_l3p4  = b_leps[2].b_p4;           b_l3pdg = b_leps[2].pdg;
        b_l4p4  = b_leps[3].b_p4;           b_l4pdg = b_leps[3].pdg;




        tree[C]->Fill();
        tree[L4]->Fill();

//      event++;

    } // END event loop




    //////////////////////
    //    WRITE FILE    //
    //////////////////////


    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();
    hTotalEvents->Write();
    hPhaseSpaceEvents->Write();
    hFiducialEvents->Write();

    outFile->Purge();
    outFile->Close();
    inFile->Close();
}
