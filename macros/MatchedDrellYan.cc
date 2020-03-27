// STL
#include <vector>
#include <cmath>
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
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

// Cuts
//#include "Cuts2018.hh"
//#include "Cuts2017.hh"
#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;



/*
**  MatchedAnalysis
**
**  Reads signal events from "Selected" level tree.  Performs lepton matching and writes out all
**  "Boosted" quantities.
*/

void MatchedDrellYan()
{

    //
    //    OPTIONS
    //

    bool debug = kFALSE;
    int analyzeEvents = 1000;
    int printEvery = 100;



    //
    //  SAMPLE INFO
    //

    const unsigned N = 2;
    unsigned                   MM = 0,  EE = 1;     // Indices
    TString selection[N]    = {"mumu",  "ee"    };
    unsigned chanIdx[N]     = {3,       4       };



    //
    //  OUTPUT FILE
    //

    TString prefix  = "matched",    suffix = "zjets_m-50";
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);



    //
    //  OUTPUT BRANCHES
    //

    // Event info
    Bool_t              isMatched,  hasTauDecay;
    Int_t               runNum,     evtNum,     lumiSec;
    UShort_t            nPV,        channel;
    Float_t             weight,     genWeight,  qtWeight,   puWeight,   ecalWeight;
    Float_t             trigWeight, idWeight;
    Bool_t              muonTrig,   siMuTrig,   diMuTrig,   elecTrig,   siElTrig,   diElTrig;


    // Lab frame objects
    TLorentzVector      z1p4,           gen_z1p4;
    UShort_t            z1pdg;

    TLorentzVector      l1p4,           l2p4;
    TLorentzVector      gen_l1p4,       gen_l2p4;
    Float_t             l1dr,           l2dr;
    Short_t             l1pdg,          l2pdg;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("isMatched",    &isMatched);

        tree[i]->Branch("runNum",       &runNum);       tree[i]->Branch("evtNum",       &evtNum);
        tree[i]->Branch("lumiSec",      &lumiSec);      tree[i]->Branch("nPV",          &nPV);
        tree[i]->Branch("weight",       &weight);       tree[i]->Branch("genWeight",    &genWeight);
        tree[i]->Branch("qtWeight",     &qtWeight);     tree[i]->Branch("puWeight",     &puWeight);
        tree[i]->Branch("ecalWeight",   &ecalWeight);   tree[i]->Branch("trigWeight",   &trigWeight);
        tree[i]->Branch("idWeight",     &idWeight);
        tree[i]->Branch("channel",      &channel);      tree[i]->Branch("hasTauDecay",  &hasTauDecay);

        tree[i]->Branch("z1p4",     &z1p4);         tree[i]->Branch("gen_z1p4",     &gen_z1p4);
        tree[i]->Branch("z1pdg",    &z1pdg);

        tree[i]->Branch("l1p4",         &l1p4);         tree[i]->Branch("gen_l1p4",     &gen_l1p4);
        tree[i]->Branch("l1dr",     &l1dr);             tree[i]->Branch("l1pdg",        &l1pdg);
        tree[i]->Branch("l2p4",         &l2p4);         tree[i]->Branch("gen_l2p4",     &gen_l2p4);
        tree[i]->Branch("l2dr",     &l2dr);             tree[i]->Branch("l2pdg",        &l2pdg);
    }



    //
    //    INPUT FILE
    //

    TString inName  = "selected_" + suffix + "_truth.root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "_update/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl;



    //
    //  HISTOGRAMS
    //

    TH1D *hTotalEvents, *hSelectedEvents, *hMatchedEvents;

    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);

    inFile->GetObject("SelectedEvents_" + suffix, hSelectedEvents);
    hSelectedEvents->SetDirectory(outFile);

    hMatchedEvents = new TH1D("MatchedEvents_" + suffix, "MatchedEvents", 10, 0.5, 10.5);
    hMatchedEvents->SetDirectory(outFile);



    //
    //  INPUT BRANCHES
    //

    for (unsigned i = 0; i < N; i++)
    {
        TTreeReader reader(selection[i] + "_" + suffix, inFile);

        TTreeReaderValue    <Int_t>             runNum_         (reader,    "runNum");
        TTreeReaderValue    <Int_t>             evtNum_         (reader,    "evtNum");
        TTreeReaderValue    <Int_t>             lumiSec_        (reader,    "lumiSec");
        TTreeReaderValue    <UShort_t>          nPV_            (reader,    "nPV");
        TTreeReaderValue    <Float_t>           weight_         (reader,    "weight");
        TTreeReaderValue    <Float_t>           genWeight_      (reader,    "genWeight");
        TTreeReaderValue    <Float_t>           qtWeight_       (reader,    "qtWeight");
        TTreeReaderValue    <Float_t>           puWeight_       (reader,    "puWeight");
        TTreeReaderValue    <Float_t>           ecalWeight_     (reader,    "ecalWeight");
        TTreeReaderValue    <Float_t>           trigWeight_     (reader,    "trigWeight");
        TTreeReaderValue    <Float_t>           idWeight_       (reader,    "idWeight");
        TTreeReaderValue    <UShort_t>          channel_        (reader,    "channel");
        TTreeReaderValue    <Bool_t>            hasTauDecay_    (reader,    "hasTauDecay");

        TTreeReaderValue    <TLorentzVector>    z1p4_           (reader,    "z1p4");
        TTreeReaderValue    <UShort_t>          z1pdg_          (reader,    "z1pdg");
        TTreeReaderValue    <TLorentzVector>    l1p4_           (reader,    "l1p4");
        TTreeReaderValue    <Short_t>           l1pdg_          (reader,    "l1pdg");
        TTreeReaderValue    <TLorentzVector>    l2p4_           (reader,    "l2p4");
        TTreeReaderValue    <Short_t>           l2pdg_          (reader,    "l2pdg");

        TTreeReaderValue    <UShort_t>          nGenMuons_      (reader,    "nDressedMuons");
        TTreeReaderValue    <UShort_t>          nGenElecs_      (reader,    "nDressedElectrons");
        TTreeReaderArray    <TLorentzVector>    genMuonP4_      (reader,    "dressedMuonP4");
        TTreeReaderValue    <vector<Short_t>>   genMuonQ_       (reader,    "dressedMuonQ");
        TTreeReaderArray    <TLorentzVector>    genElecP4_      (reader,    "dressedElectronP4");
        TTreeReaderValue    <vector<Short_t>>   genElecQ_       (reader,    "dressedElectronQ");
        TTreeReaderValue    <TLorentzVector>    genLepsP4_      (reader,    "dressedLeptonsP4");
        TTreeReaderValue    <TLorentzVector>    u_l1p4_         (reader,    "uncorr_l1p4");
        TTreeReaderValue    <TLorentzVector>    u_l2p4_         (reader,    "uncorr_l2p4");






        ////
        ////
        ////    EVENT LOOP
        ////
        ////


        int nEvents = reader.GetEntries(kTRUE);

        cout << endl;
        cout << "Running over " << nEvents << " selected " + selection[i] + " events..." << flush;
        if (debug)  cout << "Will analyze " << analyzeEvents << " events" << endl;
        cout << endl;

        long count = 0;     // track total number of analyzed events

        while (reader.Next())
        {

            //
            //  DEBUG PRINTOUT
            //

            const long currentEntry = reader.GetCurrentEntry();
            bool print = kFALSE;

            if (debug)
            {
                if (count >= analyzeEvents)
                    break;

                if (currentEntry % printEvery == 0)
                {
                    print = kTRUE;
                    cout << endl << count << endl;
                }
            }
            else if (currentEntry % 10000 == 0)
            {
                cout << "Processed " << reader.GetCurrentEntry() << " of " << nEvents << " events";
                cout << endl;
            }



            //
            //  EVENT INFO
            //                

            // Quantities copied directly to output tree
            runNum      = *runNum_;     evtNum      = *evtNum_;         lumiSec     = *lumiSec_;
            nPV         = *nPV_;        hasTauDecay = *hasTauDecay_;    weight      = *weight_;
            genWeight   = *genWeight_;  ecalWeight  = *ecalWeight_;     puWeight    = *puWeight_;
            trigWeight  = *trigWeight_; qtWeight    = *qtWeight_;       idWeight    = *idWeight_;   
            channel     = *channel_;

            z1p4    = *z1p4_;           z1pdg   = *z1pdg_;
            l1p4    = *l1p4_;           l1pdg   = *l1pdg_;
            l2p4    = *l2p4_;           l2pdg   = *l2pdg_;

            // Quantities used in analysis, but not written out
            unsigned    nGenMuons   = *nGenMuons_,          nGenElecs   = *nGenElecs_;

            // Innocent until proven guilty, unless it's a tau event
            isMatched   = !hasTauDecay;



            //
            //  RECO LEPTONS
            //

            vector<Lepton> leps(2);

            leps[0].p4  = l1p4;         leps[0].pdg = l1pdg;
            leps[1].p4  = l2p4;         leps[1].pdg = l2pdg;

            // Fill uncorrected p4
            leps[0].u_p4     = *u_l1p4_;                leps[1].u_p4     = *u_l2p4_;

            for (unsigned i = 0; i < leps.size(); i++)
                leps[i].q = -1 * copysign(1, leps[i].pdg);



            //
            //  PRESELECTION
            //

            // Make sure there are enough lepons available for matching
            // (Remember: there *can* be extras!)

            hMatchedEvents->Fill(1, weight);

            if (nGenMuons + nGenElecs != 2)
            {
                if (print)
                    cout << "Wrong number of gen leptons" << endl;

                isMatched = kFALSE;

                tree[i]->Fill();

                continue;
            }



            //
            //  GEN LEPTONS
            //

            vector<Lepton> gen_muons, gen_elecs;

            for (unsigned i = 0; i < nGenMuons; i++)
            {
                Lepton muon;

                muon.p4     = genMuonP4_.At(i);
                muon.q      = (*genMuonQ_)[i];
                muon.pdg    = -13 * muon.q;

                gen_muons.push_back(muon);
            }

            for (unsigned i = 0; i < nGenElecs; i++)
            {
                Lepton elec;

                elec.p4     = genElecP4_.At(i);
                elec.q      = (*genElecQ_)[i];
                elec.pdg    = -11 * elec.q;

                gen_elecs.push_back(elec);
            }

            vector<Lepton> gen_leps = gen_muons;
            gen_leps.insert(gen_leps.end(), gen_elecs.begin(), gen_elecs.end());






            ////
            ////
            ////    MATCHING
            ////
            ////


            bool failedDRCut = kFALSE;

            for (unsigned i = 0; i < leps.size(); i++)  // loop over reco muons
            {
                // Find DeltaR between reco lep and each gen lep
                vector<Float_t> deltaR(gen_leps.size());

                for (unsigned j = 0; j < gen_leps.size(); j++)
                {
                    deltaR[j] = leps[i].u_p4.DeltaR(gen_leps[j].p4);

                    if (print)
                        cout << deltaR[j] << "\t";
                }
                if (print)  cout << endl;


                // Find index of minimum DeltaR
                unsigned m = min_element(deltaR.begin(), deltaR.end()) - deltaR.begin();
                if (deltaR[m] > MATCH_DR_MAX)
                {
                    failedDRCut = kTRUE;
                    break;
                }

                leps[i].SetMatch(gen_leps[m]);

                if (print)
                {
                    cout << "Lepton " << i << ":\t" << leps[i].dr << "\t";
                    cout << leps[i].pdg << ", " << gen_leps[m].pdg << endl;
                }
                
                gen_leps.erase(gen_leps.begin() + m);

            }
            if (failedDRCut)
            {
                if (print)
                    cout << "Could not find a match!" << endl;

                isMatched = kFALSE;

                tree[i]->Fill();

                continue;
            }


            // Now we for sure have a match!
            hMatchedEvents->Fill(chanIdx[i], weight);



            //
            //  FILL TREE
            //

            gen_z1p4    = leps[0].m_p4 + leps[1].m_p4;
            gen_l1p4    = leps[0].m_p4;             l1dr        = leps[0].dr;
            gen_l2p4    = leps[1].m_p4;             l2dr        = leps[1].dr;

            tree[i]->Fill();

            count++;

        } // END event loop

            cout << "done!" << endl;

            cout << "Matched " << tree[i]->GetEntries("isMatched") << "/";
            cout << reader.GetCurrentEntry() << " events" << endl;
            cout << endl << endl;

    } // END tree loop


    cout << endl << endl;
    cout << "Processed all trees in " << inName << endl;



    //
    //  WRITE FILE
    //

    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();
    hTotalEvents->Write();
    hSelectedEvents->Write();
    hMatchedEvents->Write();

    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
