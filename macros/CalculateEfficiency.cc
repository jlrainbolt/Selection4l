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



void CalculateEfficiency(TString suffix)
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

    const unsigned N = 6;
    unsigned                MM = 0, EE = 1, L4 = 2, M4 = 3, ME = 4, E4 = 5;     // Indices
    TString selection[N] = {"mumu", "ee",   "4l",   "4m",   "2m2e", "4e"    };
    unsigned chanIdx[N]  = {3,      4,      5,      6,      7,      9       };
    unsigned idxChan[10]= {9,  9,  9,  0,  1,  2,  3,  4,  9,  5   };



    //
    //  OUTPUT FILE
    //

    TString prefix  = "eff";
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //    INPUT FILE
    //

    TString inName  = "eff_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "_update/" + inName;
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

    inFile->GetObject("MatchedEvents_" + suffix, hMatchedEvents);
    hMatchedEvents->SetDirectory(outFile);


    TH1D *hLeptonEff[N], *hMuonEff[N], *hElectronEff[N];
    TH1D *hLeptonDenom[N], *hMuonDenom[N], *hElectronDenom[N];
    for (unsigned i = 0; i < N; i++)
    {
        hLeptonEff[i] = new TH1D("LeptonEff_" + selection[i], selection[i] + " lepton eff", 5, -0.5, 4.5);
        hLeptonEff[i]->SetDirectory(outFile);
        hMuonEff[i] = new TH1D("MuonEff_" + selection[i], selection[i] + " muon eff", 5, -0.5, 4.5);
        hMuonEff[i]->SetDirectory(outFile);
        hElectronEff[i] = new TH1D("ElectronEff_" + selection[i], selection[i] + " electron eff", 5, -0.5, 4.5);
        hElectronEff[i]->SetDirectory(outFile);

        hLeptonDenom[i] = new TH1D("LeptonDenom_" + selection[i], "", 5, -0.5, 4.5);
        hLeptonDenom[i]->SetDirectory(outFile);
        hMuonDenom[i] = new TH1D("MuonDenom_" + selection[i], "", 5, -0.5, 4.5);
        hMuonDenom[i]->SetDirectory(outFile);
        hElectronDenom[i] = new TH1D("ElectronDenom_" + selection[i], "", 5, -0.5, 4.5);
        hElectronDenom[i]->SetDirectory(outFile);
    }





    //
    //  INPUT BRANCHES
    //

    TTreeReader reader("tree_" + suffix, inFile);

//  TTreeReaderValue    <Int_t>             runNum_         (reader,    "runNum");
//  TTreeReaderValue    <Int_t>             evtNum_         (reader,    "evtNum");
//  TTreeReaderValue    <Int_t>             lumiSec_        (reader,    "lumiSec");
//  TTreeReaderValue    <UShort_t>          nPV_            (reader,    "nPV");
    TTreeReaderValue    <Float_t>           genWeight_      (reader,    "genWeight");
    TTreeReaderValue    <Float_t>           puWeight_       (reader,    "PUWeight");
    TTreeReaderValue    <UShort_t>          channel_        (reader,    "decayChannel");
    TTreeReaderValue    <Bool_t>            isMatched_      (reader,    "isMatched");

    TTreeReaderValue    <UShort_t>          nGenLeps_       (reader,    "nDressedLeptons");
    TTreeReaderValue    <UShort_t>          nGenMuons_      (reader,    "nDressedMuons");
    TTreeReaderValue    <UShort_t>          nGenElecs_      (reader,    "nDressedElectrons");
    TTreeReaderValue    <UShort_t>          nMuons_         (reader,    "nLooseMuons");
    TTreeReaderValue    <UShort_t>          nElecs_         (reader,    "nLooseElectrons");
    TTreeReaderValue    <UShort_t>          nTightMuons_    (reader,    "nTightMuons");
    TTreeReaderValue    <UShort_t>          nTightElecs_    (reader,    "nTightElectrons");

    TTreeReaderArray    <TLorentzVector>    genMuonP4_      (reader,    "dressedMuonP4");
    TTreeReaderValue    <vector<Short_t>>   genMuonQ_       (reader,    "dressedMuonQ");
    TTreeReaderValue    <vector<Bool_t>>    genMuonMatched_ (reader,    "dressedMuonIsMatched");
    TTreeReaderValue    <vector<Short_t>>   genMuonMatchIdx_(reader,    "dressedMuonMatchIdx");

    TTreeReaderArray    <TLorentzVector>    genElecP4_      (reader,    "dressedElectronP4");
    TTreeReaderValue    <vector<Short_t>>   genElecQ_       (reader,    "dressedElectronQ");
    TTreeReaderValue    <vector<Bool_t>>    genElecMatched_ (reader,    "dressedElectronIsMatched");
    TTreeReaderValue    <vector<Short_t>>   genElecMatchIdx_(reader,    "dressedElectronMatchOdx");

    TTreeReaderArray    <TLorentzVector>    muonP4_         (reader,    "muonP4");
    TTreeReaderValue    <vector<Short_t>>   muonQ_          (reader,    "muonQ");
    TTreeReaderValue    <vector<Bool_t>>    muonIsTight_    (reader,    "muonIsTight");

    TTreeReaderArray    <TLorentzVector>    elecP4_         (reader,    "electronP4");
    TTreeReaderValue    <vector<Short_t>>   elecQ_          (reader,    "electronQ");
    TTreeReaderValue    <vector<Bool_t>>    elecIsTight_    (reader,    "electronIsTight");





    ////
    ////
    ////    EVENT LOOP
    ////
    ////


    int nEvents = reader.GetEntries(kTRUE);

    cout << endl;
    cout << "Running over " << nEvents << " events..." << flush;
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

        float       weight      = (*genWeight_) * (*puWeight_);
        bool        isMatched   = *isMatched_;
        unsigned    C           = idxChan[*channel_];

        unsigned    nGenLeps    = *nGenLeps_;
        unsigned    nGenMuons   = *nGenMuons_,          nGenElecs   = *nGenElecs_;
        unsigned    nMuons      = *nMuons_,             nElecs      = *nElecs_;
        unsigned    nTightMuons = *nTightMuons_,        nTightElecs = *nTightElecs_;





        ////
        ////
        ////    OBJECTS
        ////
        ////


        //
        //  RECO LEPTONS
        //

        vector<Lepton> muons, elecs;

        for (unsigned i = 0; i < nMuons; i++)
        {
            Lepton  muon;

            muon.p4         = muonP4_.At(i);
            muon.q          = (*muonQ_)[i];
            muon.pdg        = -13 * muon.q;
            muon.tight      = (*muonIsTight_)[i];

            muons.push_back(muon);
        }

        for (unsigned i = 0; i < nElecs; i++)
        {
            Lepton  elec;

            elec.p4         = elecP4_.At(i);
            elec.q          = (*elecQ_)[i];
            elec.pdg        = -11 * elec.q;
            elec.tight      = (*elecIsTight_)[i];

            elecs.push_back(elec);
        }



        //
        //  GEN LEPTONS
        //

        vector<Lepton> gen_muons, gen_elecs;

        for (unsigned i = 0; i < nGenMuons; i++)
        {
            Lepton muon;

            muon.p4         = genMuonP4_.At(i);
            muon.q          = (*genMuonQ_)[i];
            muon.pdg        = -13 * muon.q;
            muon.matched    = (*genMuonMatched_)[i];

            if (muon.matched)
            {
                unsigned m = (*genMuonMatchIdx_)[i];
                muon.SetMatch(muons[m]);
            }

            gen_muons.push_back(muon);
        }

        for (unsigned i = 0; i < nGenElecs; i++)
        {
            Lepton elec;

            elec.p4         = genElecP4_.At(i);
            elec.q          = (*genElecQ_)[i];
            elec.pdg        = -11 * elec.q;
            elec.matched    = (*genElecMatched_)[i];

            if (elec.matched)
            {
                unsigned m = (*genElecMatchIdx_)[i];
                elec.SetMatch(elecs[m]);
            }

            gen_elecs.push_back(elec);
        }

        vector<Lepton> gen_leps = gen_muons;
        gen_leps.insert(gen_leps.end(), gen_elecs.begin(), gen_elecs.end());
        sort(gen_leps.begin(), gen_leps.end(), DecreasingPt);





        ////
        ////
        ////    HISTOGRAMS
        ////
        ////


        for (unsigned i = 0; i < nGenLeps; i++)
        {
            double z = 0;

            hLeptonDenom[C]->Fill(z, weight);
            hLeptonDenom[C]->Fill(i + 1, weight);

            if (gen_leps[i].matched)
            {
                hLeptonEff[C]->Fill(z, weight);
                hLeptonEff[C]->Fill(i + 1, weight);
            }

            if (abs(gen_leps[i].pdg) == 13)
            {
                hMuonDenom[C]->Fill(z, weight);
                hMuonDenom[C]->Fill(i + 1, weight);

                if (gen_leps[i].matched)
                {
                    hMuonEff[C]->Fill(z, weight);
                    hMuonEff[C]->Fill(i + 1, weight);
                }
            }
            else if (abs(gen_leps[i].pdg) == 11)
            {
                hElectronDenom[C]->Fill(z, weight);
                hElectronDenom[C]->Fill(i + 1, weight);

                if (gen_leps[i].matched)
                {
                    hElectronEff[C]->Fill(z, weight);
                    hElectronEff[C]->Fill(i + 1, weight);
                }
            }
        }


        count++;

    } // END event loop



    //
    //  CALCULATE
    //

    // Add together 4l channel histograms
    hLeptonDenom[L4] = (TH1D*) hLeptonDenom[M4]->Clone();
    hLeptonDenom[L4]->Add(hLeptonDenom[ME]);
    hLeptonDenom[L4]->Add(hLeptonDenom[E4]);
    hMuonDenom[L4] = (TH1D*) hMuonDenom[M4]->Clone();
    hMuonDenom[L4]->Add(hMuonDenom[ME]);
    hElectronDenom[L4] = (TH1D*) hElectronDenom[E4]->Clone();
    hElectronDenom[L4]->Add(hElectronDenom[ME]);

    hLeptonEff[L4] = (TH1D*) hLeptonEff[M4]->Clone("LeptonEff_" + selection[L4]);
    hLeptonEff[L4]->SetTitle(selection[L4] + " lepton eff");
    hLeptonEff[L4]->Add(hLeptonEff[ME]);
    hLeptonEff[L4]->Add(hLeptonEff[E4]);
    hMuonEff[L4] = (TH1D*) hMuonEff[M4]->Clone("MuonEff_" + selection[L4]);
    hMuonEff[L4]->SetTitle(selection[L4] + " muon eff");
    hMuonEff[L4]->Add(hMuonEff[ME]);
    hElectronEff[L4] = (TH1D*) hElectronEff[E4]->Clone("ElectronEff_" + selection[L4]);
    hElectronEff[L4]->SetTitle(selection[L4] + " electron eff");
    hElectronEff[L4]->Add(hElectronEff[ME]);

    for (unsigned i = 0; i < N; i++)
    {
        hLeptonEff[i]->Divide(hLeptonDenom[i]);
        hMuonEff[i]->Divide(hMuonDenom[i]);
        hElectronEff[i]->Divide(hElectronDenom[i]);
    }


    cout << endl << endl;
    cout << "Processed all trees in " << inName << endl;



    //
    //  WRITE FILE
    //

    outFile->cd();

    hTotalEvents->Write();
    hSelectedEvents->Write();
    hMatchedEvents->Write();

    for (unsigned i = 0; i < N; i++)
    {
        if ((i == 2) || (i == 4))
            hLeptonEff[i]->Write();
    }

    for (unsigned i = 0; i < N; i++)
    {
        if ((i != 1) && (i != 5))
            hMuonEff[i]->Write();
    }
    for (unsigned i = 0; i < N; i++)
    {
        if ((i != 0) && (i != 3))
            hElectronEff[i]->Write();
    }

    outFile->Close();
    inFile->Close();

    cout << "Wrote histograms to " << outName << endl << endl << endl;
}
