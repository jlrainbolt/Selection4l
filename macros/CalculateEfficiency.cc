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
    TString selTeX[N] =    {"\\mu^{+}\\mu^{-}", "\\mbox{e}^{+}\\mbox{e}^{-}", "4l", "4\\mu", "2\\mu2\\mbox{e}", "4\\mbox{e}"};
    unsigned chanIdx[N]  = {3,      4,      5,      6,      7,      9       };
    unsigned idxChan[10]= {9,  9,  9,  0,  1,  2,  3,  4,  9,  5   };



    //
    //  OUTPUT FILE
    //

    TString prefix  = "eff";
    TString outName = prefix + "_" + suffix + "_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //    INPUT FILE
    //

    TString inName  = "eff2_" + suffix + ".root";
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
        hLeptonEff[i] = new TH1D("LeptonEff_" + selection[i], selTeX[i] + " lepton eff", 5, -0.5, 4.5);
        hLeptonEff[i]->SetDirectory(outFile);
        hLeptonEff[i]->Sumw2();
        hMuonEff[i] = new TH1D("MuonEff_" + selection[i], selTeX[i] + " muon eff", 5, -0.5, 4.5);
        hMuonEff[i]->SetDirectory(outFile);
        hMuonEff[i]->Sumw2();
        hElectronEff[i] = new TH1D("ElectronEff_" + selection[i], selTeX[i] + " electron eff", 5, -0.5, 4.5);
        hElectronEff[i]->SetDirectory(outFile);
        hElectronEff[i]->Sumw2();

        hLeptonDenom[i] = new TH1D("LeptonDenom_" + selection[i], "", 5, -0.5, 4.5);
        hLeptonDenom[i]->SetDirectory(outFile);
        hLeptonDenom[i]->Sumw2();
        hMuonDenom[i] = new TH1D("MuonDenom_" + selection[i], "", 5, -0.5, 4.5);
        hMuonDenom[i]->SetDirectory(outFile);
        hMuonDenom[i]->Sumw2();
        hElectronDenom[i] = new TH1D("ElectronDenom_" + selection[i], "", 5, -0.5, 4.5);
        hElectronDenom[i]->SetDirectory(outFile);
        hElectronDenom[i]->Sumw2();

        hLeptonEff[i]->SetXTitle("l");
        hMuonEff[i]->SetXTitle("\\mu");
        hElectronEff[i]->SetXTitle("e");
    }

    TH1D *hLeptonPtEff[N][5], *hMuonPtEff[N][5], *hElectronPtEff[N][5];
    TH1D *hLeptonPtDenom[N][5], *hMuonPtDenom[N][5], *hElectronPtDenom[N][5];
    TH1D *hLeptonEtaEff[N][5], *hMuonEtaEff[N][5], *hElectronEtaEff[N][5];
    TH1D *hLeptonEtaDenom[N][5], *hMuonEtaDenom[N][5], *hElectronEtaDenom[N][5];
    TCanvas *cMuonPtEff[N][5], *cMuonEtaEff[N][5], *cElectronPtEff[N][5], *cElectronEtaEff[N][5];
    int csize = 1000;

    for (unsigned i = 0; i < N; i++)
    {
        hLeptonPtEff[i][0] = new TH1D("hLeptonPtEff_" + selection[i], selTeX[i] + ": \\mbox{leptons}", 23, 5, 74);
        hMuonPtEff[i][0] = new TH1D("hMuonPtEff_" + selection[i], selTeX[i] + ": \\mbox{muons}", 23, 5, 74);
        hElectronPtEff[i][0] = new TH1D("hElectronPtEff_" + selection[i], selTeX[i] + ": \\mbox{electrons}", 22, 7, 73);
        hLeptonPtDenom[i][0] = new TH1D("hLeptonPtDenom_" + selection[i], "", 23, 5, 74);
        hMuonPtDenom[i][0] = new TH1D("hMuonPtDenom_" + selection[i], "", 23, 5, 74);
        hElectronPtDenom[i][0] = new TH1D("hElectronPtDenom_" + selection[i], "", 22, 7, 73);

        cMuonPtEff[i][0] = new TCanvas("MuonPtEff_" + selection[i], "", csize, csize);
        cElectronPtEff[i][0] = new TCanvas("ElectronPtEff_" + selection[i], "", csize, csize);

        hLeptonPtEff[i][1] = new TH1D("hLepton1PtEff_" + selection[i], selTeX[i] + ": l_{1}", 26, 20, 72);
        hMuonPtEff[i][1] = new TH1D("hMuon1PtEff_" + selection[i], selTeX[i] + ": \\mu_{1}", 26, 20, 72);
        hElectronPtEff[i][1] = new TH1D("hElectron1PtEff_" + selection[i], selTeX[i] + ": \\mbox{e}_{1}", 26, 20, 72);
        hLeptonPtDenom[i][1] = new TH1D("hLepton1PtDenom_" + selection[i], "", 26, 20, 72);
        hMuonPtDenom[i][1] = new TH1D("hMuon1PtDenom_" + selection[i], "", 26, 20, 72);
        hElectronPtDenom[i][1] = new TH1D("hElectron1PtDenom_" + selection[i], "", 26, 20, 72);

        cMuonPtEff[i][1] = new TCanvas("Muon1PtEff_" + selection[i], "", csize, csize);
        cElectronPtEff[i][1] = new TCanvas("Electron1PtEff_" + selection[i], "", csize, csize);

        hLeptonPtEff[i][2] = new TH1D("hLepton2PtEff_" + selection[i], selTeX[i] + ": l_{2}", 20, 10, 50);
        hMuonPtEff[i][2] = new TH1D("hMuon2PtEff_" + selection[i], selTeX[i] + ": \\mu_{2}", 20, 10, 50);
        hElectronPtEff[i][2] = new TH1D("hElectron2PtEff_" + selection[i], selTeX[i] + ": \\mbox{e}_{2}", 20, 10, 50);
        hLeptonPtDenom[i][2] = new TH1D("hLepton2PtDenom_" + selection[i], "", 20, 10, 50);
        hMuonPtDenom[i][2] = new TH1D("hMuon2PtDenom_" + selection[i], "", 20, 10, 50);
        hElectronPtDenom[i][2] = new TH1D("hElectron2PtDenom_" + selection[i], "", 20, 10, 50);

        cMuonPtEff[i][2] = new TCanvas("Muon2PtEff_" + selection[i], "", csize, csize);
        cElectronPtEff[i][2] = new TCanvas("Electron2PtEff_" + selection[i], "", csize, csize);

        hLeptonPtEff[i][3] = new TH1D("hLepton3PtEff_" + selection[i], selection[i] + ": l_{3}", 24, 5, 29);
        hMuonPtEff[i][3] = new TH1D("hMuon3PtEff_" + selection[i], selection[i] + ": \\mu_{3}", 24, 5, 29);
        hElectronPtEff[i][3] = new TH1D("hElectron3PtEff_" + selection[i], selection[i] + ": \\mbox{e}_{3}", 22, 7, 29);
        hLeptonPtDenom[i][3] = new TH1D("hLepton3PtDenom_" + selection[i], "", 24, 5, 29);
        hMuonPtDenom[i][3] = new TH1D("hMuon3PtDenom_" + selection[i], "", 24, 5, 29);
        hElectronPtDenom[i][3] = new TH1D("hElectron3PtDenom_" + selection[i], "", 22, 7, 29);

        cMuonPtEff[i][3] = new TCanvas("Muon3PtEff_" + selection[i], "", csize, csize);
        cElectronPtEff[i][3] = new TCanvas("Electron3PtEff_" + selection[i], "", csize, csize);

        hLeptonPtEff[i][4] = new TH1D("hLepton4PtEff_" + selection[i], selection[i] + ": l_{4}", 20, 5, 25);
        hMuonPtEff[i][4] = new TH1D("hMuon4PtEff_" + selection[i], selection[i] + ": \\mu_{4}", 20, 5, 25);
        hElectronPtEff[i][4] = new TH1D("hElectron4PtEff_" + selection[i], selection[i] + ": \\mbox{e}_{4}", 18, 7, 25);
        hLeptonPtDenom[i][4] = new TH1D("hLepton4PtDenom_" + selection[i], "", 20, 5, 25);
        hMuonPtDenom[i][4] = new TH1D("hMuon4PtDenom_" + selection[i], "", 20, 5, 25);
        hElectronPtDenom[i][4] = new TH1D("hElectron4PtDenom_" + selection[i], "", 18, 7, 25);

        cMuonPtEff[i][4] = new TCanvas("Muon4PtEff_" + selection[i], "", csize, csize);
        cElectronPtEff[i][4] = new TCanvas("Electron4PtEff_" + selection[i], "", csize, csize);

        for (unsigned j = 0; j < 5; j++)
        {
            TString strj = TString::Format("%i", j);

            if (j == 0)
            {
                hLeptonEtaEff[i][j] = new TH1D("hLeptonEtaEff_" + selection[i], selTeX[i] + ": leptons", 20, -2.5, 2.5);
                hMuonEtaEff[i][j] = new TH1D("hMuonEtaEff_" + selection[i], selTeX[i] + ": muons", 20, -2.4, 2.4);
                hElectronEtaEff[i][j] = new TH1D("hElectronEtaEff_" + selection[i], selTeX[i] + ": electrons", 20, -2.5, 2.5);

                cMuonEtaEff[i][0] = new TCanvas("MuonEtaEff_" + selection[i], "", csize, csize);
                cElectronEtaEff[i][0] = new TCanvas("ElectronEtaEff_" + selection[i], "", csize, csize);
            }
            else
            {
                hLeptonEtaEff[i][j] = new TH1D("hLepton" + strj + "EtaEff_" + selection[i], selTeX[i] + ": l_{" + strj + "}", 20, -2.5, 2.5);
                hMuonEtaEff[i][j] = new TH1D("hMuon" + strj + "EtaEff_" + selection[i], selTeX[i] + ": \\mu_{" + strj + "}", 20, -2.4, 2.4);
                hElectronEtaEff[i][j] = new TH1D("hElectron" + strj + "EtaEff_" + selection[i], selTeX[i] + ": \\mbox{e}_{" + strj + "}", 20, -2.5, 2.5);

                cMuonEtaEff[i][j] = new TCanvas("Muon" + strj + "EtaEff_" + selection[i], "", csize, csize);
                cElectronEtaEff[i][j] = new TCanvas("Electron" + strj + "EtaEff_" + selection[i], "", csize, csize);
            }
            hLeptonEtaDenom[i][j] = new TH1D("hLepton" + strj + "EtaDenom_" + selection[i], "", 20, -2.5, 2.5);
            hMuonEtaDenom[i][j] = new TH1D("hMuon" + strj + "EtaDenom_" + selection[i], "", 20, -2.4, 2.4);
            hElectronEtaDenom[i][j] = new TH1D("hElectron" + strj + "EtaDenom_" + selection[i], "", 20, -2.5, 2.5);

            hLeptonPtEff[i][j]->SetDirectory(outFile);
            hLeptonPtEff[i][j]->Sumw2();
            hMuonPtEff[i][j]->SetDirectory(outFile);
            hMuonPtEff[i][j]->Sumw2();
            hElectronPtEff[i][j]->SetDirectory(outFile);
            hElectronPtEff[i][j]->Sumw2();

            hLeptonPtDenom[i][j]->SetDirectory(outFile);
            hLeptonPtDenom[i][j]->Sumw2();
            hMuonPtDenom[i][j]->SetDirectory(outFile);
            hMuonPtDenom[i][j]->Sumw2();
            hElectronPtDenom[i][j]->SetDirectory(outFile);
            hElectronPtDenom[i][j]->Sumw2();

            hLeptonEtaEff[i][j]->SetDirectory(outFile);
            hLeptonEtaEff[i][j]->Sumw2();
            hMuonEtaEff[i][j]->SetDirectory(outFile);
            hMuonEtaEff[i][j]->Sumw2();
            hElectronEtaEff[i][j]->SetDirectory(outFile);
            hElectronEtaEff[i][j]->Sumw2();

            hLeptonEtaDenom[i][j]->SetDirectory(outFile);
            hLeptonEtaDenom[i][j]->Sumw2();
            hMuonEtaDenom[i][j]->SetDirectory(outFile);
            hMuonEtaDenom[i][j]->Sumw2();
            hElectronEtaDenom[i][j]->SetDirectory(outFile);
            hElectronEtaDenom[i][j]->Sumw2();
        }

        hLeptonPtEff[i][0]->SetXTitle("p_{T}^{l}");
        hMuonPtEff[i][0]->SetXTitle("p_{T}^{\\mu}");
        hElectronPtEff[i][0]->SetXTitle("p_{T}^{e}");
        hLeptonEtaEff[i][0]->SetXTitle("\\eta^{l}");
        hMuonEtaEff[i][0]->SetXTitle("\\eta^{\\mu}");
        hElectronEtaEff[i][0]->SetXTitle("\\eta^{e}");

        for (unsigned j = 1; j < 5; j++)
        {
            TString strj = TString::Format("%i", j);

            hLeptonPtEff[i][j]->SetXTitle("p_{T}^{l_{" + strj + "}}");
            hMuonPtEff[i][j]->SetXTitle("p_{T}^{\\mu_{" + strj + "}}");
            hElectronPtEff[i][j]->SetXTitle("p_{T}^{e_{" + strj + "}}");
            hLeptonEtaEff[i][j]->SetXTitle("\\eta^{l_{" + strj + "}}");
            hMuonEtaEff[i][j]->SetXTitle("\\eta^{\\mu_{" + strj + "}}");
            hElectronEtaEff[i][j]->SetXTitle("\\eta^{e_{" + strj + "}}");

            hLeptonPtEff[i][j]->SetYTitle("Efficiency");
            hMuonPtEff[i][j]->SetYTitle("Efficiency");
            hElectronPtEff[i][j]->SetYTitle("Efficiency");
            hLeptonEtaEff[i][j]->SetYTitle("Efficiency");
            hMuonEtaEff[i][j]->SetYTitle("Efficiency");
            hElectronEtaEff[i][j]->SetYTitle("Efficiency");
        }

        for (unsigned j = 0; j < 5; j++)
        {
            hLeptonPtEff[i][j]->SetMarkerStyle(20);
            hMuonPtEff[i][j]->SetMarkerStyle(20);
            hElectronPtEff[i][j]->SetMarkerStyle(20);
            hLeptonEtaEff[i][j]->SetMarkerStyle(20);
            hMuonEtaEff[i][j]->SetMarkerStyle(20);
            hElectronEtaEff[i][j]->SetMarkerStyle(20);
        }
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
    TTreeReaderValue    <vector<Short_t>>   genElecMatchIdx_(reader,    "dressedElectronMatchIdx");

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
            hLeptonPtDenom[C][0]->Fill(gen_leps[i].p4.Pt(), weight);
            hLeptonPtDenom[C][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
            hLeptonEtaDenom[C][0]->Fill(gen_leps[i].p4.Eta(), weight);
            hLeptonEtaDenom[C][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
            if (C > L4)
            {
                hLeptonDenom[L4]->Fill(z, weight);
                hLeptonDenom[L4]->Fill(i + 1, weight);
                hLeptonPtDenom[L4][0]->Fill(gen_leps[i].p4.Pt(), weight);
                hLeptonPtDenom[L4][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                hLeptonEtaDenom[L4][0]->Fill(gen_leps[i].p4.Eta(), weight);
                hLeptonEtaDenom[L4][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
            }

            if (gen_leps[i].matched)
            {
                hLeptonEff[C]->Fill(z, weight);
                hLeptonEff[C]->Fill(i + 1, weight);
                hLeptonPtEff[C][0]->Fill(gen_leps[i].p4.Pt(), weight);
                hLeptonPtEff[C][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                hLeptonEtaEff[C][0]->Fill(gen_leps[i].p4.Eta(), weight);
                hLeptonEtaEff[C][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                if (C > L4)
                {
                    hLeptonEff[L4]->Fill(z, weight);
                    hLeptonEff[L4]->Fill(i + 1, weight);
                    hLeptonPtEff[L4][0]->Fill(gen_leps[i].p4.Pt(), weight);
                    hLeptonPtEff[L4][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                    hLeptonEtaEff[L4][0]->Fill(gen_leps[i].p4.Eta(), weight);
                    hLeptonEtaEff[L4][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                }
            }

            if (abs(gen_leps[i].pdg) == 13)
            {
                hMuonDenom[C]->Fill(z, weight);
                hMuonDenom[C]->Fill(i + 1, weight);
                hMuonPtDenom[C][0]->Fill(gen_leps[i].p4.Pt(), weight);
                hMuonPtDenom[C][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                hMuonEtaDenom[C][0]->Fill(gen_leps[i].p4.Eta(), weight);
                hMuonEtaDenom[C][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                if (C > L4)
                {
                    hMuonDenom[L4]->Fill(z, weight);
                    hMuonDenom[L4]->Fill(i + 1, weight);
                    hMuonPtDenom[L4][0]->Fill(gen_leps[i].p4.Pt(), weight);
                    hMuonPtDenom[L4][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                    hMuonEtaDenom[L4][0]->Fill(gen_leps[i].p4.Eta(), weight);
                    hMuonEtaDenom[L4][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                }

                if (gen_leps[i].matched)
                {
                    hMuonEff[C]->Fill(z, weight);
                    hMuonEff[C]->Fill(i + 1, weight);
                    hMuonPtEff[C][0]->Fill(gen_leps[i].p4.Pt(), weight);
                    hMuonPtEff[C][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                    hMuonEtaEff[C][0]->Fill(gen_leps[i].p4.Eta(), weight);
                    hMuonEtaEff[C][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                    if (C > L4)
                    {
                        hMuonEff[L4]->Fill(z, weight);
                        hMuonEff[L4]->Fill(i + 1, weight);
                        hMuonPtEff[L4][0]->Fill(gen_leps[i].p4.Pt(), weight);
                        hMuonPtEff[L4][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                        hMuonEtaEff[L4][0]->Fill(gen_leps[i].p4.Eta(), weight);
                        hMuonEtaEff[L4][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                    }
                }
            }
            else if (abs(gen_leps[i].pdg) == 11)
            {
                hElectronDenom[C]->Fill(z, weight);
                hElectronDenom[C]->Fill(i + 1, weight);
                hElectronPtDenom[C][0]->Fill(gen_leps[i].p4.Pt(), weight);
                hElectronPtDenom[C][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                hElectronEtaDenom[C][0]->Fill(gen_leps[i].p4.Eta(), weight);
                hElectronEtaDenom[C][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                if (C > L4)
                {
                    hElectronDenom[L4]->Fill(z, weight);
                    hElectronDenom[L4]->Fill(i + 1, weight);
                    hElectronPtDenom[L4][0]->Fill(gen_leps[i].p4.Pt(), weight);
                    hElectronPtDenom[L4][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                    hElectronEtaDenom[L4][0]->Fill(gen_leps[i].p4.Eta(), weight);
                    hElectronEtaDenom[L4][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                }

                if (gen_leps[i].matched)
                {
                    hElectronEff[C]->Fill(z, weight);
                    hElectronEff[C]->Fill(i + 1, weight);
                    hElectronPtEff[C][0]->Fill(gen_leps[i].p4.Pt(), weight);
                    hElectronPtEff[C][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                    hElectronEtaEff[C][0]->Fill(gen_leps[i].p4.Eta(), weight);
                    hElectronEtaEff[C][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                    if (C > L4)
                    {
                        hElectronEff[L4]->Fill(z, weight);
                        hElectronEff[L4]->Fill(i + 1, weight);
                        hElectronPtEff[L4][0]->Fill(gen_leps[i].p4.Pt(), weight);
                        hElectronPtEff[L4][i+1]->Fill(gen_leps[i].p4.Pt(), weight);
                        hElectronEtaEff[L4][0]->Fill(gen_leps[i].p4.Eta(), weight);
                        hElectronEtaEff[L4][i+1]->Fill(gen_leps[i].p4.Eta(), weight);
                    }
                }
            }


        }


        count++;

    } // END event loop



    //
    //  CALCULATE
    //

    for (unsigned i = 0; i < N; i++)
    {
        hLeptonEff[i]->Divide(hLeptonEff[i], hLeptonDenom[i], 1, 1, "B");
        hMuonEff[i]->Divide(hMuonEff[i], hMuonDenom[i], 1, 1, "B");
        hElectronEff[i]->Divide(hElectronEff[i], hElectronDenom[i], 1, 1, "B");

        for (unsigned j = 0; j < 5; j++)
        {
            hLeptonPtEff[i][j]->Divide(hLeptonPtEff[i][j], hLeptonPtDenom[i][j], 1, 1, "B");
            hMuonPtEff[i][j]->Divide(hMuonPtEff[i][j], hMuonPtDenom[i][j], 1, 1, "B");
            hElectronPtEff[i][j]->Divide(hElectronPtEff[i][j], hElectronPtDenom[i][j], 1, 1, "B");

            hLeptonEtaEff[i][j]->Divide(hLeptonEtaEff[i][j], hLeptonEtaDenom[i][j], 1, 1, "B");
            hMuonEtaEff[i][j]->Divide(hMuonEtaEff[i][j], hMuonEtaDenom[i][j], 1, 1, "B");
            hElectronEtaEff[i][j]->Divide(hElectronEtaEff[i][j], hElectronEtaDenom[i][j], 1, 1, "B");
        }
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

    outFile->mkdir("muon");
    outFile->cd("muon");
    for (unsigned i = 0; i < N; i++)
    {
        if ((i != 1) && (i != 5))
        {
            hMuonEff[i]->Write();

            unsigned J = i >= L4 ? 5 : 3;

            for (unsigned j = 0; j < J; j++)
            {
                hMuonPtEff[i][j]->Write();
                hMuonEtaEff[i][j]->Write();

                cMuonPtEff[i][j]->cd();
                hMuonPtEff[i][j]->Draw();
                cMuonPtEff[i][j]->Write();
                cMuonPtEff[i][j]->SaveAs(".pdf");

                cMuonEtaEff[i][j]->cd();
                hMuonEtaEff[i][j]->Draw();
                cMuonEtaEff[i][j]->Write();
                cMuonEtaEff[i][j]->SaveAs(".pdf");
            }
        }
    }

    outFile->mkdir("electron");
    outFile->cd("electron");
    for (unsigned i = 0; i < N; i++)
    {
        if ((i != 0) && (i != 3))
        {
            hElectronEff[i]->Write();

            unsigned J = i >= L4 ? 5 : 3;

            for (unsigned j = 0; j < J; j++)
            {
                hElectronPtEff[i][j]->Write();
                hElectronEtaEff[i][j]->Write();

                cElectronPtEff[i][j]->cd();
                hElectronPtEff[i][j]->Draw();
                cElectronPtEff[i][j]->Write();
                cElectronPtEff[i][j]->SaveAs(".pdf");

                cElectronEtaEff[i][j]->cd();
                hElectronEtaEff[i][j]->Draw();
                cElectronEtaEff[i][j]->Write();
                cElectronEtaEff[i][j]->SaveAs(".pdf");
            }
        }
    }

    outFile->mkdir("lepton");
    outFile->cd("lepton");
    for (unsigned i = 0; i < N; i++)
    {
        if ((i == 2) || (i == 4))
        {
            hLeptonEff[i]->Write();

            unsigned J = i >= L4 ? 5 : 3;

            for (unsigned j = 0; j < J; j++)
            {
                hLeptonPtEff[i][j]->Write();
                hLeptonEtaEff[i][j]->Write();
            }
        }
    }

    outFile->Close();
    inFile->Close();

    cout << "Wrote histograms to " << outName << endl << endl << endl;
}
