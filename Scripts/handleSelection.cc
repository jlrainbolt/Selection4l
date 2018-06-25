#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

TH1D* Select2m(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt,
        TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2,
        TH1F *h_npv);
TH1D* Select2e(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt,
        TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2, 
        TH1F *h_npv);
TH1D* Select4m(TString rootFile, TString suffix, TH1F *h_m4l, TH1F *h_qt4l,
        TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
        TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4,
        TH1F *h_npv);
TH1D* Select4e(TString rootFile, TString suffix, TH1F *h_m4l, TH1F *h_qt4l,
        TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
        TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4,
        TH1F *h_npv);


using namespace std;

void handleSelection(const TString selection, const TString suffix, const TString id)
{
    // Choose electron or muon
    Bool_t sel2l = kFALSE, sel4l = kFALSE;
    Bool_t sel2m = kFALSE, sel2e = kFALSE;
    Bool_t sel4e = kFALSE, sel4m = kFALSE;
//  Bool_t sel2e2m = kFALSE, sel2m2e = kFALSE;
    if (selection == "mumu" || selection == "2m")
    {
        sel2l = kTRUE;      sel2m = kTRUE;
    }
    else if (selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;      sel2e = kTRUE;
    }
    else if (selection == "4mu" || selection == "4m")
    {
        sel4l = kTRUE;      sel4m = kTRUE;
    }
    else if (selection == "4e")
    {
        sel4l = kTRUE;      sel4e = kTRUE;
    }


    // Path to directory where ALL trees are stored
    TString path = "root://cmsxrootd.fnal.gov//store/user/jrainbol/Trees/2016/";


    // Names of ROOT files and trees for each sample
    TString dir     = suffix + "/";
    TString file    = suffix + "_" + id + ".root";


    // Name of output file
    TString output, lepton, Lepton;
    if (sel2m)
    {
        output = "mumu_";   lepton = "muon";        Lepton = "Muon";
    }
    else if (sel2e)
    {
        output = "ee_";     lepton = "electron";    Lepton = "Electron";
    }
    else if (sel4m)
    {
        output = "4m_";     lepton = "muon";        Lepton = "muon";
    }
    else if (sel4e)
    {
        output = "4e_";     lepton = "electron";    Lepton = "Electron";
    }
    output +=  suffix + "_" + id + ".root";


    // Histograms
    const unsigned M  = 15;   TString hname[M],     htitle[M];
    unsigned ZZ = 0;    hname[ZZ] = "4lepMass";     htitle[ZZ] = "4-lepton Mass";
    unsigned T4 = 1;    hname[T4] = "4lepPt";       htitle[T4] = "4-lepton Pt";
    unsigned Z1 = 2;    hname[Z1] = "Z1Mass";       htitle[Z1] = "Z_1 Mass";
    unsigned Z2 = 3;    hname[Z2] = "Z2Mass";       htitle[Z2] = "Z_2 Mass";
    unsigned T1 = 4;    hname[T1] = "Z1Pt";         htitle[T1] = "Z_1 Pt";
    unsigned T2 = 5;    hname[T2] = "Z2Pt";         htitle[T2] = "Z_2 Pt";
    unsigned P1 = 6;    hname[P1] = "Lep1Pt";       htitle[P1] = Lepton + " 1 Pt";
    unsigned P2 = 7;    hname[P2] = "Lep2Pt";       htitle[P2] = Lepton + " 2 Pt";
    unsigned P3 = 8;    hname[P3] = "Lep3Pt";       htitle[P3] = Lepton + " 3 Pt";
    unsigned P4 = 9;    hname[P4] = "Lep4Pt";       htitle[P4] = Lepton + " 4 Pt";
    unsigned E1 = 10;   hname[E1] = "Lep1Eta";      htitle[E1] = Lepton + " 1 Eta";
    unsigned E2 = 11;   hname[E2] = "Lep2Eta";      htitle[E2] = Lepton + " 2 Eta";
    unsigned E3 = 12;   hname[E3] = "Lep3Eta";      htitle[E3] = Lepton + " 3 Eta";
    unsigned E4 = 13;   hname[E4] = "Lep4Eta";      htitle[E4] = Lepton + " 4 Eta";
    unsigned PV = 14;   hname[PV] = "nPV";          htitle[PV] = "# of Primary Vertices";

    Int_t bins[M];      Double_t low[M],    up[M];
    bins[ZZ] = 60;      low[ZZ] = 80;       up[ZZ] = 100;
    bins[ZZ] = 30;      low[ZZ] = 60;       up[ZZ] = 120;
    bins[T4] = 30;      low[T4] = 0;        up[T4] = 120;
    bins[Z1] = 30;      low[Z1] = 0;        up[Z1] = 120;
    bins[Z2] = 30;      low[Z2] = 0;        up[Z2] = 60;
    bins[T1] = 40;      low[T1] = 0;        up[T1] = 160;
    bins[T2] = 30;      low[T2] = 0;        up[T2] = 120;
    bins[P1] = 30;      low[P1] = 0;        up[P1] = 120;
    bins[P2] = 30;      low[P2] = 0;        up[P2] = 90;
    bins[P3] = 30;      low[P3] = 0;        up[P3] = 60;
    bins[P4] = 30;      low[P4] = 0;        up[P4] = 60;
    bins[E1] = 25;      low[E1] = -2.5;     up[E1] = 2.5;
    bins[E2] = 25;      low[E2] = -2.5;     up[E2] = 2.5;
    bins[E3] = 25;      low[E3] = -2.5;     up[E3] = 2.5;
    bins[E4] = 25;      low[E4] = -2.5;     up[E4] = 2.5;
    bins[PV] = 26;      low[PV] = -1;       up[PV] = 51;

    unsigned LL, QT;
    if (sel2l)
    {
        LL = 0;     QT = 1;
        hname[LL] = "DilepMass";    htitle[LL] = "Di" + lepton + " Mass";
        hname[QT] = "DilepPt";      htitle[QT] = "Di" + lepton + " Pt";

        bins[LL] = 30;      low[LL] = 75;       up[LL] = 105;
        bins[QT] = 75;      low[QT] = 0;        up[QT] = 150;
        bins[E1] = 50;
        bins[E2] = 50;
        bins[PV] = 51;      low[PV] = -0.5;     up[PV] = 50.5;
    }
    if (sel4l)
        QT = 1;

    TH1F *h[M];
    for (unsigned j = 0; j < M; j++)
        h[j] = new TH1F(hname[j] + "_" + suffix, htitle[j], bins[j], low[j], up[j]);


    // Read trees
    TH1D *hTotalEvents;
    if (sel2m)
        hTotalEvents = Select2m(path + dir + file, suffix, h[LL], h[QT], 
                h[P1], h[P2], h[E1], h[E2],
                h[PV]);
    else if (sel2e)
        hTotalEvents = Select2e(path + dir + file, suffix, h[LL], h[QT], 
                h[P1], h[P2], h[E1], h[E2],
                h[PV]);
    else if (sel4m)
        hTotalEvents = Select4m(path + dir + file, suffix, h[ZZ], h[QT],
                h[Z1], h[Z2], h[T1], h[T2],
                h[P1], h[P2], h[P3], h[P4], h[E1], h[E2], h[E3], h[E4],
                h[PV]);
    else if (sel4e)
        hTotalEvents = Select4e(path + dir + file, suffix, h[ZZ], h[QT],
                h[Z1], h[Z2], h[T1], h[T2],
                h[P1], h[P2], h[P3], h[P4], h[E1], h[E2], h[E3], h[E4],
                h[PV]);


    // Write to file
    TFile *outFile = new TFile(output, "RECREATE");
    for (unsigned j = 0; j < M; j++)
        h[j]->Write();
    hTotalEvents->Write();
    outFile->Close();


    // Delete everything
    if (kTRUE)
    {
        delete outFile;
        for (unsigned j = 0; j < M; j++)
            delete h[j];
        delete hTotalEvents;
    }
}




TH1D* Select2m(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2,
               TH1F *h_npv)
{
    // Cuts
    Double_t MLL_MIN = 80, MLL_MAX = 100;
    Double_t PT_MIN = 25, LOOSE_PT_MIN = 10;


    // Open root file
    TFile *file = TFile::Open(rootFile);
    TH1D* h_tot;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0);
    TTreeReader reader("tree_" + suffix, file);


    // Branches
    // event
    TTreeReaderValue<Bool_t> passTrigger_(reader, "passTrigger");
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nMuons_(reader, "nMuons");
    TTreeReaderValue<UShort_t> nGoodMuons_(reader, "nStdMuons");

    // muons
    TTreeReaderArray<TLorentzVector> muonP4_(reader, "muonP4");
    TTreeReaderValue<std::vector<Short_t>> muonQ_(reader, "muonQ");
    TTreeReaderValue<std::vector<Float_t>> muonSF_(reader, "muonSF");
    TTreeReaderValue<std::vector<Float_t>> muonIso_(reader, "muonCombIso");
    TTreeReaderValue<std::vector<Bool_t>> muonIsGood_(reader, "muonPassStdCuts");
    TTreeReaderValue<std::vector<Bool_t>> muonIsPF_(reader, "muonIsPF");
    TTreeReaderValue<std::vector<Bool_t>> muonIsGLB_(reader, "muonIsGLB");

    // weighting
    TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
    TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");
    TTreeReaderValue<std::vector<Float_t>> muonIDEff_(reader, "muonIDEff");
    TTreeReaderValue<std::vector<Float_t>> muonIsoEff_(reader, "muonTightIsoEff");
    TTreeReaderValue<std::vector<Float_t>> muonTrigEffData_(reader, "muonTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> muonTrigEffMC_(reader, "muonTriggerEffMC");


    // Read in events
    while (reader.Next())
    {
        // Branch quantity vectors
        Bool_t passTrigger = *passTrigger_;
        UShort_t nPV = *nPV_;
        UShort_t nMuons = *nMuons_, nGoodMuons = *nGoodMuons_;

        vector<TLorentzVector> muonP4;
        vector<Short_t> muonQ;
        vector<Float_t> muonSF, muonIso;
        vector<Bool_t> muonIsGood, muonIsPF, muonIsGLB;

        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        vector<Float_t> muonIDEff, muonIsoEff, muonTrigEffData, muonTrigEffMC;

//      vector<Bool_t> muonIsLoose;
        UShort_t nLooseMuons = 0;


        // Reject events that do not pass trigger or have enough good leptons 
        if (!passTrigger || nGoodMuons != 2)
            continue;

        // Fill vectors
        for (const TLorentzVector& muonP4__: muonP4_)
            muonP4.push_back(muonP4__);
        for (unsigned i = 0; i < nMuons; i++)
        {
            muonQ.push_back((*muonQ_)[i]);
            muonSF.push_back((*muonSF_)[i]);
            muonIso.push_back((*muonIso_)[i]);
            muonIsGood.push_back((*muonIsGood_)[i]);
            muonIsPF.push_back((*muonIsPF_)[i]);
            muonIsGLB.push_back((*muonIsGLB_)[i]);
            muonIDEff.push_back((*muonIDEff_)[i]);
            muonIsoEff.push_back((*muonIsoEff_)[i]);
            muonTrigEffData.push_back((*muonTrigEffData_)[i]);
            muonTrigEffMC.push_back((*muonTrigEffMC_)[i]);
        }

        for (unsigned i = 0; i < nMuons; i++)
        {
            // Rochester corrections
            muonP4[i].SetPtEtaPhiM(muonP4[i].Pt() * muonSF[i], muonP4[i].Eta(), 
                                   muonP4[i].Phi(), muonP4[i].M());

            // Loose muon selection
            if (muonIsPF[i] && muonIsGLB[i] && muonP4[i].Pt() > LOOSE_PT_MIN)
            {
//              muonIsLoose.push_back(kTRUE);
                nLooseMuons++;
            }
//          else
//              muonIsLoose.push_back(kFALSE);
        }


        // Selection
        // Make sure there are no extra muons
        if (nLooseMuons != nGoodMuons)
            continue;

        // Make sure leading muon isn't crappy (change this???)
        if (!muonIsGood[0] || muonP4[0].Pt() < PT_MIN)
            continue;

        // Pair leading muon
        unsigned x = 0;
        for (unsigned j = 1; j < nMuons; j++)
        {
            if (muonQ[0] * muonQ[j] < 0 && muonIsGood[j] && muonP4[j].Pt() > PT_MIN)
            {
                // Mass window requirement
                Double_t mll = (muonP4[0] + muonP4[j]).M();
                if (mll > MLL_MIN && mll < MLL_MAX)
                {
                    x = j;      break;
                }
            }
        }
        if (x == 0)
            continue;

        eventWeight *= PUWeight;
        eventWeight *= muonIDEff[0] * muonIDEff[x];
        eventWeight *= muonIsoEff[0] * muonIsoEff[x];
        if (muonTrigEffMC[0] > 0 || muonTrigEffMC[x] > 0)
        {
            eventWeight *= 1. - (1. - muonTrigEffData[0]) * (1. - muonTrigEffData[x]);
            eventWeight /= 1. - (1. - muonTrigEffMC[0]) * (1. - muonTrigEffMC[x]);
        }

        // Fill histograms
        h_tot->Fill(6);
        h_tot->Fill(7, eventWeight);
        h_mll->Fill((muonP4[0] + muonP4[x]).M(), eventWeight);
        h_qt->Fill((muonP4[0] + muonP4[x]).Pt(), eventWeight);
        h_pt1->Fill(muonP4[0].Pt(), eventWeight);
        h_eta1->Fill(muonP4[0].Eta(), eventWeight);
        h_pt2->Fill(muonP4[x].Pt(), eventWeight);
        h_eta2->Fill(muonP4[x].Eta(), eventWeight);
        h_npv->Fill(nPV, eventWeight);
    }
    file->Close();  delete file;

    return h_tot;
}




TH1D* Select2e(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2, 
               TH1F *h_npv)
{
    // Cuts
    Double_t MLL_MIN = 80, MLL_MAX = 100;
    Double_t PT_MIN = 25;


    // Open root file
    TFile *file = TFile::Open(rootFile);
    TH1D* h_tot;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0);
    TTreeReader reader("tree_" + suffix, file);


    // Branches
    // event
    TTreeReaderValue<Bool_t> passTrigger_(reader, "passTrigger");
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nElecs_(reader, "nElectrons");
    TTreeReaderValue<UShort_t> nGoodElecs_(reader, "nStdElectrons");

    // electrons
    TTreeReaderArray<TLorentzVector> elecP4_(reader, "electronP4");
    TTreeReaderValue<std::vector<Short_t>> elecQ_(reader, "electronQ");
    TTreeReaderValue<std::vector<Float_t>> elecIso_(reader, "electronCombIso");
    TTreeReaderValue<std::vector<Float_t>> elecSF_(reader, "electronSF");
    TTreeReaderValue<std::vector<Bool_t>> elecIsGood_(reader, "electronPassStdCuts");

    // weighting
    TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
    TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");
    TTreeReaderValue<std::vector<Float_t>> elecRecoEff_(reader, "electronRecoEff");
    TTreeReaderValue<std::vector<Float_t>> elecTrigEffData_(reader, "electronTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> elecTrigEffMC_(reader, "electronTriggerEffMC");


    // Read in events
    while (reader.Next())
    {
        // Branch quantity vectors
        Bool_t passTrigger = *passTrigger_;
        UShort_t nPV = *nPV_;
        UShort_t nElecs = *nElecs_, nGoodElecs = *nGoodElecs_;

        vector<TLorentzVector> elecP4;
        vector<Short_t> elecQ;
        vector<Float_t> elecSF, elecIso;
        vector<Bool_t> elecIsGood;

        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        vector<Float_t> elecRecoEff, elecTrigEffData, elecTrigEffMC;


        // Reject events that do not pass trigger or have enough good leptons 
        if (!passTrigger || nGoodElecs != 2)
            continue;

        // Fill vectors
        for (const TLorentzVector& elecP4__: elecP4_)
            elecP4.push_back(elecP4__);
        for (unsigned i = 0; i < nElecs; i++)
        {
            elecQ.push_back((*elecQ_)[i]);
            elecSF.push_back((*elecSF_)[i]);
            elecIso.push_back((*elecIso_)[i]);
            elecIsGood.push_back((*elecIsGood_)[i]);
            elecRecoEff.push_back((*elecRecoEff_)[i]);
            elecTrigEffData.push_back((*elecTrigEffData_)[i]);
            elecTrigEffMC.push_back((*elecTrigEffMC_)[i]);
        }


        // Energy correction
        for (unsigned i = 0; i < nElecs; i++)
            elecP4[i].SetPtEtaPhiM(elecP4[i].Pt() * elecSF[i], elecP4[i].Eta(), 
                                   elecP4[i].Phi(), elecP4[i].M());


        // Selection
        // Make sure leading electron isn't crappy (change this???)
        if (!elecIsGood[0] || elecP4[0].Pt() < PT_MIN)
            continue;

        // Pair leading electron
        unsigned x = 0;
        for (unsigned j = 1; j < nElecs; j++)
        {
            if (elecQ[0] * elecQ[j] < 0 && elecIsGood[j] && elecP4[j].Pt() > PT_MIN)
            {
                // Mass window requirement
                Double_t mll = (elecP4[0] + elecP4[j]).M();
                if (mll > MLL_MIN && mll < MLL_MAX)
                {
                    x = j;      break;
                }
            }
        }
        if (x == 0)
            continue;


        // Event weighting
        eventWeight *= PUWeight;
        eventWeight *= elecRecoEff[0] * elecRecoEff[x];
        if (elecTrigEffMC[0] > 0 || elecTrigEffMC[x] > 0)
        {
            eventWeight *= 1. - (1. - elecTrigEffData[0])*(1. - elecTrigEffData[x]);
            eventWeight /= 1. - (1. - elecTrigEffMC[0])*(1. - elecTrigEffMC[x]);
        }

        // Fill histograms
        h_tot->Fill(6);
        h_tot->Fill(7, eventWeight);
        h_mll->Fill((elecP4[0] + elecP4[x]).M(), eventWeight);
        h_qt->Fill((elecP4[0] + elecP4[x]).Pt(), eventWeight);
        h_pt1->Fill(elecP4[0].Pt(), eventWeight);
        h_eta1->Fill(elecP4[0].Eta(), eventWeight);
        h_pt2->Fill(elecP4[x].Pt(), eventWeight);
        h_eta2->Fill(elecP4[x].Eta(), eventWeight);
        h_npv->Fill(nPV, eventWeight);
    }
    file->Close();  delete file;

    return h_tot;
}




TH1D* Select4m(TString rootFile, TString suffix, TH1F *h_m4l, TH1F *h_qt4l,
               TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
               TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4, 
               TH1F *h_npv)
{
    // Cuts, etc.
    Double_t MZ1_MIN = 12, MZ_MIN = 4, MZ_MAX = 120;
    Double_t M4L_MIN = 80, M4L_MAX = 100;
    Double_t PT_MIN = 5, PT1_MIN = 20, PT2_MIN = 10;
    Double_t DR_MIN = 0.02, MLL_MIN = 4;
    Double_t Z_MASS = 91.2;


    // Open root file
    TFile *file = TFile::Open(rootFile);
    TH1D* h_tot;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0);
    TTreeReader reader("tree_" + suffix, file);


    // Branches
    // event
    TTreeReaderValue<Bool_t> passTrigger_(reader, "passTrigger");
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nMuons_(reader, "nMuons");
    TTreeReaderValue<UShort_t> nGoodMuons_(reader, "nStdMuons");

    // muons
    TTreeReaderArray<TLorentzVector> muonP4_(reader, "muonP4");
    TTreeReaderValue<std::vector<Short_t>> muonQ_(reader, "muonQ");
    TTreeReaderValue<std::vector<Float_t>> muonSF_(reader, "muonSF");
    TTreeReaderValue<std::vector<Float_t>> muonIso_(reader, "muonCombIso");
    TTreeReaderValue<std::vector<Bool_t>> muonIsGood_(reader, "muonPassStdCuts");
    TTreeReaderValue<std::vector<Bool_t>> muonIsPF_(reader, "muonIsPF");
    TTreeReaderValue<std::vector<Bool_t>> muonIsGLB_(reader, "muonIsGLB");

    // weighting
    TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
    TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");
    TTreeReaderValue<std::vector<Float_t>> muonIDEff_(reader, "muonIDEff");
    TTreeReaderValue<std::vector<Float_t>> muonIsoEff_(reader, "muonTightIsoEff");
    TTreeReaderValue<std::vector<Float_t>> muonTrigEffData_(reader, "muonTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> muonTrigEffMC_(reader, "muonTriggerEffMC");


    // Read in events
    while (reader.Next())
    {
        // Branch quantity vectors
        Bool_t passTrigger = *passTrigger_;
        UShort_t nPV = *nPV_;
        UShort_t nMuons = *nMuons_, nGoodMuons = *nGoodMuons_;

        vector<TLorentzVector> muonP4;
        vector<Short_t> muonQ;
        vector<Float_t> muonSF, muonIso;
        vector<Bool_t> muonIsGood, muonIsPF, muonIsGLB;

        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        vector<Float_t> muonIDEff, muonIsoEff, muonTrigEffData, muonTrigEffMC;

        vector<Bool_t> muonIsLoose;
        UShort_t nLooseMuons = 0;


        // Veto events that are not triggered or do not have enough leptons 
        if (!passTrigger || nMuons < 4)
            continue;


        // Fill vectors
        for (const TLorentzVector& muonP4__: muonP4_)
            muonP4.push_back(muonP4__);
        for (unsigned i = 0; i < nMuons; i++)
        {
            muonQ.push_back((*muonQ_)[i]);
            muonSF.push_back((*muonSF_)[i]);
            muonIso.push_back((*muonIso_)[i]);
            muonIsGood.push_back((*muonIsGood_)[i]);
            muonIsPF.push_back((*muonIsPF_)[i]);
            muonIsGLB.push_back((*muonIsGLB_)[i]);
            muonIDEff.push_back((*muonIDEff_)[i]);
            muonIsoEff.push_back((*muonIsoEff_)[i]);
            muonTrigEffData.push_back((*muonTrigEffData_)[i]);
            muonTrigEffMC.push_back((*muonTrigEffMC_)[i]);
        }

        for (unsigned i = 0; i < nMuons; i++)
        {
            // Rochester corrections
            muonP4[i].SetPtEtaPhiM(muonP4[i].Pt() * muonSF[i], muonP4[i].Eta(), 
                                   muonP4[i].Phi(), muonP4[i].M());

            // Loose muon selection
            if (muonIsPF[i] && muonIsGLB[i])
            {
                muonIsLoose.push_back(kTRUE);
                nLooseMuons++;
            }
            else
                muonIsLoose.push_back(kFALSE);
        }


        // Make sure there are no gross muons
        if (nLooseMuons < 4)
            continue;


        // Selection
        // Find Z candidates
        vector<std::pair<unsigned, unsigned>> zCands;
        for (unsigned i = 0; i < nMuons-1; i++)
        {
            if (!muonIsLoose[i] || muonP4[i].Pt() < PT_MIN)
                continue;
            for (unsigned j = i+1; j < nMuons; j++)
            {
                if (!muonIsLoose[j] || muonP4[j].Pt() < PT_MIN || muonQ[i] * muonQ[j] > 0)
                    continue;

                Double_t mll = (muonP4[i] + muonP4[j]).M();
                if (mll > MZ_MIN && mll < MZ_MAX)
                    zCands.push_back(std::make_pair(i, j));
            }
        }
        if (zCands.size() < 2)
            continue;
  

        // Find ZZ candidates
        vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned>>> zzCands;
        for (unsigned i = 0; i < zCands.size()-1; i++)
        {
            for (unsigned j = i+1; j < zCands.size(); j++)
            {
                std::pair<unsigned, unsigned> Zi = zCands[i], Zj = zCands[j];

                // Make sure pairs do not overlap
                if (Zi.first == Zj.first || Zi.first == Zj.second
                    || Zi.second == Zj.first || Zi.second == Zj.second)
                    continue;
  
                Double_t mll_i = (muonP4[Zi.first] + muonP4[Zi.second]).M();
                Double_t mll_j = (muonP4[Zj.first] + muonP4[Zj.second]).M();

                // Choose Z1 closest to nominal Z mass
                std::pair<unsigned, unsigned> Z1, Z2;
                if (fabs(mll_i - Z_MASS) < fabs(mll_j - Z_MASS))
                {
                    Z1 = Zi;    Z2 = Zj;
                }
                else
                {
                    Z1 = Zj;    Z2 = Zi;
                }

                std::pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
                Z1_p4s = std::make_pair(muonP4[Z1.first], muonP4[Z1.second]);
                Z2_p4s = std::make_pair(muonP4[Z2.first], muonP4[Z2.second]);
                
                std::pair<Int_t, Int_t> Z1_qs, Z2_qs;
                Z1_qs = std::make_pair(muonQ[Z1.first], muonQ[Z1.second]);
                Z2_qs = std::make_pair(muonQ[Z2.first], muonQ[Z2.second]);

                TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
                Z1_p4 = Z1_p4s.first + Z1_p4s.second;
                Z2_p4 = Z2_p4s.first + Z2_p4s.second;
                ZZ_p4 = Z1_p4 + Z2_p4;


                // Mass requirements
                if (Z1_p4.M() < MZ1_MIN || ZZ_p4.M() < M4L_MIN || ZZ_p4.M() > M4L_MAX)
                    continue;


                // Smart cut
                std::pair<unsigned, unsigned> Zx, Zy, Za, Zb;
                if (Z1_qs.first * Z2_qs.first > 0)
                {
                    Zx = std::make_pair(Z1.first, Z2.second);
                    Zy = std::make_pair(Z2.first, Z1.second);
                }
                else if (Z1_qs.first * Z2_qs.second > 0)
                {
                    Zx = std::make_pair(Z1.first, Z2.first);
                    Zy = std::make_pair(Z1.second, Z2.second);
                }

                Double_t mll_x, mll_y, mll_a, mll_b;
                mll_x = (muonP4[Zx.first] + muonP4[Zx.second]).M();
                mll_y = (muonP4[Zy.first] + muonP4[Zy.second]).M();

                if (fabs(mll_x - Z_MASS) < fabs(mll_y - Z_MASS))
                {
                    Za = Zx;        Zb = Zy;
                    mll_a = mll_x;  mll_b = mll_y;
                }
                else
                {
                    Za = Zy;        Zb = Zx;
                    mll_a = mll_y;  mll_b = mll_x;
                }

                if (fabs(mll_a - Z_MASS) < fabs(Z1_p4.M() - Z_MASS) && mll_b < MLL_MIN)
                    continue;


                // Pt2 requirement
                vector<unsigned> z = {Z1.first, Z1.second, Z2.first, Z2.second};
                std::sort(z.begin(), z.end());
                vector<TLorentzVector> selMuonP4;
                vector<Short_t> selMuonQ;

                for (unsigned k = 0; k < z.size(); k++)
                {
                    selMuonP4.push_back(muonP4[z[k]]);
                    selMuonQ.push_back(muonQ[z[k]]);
                }

                if (selMuonP4[0].Pt() < PT1_MIN || selMuonP4[1].Pt() < PT2_MIN)
                    continue;


                // Ghost removal
                Bool_t foundGhost = kFALSE;
                for (unsigned x = 0; x < 4-1; x++)
                {
                    for (unsigned y = x+1; y < 4; y++)
                    {
                        if (selMuonP4[x].DeltaR(selMuonP4[y]) < DR_MIN)
                        {
                            foundGhost = kTRUE;
                            break;
                        }
                    }
                    if (foundGhost)
                        break;
                }
                if (foundGhost)
                    continue;


                // QCD suppression
                Bool_t foundQCD = kFALSE;
                for (unsigned x = 0; x < 4-1; x++)
                {
                    for (unsigned y = x+1; y < 4; y++)
                    {
                        if (selMuonQ[x] * selMuonQ[y] > 0)
                            continue;

                        if ((selMuonP4[x] + selMuonP4[y]).M() < MLL_MIN)
                        {
                            foundQCD = kTRUE;
                            break;
                        }
                    }
                    if (foundQCD)
                        break;
                }
                if (foundQCD)
                    continue;
  

                zzCands.push_back(std::make_pair(Z1, Z2));
            }
        }
        if (zzCands.size() < 1)
            continue;


        // Select best ZZ candidate


        // Fill info
        std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned>> ZZ = zzCands[0];
        std::pair<unsigned, unsigned> Z1 = ZZ.first, Z2 = ZZ.second;
        std::pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
        Z1_p4s = std::make_pair(muonP4[Z1.first], muonP4[Z1.second]);
        Z2_p4s = std::make_pair(muonP4[Z2.first], muonP4[Z2.second]);
        
        TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
        Z1_p4 = Z1_p4s.first + Z1_p4s.second;
        Z2_p4 = Z2_p4s.first + Z2_p4s.second;
        ZZ_p4 = Z1_p4 + Z2_p4;

        vector<unsigned> z = {Z1.first, Z1.second, Z2.first, Z2.second};


        // Event weighting
        unsigned x = z[0], y = z[1];
        eventWeight *= PUWeight;
        eventWeight *= muonIDEff[x] * muonIDEff[y];
        eventWeight *= muonIsoEff[x] * muonIsoEff[y];
        if (muonTrigEffMC[x] > 0 || muonTrigEffMC[y] > 0)
        {
            eventWeight *= 1. - (1. - muonTrigEffData[x]) * (1. - muonTrigEffData[y]);
            eventWeight /= 1. - (1. - muonTrigEffMC[x]) * (1. - muonTrigEffMC[y]);
        }


        // Fill histograms
        std::sort(z.begin(), z.end());

        h_tot->Fill(6);
        h_tot->Fill(7, eventWeight);
        h_m4l->Fill(ZZ_p4.M(), eventWeight);
        h_qt4l->Fill(ZZ_p4.Pt(), eventWeight);
  
        h_m12->Fill(Z1_p4.M(), eventWeight);
        h_m34->Fill(Z2_p4.M(), eventWeight);
        h_qt12->Fill(Z1_p4.Pt(), eventWeight);
        h_qt34->Fill(Z2_p4.Pt(), eventWeight);

        h_pt1->Fill(muonP4[z[0]].Pt(), eventWeight);
        h_pt2->Fill(muonP4[z[1]].Pt(), eventWeight);
        h_pt3->Fill(muonP4[z[2]].Pt(), eventWeight);
        h_pt4->Fill(muonP4[z[3]].Pt(), eventWeight);
  
        h_eta1->Fill(muonP4[z[0]].Eta(), eventWeight);
        h_eta2->Fill(muonP4[z[1]].Eta(), eventWeight);
        h_eta3->Fill(muonP4[z[2]].Eta(), eventWeight);
        h_eta4->Fill(muonP4[z[3]].Eta(), eventWeight);
    
        h_npv->Fill(nPV, eventWeight);
    }
    file->Close();  delete file;

    return h_tot;
}




TH1D* Select4e(TString rootFile, TString suffix, TH1F *h_m4l, TH1F *h_qt4l,
               TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
               TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4, 
               TH1F *h_npv)
{
    // Cuts, etc.
    Double_t MZ1_MIN = 12, MZ_MIN = 4, MZ_MAX = 120;
    Double_t M4L_MIN = 80, M4L_MAX = 100;
    Double_t PT_MIN = 5, PT1_MIN = 20, PT2_MIN = 10;
    Double_t DR_MIN = 0.02, MLL_MIN = 4;
    Double_t Z_MASS = 91.2;


    // Open root file
    TFile *file = TFile::Open(rootFile);
    TH1D* h_tot;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0);
    TTreeReader reader("tree_" + suffix, file);


    // Branches
    // event
    TTreeReaderValue<Bool_t> passTrigger_(reader, "passTrigger");
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nElecs_(reader, "nElectrons");
    TTreeReaderValue<UShort_t> nGoodElecs_(reader, "nStdElectrons");

    // electrons
    TTreeReaderArray<TLorentzVector> elecP4_(reader, "electronP4");
    TTreeReaderValue<std::vector<Short_t>> elecQ_(reader, "electronQ");
    TTreeReaderValue<std::vector<Float_t>> elecIso_(reader, "electronCombIso");
    TTreeReaderValue<std::vector<Float_t>> elecSF_(reader, "electronSF");
    TTreeReaderValue<std::vector<Bool_t>> elecIsGood_(reader, "electronPassStdCuts");

    // weighting
    TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
    TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");
    TTreeReaderValue<std::vector<Float_t>> elecRecoEff_(reader, "electronRecoEff");
    TTreeReaderValue<std::vector<Float_t>> elecTrigEffData_(reader, "electronTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> elecTrigEffMC_(reader, "electronTriggerEffMC");


    // Read in events
    Double_t nSel = 0, nSelWeight = 0;
    while (reader.Next())
    {
        // Branch quantity vectors
        Bool_t passTrigger = *passTrigger_;
        UShort_t nPV = *nPV_;
        UShort_t nElecs = *nElecs_, nGoodElecs = *nGoodElecs_;

        vector<TLorentzVector> elecP4;
        vector<Short_t> elecQ;
        vector<Float_t> elecSF, elecIso;
        vector<Bool_t> elecIsGood;

        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        vector<Float_t> elecRecoEff, elecTrigEffData, elecTrigEffMC;


        // Veto events that are not triggered or do not have enough leptons 
        if (!passTrigger || nElecs < 4)
            continue;


        // Fill vectors
        for (const TLorentzVector& elecP4__: elecP4_)
            elecP4.push_back(elecP4__);
        for (unsigned i = 0; i < nElecs; i++)
        {
            elecQ.push_back((*elecQ_)[i]);
            elecSF.push_back((*elecSF_)[i]);
            elecIso.push_back((*elecIso_)[i]);
            elecIsGood.push_back((*elecIsGood_)[i]);
            elecRecoEff.push_back((*elecRecoEff_)[i]);
            elecTrigEffData.push_back((*elecTrigEffData_)[i]);
            elecTrigEffMC.push_back((*elecTrigEffMC_)[i]);
        }


        // Energy correction
        for (unsigned i = 0; i < nElecs; i++)
            elecP4[i].SetPtEtaPhiM(elecP4[i].Pt() * elecSF[i], elecP4[i].Eta(), 
                                   elecP4[i].Phi(), elecP4[i].M());

  
        // Selection
        // Find Z candidates
        vector<std::pair<unsigned, unsigned>> zCands;
        for (unsigned i = 0; i < nElecs-1; i++)
        {
            for (unsigned j = i+1; j < nElecs; j++)
            {
                if (elecQ[i] * elecQ[j] < 0 && elecP4[i].Pt() > PT_MIN && elecP4[j].Pt() > PT_MIN)
                {
                    Double_t mll = (elecP4[0] + elecP4[j]).M();
                    if (mll > MZ_MIN && mll < MZ_MAX)
                        zCands.push_back(std::make_pair(i, j));
                }
            }
        }
        if (zCands.size() < 2)
            continue;

        // Find ZZ candidates
        vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned>>> zzCands;
        for (unsigned i = 0; i < zCands.size()-1; i++)
        {
            for (unsigned j = i+1; j < zCands.size(); j++)
            {
                std::pair<unsigned, unsigned> Zi = zCands[i], Zj = zCands[j];

                // Make sure pairs do not overlap
                if (Zi.first == Zj.first || Zi.first == Zj.second
                    || Zi.second == Zj.first || Zi.second == Zj.second)
                    continue;

                Double_t mll_i = (elecP4[Zi.first] + elecP4[Zi.second]).M();
                Double_t mll_j = (elecP4[Zj.first] + elecP4[Zj.second]).M();

                // Choose Z1 closest to nominal Z mass
                std::pair<unsigned, unsigned> Z1, Z2;
                if (fabs(mll_i - Z_MASS) < fabs(mll_j - Z_MASS))
                {
                    Z1 = Zi;    Z2 = Zj;
                }
                else
                {
                    Z1 = Zj;    Z2 = Zi;
                }
  
                std::pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
                Z1_p4s = std::make_pair(elecP4[Z1.first], elecP4[Z1.second]);
                Z2_p4s = std::make_pair(elecP4[Z2.first], elecP4[Z2.second]);
                
                std::pair<Int_t, Int_t> Z1_qs, Z2_qs;
                Z1_qs = std::make_pair(elecQ[Z1.first], elecQ[Z1.second]);
                Z2_qs = std::make_pair(elecQ[Z2.first], elecQ[Z2.second]);

                TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
                Z1_p4 = Z1_p4s.first + Z1_p4s.second;
                Z2_p4 = Z2_p4s.first + Z2_p4s.second;
                ZZ_p4 = Z1_p4 + Z2_p4;

                vector<unsigned> z = {Z1.first, Z1.second, Z2.first, Z2.second};
                std::sort(z.begin(), z.end());
  
                // Mass requirements
                if (Z1_p4.M() < MZ1_MIN || ZZ_p4.M() < M4L_MIN || ZZ_p4.M() > M4L_MAX)
                    continue;


                // Smart cut
                std::pair<unsigned, unsigned> Zx, Zy, Za, Zb;
                if (Z1_qs.first * Z2_qs.first > 0)
                {
                    Zx = std::make_pair(Z1.first, Z2.second);
                    Zy = std::make_pair(Z2.first, Z1.second);
                }
                else if (Z1_qs.first * Z2_qs.second > 0)
                {
                    Zx = std::make_pair(Z1.first, Z2.first);
                    Zy = std::make_pair(Z1.second, Z2.second);
                }

                Double_t mll_x, mll_y, mll_a, mll_b;
                mll_x = (elecP4[Zx.first] + elecP4[Zx.second]).M();
                mll_y = (elecP4[Zy.first] + elecP4[Zy.second]).M();

                if (fabs(mll_x - Z_MASS) < fabs(mll_y - Z_MASS))
                {
                    Za = Zx;        Zb = Zy;
                    mll_a = mll_x;  mll_b = mll_y;
                }
                else
                {
                    Za = Zy;        Zb = Zx;
                    mll_a = mll_y;  mll_b = mll_x;
                }

                if (fabs(mll_a - Z_MASS) < fabs(Z1_p4.M() - Z_MASS) && mll_b < MLL_MIN)
                    continue;


                // Pt1, Pt2 requirement
                vector<TLorentzVector> selElecP4;
                vector<Short_t> selElecQ;

                for (unsigned k = 0; k < z.size(); k++)
                {
                    selElecP4.push_back(elecP4[z[k]]);
                    selElecQ.push_back(elecQ[z[k]]);
                }

                if (selElecP4[0].Pt() < PT1_MIN || selElecP4[1].Pt() < PT2_MIN)
                    continue;


                // Ghost removal
                Bool_t foundGhost = kFALSE;
                for (unsigned x = 0; x < 4-1; x++)
                {
                    for (unsigned y = x+1; y < 4; y++)
                    {
                        if (selElecP4[x].DeltaR(selElecP4[y]) < DR_MIN)
                        {
                            foundGhost = kTRUE;
                            break;
                        }
                    }
                    if (foundGhost)
                        break;
                }
                if (foundGhost)
                    continue;


                // QCD suppression
                Bool_t foundQCD = kFALSE;
                for (unsigned x = 0; x < 4-1; x++)
                {
                    for (unsigned y = x+1; y < 4; y++)
                    {
                        if (selElecQ[x] * selElecQ[y] > 0)
                            continue;

                        if ((selElecP4[x] + selElecP4[y]).M() < MLL_MIN)
                        {
                            foundQCD = kTRUE;
                            break;
                        }
                    }
                    if (foundQCD)
                        break;
                }
                if (foundQCD)
                    continue;
                // Pt2 requirement, others
  
                zzCands.push_back(std::make_pair(Z1, Z2));
            }
        }
        if (zzCands.size() < 1)
            continue;


        // Choose best ZZ cand


        // Fill
        std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned>> ZZ = zzCands[0];
        std::pair<unsigned, unsigned> Z1 = ZZ.first, Z2 = ZZ.second;
        std::pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
        Z1_p4s = std::make_pair(elecP4[Z1.first], elecP4[Z1.second]);
        Z2_p4s = std::make_pair(elecP4[Z2.first], elecP4[Z2.second]);
        
        TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
        Z1_p4 = Z1_p4s.first + Z1_p4s.second;
        Z2_p4 = Z2_p4s.first + Z2_p4s.second;
        ZZ_p4 = Z1_p4 + Z2_p4;

        vector<unsigned> z = {Z1.first, Z1.second, Z2.first, Z2.second};
        std::sort(z.begin(), z.end());


        // Event weighting
        unsigned x = Z1.first, y = Z1.second;
        eventWeight *= PUWeight;
        eventWeight *= elecRecoEff[x] * elecRecoEff[y];
        if (elecTrigEffMC[x] > 0 || elecTrigEffMC[y] > 0)
        {
            eventWeight *= 1. - (1. - elecTrigEffData[x])*(1. - elecTrigEffData[y]);
            eventWeight /= 1. - (1. - elecTrigEffMC[x])*(1. - elecTrigEffMC[y]);
        }

  
        // Fill histograms
        h_tot->Fill(6);
//      nSel++;
  
//      nSelWeight += eventWeight;
        h_tot->Fill(7, eventWeight);
        h_m4l->Fill(ZZ_p4.M(), eventWeight);
        h_qt4l->Fill(ZZ_p4.Pt(), eventWeight);
  
        h_m12->Fill(Z1_p4.M(), eventWeight);
        h_m34->Fill(Z2_p4.M(), eventWeight);
        h_qt12->Fill(Z1_p4.Pt(), eventWeight);
        h_qt34->Fill(Z2_p4.Pt(), eventWeight);
  
        h_pt1->Fill(elecP4[z[0]].Pt(), eventWeight);
        h_pt2->Fill(elecP4[z[1]].Pt(), eventWeight);
        h_pt3->Fill(elecP4[z[2]].Pt(), eventWeight);
        h_pt4->Fill(elecP4[z[3]].Pt(), eventWeight);
  
        h_eta1->Fill(elecP4[z[0]].Eta(), eventWeight);
        h_eta2->Fill(elecP4[z[1]].Eta(), eventWeight);
        h_eta3->Fill(elecP4[z[2]].Eta(), eventWeight);
        h_eta4->Fill(elecP4[z[3]].Eta(), eventWeight);
  
        h_npv->Fill(nPV, eventWeight);
    }
    file->Close();  delete file;

//  h_tot->SetBinContent(6, nSel);
//  h_tot->SetBinContent(7, nSelWeight);

    TH1D *h_tot_ = new TH1D("TotalEvents_" + suffix, h_tot->GetTitle(), h_tot->GetNbinsX(), 
                            h_tot->GetXaxis()->GetXmin(), h_tot->GetXaxis()->GetXmax());
    for (unsigned i = 0; i < h_tot->GetNbinsX(); i++)
        h_tot_->SetBinContent(i, h_tot->GetBinContent(i));

    return h_tot_;
}
