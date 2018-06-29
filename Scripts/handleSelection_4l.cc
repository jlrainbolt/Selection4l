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



TH1D* Select2l(const TString rootFile, const TString suffix,
               const Bool_t mumu,
               TH1F *h_mll, TH1F *h_qt,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2,
               TH1F *h_npv);
  
TH1D* Select4l(const TString rootFile, const TString suffix,
               const Bool_t muPair1, const Bool_t muPair2,
               TH1F *h_m4l, TH1F *h_qt4l,
               TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
               TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4,
               TH1F *h_npv, TH1F* h_nzz, TH1F* h_nlep);
  

using namespace std;




///////////////////
// MAIN FUNCTION //
///////////////////

void handleSelection(const TString selection, const TString suffix, const TString id)
{
    /* OPTIONS */
    Bool_t fineBinning = kTRUE, medBinning = kFALSE;



    /* FILE INFO */
    // Parse selection
    Bool_t sel2l = kFALSE, sel4l = kFALSE;
    Bool_t muPair1, muPair2;    // TRUE for muon, FALSE for electron
    if (selection == "mumu" || selection == "2m")
    {
        sel2l = kTRUE;      muPair1 = kTRUE;
    }
    else if (selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;      muPair1 = kFALSE;
    }
    else if (selection == "4mu" || selection == "4m")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;    muPair2 = kTRUE;
    }
    else if (selection == "4e")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;   muPair2 = kFALSE;
    }

    if (sel2l)
        muPair2 = muPair1;


    // Path to directory where ALL trees are stored
    TString path = "root://cmsxrootd.fnal.gov//store/user/jrainbol/Trees/2016/";


    // Names of ROOT files and trees for each sample
    TString dir     = suffix + "/";
    TString file    = suffix + "_" + id + ".root";


    // Name of output file
    TString output = selection + "_" + suffix + "_" + id + ".root";



    /* HISTOGRAMS */
    // Lepton names 
    TString Lep12, Lep34, Lep, lep12, lep34, lep;
    Lep12   = muPair1 ? "Muon" : "Electron";        lep12   = muPair1 ? "muon" : "electron";
    Lep34   = muPair2 ? "Muon" : "Electron";        lep34   = muPair2 ? "muon" : "electron";
    Lep     = (muPair1 == muPair2) ? Lep12 : "Lepton";
    lep     = (muPair1 == muPair2) ? lep12 : "lepton";


    // Hist indices     // Names                    // Titles
    const unsigned M=17;    TString hname[M],       htitle[M];
    unsigned ZZ = 0;    hname[ZZ] = "4lepMass";     htitle[ZZ] = "4-" + lep + " Mass";
    unsigned T4 = 1;    hname[T4] = "4lepPt";       htitle[T4] = "4-" + lep + "Pt";
    unsigned Z1 = 2;    hname[Z1] = "Z1Mass";       htitle[Z1] = "Z_{1} (Di" + lep12 + ") Mass";
    unsigned Z2 = 3;    hname[Z2] = "Z2Mass";       htitle[Z2] = "Z_{2} (Di" + lep34 + ") Mass";
    unsigned T1 = 4;    hname[T1] = "Z1Pt";         htitle[T1] = "Z_{1} (Di" + lep12 + ") Pt";
    unsigned T2 = 5;    hname[T2] = "Z2Pt";         htitle[T2] = "Z_{2} (Di" + lep34 + ") Pt";
    unsigned P1 = 6;    hname[P1] = "Lep1Pt";       htitle[P1] = Lep + " 1 Pt";
    unsigned P2 = 7;    hname[P2] = "Lep2Pt";       htitle[P2] = Lep + " 2 Pt";
    unsigned P3 = 8;    hname[P3] = "Lep3Pt";       htitle[P3] = Lep + " 3 Pt";
    unsigned P4 = 9;    hname[P4] = "Lep4Pt";       htitle[P4] = Lep + " 4 Pt";
    unsigned E1 = 10;   hname[E1] = "Lep1Eta";      htitle[E1] = Lep + " 1 Eta";
    unsigned E2 = 11;   hname[E2] = "Lep2Eta";      htitle[E2] = Lep + " 2 Eta";
    unsigned E3 = 12;   hname[E3] = "Lep3Eta";      htitle[E3] = Lep + " 3 Eta";
    unsigned E4 = 13;   hname[E4] = "Lep4Eta";      htitle[E4] = Lep + " 4 Eta";
    unsigned PV = 14;   hname[PV] = "nPV";          htitle[PV] = "# of Primary Vertices";
    unsigned NC = 15;   hname[NC] = "nZZCands";     htitle[NC] = "# of ZZ Candidates";
    unsigned NL = 16;   hname[NL] = "nSelLeps";     htitle[NL] = "# of Selected " + Lep + "s";


    // Binning
    Int_t bins[M];      Double_t low[M],    up[M];
    bins[ZZ] = 30;      low[ZZ] = 60;       up[ZZ] = 180;
    bins[T4] = 15;      low[T4] = 0;        up[T4] = 120;
    bins[Z1] = 15;      low[Z1] = 0;        up[Z1] = 120;
    bins[Z2] = 15;      low[Z2] = 0;        up[Z2] = 60;
    bins[T1] = 15;      low[T1] = 0;        up[T1] = 150;
    bins[T2] = 15;      low[T2] = 0;        up[T2] = 150;
    bins[P1] = 15;      low[P1] = 0;        up[P1] = 150;
    bins[P2] = 15;      low[P2] = 0;        up[P2] = 90;
    bins[P3] = 15;      low[P3] = 0;        up[P3] = 60;
    bins[P4] = 15;      low[P4] = 0;        up[P4] = 30;
    bins[E1] = 13;      low[E1] = -2.6;     up[E1] = 2.6;
    bins[E2] = 13;      low[E2] = -2.6;     up[E2] = 2.6;
    bins[E3] = 13;      low[E3] = -2.6;     up[E3] = 2.6;
    bins[E4] = 13;      low[E4] = -2.6;     up[E4] = 2.6;
    bins[PV] = 13;      low[PV] = -1;       up[PV] = 51;
    bins[NC] = 10;      low[NC] = 0.5;      up[NC] = 10.5;
    bins[NL] = 10;      low[NL] = 0.5;      up[NL] = 10.5;

    if (fineBinning)
    {
        bins[ZZ] = 235;     low[ZZ] = 60;       up[ZZ] = 1000;
        bins[T4] = 250;     low[T4] = 0;        up[T4] = 1000;
        bins[Z1] = 150;     bins[Z2] = 150;     bins[T1] = 150;     bins[T2] = 150;
        bins[P1] = 150;     bins[P2] = 150;     bins[P3] = 150;     bins[P4] = 150;
        bins[E1] = 130;     bins[E2] = 130;     bins[E3] = 130;     bins[E4] = 130;
        bins[PV] = 130;
    }
    else if (medBinning)
    {
        bins[ZZ] = 60;      low[ZZ] = 60;       up[ZZ] = 180;
        bins[T4] = 60;      low[T4] = 0;        up[T4] = 180;
        bins[Z1] = 60;      bins[Z2] = 60;      bins[T1] = 60;      bins[T2] = 60;
        bins[P1] = 60;      bins[P2] = 60;      bins[P3] = 60;      bins[P4] = 60;
        bins[E1] = 52;      bins[E2] = 52;      bins[E3] = 52;      bins[E4] = 52;
        bins[PV] = 52;
    }

    unsigned LL, QT;
    if (sel2l)
    {
        LL = 0;     QT = 1;
        hname[LL] = "DilepMass";    htitle[LL] = "Di" + lep12 + " Mass";
        hname[QT] = "DilepPt";      htitle[QT] = "Di" + lep12 + " Pt";

        bins[LL] = 60;      low[LL] = 75;       up[LL] = 105;
        bins[QT] = 75;      low[QT] = 0;        up[QT] = 150;
        bins[P1] = 100;     low[P1] = 20;       up[P1] = 120;
        bins[P2] = 100;     low[P2] = 20;       up[P2] = 120;
        bins[E1] = 50;      low[E1] = -2.5;     up[E1] = 2.5;
        bins[E2] = 50;      low[E2] = -2.5;     up[E2] = 2.5;
        bins[PV] = 51;      low[PV] = -0.5;     up[PV] = 50.5;
    }
    else if (sel4l)
        QT = 1;


    // Make histograms
    TH1F *h[M];
    for (unsigned j = 0; j < M; j++)
        h[j] = new TH1F(hname[j] + "_" + suffix, htitle[j], bins[j], low[j], up[j]);



    /* SELECTION */
    // Read trees
    TH1D *hTotalEvents;
    if (sel2l)
        hTotalEvents = Select2l(path + dir + file, suffix, muPair1,
                                h[LL], h[QT], 
                                h[P1], h[P2], h[E1], h[E2],
                                h[PV]);
    else if (sel4l)
        hTotalEvents = Select4l(path + dir + file, suffix, muPair1, muPair2,
                                h[ZZ], h[QT],
                                h[Z1], h[Z2], h[T1], h[T2],
                                h[P1], h[P2], h[P3], h[P4], h[E1], h[E2], h[E3], h[E4],
                                h[PV], h[NC], h[NL]);
 

    // Write to file
    TFile *outFile = new TFile(output, "RECREATE");
    for (unsigned j = 0; j < M; j++)
        h[j]->Write();
    hTotalEvents->Write();
    outFile->Close();



    /* CLEANUP */
    delete outFile;
    for (unsigned j = 0; j < M; j++)
        delete h[j];
    delete hTotalEvents;

    return;
}




/////////////////////////
// SELECTION FUNCTIONS //
/////////////////////////



/* DILEPTON */
TH1D* Select2l(const TString rootFile, const TString suffix, const Bool_t mumu,
               TH1F *h_mll, TH1F *h_qt,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2,
               TH1F *h_npv)
{
    ////////////////////
    // INITIALIZATION //
    ////////////////////

    /* CUTS & CONSTANTS */
    Float_t MLL_MIN = 80, MLL_MAX = 100;
    Float_t TIGHT_PT_MIN = 25, LOOSE_PT_MIN = 10;

    Double_t LEP_MASS = mumu ? 0.105658369 : 0.000511;



    /* ROOT FILE */
    TFile *file = TFile::Open(rootFile);
    TH1D* h_tot;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0);
    TTreeReader reader("tree_" + suffix, file);



    /* BRANCHES */
    TString Lep = mumu ? "Muon" : "Electron";
    TString lep = mumu ? "muon" : "electron";
    TString BadLep = mumu ? "Electron" : "Muon";

    // Event
    TTreeReaderValue<Bool_t> evtGoodTrigger_(reader, "evt"+Lep+"Triggered");
    TTreeReaderValue<Bool_t> evtBadTrigger_(reader, "evt"+BadLep+"Triggered");
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nLeps_(reader, "n"+Lep+"s");
    TTreeReaderValue<UShort_t> nTightLeps_(reader, "nTight"+Lep+"s");
    TTreeReaderValue<UShort_t> nLooseLeps_(reader, "nLoose"+Lep+"s");


    // Leptons
    TTreeReaderArray<TLorentzVector> lepP4_(reader, lep+"P4");
    TTreeReaderValue<vector<Short_t>> lepQ_(reader, lep+"Q");
    TTreeReaderValue<vector<Bool_t>> lepIsTight_(reader, lep+"IsTight");
    TTreeReaderValue<vector<Bool_t>> lepIsLoose_(reader, lep+"IsLoose");


    // Scaling & weighting
    TTreeReaderValue<vector<Float_t>> lepSF_(reader, lep+"SF");
    TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
    TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");

    TTreeReaderValue<vector<Float_t>> lepIDEff_(reader, "muonTightIDEff");
    TTreeReaderValue<vector<Float_t>> lepIsoEff_(reader, "muonTightIsoEff");
    TTreeReaderValue<vector<Float_t>> lepRecoEff_(reader, "electronRecoEff");

    TTreeReaderValue<vector<Float_t>> lepTrigEffData_(reader, lep+"TriggerEffData");
    TTreeReaderValue<vector<Float_t>> lepTrigEffMC_(reader, lep+"TriggerEffMC");
    TTreeReaderValue<vector<Bool_t>> lepTriggered_(reader, lep+"Triggered");




    ////////////////
    // EVENT LOOP //
    ////////////////

    while (reader.Next())
    {
        /* BRANCHES */
        Bool_t passTrigger = (*evtGoodTrigger_ && !*evtBadTrigger_);
        UShort_t nPV = *nPV_;
        UShort_t nLeps = *nLeps_, nTightLeps = *nTightLeps_, nLooseLeps = *nLooseLeps_;

        vector<TLorentzVector> lepP4;
        vector<Short_t> lepQ;
        vector<Bool_t> lepIsTight, lepIsLoose;

        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        vector<Float_t> lepSF, lepRecoEff, lepTrigEffData, lepTrigEffMC;
        vector<Bool_t> lepTriggered;



        /* PRESELECTION */
        if (!passTrigger || nTightLeps < 2)
            continue;



        /* LEPTONS */
        for (const TLorentzVector& lepP4__: lepP4_)
            lepP4.push_back(lepP4__);
        for (unsigned i = 0; i < nLeps; i++)
        {
            lepQ.push_back((*lepQ_)[i]);
            lepIsTight.push_back((*lepIsTight_)[i]);
            lepIsLoose.push_back((*lepIsLoose_)[i]);

            lepSF.push_back((*lepSF_)[i]);
            if (mumu)
                lepRecoEff.push_back((*lepIDEff_)[i] * (*lepIsoEff_)[i]);
            else
                lepRecoEff.push_back((*lepRecoEff_)[i]);
            lepTrigEffData.push_back((*lepTrigEffData_)[i]);
            lepTrigEffMC.push_back((*lepTrigEffMC_)[i]);
            lepTriggered.push_back((*lepTriggered_)[i]);
        }


        // Adjust lepton requirements
        for (unsigned i = 0; i < nLeps; i++)
        {
            // Apply Pt scaling
            lepP4[i].SetPtEtaPhiM(lepP4[i].Pt()*lepSF[i], lepP4[i].Eta(), lepP4[i].Phi(), LEP_MASS);

            // Remove leptons with Pt below threshold
            if (lepIsTight[i] && lepP4[i].Pt() < TIGHT_PT_MIN)
            {
                lepIsTight[i] = kFALSE;     nTightLeps--;
            }
            if (lepIsLoose[i] && lepP4[i].Pt() < LOOSE_PT_MIN)
            {
                lepIsLoose[i] = kFALSE;     nLooseLeps--;
            }
        }



        /* SELECTION */
        // Make sure we have the right number of leps
        if (nTightLeps != 2 || nLooseLeps != nTightLeps)
            continue;


        // Find indices of pair
        unsigned x = 0, y = 0;

        if (!lepIsTight[x])
            continue;

        bool foundPair = kFALSE;
        for (unsigned j = 1; j < nLeps; j++)
        {
            if (!lepIsTight[j])
                continue;
            if (lepQ[x] == lepQ[j])
                continue;

            Double_t mll = (lepP4[x] + lepP4[j]).M();
            if (mll > MLL_MIN && mll < MLL_MAX)
            {
                y = j;  break;
            }
        }
        if (x == y)
            continue;



        /* WEIGHTING */
        // Calculate trigger efficiency scale factor
        float triggerEff = 1.;
        if (lepTriggered[x] && lepTriggered[y])
        {
            if (lepTrigEffMC[x] > 0 || lepTrigEffMC[y] > 0)
            {
                triggerEff *= 1. - (1. - lepTrigEffData[x]) * (1. - lepTrigEffData[y]);
                triggerEff /= 1. - (1. - lepTrigEffMC[x]) * (1. - lepTrigEffMC[y]);
            }
        }
        else if (lepTriggered[x])
            triggerEff = lepTrigEffData[x] / lepTrigEffMC[x];
        else if (lepTriggered[y])
            triggerEff = lepTrigEffData[y] / lepTrigEffMC[y];


        // Calculate reco efficiency scale factor
        float recoEff = lepRecoEff[x] * lepRecoEff[y];


        // Apply weights
        eventWeight *= PUWeight * recoEff * triggerEff;



        /* HISTOGRAMS */
        h_tot->Fill(6);
        h_tot->Fill(7, eventWeight);
        h_mll->Fill((lepP4[x] + lepP4[y]).M(), eventWeight);
        h_qt->Fill((lepP4[x] + lepP4[y]).Pt(), eventWeight);
        h_pt1->Fill(lepP4[x].Pt(), eventWeight);
        h_eta1->Fill(lepP4[x].Eta(), eventWeight);
        h_pt2->Fill(lepP4[y].Pt(), eventWeight);
        h_eta2->Fill(lepP4[y].Eta(), eventWeight);
        h_npv->Fill(nPV, eventWeight);
    }
    file->Close();  delete file;

    return h_tot;
}



/* 4-LEPTON */
TH1D* Select4l(const TString rootFile, const TString suffix,
               const Bool_t muPair1, const Bool_t muPair2,
               TH1F *h_m4l, TH1F *h_qt4l,
               TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
               TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4,
               TH1F *h_npv, TH1F* h_nzz, TH1F* h_nlep);
{
    ////////////////////
    // INITIALIZATION //
    ////////////////////

    /* CUTS & CONSTANTS */
    Float_t PT1_MIN = 20, PT2_MIN = 10;
    Float_t M4L_MIN = 70, M4L_MAX = 1000;
//  Float_t M4L_MIN = 60, M4L_MAX = 180;
    Float_t MZ1_MIN = 40, MZ_MIN = 12, MZ_MAX = 120;
//  Float_t MZ1_MIN = 12, MZ_MIN = 4, MZ_MAX = 120;
    Float_t MLL_MIN = 4, DR_MIN = 0.02;

    Double_t Z_MASS = 91.2;
    Double_t LEP_MASS = mus ? 0.105658369 : 0.000511;



    /* ROOT FILE */
    TFile *file = TFile::Open(rootFile);
    TH1D* h_tot;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0);
    TTreeReader reader("tree_" + suffix, file);




    /////////////////
    // SAME FLAVOR //
    /////////////////

    if (muPair1 == muPair2)
    {
        Bool_t mus = muPair1;

        /* BRANCHES */
        TString Lep = mus ? "Muon" : "Electron";
        TString lep = mus ? "muon" : "electron";
        TString BadLep = mus ? "Electron" : "Muon";

        // Event
        TTreeReaderValue<Bool_t> evtGoodTrigger_(reader, "evt"+Lep+"Triggered");
        TTreeReaderValue<Bool_t> evtBadTrigger_(reader, "evt"+BadLep+"Triggered");
        TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
        TTreeReaderValue<UShort_t> nLeps_(reader, "n"+Lep+"s");
        TTreeReaderValue<UShort_t> nTightLeps_(reader, "nTight"+Lep+"s");
        TTreeReaderValue<UShort_t> nLooseLeps_(reader, "nLoose"+Lep+"s");


        // Leptons
        TTreeReaderArray<TLorentzVector> lepP4_(reader, lep+"P4");
        TTreeReaderValue<vector<Short_t>> lepQ_(reader, lep+"Q");
        TTreeReaderValue<vector<Bool_t>> lepIsHZZ_(reader, lep+"IsHZZ");


        // Scaling & weighting
        TTreeReaderValue<vector<Float_t>> lepSF_(reader, lep+"SF");
        TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
        TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");

        TString recoBranchName = mus ? "muonHZZIDEff" : "electronHZZRecoEff";
        TTreeReaderValue<vector<Float_t>> lepRecoEff_(reader, recoBranchName);

        TTreeReaderValue<vector<Float_t>> lepTrigEffData_(reader, lep+"TriggerEffData");
        TTreeReaderValue<vector<Float_t>> lepTrigEffMC_(reader, lep+"TriggerEffMC");
        TTreeReaderValue<vector<Bool_t>> lepTriggered_(reader, lep+"Triggered");




        ////////////////
        // EVENT LOOP //
        ////////////////

        while (reader.Next())
        {
            /* BRANCHES */
            Bool_t passTrigger = (*evtGoodTrigger_ && !*evtBadTrigger_);
            UShort_t nPV = *nPV_;
            UShort_t nLeps = *nLeps_, nHZZLeps = *nHZZLeps_;

            vector<TLorentzVector> lepP4;
            vector<Short_t> lepQ;
            vector<Bool_t> lepIsHZZ;

            Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
            vector<Float_t> lepSF, lepRecoEff, lepTrigEffData, lepTrigEffMC;
            vector<Bool_t> lepTriggered;



            /* PRESELECTION */
            if (!passTrigger || nHZZLeps < 4)
                continue;



            /* LEPTONS */
            for (const TLorentzVector& lepP4__: lepP4_)
                lepP4.push_back(lepP4__);
            for (unsigned i = 0; i < nLeps; i++)
            {
                lepQ.push_back((*lepQ_)[i]);
                lepIsHZZ.push_back((*lepIsHZZ_)[i]);

                lepSF.push_back((*lepSF_)[i]);
                lepRecoEff.push_back((*lepRecoEff_)[i]);
                lepTrigEffData.push_back((*lepTrigEffData_)[i]);
                lepTrigEffMC.push_back((*lepTrigEffMC_)[i]);
                lepTriggered.push_back((*lepTriggered_)[i]);
            }


            // Apply Pt scaling
            for (unsigned i = 0; i < nLeps; i++)
                lepP4[i].SetPtEtaPhiM(lepP4[i].Pt()*lepSF[i], lepP4[i].Eta(), lepP4[i].Phi(), LEP_MASS);



            /* Z SELECTION */
            vector<pair<unsigned, unsigned>> zCand;
            for (unsigned j = 1; j < nLeps; j++)
            {
                if (!lepIsHZZ[j])
                    continue;

                for (unsigned i = 0; i < j; i++)
                {
                    if (!lepIsHZZ[i])
                        continue;
                    if (lepQ[i] == lepQ[j])
                        continue;

                    Double_t mll = (lepP4[i] + lepP4[j]).M();
                    if (mll > MZ_MIN && mll < MZ_MAX)
                        zCand.push_back(make_pair(i, j));
                }
            }
            if (zCand.size() < 2)
                continue;



            /* ZZ SELECTION */
            vector<pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>>> zzCand;
            vector<Double_t> zzCand_m4l, zzCand_mZ1;
            for (unsigned j = 1; j < zCand.size(); j++)
            {
                for (unsigned i = 0; i < j; i++)
                {
                    pair<unsigned, unsigned> Zi = zCand[i], Zj = zCand[j];


                    // Make sure pairs do not overlap
                    if (Zi.first == Zj.first || Zi.first == Zj.second
                            || Zi.second == Zj.first || Zi.second == Zj.second)
                        continue;

                    Double_t mll_i = (lepP4[Zi.first] + lepP4[Zi.second]).M();
                    Double_t mll_j = (lepP4[Zj.first] + lepP4[Zj.second]).M();


                    // Choose Z1 closest to nominal Z mass
                    pair<unsigned, unsigned> Z1, Z2;
                    if (fabs(mll_i - Z_MASS) < fabs(mll_j - Z_MASS))
                    {
                        Z1 = Zi;    Z2 = Zj;
                    }
                    else
                    {
                        Z1 = Zj;    Z2 = Zi;
                    }

                    pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
                    Z1_p4s = make_pair(lepP4[Z1.first], lepP4[Z1.second]);
                    Z2_p4s = make_pair(lepP4[Z2.first], lepP4[Z2.second]);

                    pair<Int_t, Int_t> Z1_qs, Z2_qs;
                    Z1_qs = make_pair(lepQ[Z1.first], lepQ[Z1.second]);
                    Z2_qs = make_pair(lepQ[Z2.first], lepQ[Z2.second]);

                    TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
                    Z1_p4 = Z1_p4s.first + Z1_p4s.second;
                    Z2_p4 = Z2_p4s.first + Z2_p4s.second;
                    ZZ_p4 = Z1_p4 + Z2_p4;


                    // Mass requirements
                    if (Z1_p4.M() < MZ1_MIN || ZZ_p4.M() < M4L_MIN || ZZ_p4.M() > M4L_MAX)
                        continue;


                    // Smart cut
                    pair<unsigned, unsigned> Zx, Zy, Za, Zb;
                    if (Z1_qs.first == Z2_qs.first)
                    {
                        Zx = make_pair(Z1.first, Z2.second);
                        Zy = make_pair(Z2.first, Z1.second);
                    }
                    else if (Z1_qs.first == Z2_qs.second)
                    {
                        Zx = make_pair(Z1.first, Z2.first);
                        Zy = make_pair(Z1.second, Z2.second);
                    }

                    Double_t mll_x, mll_y, mll_a, mll_b;
                    mll_x = (lepP4[Zx.first] + lepP4[Zx.second]).M();
                    mll_y = (lepP4[Zy.first] + lepP4[Zy.second]).M();

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

                    if (fabs(mll_a - Z_MASS) < fabs(Z1_p4.M() - Z_MASS) && mll_b < MZ_MIN)
                        continue;


                    // Pt2 requirement
                    vector<unsigned> z = {Z1.first, Z1.second, Z2.first, Z2.second};
                    sort(z.begin(), z.end());
                    vector<TLorentzVector> selLepP4;
                    vector<Short_t> selLepQ;

                    for (unsigned k = 0; k < z.size(); k++)
                    {
                        selLepP4.push_back(lepP4[z[k]]);
                        selLepQ.push_back(lepQ[z[k]]);
                    }

                    if (selLepP4[0].Pt() < PT1_MIN || selLepP4[1].Pt() < PT2_MIN)
                        continue;


                    // Ghost removal
                    Bool_t foundGhost = kFALSE;
                    for (unsigned y = 1; y < selLepP4.size(); y++)
                    {
                        for (unsigned x = 0; x < y; x++)
                        {
                            if (selLepP4[x].DeltaR(selLepP4[y]) < DR_MIN)
                            {
                                foundGhost = kTRUE; break;
                            }
                        }
                        if (foundGhost)
                            break;
                    }
                    if (foundGhost)
                        continue;


                    // QCD suppression
                    Bool_t foundQCD = kFALSE;
                    for (unsigned y = 1; y < selLepP4.size(); y++)
                    {
                        for (unsigned x = 0; x < y; x++)
                        {
                            if (selLepQ[x] == selLepQ[y])
                                continue;

                            if ((selLepP4[x] + selLepP4[y]).M() < MLL_MIN)
                            {
                                foundQCD = kTRUE;   break;
                            }
                        }
                        if (foundQCD)
                            break;
                    }
                    if (foundQCD)
                        continue;


                    zzCand.push_back(make_pair(Z1, Z2));
                    zzCand_mZ1.push_back(Z1_p4.M());
                    zzCand_m4l.push_back(ZZ_p4.M());
                }
            }
            if (zzCand.size() < 1)
                continue;



            /* BEST ZZ SELECTION */
            unsigned zz = 0;
            Float_t massDiff = 1000.;
            for (unsigned m = 0; m < zzCand.size(); m++)
            {
                Float_t massDiff_ = fabs(zzCand_mZ1[m] - Z_MASS);
//              Float_t massDiff_ = fabs(zzCand_m4l[m] - Z_MASS);
                if (massDiff_ < massDiff)
                {
                    massDiff = massDiff_;
                    zz = m;
                }
            }
            cout << endl;


            // Fill info
            pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>> ZZ = zzCand[zz];
            pair<unsigned, unsigned> Z1 = ZZ.first, Z2 = ZZ.second;
            pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
            Z1_p4s = make_pair(lepP4[Z1.first], lepP4[Z1.second]);
            Z2_p4s = make_pair(lepP4[Z2.first], lepP4[Z2.second]);

            TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
            Z1_p4 = Z1_p4s.first + Z1_p4s.second;
            Z2_p4 = Z2_p4s.first + Z2_p4s.second;
            ZZ_p4 = Z1_p4 + Z2_p4;

            vector<unsigned> z = {Z1.first, Z1.second, Z2.first, Z2.second};
            sort(z.begin(), z.end());



            /* WEIGHTING */
            // Calculate reco & trigger efficiency scale factors
            float recoEff = 1., triggerEff = 1.;
            float num = 1., denom = 1.;
            for (unsigned k = 0; k < z.size(); k_++)
            {
                if (lepIsHZZ[k])    // Should always be true
                    recoEff *= lepRecoEff[k];

                if (lepTriggered[k])
                {
                    num *= 1. - lepTrigEffData[k];      denom *= 1. - lepTrigEffMC[k];
                }
            }
            if (denom != 1.)
                triggerEff = (1. - num) / (1. - denom);


            // Apply weights
            eventWeight *= PUWeight * recoEff * triggerEff;



            /* HISTOGRAMS */
            h_tot->Fill(6);
            h_tot->Fill(7, eventWeight);

            h_m4l->Fill(ZZ_p4.M(), eventWeight);
            h_qt4l->Fill(ZZ_p4.Pt(), eventWeight);

            h_m12->Fill(Z1_p4.M(), eventWeight);
            h_m34->Fill(Z2_p4.M(), eventWeight);
            h_qt12->Fill(Z1_p4.Pt(), eventWeight);
            h_qt34->Fill(Z2_p4.Pt(), eventWeight);

            h_pt1->Fill(lepP4[z[0]].Pt(), eventWeight);
            h_pt2->Fill(lepP4[z[1]].Pt(), eventWeight);
            h_pt3->Fill(lepP4[z[2]].Pt(), eventWeight);
            h_pt4->Fill(lepP4[z[3]].Pt(), eventWeight);

            h_eta1->Fill(lepP4[z[0]].Eta(), eventWeight);
            h_eta2->Fill(lepP4[z[1]].Eta(), eventWeight);
            h_eta3->Fill(lepP4[z[2]].Eta(), eventWeight);
            h_eta4->Fill(lepP4[z[3]].Eta(), eventWeight);

            h_npv->Fill(nPV, eventWeight);
        }
    }
    file->Close();  delete file;

    return h_tot;
}
