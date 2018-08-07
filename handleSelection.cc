#include <vector>
#include <iostream>
#include <utility>
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




///////////////////
//  DECLARATION  //
///////////////////


//--- FUNCTIONS ---//

bool sort_dec_pt (const TLorentzVector &i, const TLorentzVector &j) { return (i.Pt() > j.Pt()); }

double GetBinContentPtEta(const TH2 *hist, const TLorentzVector &p4);
int GetXbin(const TH2 *hist, const double xval);
int GetYbin(const TH2 *hist, const double yval);

TH1D* Select2l(const TString rootFile, const TString suffix,
               const Bool_t mumu,
               TH1D *h_acc, TH1F *h_mll, TH1F *h_qt,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2,
               TH1F *h_npv, TH1F *h_sf);
  
TH1D* Select4l(const TString rootFile, const TString suffix,
               const Bool_t muPair1, const Bool_t muPair2,
               TH1D *h_acc, TH1F *h_m4l, TH1F *h_qt4l,
               TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
               TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4,
               TH1F *h_npv, TH1F* h_m4l2, TH1F* h_nlep, TH1F* h_sf);



//--- GLOBAL VARIABLES ---//

bool systOn     = kFALSE;     // TString SFType  = "Total";
// if this is kFALSE, the below don't matter

bool smearID    = kTRUE;     TString SFType  = "ID/Reco";
bool smearPtMC  = kFALSE;    // TString SFType  = "Energy";
bool smearPtData= kFALSE;    // TString SFType  = "Energy";



//--- NAMESPACE ---//

using namespace std;




///////////////////
// MAIN FUNCTION //
///////////////////


void handleSelection(const TString selection, const TString suffix, const TString id)
{
    //--- OPTIONS ---//

    Bool_t fineBinning = kTRUE, medBinning = kFALSE;



    //--- FILE INFO ---//

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
    else if (selection == "2mu2e" || selection == "2m2e")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;   muPair2 = kFALSE;
    }
    else if (selection == "2e2mu" || selection == "2e2m")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;    muPair2 = kTRUE;
    }

    if (sel2l)
        muPair2 = muPair1;


    // Path to directory where ALL trees are stored
    TString path = "root://cmsxrootd.fnal.gov//store/user/jrainbol/Trees/2016_12a/";


    // Names of ROOT files and trees for each sample
    TString dir     = suffix + "/";
    TString file    = suffix + "_" + id + ".root";


    // Name of output file
    TString output = selection + "_" + suffix + "_" + id + ".root";



    //--- HISTOGRAMS ---//

    // Lepton names 
    TString Lep12, Lep34, Lep, lep12, lep34, lep;
    Lep12   = muPair1 ? "Muon" : "Electron";
    lep12   = muPair1 ? "muon" : "electron";
    Lep34   = muPair2 ? "Muon" : "Electron";
    lep34   = muPair2 ? "muon" : "electron";
    Lep     = (muPair1 == muPair2) ? Lep12 : "Lepton";
    lep     = (muPair1 == muPair2) ? lep12 : "lepton";


    // Hist indices         // Names                    
    const unsigned M=17;    TString hname[M];
    unsigned ZZ = 0;        hname[ZZ] = "4lepMass";
    unsigned T4 = 1;        hname[T4] = "4lepPt";
    unsigned Z1 = 2;        hname[Z1] = "Z1Mass";
    unsigned Z2 = 3;        hname[Z2] = "Z2Mass";
    unsigned T1 = 4;        hname[T1] = "Z1Pt";
    unsigned T2 = 5;        hname[T2] = "Z2Pt";
    unsigned P1 = 6;        hname[P1] = "Lep1Pt";
    unsigned P2 = 7;        hname[P2] = "Lep2Pt";
    unsigned P3 = 8;        hname[P3] = "Lep3Pt";
    unsigned P4 = 9;        hname[P4] = "Lep4Pt";
    unsigned E1 = 10;       hname[E1] = "Lep1Eta";
    unsigned E2 = 11;       hname[E2] = "Lep2Eta";
    unsigned E3 = 12;       hname[E3] = "Lep3Eta";
    unsigned E4 = 13;       hname[E4] = "Lep4Eta";
    unsigned PV = 14;       hname[PV] = "nPV";
    unsigned NL = 15;       hname[NL] = "nSelLeps";
    unsigned SF = 16;       hname[SF] = "ScaleFactor";

    // Titles
    TString htitle[M];
    htitle[ZZ] = "4-" + lep + " Mass";
    htitle[T4] = "4-" + lep + "Pt";
    htitle[Z1] = "Z_{1} (Di" + lep12 + ") Mass";
    htitle[Z2] = "Z_{2} (Di" + lep34 + ") Mass";
    htitle[T1] = "Z_{1} (Di" + lep12 + ") Pt";
    htitle[T2] = "Z_{2} (Di" + lep34 + ") Pt";
    htitle[P1] = Lep + " 1 Pt";
    htitle[P2] = Lep + " 2 Pt";
    htitle[P3] = Lep + " 3 Pt";
    htitle[P4] = Lep + " 4 Pt";
    htitle[E1] = Lep + " 1 Eta";
    htitle[E2] = Lep + " 2 Eta";
    htitle[E3] = Lep + " 3 Eta";
    htitle[E4] = Lep + " 4 Eta";
    htitle[PV] = "# of Primary Vertices";
    htitle[NL] = "# of Selected " + Lep + "s";
    htitle[SF] = "SF (" + SFType + ")";


    // Binning
    Int_t bins[M];      Double_t low[M],    up[M];
    bins[ZZ] = 20;      low[ZZ] = 70;       up[ZZ] = 170;
    bins[T4] = 15;      low[T4] = 0;        up[T4] = 120;
    bins[Z1] = 15;      low[Z1] = 0;        up[Z1] = 120;
    bins[Z2] = 15;      low[Z2] = 0;        up[Z2] = 120;
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
    bins[NL] = 10;      low[NL] = 0.5;      up[NL] = 10.5;
    bins[SF] = 100;     low[SF] = 0;        up[SF] = 2;

    if (fineBinning)
    {
        bins[ZZ] = 50;
        bins[T4] = 250;     low[T4] = 0;        up[T4] = 1000;
        bins[Z1] = 150;
        bins[Z2] = 150;
        bins[T1] = 150;
        bins[T2] = 150;
        bins[P1] = 75;                          
        bins[P2] = 75;                          up[P2] = 150;     
        bins[P3] = 75;                          up[P3] = 150;
        bins[P4] = 75;                          up[P4] = 150;
        bins[E1] = 50;      low[E1] = -2.5;     up[E1] = 2.5;
        bins[E2] = 50;      low[E2] = -2.5;     up[E2] = 2.5;
        bins[E3] = 50;      low[E3] = -2.5;     up[E3] = 2.5;
        bins[E4] = 50;      low[E4] = -2.5;     up[E4] = 2.5;
        bins[PV] = 51;      low[PV] = -0.5;     up[PV] = 50.5;
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
        bins[P1] = 100;     low[P1] = 25;       up[P1] = 125;
        bins[P2] = 100;     low[P2] = 10;       up[P2] = 110;
        bins[E1] = 50;      low[E1] = -2.5;     up[E1] = 2.5;
        bins[E2] = 50;      low[E2] = -2.5;     up[E2] = 2.5;
        bins[PV] = 51;      low[PV] = -0.5;     up[PV] = 50.5;
    }
    else if (sel4l)
        QT = 1;


    // Make histograms
    TH1F *h[M];
    for (unsigned j = 0; j < M; j++)
    {
        h[j] = new TH1F(hname[j] + "_" + suffix, htitle[j], bins[j], low[j], up[j]);
        h[j]->Sumw2(kTRUE);
    }
    TH1D *hTotalEvents;
    TH1D *hAcceptedEvents = new TH1D("AcceptedEvents_" + suffix, "AcceptedEvents", 10, 0.5, 10.5);
    hAcceptedEvents->Sumw2(kTRUE);



    //--- SELECTION ---//

    // Read trees
    if (sel2l)
        hTotalEvents = Select2l(path + dir + file, suffix, muPair1,
                                hAcceptedEvents, h[LL], h[QT], 
                                h[P1], h[P2], h[E1], h[E2],
                                h[PV], h[SF]);
    else if (sel4l)
        hTotalEvents = Select4l(path + dir + file, suffix, muPair1, muPair2,
                                hAcceptedEvents, h[ZZ], h[QT],
                                h[Z1], h[Z2], h[T1], h[T2],
                                h[P1], h[P2], h[P3], h[P4], h[E1], h[E2], h[E3], h[E4],
                                h[PV], h[ZM], h[NL], h[SF]);


    // Write to file
    TFile *outFile = new TFile(output, "RECREATE");
    for (unsigned j = 0; j < M; j++)
        h[j]->Write();
    hTotalEvents->Write();
    hAcceptedEvents->Write();
    outFile->Close();



    //--- CLEANUP ---//

    delete outFile;
    for (unsigned j = 0; j < M; j++)
        delete h[j];
    delete hTotalEvents;
    delete hAcceptedEvents;

    return;
}




/////////////////////////
// SELECTION FUNCTIONS //
/////////////////////////


//--- DILEPTON ---//

TH1D* Select2l(const TString rootFile, const TString suffix, const Bool_t mumu,
               TH1D *h_acc, TH1F *h_mll, TH1F *h_qt,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2,
               TH1F *h_npv, TH1F *h_sf)
{



    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //--- CUTS & CONSTS ---//

    Float_t MLL_MIN = 80, MLL_MAX = 100;
    Float_t PT1_MIN = 25, PT2_MIN = 10;
//  Float_t ID_PT_MIN = mumu ? 7 : 5;
    Double_t LEP_MASS = mumu ? 0.105658369 : 0.000511;

    TString Lep = mumu ? "Muon" : "Electron";
    TString lep = mumu ? "muon" : "electron";
//  TString ID = mumu ? "Loose" : "Medium";
    TString Reco = mumu ? "ID" : "Reco";

    bool isData = suffix.Contains("2016");


    
    //--- SYSTEMATICS ---//

    // Lepton ID SF (temp fix)
    TFile *file_id = TFile::Open("hzz_muon_id_sf.root");
    TH2 *h_id;
    file_id->GetObject("FINAL", h_id);
    h_id->SetDirectory(0);
    file_id->Close();


    // Lepton ID smear
    TFile *file_unc = TFile::Open("hzz_muon_id_smear.root");
    TH2 *h_unc;
    file_unc->GetObject("SMEAR", h_unc);
//  file_id->GetObject("ERROR", h_unc);
    h_unc->SetDirectory(0);
    file_unc->Close();
//  file_id->Close();


    // Lepton Pt smear
    Double_t PT_UNC = 0.002;



    //--- TREE FILE ---//

    TFile *file = TFile::Open(rootFile);

    // Branches
    TTreeReader reader("tree_" + suffix, file);


    // Histograms
    TH1D *h_tot, *h_acc_;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0); h_tot->Sumw2(kTRUE);
    file->GetObject("AcceptedEvents_" + suffix, h_acc_);

    // Copy applicable result from accepted events histogram
    if (h_acc_)
    {
        unsigned binIdx = 1;
        if (mumu)       // mumu = 3
            binIdx = 3;
        else            // ee   = 4
            binIdx = 4;

        h_acc->SetBinContent(1, h_acc_->GetBinContent(1));
        h_acc->SetBinContent(2, h_acc_->GetBinContent(binIdx));
        delete h_acc_;
    }



    //--- BRANCHES ---//

    // Event
    TTreeReaderValue<Bool_t> passTrigger_(reader, "evt"+Lep+"Triggered");
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nPartons_(reader, "nPartons");
    TTreeReaderValue<UShort_t> nLeps_(reader, "n"+Lep+"s");
    TTreeReaderValue<UShort_t> nHZZLeps_(reader, "nHZZ"+Lep+"s");
//  TTreeReaderValue<UShort_t> nIDLeps_(reader, "n"+ID+Lep+"s");


    // Leptons
    TTreeReaderArray<TLorentzVector> lepP4_(reader, lep+"P4");
    TTreeReaderValue<vector<Short_t>> lepQ_(reader, lep+"Q");
    TTreeReaderValue<vector<Bool_t>> lepIsHZZ_(reader, lep+"IsHZZ");
//  TTreeReaderValue<vector<Bool_t>> lepIsID_(reader, lep+"Is"+ID);


    // Scaling & weighting
    TTreeReaderValue<vector<Float_t>> lepEnergySF_(reader, lep+"SF");
    TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
    TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");

    TTreeReaderValue<vector<Float_t>> lepRecoSF_(reader, lep+"HZZ"+Reco+"Weight");
    TTreeReaderValue<vector<Float_t>> lepTrigEffData_(reader, lep+"TriggerEffData");
    TTreeReaderValue<vector<Float_t>> lepTrigEffMC_(reader, lep+"TriggerEffMC");
    TTreeReaderValue<vector<Bool_t>> lepTriggered_(reader, lep+"Triggered");




    ////////////////
    // EVENT LOOP //
    ////////////////


    while (reader.Next())
    {
        //--- BRANCHES ---//

        Bool_t passTrigger = (*passTrigger_);
        UShort_t nPV = *nPV_, nPartons = *nPartons_;
        UShort_t nLeps = *nLeps_, nHZZLeps = *nHZZLeps_;    // , nIDLeps = *nIDLeps_;

        vector<TLorentzVector> lepP4;
        vector<Short_t> lepQ;
        vector<Bool_t> lepIsHZZ;    // , lepIsID;

        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        vector<Float_t> lepEnergySF, lepRecoSF, lepTrigEffData, lepTrigEffMC;
        vector<Bool_t> lepTriggered;



        //--- PRESELECTION ---//

        if (!passTrigger || nHZZLeps != 2)
            continue;



        //--- LEPTONS ---//

        for (const TLorentzVector& lepP4__: lepP4_)
            lepP4.push_back(lepP4__);
        for (unsigned i = 0; i < nLeps; i++)
        {
            lepQ.push_back((*lepQ_)[i]);
            lepIsHZZ.push_back((*lepIsHZZ_)[i]);
//          lepIsID.push_back((*lepIsID_)[i]);

            lepEnergySF.push_back((*lepEnergySF_)[i]);
            lepRecoSF.push_back((*lepRecoSF_)[i]);
            lepTrigEffData.push_back((*lepTrigEffData_)[i]);
            lepTrigEffMC.push_back((*lepTrigEffMC_)[i]);
            lepTriggered.push_back((*lepTriggered_)[i]);
        }


        // Adjust lepton requirements
        for (unsigned i = 0; i < nLeps; i++)
        {
            // Apply Pt smear
            if (systOn)
            {
                if (smearPtMC && !isData)
                    lepEnergySF[i] += PT_UNC * lepEnergySF[i];
                else if (smearPtData && isData)
                    lepEnergySF[i] += PT_UNC * lepEnergySF[i];
            }


            // Apply Pt scaling
            lepP4[i].SetPtEtaPhiM(lepP4[i].Pt()*lepEnergySF[i],
                                  lepP4[i].Eta(), lepP4[i].Phi(), LEP_MASS);


            // Remove leptons with Pt below threshold
            if (lepIsHZZ[i] && lepP4[i].Pt() < PT2_MIN)
            {
                lepIsHZZ[i] = kFALSE;       nHZZLeps--;
            }


            // Fix reco weight (temp)
            if (!isData)
                lepRecoSF[i] = GetBinContentPtEta(h_id, lepP4[i]);


            // Apply ID smear
            if (systOn)
            {
                if (smearID && !isData)
                    lepRecoSF[i] += GetBinContentPtEta(h_unc, lepP4[i]);
            }
        }



        //--- SELECTION ---//

        // Make sure we have the right number of leps
//      if (nHZZLeps != 2 || nIDLeps > nHZZLeps)
        if (nHZZLeps != 2)
            continue;


        // Make sure leading lepton is good
        if (!lepIsHZZ[0] || lepP4[0].Pt() < PT1_MIN)
            continue;


        // Pair leading lepton
        unsigned x = 0, y = 0;
        bool foundPair = kFALSE;
        for (unsigned j = 1; j < nLeps; j++)
        {
            if (!lepIsHZZ[j])
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



        //--- WEIGHTING ---//

        // Calculate trigger efficiency scale factor
        float triggerWeight = 1.;
        if (lepTriggered[x] && lepTriggered[y])
        {
            if (lepTrigEffMC[x] > 0 || lepTrigEffMC[y] > 0)
            {
                triggerWeight *= 1. - (1. - lepTrigEffData[x]) * (1. - lepTrigEffData[y]);
                triggerWeight /= 1. - (1. - lepTrigEffMC[x]) * (1. - lepTrigEffMC[y]);
            }
        }
        else if (lepTriggered[x])
            triggerWeight = lepTrigEffData[x] / lepTrigEffMC[x];
        else if (lepTriggered[y])
            triggerWeight = lepTrigEffData[y] / lepTrigEffMC[y];


        // Calculate reco efficiency scale factor
        float recoWeight = lepRecoSF[x] * lepRecoSF[y];


        // Apply weights
        eventWeight *= PUWeight * recoWeight * triggerWeight;



        //--- HISTOGRAMS ---//

        // Record MC selection
        h_acc->Fill(3, eventWeight);

        if (suffix.Contains("dy") && nPartons != 0)
            continue;


        // Fill other histograms
        h_tot->Fill(6);
        h_tot->Fill(7, eventWeight);
        h_mll->Fill((lepP4[x] + lepP4[y]).M(), eventWeight);
        h_qt->Fill((lepP4[x] + lepP4[y]).Pt(), eventWeight);
        h_pt1->Fill(lepP4[x].Pt(), eventWeight);
        h_eta1->Fill(lepP4[x].Eta(), eventWeight);
        h_pt2->Fill(lepP4[y].Pt(), eventWeight);
        h_eta2->Fill(lepP4[y].Eta(), eventWeight);
        h_npv->Fill(nPV, eventWeight);

        if (smearID)
            h_sf->Fill(recoWeight);
        else if (smearPtData || smearPtMC)
        {
            h_sf->Fill(lepEnergySF[x]);
            h_sf->Fill(lepEnergySF[y]);
        }

    }
    file->Close();  delete file;

    return h_tot;
}



//--- 4-LEPTON ---//

TH1D* Select4l(const TString rootFile, const TString suffix,
               const Bool_t muPair1, const Bool_t muPair2,
               TH1D *h_acc, TH1F *h_m4l, TH1F *h_qt4l,
               TH1F *h_m12, TH1F *h_m34, TH1F *h_qt12, TH1F *h_qt34,
               TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
               TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4,
               TH1F *h_npv, TH1F* h_m4l2, TH1F *h_nlep, TH1F *h_sf)
{



    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //--- CUTS & CONSTS ---//

    Float_t PT1_MIN = 25, PT2_MIN = 10;
//  Float_t M4L_MIN = 70, M4L_MAX = 1000;
    Float_t M4L_MIN = 80, M4L_MAX = 100;
//  Float_t MZ1_MIN = 40, MZ_MIN = 12, MZ_MAX = 120;
    Float_t MZ1_MIN = 12, MZ_MIN = 4, MZ_MAX = 120;
    Float_t MLL_MIN = 4, DR_MIN = 0.02;

    Double_t Z_MASS = 91.2;
    Double_t MU_MASS = 0.105658369, ELE_MASS = 0.000511;

    bool isData = suffix.Contains("2016");


    
    //--- SYSTEMATICS ---//

    // Muon ID SF (temp fix)
    TFile *file_id = TFile::Open("hzz_muon_id_sf.root");
    TH2 *h_id;
    file_id->GetObject("FINAL", h_id);
    h_id->SetDirectory(0);
    file_id->Close();


    // Lepton ID smear
    TFile *file_unc = TFile::Open("hzz_muon_id_smear.root");
    TH2 *h_unc;
    file_unc->GetObject("SMEAR", h_unc);
//  file_id->GetObject("ERROR", h_unc);
    h_unc->SetDirectory(0);
    file_unc->Close();
//  file_id->Close();


    // Lepton Pt smear
    Double_t PT_UNC = 0.002;



    //--- TREE FILE ---//

    TFile *file = TFile::Open(rootFile);

    // Branches
    TTreeReader reader("tree_" + suffix, file);


    // Histograms
    TH1D *h_tot, *h_acc_;
    file->GetObject("TotalEvents_" + suffix, h_tot);
    h_tot->SetDirectory(0); h_tot->Sumw2(kTRUE);
    file->GetObject("AcceptedEvents_" + suffix, h_acc_);


    // Copy applicable result from accepted events histogram
    if (h_acc_)
    {
        unsigned binIdx = 1;
        if (muPair1 && muPair2)         // 4m   = 6
            binIdx = 6;
        else if (muPair1 && !muPair2)   // 2m2e = 7
            binIdx = 7;
        else if (!muPair1 && muPair2)   // 2e2m = 8
            binIdx = 8;
        else if (!muPair1 && !muPair2)  // 4e   = 9
            binIdx = 9;

        h_acc->SetBinContent(1, h_acc_->GetBinContent(1));
        h_acc->SetBinContent(2, h_acc_->GetBinContent(binIdx));
        delete h_acc_;
    }




    /////////////////
    // SAME FLAVOR //
    /////////////////


    if (muPair1 == muPair2)
    {
        Bool_t mus = muPair1;
        Double_t LEP_MASS = mus ? MU_MASS : ELE_MASS;



        //--- BRANCHES ---//

        TString Lep = mus ? "Muon" : "Electron";
        TString lep = mus ? "muon" : "electron";
//      TString ID = mus ? "Loose" : "Loose";
        TString Reco = mus ? "ID" : "Reco";


        // Event
        TTreeReaderValue<Bool_t> passTrigger_(reader, "evt"+Lep+"Triggered");
        TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
        TTreeReaderValue<UShort_t> nPartons_(reader, "nPartons");
        TTreeReaderValue<UShort_t> nLeps_(reader, "n"+Lep+"s");
        TTreeReaderValue<UShort_t> nHZZLeps_(reader, "nHZZ"+Lep+"s");


        // Leptons
        TTreeReaderArray<TLorentzVector> lepP4_(reader, lep+"P4");
        TTreeReaderValue<vector<Short_t>> lepQ_(reader, lep+"Q");
//      TTreeReaderValue<vector<Bool_t>> lepIsID_(reader, lep+"Is"+ID);
        TTreeReaderValue<vector<Bool_t>> lepIsHZZ_(reader, lep+"IsHZZ");


        // Scaling & weighting
        TTreeReaderValue<vector<Float_t>> lepEnergySF_(reader, lep+"SF");
        TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
        TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");

        TTreeReaderValue<vector<Float_t>> lepRecoSF_(reader, lep+"HZZ"+Reco+"Weight");
        TTreeReaderValue<vector<Float_t>> lepTrigEffData_(reader, lep+"TriggerEffData");
        TTreeReaderValue<vector<Float_t>> lepTrigEffMC_(reader, lep+"TriggerEffMC");
        TTreeReaderValue<vector<Bool_t>> lepTriggered_(reader, lep+"Triggered");




        ////////////////
        // EVENT LOOP //
        ////////////////


        unsigned count = 0, max = 100;
//      while (reader.Next() && count < max)
        while (reader.Next())
        {
            //--- BRANCHES ---//

            Bool_t passTrigger = (*passTrigger_);
            UShort_t nPV = *nPV_, nPartons = *nPartons_;
            UShort_t nLeps = *nLeps_, nHZZLeps = *nHZZLeps_;

            vector<TLorentzVector> lepP4;
            vector<Short_t> lepQ;
            vector<Bool_t> lepIsHZZ;    // , lepIsID;

            Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
            vector<Float_t> lepEnergySF, lepRecoSF, lepTrigEffData, lepTrigEffMC;
            vector<Bool_t> lepTriggered;

            vector<Bool_t> elecPassID;



            //--- PRESELECTION ---//

            if (!passTrigger || nHZZLeps < 4)
                continue;



            //--- LEPTONS ---//

            for (const TLorentzVector& lepP4__: lepP4_)
                lepP4.push_back(lepP4__);
            for (unsigned i = 0; i < nLeps; i++)
            {
                lepQ.push_back((*lepQ_)[i]);
//              lepIsID.push_back((*lepIsID_)[i]);
                lepIsHZZ.push_back((*lepIsHZZ_)[i]);

                lepEnergySF.push_back((*lepEnergySF_)[i]);
                lepRecoSF.push_back((*lepRecoSF_)[i]);
                lepTrigEffData.push_back((*lepTrigEffData_)[i]);
                lepTrigEffMC.push_back((*lepTrigEffMC_)[i]);
                lepTriggered.push_back((*lepTriggered_)[i]);
            }


            // Adjust lepton requirements
            for (unsigned i = 0; i < nLeps; i++)
            {
                // Apply Pt smear
                if (systOn)
                {
                    if (smearPtMC && !isData)
                        lepEnergySF[i] += PT_UNC * lepEnergySF[i];
                    else if (smearPtData && isData)
                        lepEnergySF[i] += PT_UNC * lepEnergySF[i];
                }


                // Apply Pt scaling
                lepP4[i].SetPtEtaPhiM(lepP4[i].Pt() * lepEnergySF[i],
                        lepP4[i].Eta(), lepP4[i].Phi(), LEP_MASS);


                if (mus)
                {
                    // Fix reco weight (temp)
                    if (!isData)
                        lepRecoSF[i] = GetBinContentPtEta(h_id, lepP4[i]);


                    // Apply systematics
                    if (systOn)
                    {
                        if (smearID && !isData)
                            lepRecoSF[i] += GetBinContentPtEta(h_unc, lepP4[i]);
                    }
                }
            }



            //--- Z SELECTION ---//

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



            //--- ZZ SELECTION ---//

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



            //--- BEST ZZ SELECTION ---//

            unsigned zz = 0;
            Float_t massDiff = 1000.;
            for (unsigned m = 0; m < zzCand.size(); m++)
            {
                Float_t massDiff_ = fabs(zzCand_mZ1[m] - Z_MASS); //= fabs(zzCand_m4l[m] - Z_MASS);
                if (massDiff_ < massDiff)
                {
                    massDiff = massDiff_;
                    zz = m;
                }
            }


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



            //--- WEIGHTING ---//

            // Calculate trigger efficiency scale factor using leading leptons
            unsigned x = z[0], y = z[1];
            float triggerWeight = 1.;
            if (lepTriggered[x] && lepTriggered[y])
            {
                if (lepTrigEffMC[x] > 0 || lepTrigEffMC[y] > 0)
                {
                    triggerWeight *= 1. - (1. - lepTrigEffData[x]) * (1. - lepTrigEffData[y]);
                    triggerWeight /= 1. - (1. - lepTrigEffMC[x]) * (1. - lepTrigEffMC[y]);
                }
            }
            else if (lepTriggered[x])
                triggerWeight = lepTrigEffData[x] / lepTrigEffMC[x];
            else if (lepTriggered[y])
                triggerWeight = lepTrigEffData[y] / lepTrigEffMC[y];


            // Calculate reco efficiency scale factor using all leptons
            float recoWeight = 1.;
            for (unsigned k_ = 0; k_ < z.size(); k_++)
            {
                unsigned k = z[k_];
                recoWeight *= lepRecoSF[k];
            }
/*
            for (unsigned k_ = 0; k_ < z.size(); k_++)
            {
                unsigned k = z[k_];
                cout << lepRecoSF[k] << " = ";
                cout << GetBinContentPtEta(h_id, lepP4[k]) << " + ";
                cout << GetBinContentPtEta(h_unc, lepP4[k]) << endl;
            }
            cout << endl;
*/

            // Apply weights
            eventWeight *= PUWeight * recoWeight * triggerWeight;



            //--- HISTOGRAMS ---//

            // Record MC selection
            h_acc->Fill(3, eventWeight);
            if (suffix.Contains("dy") && nPartons != 0)
                continue;
            count++;


            // Fill other histograms
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
            h_m4l2->Fill(ZZ_p4.M(), eventWeight);
            h_nlep->Fill(nHZZLeps, eventWeight);

            if (smearID)
                h_sf->Fill(recoWeight);
            else if (smearPtData || smearPtMC)
            {
                for (unsigned k_ = 0; k_ < z.size(); k_++)
                {
                    unsigned k = z[k_];
                    h_sf->Fill(lepEnergySF[k]);
                }
            }
        }
    }




    /////////////////////
    // OPPOSITE FLAVOR //
    /////////////////////


    else
    {
        Double_t LEP1_MASS = muPair1 ? MU_MASS : ELE_MASS;
        Double_t LEP2_MASS = muPair2 ? MU_MASS : ELE_MASS;



        //--- BRANCHES ---//

        TString Lep1 = muPair1 ? "Muon" : "Electron",   Lep2 = muPair2 ? "Muon" : "Electron";
        TString lep1 = muPair1 ? "muon" : "electron",   lep2 = muPair2 ? "muon" : "electron";
//      TString ID1 = muPair1 ? "Loose" : "Medium",     ID2 = muPair2 ? "Loose" : "Medium";
        TString Reco1 = muPair1 ?  "ID" : "Reco",       Reco2 = muPair2 ?  "ID" : "Reco";


        // Event
        TTreeReaderValue<Bool_t> passTrigger_(reader, "evt"+Lep1+"Triggered");
        TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
        TTreeReaderValue<UShort_t> nPartons_(reader, "nPartons");

        TTreeReaderValue<UShort_t> nLeps_(reader, "nLeptons");
        TTreeReaderValue<UShort_t> nHZZLeps_(reader, "nHZZLeptons");

        TTreeReaderValue<UShort_t> nLep1s_(reader, "n"+Lep1+"s");
        TTreeReaderValue<UShort_t> nHZZLep1s_(reader, "nHZZ"+Lep1+"s");

        TTreeReaderValue<UShort_t> nLep2s_(reader, "n"+Lep2+"s");
        TTreeReaderValue<UShort_t> nHZZLep2s_(reader, "nHZZ"+Lep2+"s");


        // Leptons
        TTreeReaderArray<TLorentzVector> lep1P4_(reader, lep1+"P4");
        TTreeReaderValue<vector<Short_t>> lep1Q_(reader, lep1+"Q");
//      TTreeReaderValue<vector<Bool_t>> lep1IsID_(reader, lep1+"Is"+ID1);
        TTreeReaderValue<vector<Bool_t>> lep1IsHZZ_(reader, lep1+"IsHZZ");

        TTreeReaderArray<TLorentzVector> lep2P4_(reader, lep2+"P4");
        TTreeReaderValue<vector<Short_t>> lep2Q_(reader, lep2+"Q");
//      TTreeReaderValue<vector<Bool_t>> lep2IsID_(reader, lep2+"Is"+ID2);
        TTreeReaderValue<vector<Bool_t>> lep2IsHZZ_(reader, lep2+"IsHZZ");


        // Scaling & weighting
        TTreeReaderValue<vector<Float_t>> lep1EnergySF_(reader, lep1+"SF");
        TTreeReaderValue<vector<Float_t>> lep2EnergySF_(reader, lep2+"SF");

        TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
        TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");

        TString lep1RecoBranchName = muPair1 ? "muonHZZIDWeight" : "electronHZZRecoWeight";
        TString lep2RecoBranchName = muPair2 ? "muonHZZIDWeight" : "electronHZZRecoWeight";
        TTreeReaderValue<vector<Float_t>> lep1RecoSF_(reader, lep1+"HZZ"+Reco1+"Weight");
        TTreeReaderValue<vector<Float_t>> lep2RecoSF_(reader, lep2+"HZZ"+Reco2+"Weight");

        TTreeReaderValue<vector<Float_t>> lep1TrigEffData_(reader, lep1+"TriggerEffData");
        TTreeReaderValue<vector<Float_t>> lep1TrigEffMC_(reader, lep1+"TriggerEffMC");
        TTreeReaderValue<vector<Bool_t>> lep1Triggered_(reader, lep1+"Triggered");


        ////////////////
        // EVENT LOOP //
        ////////////////


        while (reader.Next())
        {
            //--- BRANCHES ---//

            Bool_t passTrigger = (*passTrigger_);
            UShort_t nPV = *nPV_, nPartons = *nPartons_;
            UShort_t nLeps = *nLeps_, nHZZLeps = *nHZZLeps_;
            UShort_t nLep1s = *nLep1s_, nHZZLep1s = *nHZZLep1s_;
            UShort_t nLep2s = *nLep2s_, nHZZLep2s = *nHZZLep2s_;

            vector<TLorentzVector> lep1P4, lep2P4;
            vector<Short_t> lep1Q, lep2Q;
            vector<Bool_t> lep1IsHZZ, lep2IsHZZ;    // lep1IsID, lep2IsID; 

            Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
            vector<Float_t> lep1EnergySF, lep2EnergySF, lep1RecoSF, lep2RecoSF;
            vector<Float_t> lep1TrigEffData, lep1TrigEffMC;
            vector<Bool_t> lep1Triggered;



            //--- PRESELECTION ---//

            if (!passTrigger || nHZZLep1s < 2 || nHZZLep2s < 2)
                continue;



            //--- LEPTONS ---//

            // Leading flavor
            for (const TLorentzVector& lep1P4__: lep1P4_)
                lep1P4.push_back(lep1P4__);
            for (unsigned i = 0; i < nLep1s; i++)
            {
                lep1Q.push_back((*lep1Q_)[i]);
//              lep1IsID.push_back((*lep1IsID_)[i]);
                lep1IsHZZ.push_back((*lep1IsHZZ_)[i]);

                lep1EnergySF.push_back((*lep1EnergySF_)[i]);
                lep1RecoSF.push_back((*lep1RecoSF_)[i]);
                lep1TrigEffData.push_back((*lep1TrigEffData_)[i]);
                lep1TrigEffMC.push_back((*lep1TrigEffMC_)[i]);
                lep1Triggered.push_back((*lep1Triggered_)[i]);
            }

            // Adjust lepton requirements
            for (unsigned i = 0; i < nLep1s; i++)
            {
                // Apply Pt smear
                if (systOn)
                {
                    if (smearPtMC && !isData)
                        lep1EnergySF[i] += PT_UNC * lep1EnergySF[i];
                    else if (smearPtData && isData)
                        lep1EnergySF[i] += PT_UNC * lep1EnergySF[i];
                }


                // Apply Pt scaling
                lep1P4[i].SetPtEtaPhiM(lep1P4[i].Pt() * lep1EnergySF[i],
                                       lep1P4[i].Eta(), lep1P4[i].Phi(), LEP1_MASS);


                if (muPair1)
                {
                    // Fix reco weight (temp)
                    if (!isData)
                        lep1RecoSF[i] = GetBinContentPtEta(h_id, lep1P4[i]);


                    // Apply systematics
                    if (systOn)
                    {
                        if (smearID)
                            lep1RecoSF[i] += GetBinContentPtEta(h_unc, lep1P4[i]);
                    }
                }
            }


            // Subleading flavor
            for (const TLorentzVector& lep2P4__: lep2P4_)
                lep2P4.push_back(lep2P4__);
            for (unsigned i = 0; i < nLep2s; i++)
            {
                lep2Q.push_back((*lep2Q_)[i]);
//              lep2IsID.push_back((*lep2IsID_)[i]);
                lep2IsHZZ.push_back((*lep2IsHZZ_)[i]);

                lep2EnergySF.push_back((*lep2EnergySF_)[i]);
                lep2RecoSF.push_back((*lep2RecoSF_)[i]);
            }

            // Adjust lepton requirements
            for (unsigned i = 0; i < nLep2s; i++)
            {
                // Apply Pt smear
                if (systOn)
                {
                    if (smearPtMC && !isData)
                        lep2EnergySF[i] += PT_UNC * lep2EnergySF[i];
                    else if (smearPtData && isData)
                        lep2EnergySF[i] += PT_UNC * lep2EnergySF[i];
                }


                // Apply Pt scaling
                lep2P4[i].SetPtEtaPhiM(lep2P4[i].Pt() * lep2EnergySF[i],
                                       lep2P4[i].Eta(), lep2P4[i].Phi(), LEP2_MASS);


                if (muPair2)
                {
                    // Fix reco weight (temp)
                    if (!isData)
                        lep2RecoSF[i] = GetBinContentPtEta(h_id, lep2P4[i]);


                    // Apply systematics
                    if (systOn && !isData)
                    {
                        if (smearID)
                            lep2RecoSF[i] += GetBinContentPtEta(h_unc, lep2P4[i]);
                    }
                }
            }



            //--- Z SELECTION ---//

            vector<pair<unsigned, unsigned>> z1Cand, z2Cand;


            // Z1 candidates
            for (unsigned j = 1; j < nLep1s; j++)
            {
                if (!lep1IsHZZ[j])
                    continue;

                for (unsigned i = 0; i < j; i++)
                {
                    if (!lep1IsHZZ[i] || lep1P4[i].Pt() < PT1_MIN)
                        continue;
                    if (lep1Q[i] == lep1Q[j])
                        continue;

                    Double_t mll = (lep1P4[i] + lep1P4[j]).M();
                    if (mll > MZ1_MIN && mll < MZ_MAX)
                        z1Cand.push_back(make_pair(i, j));
                }
            }
            if (z1Cand.size() < 1)
                continue;


            // Z2 candidates
            for (unsigned j = 1; j < nLep2s; j++)
            {
                if (!lep2IsHZZ[j])
                    continue;

                for (unsigned i = 0; i < j; i++)
                {
                    if (!lep2IsHZZ[i])
                        continue;
                    if (lep2Q[i] == lep2Q[j])
                        continue;

                    Double_t mll = (lep2P4[i] + lep2P4[j]).M();
                    if (mll > MZ_MIN && mll < MZ_MAX)
                        z2Cand.push_back(make_pair(i, j));
                }
            }
            if (z2Cand.size() < 1)
                continue;



            //--- ZZ SELECTION ---//

            vector<pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>>> zzCand;
            vector<Double_t> zzCand_m4l, zzCand_mZ1;
            for (unsigned i = 0; i < z1Cand.size(); i++)
            {
                for (unsigned j = 0; j < z2Cand.size(); j++)
                {
                    pair<unsigned, unsigned> Z1 = z1Cand[i], Z2 = z2Cand[j];

                    pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
                    Z1_p4s = make_pair(lep1P4[Z1.first], lep1P4[Z1.second]);
                    Z2_p4s = make_pair(lep2P4[Z2.first], lep2P4[Z2.second]);

                    pair<Short_t, Short_t> Z1_qs, Z2_qs;
                    Z1_qs = make_pair(lep1Q[Z1.first], lep1Q[Z1.second]);
                    Z2_qs = make_pair(lep2Q[Z2.first], lep2Q[Z2.second]);

                    TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
                    Z1_p4 = Z1_p4s.first + Z1_p4s.second;
                    Z2_p4 = Z2_p4s.first + Z2_p4s.second;
                    ZZ_p4 = Z1_p4 + Z2_p4;


                    // Ensure Z1 is closest to nominal Z mass
                    if (fabs(Z1_p4.M() - Z_MASS) > fabs(Z2_p4.M() - Z_MASS))
                        continue;


                    // Mass requirement
                    if (ZZ_p4.M() < M4L_MIN || ZZ_p4.M() > M4L_MAX)
                        continue;


                    // Pt2 requirement
                    if (Z1_p4s.second.Pt() < PT2_MIN && Z2_p4s.first.Pt() < PT2_MIN)
                        continue;


                    // Ghost removal
                    vector<TLorentzVector> selLepP4 = {Z1_p4s.first, Z1_p4s.second,
                                                       Z2_p4s.first, Z2_p4s.second};
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
                    vector<Short_t> selLepQ = {Z1_qs.first, Z1_qs.second,
                                               Z2_qs.first, Z2_qs.second};
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



            //--- BEST ZZ SELECTION ---//

            unsigned zz = 0;
            Float_t massDiff = 1000.;
            for (unsigned m = 0; m < zzCand.size(); m++)
            {
                Float_t massDiff_ = fabs(zzCand_mZ1[m] - Z_MASS); //= fabs(zzCand_m4l[m] - Z_MASS);
                if (massDiff_ < massDiff)
                {
                    massDiff = massDiff_;
                    zz = m;
                }
            }


            // Fill info
            pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>> ZZ = zzCand[zz];
            pair<unsigned, unsigned> Z1 = ZZ.first, Z2 = ZZ.second;
            pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
            Z1_p4s = make_pair(lep1P4[Z1.first], lep1P4[Z1.second]);
            Z2_p4s = make_pair(lep2P4[Z2.first], lep2P4[Z2.second]);

            TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
            Z1_p4 = Z1_p4s.first + Z1_p4s.second;
            Z2_p4 = Z2_p4s.first + Z2_p4s.second;
            ZZ_p4 = Z1_p4 + Z2_p4;

            vector<TLorentzVector> z_p4s = {Z1_p4s.first, Z1_p4s.second,
                                            Z2_p4s.first, Z2_p4s.second};
            sort(z_p4s.begin(), z_p4s.end(), sort_dec_pt);



            //--- WEIGHTING ---//

            // Calculate trigger efficiency scale factor using leading leptons
            unsigned x = Z1.first, y = Z1.second;
            float triggerWeight = 1.;
            if (lep1Triggered[x] && lep1Triggered[y])
            {
                if (lep1TrigEffMC[x] > 0 || lep1TrigEffMC[y] > 0)
                {
                    triggerWeight *= 1. - (1. - lep1TrigEffData[x]) * (1. - lep1TrigEffData[y]);
                    triggerWeight /= 1. - (1. - lep1TrigEffMC[x]) * (1. - lep1TrigEffMC[y]);
                }
            }
            else if (lep1Triggered[x])
                triggerWeight = lep1TrigEffData[x] / lep1TrigEffMC[x];
            else if (lep1Triggered[y])
                triggerWeight = lep1TrigEffData[y] / lep1TrigEffMC[y];


            // Calculate reco efficiency scale factor using all leptons
            float recoWeight = 1.;
            recoWeight *= lep1RecoSF[x] * lep1RecoSF[y];
            recoWeight *= lep2RecoSF[Z2.first] * lep2RecoSF[Z2.second];


            // Apply weights
            eventWeight *= PUWeight * recoWeight * triggerWeight;



            //--- HISTOGRAMS ---//

            // Record MC selection
            h_acc->Fill(3, eventWeight);

            if (suffix.Contains("dy") && nPartons != 0)
                continue;


            // Fill other histograms
            h_tot->Fill(6);
            h_tot->Fill(7, eventWeight);

            h_m4l->Fill(ZZ_p4.M(), eventWeight);
            h_qt4l->Fill(ZZ_p4.Pt(), eventWeight);

            h_m12->Fill(Z1_p4.M(), eventWeight);
            h_m34->Fill(Z2_p4.M(), eventWeight);
            h_qt12->Fill(Z1_p4.Pt(), eventWeight);
            h_qt34->Fill(Z2_p4.Pt(), eventWeight);

            h_pt1->Fill(z_p4s[0].Pt(), eventWeight);
            h_pt2->Fill(z_p4s[1].Pt(), eventWeight);
            h_pt3->Fill(z_p4s[2].Pt(), eventWeight);
            h_pt4->Fill(z_p4s[3].Pt(), eventWeight);

            h_eta1->Fill(z_p4s[0].Eta(), eventWeight);
            h_eta2->Fill(z_p4s[1].Eta(), eventWeight);
            h_eta3->Fill(z_p4s[2].Eta(), eventWeight);
            h_eta4->Fill(z_p4s[3].Eta(), eventWeight);

            h_npv->Fill(nPV, eventWeight);
            h_m4l2->Fill(ZZ_p4.M(), eventWeight);
            h_nlep->Fill(nHZZLeps, eventWeight);

            if (systOn)
            {
                if (smearID)
                    h_sf->Fill(recoWeight);
                else if (smearPtData || smearPtMC)
                {
                    h_sf->Fill(lep1EnergySF[x]);
                    h_sf->Fill(lep1EnergySF[y]);
                    h_sf->Fill(lep2EnergySF[Z2.first]);
                    h_sf->Fill(lep2EnergySF[Z2.second]);
                }
            }
            else
                h_sf->Fill(eventWeight);
        }
    }
    file->Close();  delete file;

    return h_tot;
}




//////////////////////
// HELPER FUNCTIONS //
//////////////////////


double GetBinContentPtEta(const TH2 *hist, const TLorentzVector &p4)
{
    int xbin = GetXbin(hist, p4.Eta());
    int ybin = GetYbin(hist, p4.Pt());

    return hist->GetBinContent(xbin, ybin);
}

int GetXbin(const TH2 *hist, const double xval)
{
    int bin;
    int nbins = hist->GetNbinsX();

    if (xval >= hist->GetXaxis()->GetBinLowEdge(nbins))
        bin = nbins;
    else
    {
        for (int i = 1; i < nbins; i++)
        {
            if (xval >= hist->GetXaxis()->GetBinLowEdge(i)
                && xval < hist->GetXaxis()->GetBinLowEdge(i+1))
            {
                bin = i;
                break;
            }
        }
    }
    
    return bin;
}

int GetYbin(const TH2 *hist, const double yval)
{
    int bin;
    int nbins = hist->GetNbinsY();

    if (yval >= hist->GetYaxis()->GetBinLowEdge(nbins))
        bin = nbins;
    else
    {
        for (int i = 1; i < nbins; i++)
        {
            if (yval >= hist->GetYaxis()->GetBinLowEdge(i)
                && yval < hist->GetYaxis()->GetBinLowEdge(i+1))
            {
                bin = i;
                break;
            }
        }
    }

    return bin;
}
