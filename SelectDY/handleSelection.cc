#include <vector>
#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

TH1D* select_mumu(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_pt1, TH1F *h_pt2,
        TH1F *h_eta1, TH1F *h_eta2);
TH1D* select_ee(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_pt1, TH1F *h_pt2,
        TH1F *h_eta1, TH1F *h_eta2);


using namespace std;

void handleSelection(const TString selection, const TString suffix, const TString id)
{
    // Choose electron or muon
    Bool_t selMuMu;
    if (selection == "mumu")
        selMuMu = kTRUE;
    else if (selection == "ee")
        selMuMu = kFALSE;


    // Path to directory where ALL trees are stored
    TString path = "root://cmsxrootd.fnal.gov//store/user/jrainbol/Trees/2016/";


    // Names of ROOT files and trees for each sample
    TString dir     = suffix + "/";
    TString file    = suffix + "_" + id + ".root";


    // Name of output file
    TString output, lepton, Lepton;
    if (selMuMu)
    {
        output = "mumu_";   lepton = "muon";        Lepton = "Muon";
    }
    else
    {
        output = "ee_";     lepton = "electron";    Lepton = "Electron";
    }
    output +=  suffix + "_" + id + ".root";


    // Histograms
    const unsigned M  = 5;    TString hname[M],               htitle[M];
    const unsigned LL = 0;    hname[LL] = "DileptonMass";     htitle[LL] = "Di" + lepton + " Mass";
    const unsigned P1 = 1;    hname[P1] = "Lepton1Pt";        htitle[P1] = Lepton + " 1 Pt";
    const unsigned E1 = 2;    hname[E1] = "Lepton1Eta";       htitle[E1] = Lepton + " 1 Eta";
    const unsigned P2 = 3;    hname[P2] = "Lepton2Pt";        htitle[P2] = Lepton + " 2 Pt";
    const unsigned E2 = 4;    hname[E2] = "Lepton2Eta";       htitle[E2] = Lepton + " 2 Eta";


    Int_t bins[M];      Double_t low[M],    up[M];
    bins[LL] = 30;      low[LL] = 75;       up[LL] = 105;
    bins[P1] = 75;      low[P1] = 0;        up[P1] = 150;
    bins[P2] = 50;      low[P2] = 0;        up[P2] = 100;
    bins[E1] = 50;      low[E1] = -2.5;     up[E1] = 2.5;
    bins[E2] = 50;      low[E2] = -2.5;     up[E2] = 2.5;

    TH1F *h[M];
    for (unsigned j = 0; j < M; j++)
        h[j] = new TH1F(hname[j] + "_" + suffix, htitle[j], bins[j], low[j], up[j]);


    // Read trees
    TH1D *hTotalEvents;
    if (selMuMu)
        hTotalEvents = select_mumu(path + dir + file, suffix, h[LL], h[P1], h[P2], h[E1], h[E2]);
    else
        hTotalEvents = select_ee(path + dir + file, suffix, h[LL], h[P1], h[P2], h[E1], h[E2]);


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




TH1D* select_mumu(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_pt1, TH1F *h_pt2,
        TH1F *h_eta1, TH1F *h_eta2)
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
    TTreeReaderValue<Bool_t> passTrigger_(reader, "passTrigger");
    TTreeReaderValue<Float_t> eventWeight_(reader, "eventWeight");
    TTreeReaderValue<Float_t> PUWeight_(reader, "PUWeight");
    TTreeReaderArray<TLorentzVector> muonP4_(reader, "muonP4");
    TTreeReaderValue<std::vector<Short_t>> muonQ_(reader, "muonQ");
    TTreeReaderValue<std::vector<Float_t>> muonSF_(reader, "muonSF");
    TTreeReaderValue<std::vector<Float_t>> muonIDEff_(reader, "muonIDEff");
    TTreeReaderValue<std::vector<Float_t>> muonIsoEff_(reader, "muonIsoEff");
    TTreeReaderValue<std::vector<Float_t>> muonTriggerEff_(reader, "muonTriggerEff");
//  TTreeReaderValue<std::vector<Bool_t>> muonIsStd_(reader, "muonPassStdCuts");
    TTreeReaderValue<UShort_t> nMuons_(reader, "nMuons");
    TTreeReaderValue<UShort_t> nStdMuons_(reader, "nStdMuons");

//  reader.GetTree()->Print();

    while (reader.Next())
    {
        Bool_t passTrigger = *passTrigger_;
        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        vector<TLorentzVector> muonP4;
        vector<Int_t> muonQ;
        vector<Float_t> muonSF, muonIDEff, muonIsoEff, muonTriggerEff;
//      vector<Bool_t> muonIsStd;
        UInt_t nMuons = *nMuons_, nStdMuons = *nStdMuons_;

        if (passTrigger && nStdMuons == 2 && nMuons == 2)
        {
            for (const TLorentzVector& muonP4__: muonP4_)
                muonP4.push_back(muonP4__);

            for (unsigned i = 0; i < nMuons; i++)
            {
                muonQ.push_back((*muonQ_)[i]);
                muonSF.push_back((*muonSF_)[i]);
                muonIDEff.push_back((*muonIDEff_)[i]);
                muonIsoEff.push_back((*muonIsoEff_)[i]);
                muonTriggerEff.push_back((*muonTriggerEff_)[i]);
//              muonIsStd.push_back((*muonIsStd_)[i]);
            }

            // Apply momentum correction
            for (unsigned i = 0; i < nMuons; i++)
                muonP4[i].SetPtEtaPhiM(muonP4[i].Pt() * muonSF[i], muonP4[i].Eta(), muonP4[i].Phi(),
                        muonP4[i].M());

            for (unsigned j = 1; j < nMuons; j++)
            {
                if (muonQ[0] * muonQ[j] < 1
                        && muonP4[0].Pt() > PT_MIN
                        && muonP4[j].Pt() > PT_MIN) 
                {
                    Double_t mll = (muonP4[0] + muonP4[j]).M();
                    if (mll > MLL_MIN && mll < MLL_MAX)
                    {
                        eventWeight *= PUWeight;
                        eventWeight *= muonIDEff[0] * muonIDEff[j];
                        eventWeight *= muonIsoEff[0] * muonIsoEff[j];
                        eventWeight *= 1 - (1 - muonTriggerEff[0]) * (1 - muonTriggerEff[j]);

                        h_tot->Fill(6);
                        h_tot->Fill(7, eventWeight);
                        h_mll->Fill(mll, eventWeight);
                        h_pt1->Fill(muonP4[0].Pt(), eventWeight);
                        h_eta1->Fill(muonP4[0].Eta(), eventWeight);
                        h_pt2->Fill(muonP4[j].Pt(), eventWeight);
                        h_eta2->Fill(muonP4[j].Eta(), eventWeight);
                        break;
                    }
                }
            }
        }
    }
    file->Close();  delete file;

    return h_tot;
}




TH1D* select_ee(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_pt1, TH1F *h_pt2,
        TH1F *h_eta1, TH1F *h_eta2)
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
    TTreeReaderValue<Bool_t> passTrigger_(reader, "passTrigger");
    TTreeReaderArray<TLorentzVector> elecP4_(reader, "electronP4");
    TTreeReaderValue<std::vector<Short_t>> elecQ_(reader, "electronQ");
//  TTreeReaderValue<std::vector<Bool_t>> elecIsStd_(reader, "elecPassStdCuts");
    TTreeReaderValue<UShort_t> nElecs_(reader, "nElectrons");
    TTreeReaderValue<UShort_t> nStdElecs_(reader, "nStdElectrons");

//  reader.GetTree()->Print();

    while (reader.Next())
    {
        Bool_t passTrigger = *passTrigger_;
        vector<TLorentzVector> elecP4;
        vector<Int_t> elecQ;
//      vector<Bool_t> elecIsStd;
        UInt_t nElecs = *nElecs_, nStdElecs = *nStdElecs_;

        if (passTrigger && nStdElecs == 2 && nElecs == 2)
        {
            for (const TLorentzVector& elecP4__: elecP4_)
                elecP4.push_back(elecP4__);

            for (unsigned i = 0; i < nElecs; i++)
            {
                elecQ.push_back((*elecQ_)[i]);
//              elecIsStd.push_back((*elecIsStd_)[i]);
            }

            for (unsigned j = 1; j < nElecs; j++)
            {
                if (elecQ[0] * elecQ[j] < 1
                        && elecP4[0].Pt() > PT_MIN
                        && elecP4[j].Pt() > PT_MIN) 
                {
                    Double_t mll = (elecP4[0] + elecP4[j]).M();
                    if (mll > MLL_MIN && mll < MLL_MAX)
                    {
                        h_tot->Fill(6);
                        h_mll->Fill(mll);
                        h_pt1->Fill(elecP4[0].Pt());
                        h_eta1->Fill(elecP4[0].Eta());
                        h_pt2->Fill(elecP4[j].Pt());
                        h_eta2->Fill(elecP4[j].Eta());
                        break;
                    }
                }
            }
        }
    }
    file->Close();  delete file;

    return h_tot;
}


// Extremely deprecated
/*
# ifndef __CINT__
int main(int argc, char *argv[])
{
    if (argc == 2)
    {
        testTree(argv[1]);

        return 0;
    }
    else
        return 1;
}
# endif
*/
