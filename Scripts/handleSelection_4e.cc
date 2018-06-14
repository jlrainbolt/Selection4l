#include <vector>
#include <iostream>
#include <utility>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

TH1D* select_mumu(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_pt1, TH1F *h_pt2,
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv);
TH1D* select_ee(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_pt1, TH1F *h_pt2,
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv);
TH1D* select_4e(TString rootFile, TString suffix, TH1F *h_m4l, TH1F *h_m12, TH1F *h_m34,
        TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4, TH1F *h_npv);


using namespace std;

void handleSelection(const TString selection, const TString suffix, const TString id)
{
    // Choose electron or muon
    Bool_t selDilep = kFALSE, sel4lep = kFALSE;
    Bool_t selMuMu = kFALSE, selEE = kFALSE;
    Bool_t sel4e = kFALSE;
    if (selection == "mumu")
    {
        selDilep = kTRUE;
        selMuMu = kTRUE;
    }
    else if (selection == "ee")
    {
        selDilep = kTRUE;
        selEE = kTRUE;
    }
    else if (selection == "4e")
    {
        sel4lep = kTRUE;
        sel4e = kTRUE;
    }


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
    else if (selEE)
    {
        output = "ee_";     lepton = "electron";    Lepton = "Electron";
    }
    else if (sel4e)
    {
        output = "4e_";     lepton = "electron";    Lepton = "Electron";
    }
    output +=  suffix + "_" + id + ".root";


    // Histograms
    const unsigned M  = 13;   TString hname[M],             htitle[M];
    const unsigned LL = 0;    hname[LL] = "DilepMass";      htitle[LL] = "Di" + lepton + " Mass";
    const unsigned ZZ = 1;    hname[ZZ] = "4lepMass";       htitle[ZZ] = "4-" + lepton + " Mass";
    const unsigned Z1 = 2;    hname[Z1] = "Z1Mass";         htitle[Z1] = "Z_1 Mass";
    const unsigned Z2 = 3;    hname[Z2] = "Z2Mass";         htitle[Z2] = "Z_2 Mass";
    const unsigned P1 = 4;    hname[P1] = "Lep1Pt";         htitle[P1] = Lepton + " 1 Pt";
    const unsigned P2 = 5;    hname[P2] = "Lep2Pt";         htitle[P2] = Lepton + " 2 Pt";
    const unsigned P3 = 6;    hname[P3] = "Lep3Pt";         htitle[P3] = Lepton + " 3 Pt";
    const unsigned P4 = 7;    hname[P4] = "Lep4Pt";         htitle[P4] = Lepton + " 4 Pt";
    const unsigned E1 = 8;    hname[E1] = "Lep1Eta";        htitle[E1] = Lepton + " 1 Eta";
    const unsigned E2 = 9;    hname[E2] = "Lep2Eta";        htitle[E2] = Lepton + " 2 Eta";
    const unsigned E3 = 10;   hname[E3] = "Lep3Eta";        htitle[E3] = Lepton + " 3 Eta";
    const unsigned E4 = 11;   hname[E4] = "Lep4Eta";        htitle[E4] = Lepton + " 4 Eta";
    const unsigned PV = 12;   hname[PV] = "nPV";            htitle[PV] = "# of Primary Vertices";


    Int_t bins[M];      Double_t low[M],    up[M];
    bins[LL] = 30;      low[LL] = 75;       up[LL] = 105;
    bins[ZZ] = 30;      low[ZZ] = 75;       up[ZZ] = 105;
    bins[Z1] = 40;      low[Z1] = 40;       up[Z1] = 120;
    bins[Z2] = 60;      low[Z2] = 0;        up[Z2] = 120;
    bins[P1] = 75;      low[P1] = 0;        up[P1] = 150;
    bins[P2] = 50;      low[P2] = 0;        up[P2] = 100;
    bins[P3] = 50;      low[P3] = 0;        up[P3] = 100;
    bins[P4] = 50;      low[P4] = 0;        up[P4] = 100;
    bins[E1] = 50;      low[E1] = -2.5;     up[E1] = 2.5;
    bins[E2] = 50;      low[E2] = -2.5;     up[E2] = 2.5;
    bins[E3] = 50;      low[E3] = -2.5;     up[E3] = 2.5;
    bins[E4] = 50;      low[E4] = -2.5;     up[E4] = 2.5;
    bins[PV] = 51;      low[PV] = -0.5;     up[PV] = 50.5;

    TH1F *h[M];
    for (unsigned j = 0; j < M; j++)
        h[j] = new TH1F(hname[j] + "_" + suffix, htitle[j], bins[j], low[j], up[j]);


    // Read trees
    TH1D *hTotalEvents;
    if (selMuMu)
        hTotalEvents = select_mumu(path + dir + file, suffix, h[LL], h[P1], h[P2], h[E1], h[E2],
                h[PV]);
    else if (selEE)
        hTotalEvents = select_ee(path + dir + file, suffix, h[LL], h[P1], h[P2], h[E1], h[E2],
                h[PV]);
    else if (sel4e)
        hTotalEvents = select_4e(path + dir + file, suffix, h[ZZ], h[Z1], h[Z2], h[P1], h[P2],
                h[P3], h[P4], h[E1], h[E2], h[E3], h[E4], h[PV]);


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
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv)
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
    TTreeReaderValue<Int_t> nPV_(reader, "nPV");
    TTreeReaderArray<TLorentzVector> muonP4_(reader, "muonP4");
    TTreeReaderValue<std::vector<Short_t>> muonQ_(reader, "muonQ");
    TTreeReaderValue<std::vector<Float_t>> muonSF_(reader, "muonSF");
    TTreeReaderValue<std::vector<Float_t>> muonIDEff_(reader, "muonIDEff");
    TTreeReaderValue<std::vector<Float_t>> muonIsoEff_(reader, "muonTightIsoEff");
    TTreeReaderValue<std::vector<Float_t>> muonTriggerEffData_(reader, "muonTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> muonTriggerEffMC_(reader, "muonTriggerEffMC");
    TTreeReaderValue<std::vector<Bool_t>> muonIsStd_(reader, "muonPassStdCuts");
    TTreeReaderValue<UShort_t> nMuons_(reader, "nMuons");
    TTreeReaderValue<UShort_t> nStdMuons_(reader, "nStdMuons");

//  reader.GetTree()->Print();

    while (reader.Next())
    {
        Bool_t passTrigger = *passTrigger_;
        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        Int_t nPV = *nPV_;
        vector<TLorentzVector> muonP4;
        vector<Int_t> muonQ;
        vector<Float_t> muonSF, muonIDEff, muonIsoEff, muonTrigEffData, muonTrigEffMC;
        vector<Bool_t> muonIsStd;
        UInt_t nMuons = *nMuons_, nStdMuons = *nStdMuons_;

        if (passTrigger && nStdMuons == 2)// && nMuons == 2)
        {
            for (const TLorentzVector& muonP4__: muonP4_)
                muonP4.push_back(muonP4__);

            for (unsigned i = 0; i < nMuons; i++)
            {
                muonQ.push_back((*muonQ_)[i]);
                muonSF.push_back((*muonSF_)[i]);
                muonIDEff.push_back((*muonIDEff_)[i]);
                muonIsoEff.push_back((*muonIsoEff_)[i]);
                muonTrigEffData.push_back((*muonTriggerEffData_)[i]);
                muonTrigEffMC.push_back((*muonTriggerEffMC_)[i]);
                muonIsStd.push_back((*muonIsStd_)[i]);
            }

            // Apply momentum correction
            for (unsigned i = 0; i < nMuons; i++)
                muonP4[i].SetPtEtaPhiM(muonP4[i].Pt() * muonSF[i], muonP4[i].Eta(), muonP4[i].Phi(),
                        muonP4[i].M());

            for (unsigned j = 1; j < nMuons; j++)
            {
                if (muonQ[0] * muonQ[j] < 1
                        && muonP4[0].Pt() > PT_MIN && muonP4[j].Pt() > PT_MIN
                        && muonIsStd[0] && muonIsStd[j])    // remove this?
                {
                    Double_t mll = (muonP4[0] + muonP4[j]).M();
                    if (mll > MLL_MIN && mll < MLL_MAX)
                    {
                        eventWeight *= PUWeight;
                        eventWeight *= muonIDEff[0] * muonIDEff[j];
                        eventWeight *= muonIsoEff[0] * muonIsoEff[j];
                        if (muonTrigEffMC[0] > 0 || muonTrigEffMC[j] > 0)
                        {
                            eventWeight *= 1. - (1. - muonTrigEffData[0])*(1. - muonTrigEffData[j]);
                            eventWeight /= 1. - (1. - muonTrigEffMC[0])*(1. - muonTrigEffMC[j]);
                        }

                        h_tot->Fill(6);
                        h_tot->Fill(7, eventWeight);
                        h_mll->Fill(mll, eventWeight);
                        h_pt1->Fill(muonP4[0].Pt(), eventWeight);
                        h_eta1->Fill(muonP4[0].Eta(), eventWeight);
                        h_pt2->Fill(muonP4[j].Pt(), eventWeight);
                        h_eta2->Fill(muonP4[j].Eta(), eventWeight);
                        h_npv->Fill(nPV);
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
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv)
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
    TTreeReaderValue<Int_t> nPV_(reader, "nPV");
    TTreeReaderArray<TLorentzVector> elecP4_(reader, "electronP4");
    TTreeReaderValue<std::vector<Short_t>> elecQ_(reader, "electronQ");
    TTreeReaderValue<std::vector<Float_t>> elecSF_(reader, "electronSF");
    TTreeReaderValue<std::vector<Float_t>> elecRecoEff_(reader, "electronRecoEff");
    TTreeReaderValue<std::vector<Float_t>> elecTriggerEffData_(reader, "electronTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> elecTriggerEffMC_(reader, "electronTriggerEffMC");
    TTreeReaderValue<std::vector<Bool_t>> elecIsStd_(reader, "elecPassStdCuts");
    TTreeReaderValue<UShort_t> nElecs_(reader, "nElectrons");
    TTreeReaderValue<UShort_t> nStdElecs_(reader, "nStdElectrons");

//  reader.GetTree()->Print();

    while (reader.Next())
    {
        Bool_t passTrigger = *passTrigger_;
        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        Int_t nPV = *nPV_;
        vector<TLorentzVector> elecP4;
        vector<Int_t> elecQ;
        vector<Float_t> elecSF, elecRecoEff, elecTrigEffData, elecTrigEffMC;
        vector<Bool_t> elecIsStd;
        UInt_t nElecs = *nElecs_, nStdElecs = *nStdElecs_;

        if (passTrigger && nStdElecs == 2) // && nElecs == 2)
        {
            for (const TLorentzVector& elecP4__: elecP4_)
                elecP4.push_back(elecP4__);

            for (unsigned i = 0; i < nElecs; i++)
            {
                elecQ.push_back((*elecQ_)[i]);
                elecSF.push_back((*elecSF_)[i]);
                elecRecoEff.push_back((*elecRecoEff_)[i]);
                elecTrigEffData.push_back((*elecTriggerEffData_)[i]);
                elecTrigEffMC.push_back((*elecTriggerEffMC_)[i]);
                elecIsStd.push_back((*elecIsStd_)[i]);
            }

            // Apply energy correction
            for (unsigned i = 0; i < nElecs; i++)
                elecP4[i].SetPtEtaPhiM(elecP4[i].Pt() * elecSF[i], elecP4[i].Eta(), elecP4[i].Phi(),
                        elecP4[i].M());

            for (unsigned j = 1; j < nElecs; j++)
            {
                if (elecQ[0] * elecQ[j] < 1
                        && elecP4[0].Pt() > PT_MIN && elecP4[j].Pt() > PT_MIN
                        && elecIsStd[0] && elecIsStd[j]) 
                {
                    Double_t mll = (elecP4[0] + elecP4[j]).M();
                    if (mll > MLL_MIN && mll < MLL_MAX)
                    {
                        eventWeight *= PUWeight;
                        eventWeight *= elecRecoEff[0] * elecRecoEff[j];
                        if (elecTrigEffMC[0] > 0 || elecTrigEffMC[j] > 0)
                        {
                            eventWeight *= 1. - (1. - elecTrigEffData[0])*(1. - elecTrigEffData[j]);
                            eventWeight /= 1. - (1. - elecTrigEffMC[0])*(1. - elecTrigEffMC[j]);
                        }

                        h_tot->Fill(6);
                        h_tot->Fill(7, eventWeight);
                        h_mll->Fill(mll, eventWeight);
                        h_pt1->Fill(elecP4[0].Pt(), eventWeight);
                        h_eta1->Fill(elecP4[0].Eta(), eventWeight);
                        h_pt2->Fill(elecP4[j].Pt(), eventWeight);
                        h_eta2->Fill(elecP4[j].Eta(), eventWeight);
                        h_npv->Fill(nPV);
                        break;
                    }
                }
            }
        }
    }
    file->Close();  delete file;

    return h_tot;
}




TH1D* select_4e(TString rootFile, TString suffix, TH1F *h_m4l, TH1F *h_m12, TH1F *h_m34,
        TH1F *h_pt1, TH1F *h_pt2, TH1F *h_pt3, TH1F *h_pt4,
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_eta3, TH1F *h_eta4, TH1F *h_npv)
{
    // Cuts, etc.
    Double_t MZ1_MIN = 40, MZ_MIN = 12, MZ_MAX = 120, M4L_MIN = 70;
    Double_t PT_MIN = 7, PT1_MIN = 20, PT2_MIN = 10;
    Double_t Z_MASS = 91.2;

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
    TTreeReaderValue<Int_t> nPV_(reader, "nPV");
    TTreeReaderArray<TLorentzVector> elecP4_(reader, "electronP4");
    TTreeReaderValue<std::vector<Short_t>> elecQ_(reader, "electronQ");
    TTreeReaderValue<std::vector<Float_t>> elecSF_(reader, "electronSF");
    TTreeReaderValue<std::vector<Float_t>> elecRecoEff_(reader, "electronRecoEff");
    TTreeReaderValue<std::vector<Float_t>> elecTriggerEffData_(reader, "electronTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> elecTriggerEffMC_(reader, "electronTriggerEffMC");
    TTreeReaderValue<std::vector<Bool_t>> elecIsStd_(reader, "elecPassStdCuts");
    TTreeReaderValue<UShort_t> nElecs_(reader, "nElectrons");
    TTreeReaderValue<UShort_t> nStdElecs_(reader, "nStdElectrons");

//  reader.GetTree()->Print();

    while (reader.Next())
    {
        // Stuff from tree
        Bool_t passTrigger = *passTrigger_;
        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        Int_t nPV = *nPV_;
        vector<TLorentzVector> elecP4;
        vector<Int_t> elecQ;
        vector<Float_t> elecSF, elecRecoEff, elecTrigEffData, elecTrigEffMC;
        vector<Bool_t> elecIsStd;
        UInt_t nElecs = *nElecs_, nStdElecs = *nStdElecs_;

        // Stuff for selection
        vector<pair<unsigned, unsigned>> zCands;

        if (passTrigger && nElecs >= 4)
        {
            for (const TLorentzVector& elecP4__: elecP4_)
                elecP4.push_back(elecP4__);

            for (unsigned i = 0; i < nElecs; i++)
            {
                elecQ.push_back((*elecQ_)[i]);
                elecSF.push_back((*elecSF_)[i]);
                elecRecoEff.push_back((*elecRecoEff_)[i]);
                elecTrigEffData.push_back((*elecTriggerEffData_)[i]);
                elecTrigEffMC.push_back((*elecTriggerEffMC_)[i]);
                elecIsStd.push_back((*elecIsStd_)[i]);
            }

            // Apply energy correction
            for (unsigned i = 0; i < nElecs; i++)
                elecP4[i].SetPtEtaPhiM(elecP4[i].Pt() * elecSF[i], elecP4[i].Eta(), elecP4[i].Phi(),
                        elecP4[i].M());

            // Find Z candidates
            for (unsigned i = 0; i < nElecs-1; i++)
            {
                for (unsigned j = i+1; j < nElecs; j++)
                {
                    if (elecQ[i] * elecQ[j] < 1
                            && elecP4[i].Pt() > PT_MIN && elecP4[j].Pt() > PT_MIN)
                    {
                        Double_t mll = (elecP4[0] + elecP4[j]).M();
                        if (mll > MZ_MIN && mll < MZ_MAX)
                            zCands.push_back(make_pair(i, j));
                    }
                }
            }

            for (unsigned i = 0; i < zCands.size()-1; i++)
            {
                for (unsigned j = j+1; j < zCands.size(); j++)
                {
                    // Make sure pairs do not overlap
                    if (zCands[i].first != zCands[j].first
                            && zCands[i].first != zCands[j].second
                            && zCands[i].second != zCands[j].first
                            && zCands[i].second != zCands[j].second)
                    {
                        Double_t mll_i = (elecP4[zCands[i].first] + elecP4[zCands[i].second]).M();
                        Double_t mll_j = (elecP4[zCands[j].first] + elecP4[zCands[j].second]).M();

                        // Choose Z1 closest to nominal Z mass
                        std::pair<unsigned, unsigned> Z1, Z2;
                        std::pair<TLorentzVector, TLorentzVector> Z1_p4, Z2_p4;
                        std::pair<Int_t, Int_t> Z1_q, Z2_q;
                        Double_t Z1_m, Z2_m, m4l;
                        if (abs(mll_i - Z_MASS) < abs(mll_j - Z_MASS))
                        {
                            Z1 = zCands[i]; Z2 = zCands[j];
                        }
                        else
                        {
                            Z1 = zCands[j]; Z2 = zCands[i];
                        }
                        Z1_p4 = std::make_pair(elecP4[Z1.first], elecP4[Z1.second]);
                        Z2_p4 = std::make_pair(elecP4[Z2.first], elecP4[Z2.second]);
                        Z1_q = std::make_pair(elecQ[Z1.first], elecQ[Z1.second]);
                        Z2_q = std::make_pair(elecQ[Z2.first], elecQ[Z2.second]);
                        Z1_m = (Z1_p4.first + Z1_p4.second).M();
                        Z2_m = (Z2_p4.first + Z2_p4.second).M();
                        m4l = (Z1_p4.first + Z1_p4.second + Z2_p4.first + Z2_p4.second).M();

                        // Mass requirements
                        if (Z1_m > MZ1_MIN && m4l > M4L_MIN)
                        {
                            // Pt2 requirement
                            if ((Z1_p4.first.Pt() > PT2_MIN || Z1_p4.second.Pt() > PT2_MIN
                                        || Z2_p4.first.Pt() > PT2_MIN || Z2_p4.first.Pt() > PT2_MIN)
                                    && 
                        }
                    }

                }

            }


        }
    }
    file->Close();  delete file;

    return h_tot;
}
