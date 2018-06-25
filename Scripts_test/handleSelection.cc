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

TH1D* select_mumu(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt, TH1F *h_pt1,
        TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv, TH1F *h_njets, TH1F *h_n_bad,
        TH1F *x_pt, TH1F *x_eta, TH1F *x_q, TH1F *x_iso, TH1F *x_tiso, TH1F *x_d0, TH1F *x_dz, 
        TH1F *x_trig, TH1F *x_gbl, TH1F *x_nstn, TH1F *x_npix, TH1F *x_nlyr, TH1F *x_nhits, 
        TH1F *x_chi2);
TH1D* select_ee(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt, TH1F *h_pt1,
        TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv, TH1F *h_njets, TH1F *h_n_bad,
        TH1F *x_pt, TH1F *x_eta, TH1F *x_q, TH1F *x_iso, TH1F *x_tiso, TH1F *x_d0, TH1F *x_dz,
        TH1F *x_trig);
        //TH1F *x_sceta, TH1F *x_sieie, TH1F *x_hovere, TH1F *x_einv, TH1F *x_nmiss, TH1F *x_conv);


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
    const unsigned M  = 23;   TString hname[M],             htitle[M];
    const unsigned LL = 0;    hname[LL] = "DileptonMass";   htitle[LL] = "Di" + lepton + " Mass";
    const unsigned QT = 1;    hname[QT] = "DileptonPt";     htitle[QT] = "Di" + lepton + " Pt";
    const unsigned P1 = 2;    hname[P1] = "Lepton1Pt";      htitle[P1] = Lepton + " 1 Pt";
    const unsigned E1 = 3;    hname[E1] = "Lepton1Eta";     htitle[E1] = Lepton + " 1 Eta";
    const unsigned P2 = 4;    hname[P2] = "Lepton2Pt";      htitle[P2] = Lepton + " 2 Pt";
    const unsigned E2 = 5;    hname[E2] = "Lepton2Eta";     htitle[E2] = Lepton + " 2 Eta";
    const unsigned PV = 6;    hname[PV] = "nPV";            htitle[PV] = "# Primary Vertices";
    const unsigned JT = 7;    hname[JT] = "nJets";          htitle[JT] = "# Jets";
    const unsigned NX = 8;    hname[NX] = "nExtraLeps";     htitle[NX] = "# extra " + lepton+"s";
    const unsigned PX = 9;    hname[PX] = "xLepPt";         htitle[PX] = "Extra " + lepton+"s Pt";
    const unsigned EX = 10;   hname[EX] = "xLepEta";        htitle[EX] = "Extra " + lepton+"s Eta";
    const unsigned QX = 11;   hname[QX] = "xLepQ";          htitle[QX] = "Extra " + lepton+"s Q";
    const unsigned CI = 12;   hname[CI] = "xLepIso";        htitle[CI] = "Extra " + lepton+"s Iso/Pt";
    const unsigned TI = 13;   hname[TI] = "xLepTrkIso";     htitle[TI] = "Extra " + lepton+"s Track Iso/Pt";
    const unsigned D0 = 14;   hname[D0] = "xLepD0";         htitle[D0] = "Extra " + lepton+"s d_0";
    const unsigned DZ = 15;   hname[DZ] = "xLepDz";         htitle[DZ] = "Extra " + lepton+"s d_z";
    const unsigned TG = 16;   hname[TG] = "xLepIsTrig";     htitle[TG] = "Extra " + lepton+"s Trig";

    // these are different for muons/electrons
    const unsigned GB = 17,   SE = 17;
    const unsigned MS = 18,   SI = 18;
    const unsigned PH = 19,   HE = 19;
    const unsigned TL = 20,   EI = 20;
    const unsigned VH = 21,   MH = 21;
    const unsigned X2 = 22,   IC = 22;
    if (selMuMu)
    {
        hname[GB] = "xLepIsGLB";            htitle[GB] = "Extra muons GLB";
        hname[MS] = "xLepNMatchStn";        htitle[MS] = "Extra muons NMatchStn";
        hname[PH] = "xLepNPixHits";         htitle[PH] = "Extra muons NPixHits";
        hname[TL] = "xLepNTkLayers";        htitle[TL] = "Extra muons NTkLayers";
        hname[VH] = "xLepNValidHits";       htitle[VH] = "Extra muons NValidHits";
        hname[X2] = "xLepMuNChi2";          htitle[X2] = "Extra muons muNChi2";
    }
    else
    {
        hname[SE] = "xLepScEta";            htitle[SE] = "Extra electrons sc eta";
        hname[EI] = "xLepSieie";            htitle[EI] = "Extra electrons sieie";
        hname[HE] = "xLepHOverE";           htitle[HE] = "Extra electrons H over E";
        hname[EI] = "xLepEnergyInv";        htitle[EI] = "Extra electrons inverse energy";
        hname[MH] = "xLepNMissHits";        htitle[MH] = "Extra electrons missing hits";
        hname[IC] = "xLepIsConv";           htitle[IC] = "Extra electrons converted";
    }

    // electrons only
//  const unsigned DE = 23;   hname[DE] = "xLepDEtaIn";     htitle[DE] = "Extra electrons dEtaIn";
//  const unsigned DP = 24;   hname[DE] = "xLepDPhiIn";     htitle[DE] = "Extra electrons dPhiIn";

    Int_t bins[M];      Double_t low[M],    up[M];
    bins[LL] = 30;      low[LL] = 75;       up[LL] = 105;
    bins[QT] = 75;      low[QT] = 0;        up[QT] = 150;
    bins[P1] = 75;      low[P1] = 0;        up[P1] = 150;
    bins[P2] = 50;      low[P2] = 0;        up[P2] = 100;
    bins[E1] = 50;      low[E1] = -2.5;     up[E1] = 2.5;
    bins[E2] = 50;      low[E2] = -2.5;     up[E2] = 2.5;
    bins[PV] = 51;      low[PV] = -0.5;     up[PV] = 50.5;
    bins[JT] = 11;      low[JT] = -0.5;     up[JT] = 10.5;
    bins[NX] = 7;       low[NX] = -0.5;     up[NX] = 6.5;
    if (selMuMu)
    {
//      bins[PX] = 35;      low[PX] = 15;       up[PX] = 50;
        bins[PX] = 75;      low[PX] = 0;        up[PX] = 75;
        bins[CI] = 30;      low[CI] = 0;        up[CI] = 0.15;
        bins[TI] = 30;      low[TI] = 0;        up[TI] = 0.15;
    }
    else
    {
        bins[PX] = 75;      low[PX] = 0;        up[PX] = 75;
        bins[CI] = 75;      low[CI] = 0;        up[CI] = 1.5;
        bins[TI] = 75;      low[TI] = 0;        up[TI] = 1.5;
    }
    bins[EX] = 50;      low[EX] = -2.5;     up[EX] = 2.5;
    bins[QX] = 3;       low[QX] = -1.5;     up[QX] = 1.5;
    bins[D0] = 50;      low[D0] = -2.5;     up[D0] = 2.5;
    bins[DZ] = 50;      low[DZ] = -50;      up[DZ] = 50;
    if (selMuMu)
    {
        bins[GB] = 2;       low[GB] = -0.5;     up[GB] = 1.5;
        bins[TG] = 2;       low[TG] = -0.5;     up[TG] = 1.5;
        bins[MS] = 7;       low[MS] = -0.5;     up[MS] = 6.5;
        bins[PH] = 7;       low[PH] = -0.5;     up[PH] = 6.5;
        bins[TL] = 19;      low[TL] = -0.5;     up[TL] = 18.5;
        bins[VH] = 50;      low[VH] = 0;        up[VH] = 50;
        bins[X2] = 31;      low[X2] = -0.5;     up[X2] = 30.5;
    }
    else
    {
        bins[SE] = 50;      low[SE] = -2.5;     up[SE] = 2.5;
        bins[EI] = 50;      low[EI] = -0.05;    up[EI] = 0.05;
        bins[HE] = 50;      low[HE] = -0.1;     up[HE] = 0.1;
        bins[HE] = 50;      low[HE] = -0.1;     up[HE] = 0.1;
    }
//  bins[DE]
//  bins[DP]

    TH1F *h[M];
    for (unsigned j = 0; j < M; j++)
        h[j] = new TH1F(hname[j] + "_" + suffix, htitle[j], bins[j], low[j], up[j]);


    // Read trees
    TH1D *hTotalEvents;
    if (selMuMu)
        hTotalEvents = select_mumu(path + dir + file, suffix, h[LL], h[QT], h[P1], h[P2], h[E1],
                h[E2], h[PV], h[JT], h[NX], h[PX], h[EX], h[QX], h[CI], h[TI], h[D0], h[DZ], h[TG], 
                h[GB], h[MS], h[PH], h[TL], h[VH], h[X2]);
    else
        hTotalEvents = select_ee(path + dir + file, suffix, h[LL], h[QT], h[P1], h[P2], h[E1],
                h[E2], h[PV], h[JT], h[NX], h[PX], h[EX], h[QX], h[CI], h[TI], h[D0], h[DZ], h[TG]);


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




TH1D* select_mumu(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt, TH1F *h_pt1,
        TH1F *h_pt2, TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv, TH1F *h_njets, TH1F *h_n_bad, 
        TH1F *x_pt, TH1F *x_eta, TH1F *x_q, TH1F *x_iso, TH1F *x_tiso, TH1F *x_d0, TH1F *x_dz,
        TH1F *x_trig, TH1F *x_gbl, TH1F *x_nstn, TH1F *x_npix, TH1F *x_nlyr, TH1F *x_nhits,
        TH1F *x_chi2)
{
    // Cuts
    Double_t MLL_MIN = 80, MLL_MAX = 100;
    Double_t PT_MIN = 25, ISO_MAX = 0.15;
    Double_t E_PT_MIN = 15, E_PT_MAX = 50, LOOSE_PT_MIN = 10;

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
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nJets_(reader, "nJets");
    TTreeReaderArray<TLorentzVector> muonP4_(reader, "muonP4");
    TTreeReaderValue<std::vector<Short_t>> muonQ_(reader, "muonQ");
    TTreeReaderValue<std::vector<Float_t>> muonIso_(reader, "muonCombIso");
    TTreeReaderValue<std::vector<Float_t>> muonTrkIso_(reader, "muonTrkIso");
    TTreeReaderValue<std::vector<Float_t>> muonD0_(reader, "muonD0");
    TTreeReaderValue<std::vector<Float_t>> muonDz_(reader, "muonDz");
    TTreeReaderValue<std::vector<Float_t>> muonSF_(reader, "muonSF");
    TTreeReaderValue<std::vector<Float_t>> muonIDEff_(reader, "muonIDEff");
    TTreeReaderValue<std::vector<Float_t>> muonIsoEff_(reader, "muonTightIsoEff");
    TTreeReaderValue<std::vector<Float_t>> muonTriggerEffData_(reader, "muonTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> muonTriggerEffMC_(reader, "muonTriggerEffMC");
    TTreeReaderValue<std::vector<Bool_t>> muonIsStd_(reader, "muonPassStdCuts");
    TTreeReaderValue<std::vector<Bool_t>> muonIsPF_(reader, "muonIsPF");
    TTreeReaderValue<std::vector<Bool_t>> muonIsGLB_(reader, "muonIsGLB");
    TTreeReaderValue<std::vector<Bool_t>> muonTrigger_(reader, "muonPassTrigger");
    TTreeReaderValue<std::vector<UShort_t>> muonNMatchStn_(reader, "muonNMatchStn");
    TTreeReaderValue<std::vector<UShort_t>> muonNPixHits_(reader, "muonNPixHits");
    TTreeReaderValue<std::vector<UShort_t>> muonNTkLayers_(reader, "muonNTkLayers");
    TTreeReaderValue<std::vector<UShort_t>> muonNValidHits_(reader, "muonNValidHits");
    TTreeReaderValue<std::vector<Float_t>> muonMuNChi2_(reader, "muonMuNChi2");
    TTreeReaderValue<UShort_t> nMuons_(reader, "nMuons");
    TTreeReaderValue<UShort_t> nStdMuons_(reader, "nStdMuons");

//  reader.GetTree()->Print();

    while (reader.Next())
    {
        // From tree
        Bool_t passTrigger = *passTrigger_;
        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        UShort_t nPV = *nPV_, nJets = *nJets_;
        vector<TLorentzVector> muonP4;
        vector<Short_t> muonQ;
        vector<Float_t> muonIso, muonTrkIso, muonSF, muonD0, muonDz;
        vector<Float_t> muonIDEff, muonIsoEff, muonTrigEffData, muonTrigEffMC, muonMuNChi2;
        vector<Bool_t> muonIsStd, muonIsPF, muonIsGLB, muonTrigger;
        vector<UShort_t> muonNMatchStn, muonNPixHits, muonNTkLayers, muonNValidHits;
        UShort_t nMuons = *nMuons_, nStdMuons = *nStdMuons_;

        // For selection
        vector<Bool_t> muonIsIso, muonIsLoose;
        Short_t nGoodMuons = nMuons, nIsoMuons = 0, nLooseMuons = 0;


        if (passTrigger && nStdMuons == 2)  // && nMuons > 2)
        {
            for (const TLorentzVector& muonP4__: muonP4_)
                muonP4.push_back(muonP4__);

            for (unsigned i = 0; i < nMuons; i++)
            {
                muonQ.push_back((*muonQ_)[i]);
                muonIso.push_back((*muonIso_)[i]);
                muonTrkIso.push_back((*muonTrkIso_)[i]);
                muonSF.push_back((*muonSF_)[i]);
                muonD0.push_back((*muonD0_)[i]);
                muonDz.push_back((*muonDz_)[i]);
                muonIDEff.push_back((*muonIDEff_)[i]);
                muonIsoEff.push_back((*muonIsoEff_)[i]);
                muonTrigEffData.push_back((*muonTriggerEffData_)[i]);
                muonTrigEffMC.push_back((*muonTriggerEffMC_)[i]);
                muonIsStd.push_back((*muonIsStd_)[i]);
                muonIsPF.push_back((*muonIsPF_)[i]);
                muonIsGLB.push_back((*muonIsGLB_)[i]);
                muonTrigger.push_back((*muonTrigger_)[i]);
                muonNMatchStn.push_back((*muonNMatchStn_)[i]);
                muonNPixHits.push_back((*muonNPixHits_)[i]);
                muonNTkLayers.push_back((*muonNTkLayers_)[i]);
                muonNValidHits.push_back((*muonNValidHits_)[i]);
                muonMuNChi2.push_back((*muonMuNChi2_)[i]);
            }

            // Apply momentum correction
            for (unsigned i = 0; i < nMuons; i++)
                muonP4[i].SetPtEtaPhiM(muonP4[i].Pt() * muonSF[i], muonP4[i].Eta(), muonP4[i].Phi(),
                        muonP4[i].M());

            // Find isolated and loose muons
            for (unsigned i = 0; i < nMuons; i++)
            {
//              if (muonIso[i] / muonP4[i].Pt() < ISO_MAX)
//              {
//                  muonIsIso.push_back(kTRUE);
//                  nIsoMuons++;
//              }
//              else
//                  muonIsIso.push_back(kFALSE);
  
                if (muonIsPF[i] && muonIsGLB[i] 
                        && muonP4[i].Pt() > LOOSE_PT_MIN)
                {
                    muonIsLoose.push_back(kTRUE);
                    nLooseMuons++;
                }
                else
                    muonIsLoose.push_back(kFALSE);
            }

//          if (nIsoMuons == nStdMuons)
            if (nLooseMuons == nStdMuons)
            {
                for (unsigned j = 1; j < nMuons; j++)
                {
                    if (muonQ[0] * muonQ[j] < 1
                            && muonP4[0].Pt() > PT_MIN && muonP4[j].Pt() > PT_MIN
                            && muonIsStd[0] && muonIsStd[j]
                            && muonIsPF[0] && muonIsPF[j])
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
                            h_qt->Fill((muonP4[0] + muonP4[j]).Pt(), eventWeight);
                            h_pt1->Fill(muonP4[0].Pt(), eventWeight);
                            h_eta1->Fill(muonP4[0].Eta(), eventWeight);
                            h_pt2->Fill(muonP4[j].Pt(), eventWeight);
                            h_eta2->Fill(muonP4[j].Eta(), eventWeight);
                            h_npv->Fill(nPV, eventWeight);
                            h_njets->Fill(nJets, eventWeight);

                            // Fill bad muon histograms
                            for (unsigned k = 0; k < nMuons; k++)
                            {
                                if (!muonIsLoose[k] && !muonIsStd[k])
//                              if (muonIsIso[k] && !muonIsStd[k]
//                                      && muonP4[k].Pt() > E_PT_MIN && muonP4[k].Pt() < E_PT_MAX)
                                {
                                    x_pt->Fill(muonP4[k].Pt(), eventWeight);
                                    x_eta->Fill(muonP4[k].Eta(), eventWeight);
                                    x_q->Fill(muonQ[k], eventWeight);
                                    x_iso->Fill(muonIso[k] / muonP4[k].Pt(), eventWeight);
                                    x_tiso->Fill(muonTrkIso[k] / muonP4[k].Pt(), eventWeight);
                                    x_d0->Fill(muonD0[k], eventWeight);
                                    x_dz->Fill(muonDz[k], eventWeight);
                                    x_gbl->Fill(muonIsGLB[k], eventWeight);
                                    x_trig->Fill(muonTrigger[k], eventWeight);
                                    x_nstn->Fill(muonNMatchStn[k], eventWeight);
                                    x_npix->Fill(muonNPixHits[k], eventWeight);
                                    x_nlyr->Fill(muonNTkLayers[k], eventWeight);
                                    x_nhits->Fill(muonNValidHits[k], eventWeight);
                                    x_chi2->Fill(muonMuNChi2[k], eventWeight);
                                }
                            }
//                          h_n_bad->Fill(nIsoMuons - nStdMuons, eventWeight);
                            h_n_bad->Fill(nLooseMuons - nStdMuons, eventWeight);

                            break;
                        }
                    }
                }
            }
        }
    }
    file->Close();  delete file;

    return h_tot;
}




TH1D* select_ee(TString rootFile, TString suffix, TH1F *h_mll, TH1F *h_qt, TH1F *h_pt1, TH1F *h_pt2,
        TH1F *h_eta1, TH1F *h_eta2, TH1F *h_npv, TH1F *h_njets, TH1F *h_n_bad, TH1F *x_pt, 
        TH1F *x_eta, TH1F *x_q, TH1F *x_iso, TH1F *x_tiso, TH1F *x_d0, TH1F *x_dz, TH1F *x_trig)
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
    TTreeReaderValue<UShort_t> nPV_(reader, "nPV");
    TTreeReaderValue<UShort_t> nJets_(reader, "nJets");
    TTreeReaderArray<TLorentzVector> elecP4_(reader, "electronP4");
    TTreeReaderValue<std::vector<Short_t>> elecQ_(reader, "electronQ");
    TTreeReaderValue<std::vector<Float_t>> elecIso_(reader, "electronCombIso");
    TTreeReaderValue<std::vector<Float_t>> elecTrkIso_(reader, "electronTrkIso");
    TTreeReaderValue<std::vector<Float_t>> elecD0_(reader, "electronD0");
    TTreeReaderValue<std::vector<Float_t>> elecDz_(reader, "electronDz");
    TTreeReaderValue<std::vector<Float_t>> elecSF_(reader, "electronSF");
    TTreeReaderValue<std::vector<Float_t>> elecRecoEff_(reader, "electronRecoEff");
    TTreeReaderValue<std::vector<Float_t>> elecTriggerEffData_(reader, "electronTriggerEffData");
    TTreeReaderValue<std::vector<Float_t>> elecTriggerEffMC_(reader, "electronTriggerEffMC");
    TTreeReaderValue<std::vector<Bool_t>> elecIsStd_(reader, "electronPassStdCuts");
    TTreeReaderValue<std::vector<Bool_t>> elecTrigger_(reader, "electronPassTrigger");
    TTreeReaderValue<UShort_t> nElecs_(reader, "nElectrons");
    TTreeReaderValue<UShort_t> nStdElecs_(reader, "nStdElectrons");

//  reader.GetTree()->Print();

    while (reader.Next())
    {
        Bool_t passTrigger = *passTrigger_;
        Float_t eventWeight = *eventWeight_, PUWeight = *PUWeight_;
        UShort_t nPV = *nPV_, nJets = *nJets_;
        vector<TLorentzVector> elecP4;
        vector<Short_t> elecQ;
        vector<Float_t> elecIso, elecTrkIso, elecD0, elecDz;
        vector<Float_t> elecSF, elecRecoEff, elecTrigEffData, elecTrigEffMC;
        vector<Bool_t> elecIsStd, elecTrigger;
        UInt_t nElecs = *nElecs_, nStdElecs = *nStdElecs_;

        if (passTrigger && nStdElecs == 2) // && nElecs == 2)
        {
            for (const TLorentzVector& elecP4__: elecP4_)
                elecP4.push_back(elecP4__);

            for (unsigned i = 0; i < nElecs; i++)
            {
                elecQ.push_back((*elecQ_)[i]);
                elecIso.push_back((*elecIso_)[i]);
                elecTrkIso.push_back((*elecTrkIso_)[i]);
                elecSF.push_back((*elecSF_)[i]);
                elecD0.push_back((*elecD0_)[i]);
                elecDz.push_back((*elecDz_)[i]);
                elecRecoEff.push_back((*elecRecoEff_)[i]);
                elecTrigEffData.push_back((*elecTriggerEffData_)[i]);
                elecTrigEffMC.push_back((*elecTriggerEffMC_)[i]);
                elecIsStd.push_back((*elecIsStd_)[i]);
                elecTrigger.push_back((*elecTrigger_)[i]);
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
                        h_qt->Fill((elecP4[0] + elecP4[j]).Pt(), eventWeight);
                        h_pt1->Fill(elecP4[0].Pt(), eventWeight);
                        h_eta1->Fill(elecP4[0].Eta(), eventWeight);
                        h_pt2->Fill(elecP4[j].Pt(), eventWeight);
                        h_eta2->Fill(elecP4[j].Eta(), eventWeight);
                        h_npv->Fill(nPV, eventWeight);
                        h_njets->Fill(nJets, eventWeight);

                        // Fill bad elec histograms
                        for (unsigned k = 0; k < nElecs; k++)
                        {
                            if (!elecIsStd[k])
                            {
                                x_pt->Fill(elecP4[k].Pt(), eventWeight);
                                x_eta->Fill(elecP4[k].Eta(), eventWeight);
                                x_q->Fill(elecQ[k], eventWeight);
                                x_iso->Fill(elecIso[k] / elecP4[k].Pt(), eventWeight);
                                x_tiso->Fill(elecTrkIso[k] / elecP4[k].Pt(), eventWeight);
                                x_d0->Fill(elecD0[k], eventWeight);
                                x_dz->Fill(elecDz[k], eventWeight);
                                x_trig->Fill(elecTrigger[k], eventWeight);
                            }
                        }
                        h_n_bad->Fill(nElecs - nStdElecs, eventWeight);

                        break;
                    }
                }
            }
        }
    }
    file->Close();  delete file;

    return h_tot;
}
