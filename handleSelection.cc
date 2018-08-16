#include <vector>
#include <iostream>
#include <utility>
#include <tuple>
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

using namespace std;




//--- HELPERS ---//

bool SortDecPt(const tuple<TLorentzVector, Short_t, Float_t> &i_,
                const tuple<TLorentzVector, Short_t, Float_t> &j_);
double GetBinContentPtEta(const TH2 *hist, const TLorentzVector &p4);
int GetXbin(const TH2 *hist, const double xval);
int GetYbin(const TH2 *hist, const double yval);




void handleSelection(const TString selection, const TString suffix, const TString id)
{


    ///////////////////////
    //      OPTIONS      //
    ///////////////////////


    //--- SYSTEMATICS ---//

    bool systOn = kFALSE, smearID = kFALSE, smearPtMC = kFALSE, smearPtData = kFALSE;
    TString SFType  = "Total";

    //  smearID    = kTRUE;     systOn = kTRUE;     SFType  = "ID/Reco";
    //  smearPtMC  = kTRUE;     systOn = kTRUE;     SFType  = "Energy";
    //  smearPtData= kTRUE;     systOn = kTRUE;     SFType  = "Energy";


    // Lepton Pt smear
    Float_t PT_UNC = 0.002;



    //--- PARSING ---//

    // Selection
    Bool_t sel2l = kFALSE, sel4l = kFALSE;
    Bool_t muPair1, muPair2;    // TRUE for muon, FALSE for electron
    unsigned binIdx = 0;
    if (selection == "mumu" || selection == "2mu" || selection == "2m")
    {
        sel2l = kTRUE;      muPair1 = kTRUE;    muPair2 = muPair1;      binIdx = 3;
    }
    else if (selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;      muPair1 = kFALSE;   muPair2 = muPair1;      binIdx = 4;
    }
    else if (selection == "4mu" || selection == "4m")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;    muPair2 = kTRUE;        binIdx = 6;
    }
    else if (selection == "2mu2e" || selection == "2m2e")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;    muPair2 = kFALSE;       binIdx = 7;
    }
    else if (selection == "2e2mu" || selection == "2e2m")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;   muPair2 = kTRUE;        binIdx = 8;
    }
    else if (selection == "4e")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;   muPair2 = kFALSE;       binIdx = 9;
    }


    // Dataset
    int year = 2016;
    bool trig1l = kTRUE, trig2l = kFALSE;

    TString yearStr;        yearStr.Form("%i", year);
    bool isData = suffix.Contains(yearStr);

    TString dir, suffIncDY;
    if (trig1l)
    {
        dir = yearStr + "_12a";
        suffIncDY = "dy";
    }
    else if (trig2l)
    {
        dir = yearStr + "_12d";
        suffIncDY = "zjets";
    }
    bool isIncDY = suffix.Contains(suffIncDY);
    

    // Branches
    TString genWeightStr = trig1l ? "eventWeight" : "genWeight";
    TString FiredLeg1Str = trig1l ? "Triggered" : "FiredLeg1";
    TString FiredLeg2Str = trig1l ? "Triggered" : "FiredLeg2";
    TString EnergySFStr = trig1l ? "SF" : "EnergySF";
    TString muonIDSFStr = trig1l? "muonHZZIDWeight" : "muonHZZIDSF";
    TString elecIDSFStr = trig1l? "electronHZZRecoWeight" : "electronHZZIDSF";
    TString TrigEffLeg1DataStr = trig1l ? "TriggerEffData" : "TrigEffLeg1Data";
    TString TrigEffLeg1MCStr = trig1l ? "TriggerEffMC" : "TrigEffLeg1MC";
    TString TrigEffLeg2DataStr = trig1l ? "TriggerEffData" : "TrigEffLeg2Data";
    TString TrigEffLeg2MCStr = trig1l ? "TriggerEffMC" : "TrigEffLeg2MC";


    // Cuts & consts
    Float_t M_MIN = 80,     M_MAX = 100;
    Float_t PT1_MIN = trig1l ? 25 : 20,     PT2_MIN = 10;
    Float_t MU_PT_MIN = 5,                  ELE_PT_MIN = 7;
    Float_t MLL_MIN = 4,    MZ1_MIN = 12,   MZ_MIN = 4,     MZ_MAX = 120;
    Float_t DR_MIN = 0.02;
    Float_t Z_MASS = 91.2,  MU_MASS = 0.105658369,  ELE_MASS = 0.000511;



    //--- OUTPUT ---//

    // File
    TString output  = selection + "_" + suffix + "_" + id + ".root";
    TFile *outFile  = new TFile(output, "RECREATE");
    TTree *tree     = new TTree("tree_" + suffix, suffix);


    // Initialize branch variables
    UInt_t          evtNum;
    UShort_t        nPV;
    Float_t         met, weight;
    TLorentzVector  z1p4, z2p4, zzp4, l1p4, l2p4, l3p4, l4p4;
    Short_t         z1pdg, z2pdg, l1pdg, l2pdg, l3pdg, l4pdg;
    Float_t         l1iso, l2iso, l3iso, l4iso;

    TLorentzVector  z1l1p4, z1l2p4, z2l1p4, z2l2p4, tlp4;
    Short_t         z1l1pdg, z1l2pdg, z2l1pdg, z2l2pdg;
    Float_t         z1l1iso, z1l2iso, z2l1iso, z2l2iso, angle;


    // Branches
    tree->Branch("evtNum", &evtNum);
    tree->Branch("nPV", &nPV);      tree->Branch("met", &met);      tree->Branch("weight", &weight);

    if (sel4l) {tree->Branch("zzp4", &zzp4);}
    tree->Branch("z1p4", &z1p4);    if (sel4l) {tree->Branch("z2p4", &z2p4);}

    tree->Branch("l1p4", &l1p4);    tree->Branch("l1pdg", &l1pdg);  tree->Branch("l1iso", &l1iso);
    tree->Branch("l2p4", &l2p4);    tree->Branch("l2pdg", &l2pdg);  tree->Branch("l2iso", &l2iso);
    if (sel4l)
    {
        tree->Branch("l3p4", &l3p4);  tree->Branch("l3pdg", &l3pdg);  tree->Branch("l3iso", &l3iso);
        tree->Branch("l4p4", &l4p4);  tree->Branch("l4pdg", &l4pdg);  tree->Branch("l4iso", &l4iso);

        tree->Branch("z1l1p4", &z1l1p4);   // tree->Branch("z1l1pdg", &z1l1pdg);  tree->Branch*"z1l1iso
        tree->Branch("z1l2p4", &z1l2p4);
        tree->Branch("z2l1p4", &z2l1p4);   // tree->Branch("z2l2p4", &z2l2p4);
        tree->Branch("z2l2p4", &z2l2p4);
        tree->Branch("tlp4", &tlp4);
        tree->Branch("angle", &angle);
    }




    //////////////////////
    //    OPEN TFILE    //
    //////////////////////


    // Path to directory
    TString path = "root://cmsxrootd.fnal.gov//store/user/jrainbol/Trees/";


    // Input ROOT file
    TString input = path + dir + "/" + suffix + "/" + suffix + "_" + id + ".root";
    TFile *file = TFile::Open(input);
    TTreeReader reader("tree_" + suffix, file);



    //--- HISTOGRAMS ---//
    
    TH1D *hTotalEvents, *hAcceptedEvents;
    file->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);        hTotalEvents->Sumw2(kTRUE);
    file->GetObject("AcceptedEvents_" + suffix, hAcceptedEvents);
    hAcceptedEvents->SetDirectory(outFile);     //  hAcceptedEvents->Sumw2(kTRUE);


    // Fill hAcceptedEvents with correct info
    Double_t bin1 = hAcceptedEvents->GetBinContent(1); 
    Double_t bin2 = hAcceptedEvents->GetBinContent(binIdx); 

    hAcceptedEvents->Reset();
    hAcceptedEvents->SetBinContent(1, bin1);
    hAcceptedEvents->SetBinContent(2, bin2);




    /////////////////////////
    //    SINGLE FLAVOR    //
    /////////////////////////


    if (muPair1 == muPair2)
    {


        //--- BRANCHES ---//

        TString Lep = muPair1 ? "Muon" : "Electron",   lep = muPair1 ? "muon" : "electron";
        Short_t pdg = muPair1 ? 13 : 11;
        Float_t LEP_MASS = muPair1 ? MU_MASS : ELE_MASS;
        Float_t PT_MIN = muPair1 ? MU_PT_MIN : ELE_PT_MIN;
        TString lepIDSFStr = muPair1 ? muonIDSFStr : elecIDSFStr;
        TString th2Name = muPair1 ? "FINAL" : "EGamma_SF2D";


        // Event
        TTreeReaderValue<UInt_t>    evtNum_(reader, "evtNumber.eventNumber");
        TTreeReaderValue<UShort_t>  nPV_(reader, "nPV");
        TTreeReaderValue<Float_t>   met_(reader, "met");

        TTreeReaderValue<Bool_t>    passMuonTrig_(reader, "evtMuonTriggered");
        TTreeReaderValue<Bool_t>    passElecTrig_(reader, "evtElectronTriggered");
        TTreeReaderValue<UShort_t>  nLeps_(reader, "n"+Lep+"s");
        TTreeReaderValue<UShort_t>  nHZZLeps_(reader, "nHZZ"+Lep+"s");
        TTreeReaderValue<UShort_t>  nHZZAll_(reader, "nHZZLeptons");
        TTreeReaderValue<UShort_t>  nPartons_(reader, "nPartons");

        TTreeReaderValue<Float_t>   genWeight_(reader, genWeightStr);
        TTreeReaderValue<Float_t>   PUWeight_(reader, "PUWeight");


        // Leptons
        TTreeReaderArray<TLorentzVector>    lepP4_(reader, lep+"P4");
        TTreeReaderValue<vector<Short_t>>   lepQ_(reader, lep+"Q");
        TTreeReaderValue<vector<Bool_t>>    lepFiredLeg1_(reader, lep+FiredLeg1Str);
        TTreeReaderValue<vector<Bool_t>>    lepFiredLeg2_(reader, lep+FiredLeg2Str);
        TTreeReaderValue<vector<Bool_t>>    lepIsHZZ_(reader, lep+"IsHZZ");

        TTreeReaderValue<vector<Float_t>>   lepEnergySF_(reader, lep+EnergySFStr);
        TTreeReaderValue<vector<Float_t>>   lepIDSF_(reader, lepIDSFStr);
        TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg1Data_(reader, lep+TrigEffLeg1DataStr);
        TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg1MC_(reader, lep+TrigEffLeg1MCStr);
        TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg2Data_(reader, lep+TrigEffLeg2DataStr);
        TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg2MC_(reader, lep+TrigEffLeg2MCStr);

        TTreeReaderValue<vector<Float_t>>   lepIso_(reader, lep+"CombIso");
/*
        if (trig2l)
        {
            // Event
            TTreeReaderValue<UInt_t>    evtNum_(reader, "evtNumber.eventNumber");
            TTreeReaderValue<UShort_t>  nPV_(reader, "nPV");
            TTreeReaderValue<Float_t>   met_(reader, "met");

            TTreeReaderValue<Bool_t>    passTrigger_(reader, "evt"+Lep+"Triggered");
            TTreeReaderValue<UShort_t>  nLeps_(reader, "n"+Lep+"s");
            TTreeReaderValue<UShort_t>  nHZZLeps_(reader, "nHZZ"+Lep+"s");
            TTreeReaderValue<UShort_t>  nPartons_(reader, "nPartons");

            TTreeReaderValue<Float_t>   genWeight_(reader, "genWeight");
            TTreeReaderValue<Float_t>   PUWeight_(reader, "PUWeight");


            // Leptons
            TTreeReaderArray<TLorentzVector>    lepP4_(reader, lep+"P4");
            TTreeReaderValue<vector<Short_t>>   lepQ_(reader, lep+"Q");
            TTreeReaderValue<vector<Bool_t>>    lepFiredLeg1_(reader, lep+"FiredLeg1");
            TTreeReaderValue<vector<Bool_t>>    lepFiredLeg2_(reader, lep+"FiredLeg2");
            TTreeReaderValue<vector<Bool_t>>    lepIsHZZ_(reader, lep+"IsHZZ");

            TTreeReaderValue<vector<Float_t>>   lepEnergySF_(reader, lep+"EnergySF");
            TTreeReaderValue<vector<Float_t>>   lepIDSF_(reader, lep+"HZZIDSF");
            TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg1Data_(reader, lep+"TrigEffLeg1Data");
            TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg1MC_(reader, lep+"TrigEffLeg1MC");
            TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg2Data_(reader, lep+"TrigEffLeg2Data");
            TTreeReaderValue<vector<Float_t>>   lepTrigEffLeg2MC_(reader, lep+"TrigEffLeg2MC");

            TTreeReaderValue<vector<Float_t>>   lepIso_(reader, lep+"CombIso");
        }
*/


        //--- HISTOGRAMS ---//

        // Fix for single-lepton trigger files
        TFile *file_id = TFile::Open("hzz_"+lep+"_id_sf.root");
        TH2 *h_id, *h_unc;
        file_id->GetObject(th2Name, h_id);
        h_id->SetDirectory(0);
        file_id->Close();
        

        // Lepton ID smear
        // FIXME
        if (systOn && smearID)
        {
            TFile *file_unc = TFile::Open("hzz_"+lep+"_id_smear.root");
            file_unc->GetObject("SMEAR", h_unc);
            h_unc->SetDirectory(0);
            file_unc->Close();
        }




        //////////////////////
        //    EVENT LOOP    //
        //////////////////////


        unsigned event = 0;
        while (reader.Next() && event < 100)
        {
            //  event++;

            Bool_t   passTrigger;
            if (!isData)
            {
                if (muPair1)
                    passTrigger = *passMuonTrig_ && !*passElecTrig_;
                else
                    passTrigger = *passElecTrig_;
            }
            else
                passTrigger = kTRUE;

            UShort_t nLeps = *nLeps_,           nHZZLeps = *nHZZLeps_;
            UShort_t nHZZAll = *nHZZAll_,       nPartons = *nPartons_;
            Float_t  genWeight = *genWeight_,   PUWeight = *PUWeight_;

            vector<TLorentzVector> lepP4;
            vector<Short_t> lepQ;
            vector<Bool_t> lepFiredLeg1, lepFiredLeg2, lepIsHZZ;
            vector<Float_t> lepEnergySF, lepIDSF;
            vector<Float_t> lepTrigEffLeg1Data, lepTrigEffLeg1MC;
            vector<Float_t> lepTrigEffLeg2Data, lepTrigEffLeg2MC;



            //--- PRESELECTION ---//

            if (!passTrigger)
                continue;

            if (sel2l && (nHZZLeps < 2 || nLeps < 2))
                continue;

            if (sel4l && (nHZZLeps < 4 || nLeps < 4))
                continue;



            //--- LEPTONS ---//

            for (const TLorentzVector& lepP4__: lepP4_)
                lepP4.push_back(lepP4__);
            for (unsigned i = 0; i < nLeps; i++)
            {
                lepQ.push_back((*lepQ_)[i]);
                lepFiredLeg1.push_back((*lepFiredLeg1_)[i]);
                lepFiredLeg2.push_back((*lepFiredLeg2_)[i]);
                lepIsHZZ.push_back((*lepIsHZZ_)[i]);

                lepEnergySF.push_back((*lepEnergySF_)[i]);
                lepIDSF.push_back((*lepIDSF_)[i]);
                lepTrigEffLeg1Data.push_back((*lepTrigEffLeg1Data_)[i]);
                lepTrigEffLeg1MC.push_back((*lepTrigEffLeg1MC_)[i]);
                lepTrigEffLeg2Data.push_back((*lepTrigEffLeg2Data_)[i]);
                lepTrigEffLeg2MC.push_back((*lepTrigEffLeg2MC_)[i]);
            }



            //--- CORRECTIONS ---//

            for (unsigned i = 0; i < nLeps; i++)
            {
                // Energy correction
                if (systOn)
                {
                    if (smearPtMC && !isData)
                        lepEnergySF[i] += PT_UNC * lepEnergySF[i];
                    else if (smearPtData && isData)
                        lepEnergySF[i] += PT_UNC * lepEnergySF[i];
                }
                lepP4[i].SetPtEtaPhiM(lepP4[i].Pt() * lepEnergySF[i], lepP4[i].Eta(),
                                        lepP4[i].Phi(), LEP_MASS);


                // Remove leptons with Pt below threshold
                if (lepIsHZZ[i] && lepP4[i].Pt() < PT_MIN)
                {
                    lepIsHZZ[i] = kFALSE;       nHZZLeps--;
                }


                // ID scale factor
                if (!isData && trig1l)
                    lepIDSF[i] = GetBinContentPtEta(h_id, lepP4[i]);
                if (systOn)
                {
                    if (smearID && !isData)
                        lepIDSF[i] += GetBinContentPtEta(h_unc, lepP4[i]);
                }
            }



            //--- SELECTION ---//


            ////////////////////
            //    DILEPTON    //
            ////////////////////


            if (sel2l)
            {
                if (nHZZLeps != 2)
                    continue;

                // Make sure leading lepton is good
                if (!lepIsHZZ[0] || lepP4[0].Pt() < PT1_MIN || !lepFiredLeg1[0])
                    continue;


                // Pair leading lepton
                unsigned x = 0, y = 0;
                bool foundPair = kFALSE;
                for (unsigned j = 1; j < nLeps; j++)
                {
                    if (!lepIsHZZ[j] || lepP4[j].Pt() < PT2_MIN)
                        continue;
                    if (lepQ[x] == lepQ[j])
                        continue;
                    //  if (trig2l && !lepFiredLeg2[j])
                    //      continue;

                    Float_t mll = (lepP4[x] + lepP4[j]).M();
                    if (mll > M_MIN && mll < M_MAX)
                    {
                        y = j;  foundPair = kTRUE;  break;
                    }
                    if (foundPair)
                        break;
                }
                if (!foundPair)
                    continue;



                //--- WEIGHTING ---//

                float triggerWeight = 1, idWeight = 1;

                if (!isData)
                {
                    // Trigger efficiency
                    if (trig1l)
                    {
                        if (lepFiredLeg1[x] && lepFiredLeg1[y])
                        {
                            if (lepTrigEffLeg1MC[x] > 0 || lepTrigEffLeg1MC[y] > 0)
                            {
                                triggerWeight *= 1. - (1. - lepTrigEffLeg1Data[x]) * (1. - lepTrigEffLeg1Data[y]);
                                triggerWeight /= 1. - (1. - lepTrigEffLeg1MC[x]) * (1. - lepTrigEffLeg1MC[y]);
                            }
                        }
                        else if (lepFiredLeg1[x])
                            triggerWeight = lepTrigEffLeg1Data[x] / lepTrigEffLeg1MC[x];
                        else if (lepFiredLeg1[y])
                            triggerWeight = lepTrigEffLeg1Data[y] / lepTrigEffLeg1MC[y];
                    }
                    else if (trig2l)
                    {
                        if ((lepFiredLeg1[x] || lepFiredLeg2[x]) && (lepFiredLeg1[y] || lepFiredLeg2[y]))
                        {
                            float trigEffData, trigEffMC;

                            trigEffData  = lepTrigEffLeg1Data[x] * lepTrigEffLeg2Data[y];
                            trigEffData += lepTrigEffLeg1Data[y] * lepTrigEffLeg2Data[x];
                            trigEffData -= lepTrigEffLeg1Data[x] * lepTrigEffLeg1Data[y];

                            trigEffMC  = lepTrigEffLeg1MC[x] * lepTrigEffLeg2MC[y];
                            trigEffMC += lepTrigEffLeg1MC[y] * lepTrigEffLeg2MC[x];
                            trigEffMC -= lepTrigEffLeg1MC[x] * lepTrigEffLeg1MC[y];
                            if (trigEffMC == 0)
                                trigEffMC = 1;

                            triggerWeight = trigEffData / trigEffMC;
                        }
                    }


                    // ID efficiency
                    idWeight = lepIDSF[x] * lepIDSF[y];
                }
                float eventWeight = genWeight * PUWeight * triggerWeight * idWeight;



                //--- HISTOGRAMS ---//

                hAcceptedEvents->Fill(3, eventWeight);

                if (isIncDY && nPartons != 0)
                    continue;

                hTotalEvents->Fill(6);
                hTotalEvents->Fill(7, eventWeight);



                //--- CONTAINERS ---//

                evtNum = *evtNum_;                          nPV = *nPV_;        
                met = *met_;                                weight = eventWeight;
                l1p4 = lepP4[x];    l2p4 = lepP4[y];        z1p4 = l1p4 + l2p4;
                l1pdg = lepQ[x] * pdg;                      l2pdg = lepQ[y] * pdg;
                l1iso = (*lepIso_)[x] / l1p4.Pt();          l2iso = (*lepIso_)[y] / l2p4.Pt();
            
                tree->Fill();
            }




            ////////////////////
            //    4-LEPTON    //
            ////////////////////


            else if (sel4l)
            {
                // Remove parton events from inclusive DY sample
                if (isIncDY && nPartons != 0)
                    continue;

                if (nHZZLeps < 4)
                    continue;



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

                        Float_t mll = (lepP4[i] + lepP4[j]).M();
                        if (mll > MZ_MIN && mll < MZ_MAX)
                            zCand.push_back(make_pair(i, j));
                    }
                }
                if (zCand.size() < 2)
                    continue;



                //--- ZZ SELECTION ---//

                vector<pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>>> zzCand;
                vector<Float_t> zzCand_m4l, zzCand_mZ1;
                for (unsigned j = 1; j < zCand.size(); j++)
                {
                    for (unsigned i = 0; i < j; i++)
                    {
                        pair<unsigned, unsigned> Zi = zCand[i], Zj = zCand[j];


                        // Make sure pairs do not overlap
                        if (Zi.first == Zj.first || Zi.first == Zj.second
                                || Zi.second == Zj.first || Zi.second == Zj.second)
                            continue;


                        // Choose Z1 closest to nominal Z mass
                        Float_t mll_i = (lepP4[Zi.first] + lepP4[Zi.second]).M();
                        Float_t mll_j = (lepP4[Zj.first] + lepP4[Zj.second]).M();
                        pair<unsigned, unsigned> Z1, Z2;

                        if (fabs(mll_i - Z_MASS) < fabs(mll_j - Z_MASS))
                        {
                            Z1 = Zi;        Z2 = Zj;
                        }
                        else
                        {
                            Z1 = Zj;        Z2 = Zi;
                        }


                        // Pair Pt, Q for easier manipulation
                        pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
                        Z1_p4s = make_pair(lepP4[Z1.first], lepP4[Z1.second]);
                        Z2_p4s = make_pair(lepP4[Z2.first], lepP4[Z2.second]);

                        pair<Short_t, Short_t> Z1_qs, Z2_qs;
                        Z1_qs = make_pair(lepQ[Z1.first], lepQ[Z1.second]);
                        Z2_qs = make_pair(lepQ[Z2.first], lepQ[Z2.second]);

                        TLorentzVector ZZ_p4, Z1_p4, Z2_p4;
                        Z1_p4 = Z1_p4s.first + Z1_p4s.second;
                        Z2_p4 = Z2_p4s.first + Z2_p4s.second;
                        ZZ_p4 = Z1_p4 + Z2_p4;


                        // Mass requirements
                        if (Z1_p4.M() < MZ1_MIN)
                            continue;
                        if (ZZ_p4.M() < M_MIN || ZZ_p4.M() > M_MAX)
                            continue;


                        // "Smart cut"
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

                        Float_t mll_x, mll_y, mll_a, mll_b;
                        mll_x = (lepP4[Zx.first] + lepP4[Zx.second]).M();
                        mll_y = (lepP4[Zy.first] + lepP4[Zy.second]).M();

                        if (fabs(mll_x - Z_MASS) < fabs(mll_y - Z_MASS))
                        {
                            Za = Zx;        Zb = Zy;        mll_a = mll_x;  mll_b = mll_y;
                        }
                        else
                        {
                            Za = Zy;        Zb = Zx;        mll_a = mll_y;  mll_b = mll_x;
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
                                    foundGhost = kTRUE;     break;
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
                                    foundQCD = kTRUE;       break;
                                }
                            }
                            if (foundQCD)
                                break;
                        }
                        if (foundQCD)
                            continue;


                        // Assemble qualifying ZZ candidate
                        zzCand.push_back(make_pair(Z1, Z2));
                        zzCand_mZ1.push_back(Z1_p4.M());
                        zzCand_m4l.push_back(ZZ_p4.M());
                    }
                }
                if (zzCand.size() < 1)
                    continue;



                //--- BEST ZZ ---//

                unsigned zz = 0;
                Float_t massDiff = 1000.;
                for (unsigned m = 0; m < zzCand.size(); m++)
                {
                    Float_t massDiff_ = fabs(zzCand_mZ1[m] - Z_MASS); //= fabs(zzCand_m4l[m] - Z_MASS);
                    if (massDiff_ < massDiff)
                    {
                        massDiff = massDiff_;       zz = m;
                    }
                }

                // Assemble indices of selected pairs
                pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>> ZZ = zzCand[zz];
                pair<unsigned, unsigned> Z1 = ZZ.first, Z2 = ZZ.second;
                unsigned x = Z1.first, y = Z1.second, xx = Z2.first, yy = Z2.second;


                // Sort lepton PDG ID and iso by Pt
                vector<tuple<TLorentzVector, Short_t, Float_t>> lepTuple;
                lepTuple.push_back(make_tuple(lepP4[x], lepQ[x]*pdg, (*lepIso_)[x]/lepP4[x].Pt()));
                lepTuple.push_back(make_tuple(lepP4[y], lepQ[y]*pdg, (*lepIso_)[y]/lepP4[y].Pt()));
                lepTuple.push_back(make_tuple(lepP4[xx], lepQ[xx]*pdg, (*lepIso_)[xx]/lepP4[xx].Pt()));
                lepTuple.push_back(make_tuple(lepP4[yy], lepQ[yy]*pdg, (*lepIso_)[yy]/lepP4[yy].Pt()));
                sort(lepTuple.begin(), lepTuple.end(), SortDecPt);



                //--- WEIGHTING ---//

                float triggerWeight = 1, idWeight = 1;

                if (!isData)
                {
                    // Trigger efficiency
                    // Uses only leading leptons, FIXME
                    if (trig1l)
                    {
                        if (lepFiredLeg1[x] && lepFiredLeg1[y])
                        {
                            if (lepTrigEffLeg1MC[x] > 0 || lepTrigEffLeg1MC[y] > 0)
                            {
                                triggerWeight *= 1. - (1. - lepTrigEffLeg1Data[x]) * (1. - lepTrigEffLeg1Data[y]);
                                triggerWeight /= 1. - (1. - lepTrigEffLeg1MC[x]) * (1. - lepTrigEffLeg1MC[y]);
                            }
                        }
                        else if (lepFiredLeg1[x])
                            triggerWeight = lepTrigEffLeg1Data[x] / lepTrigEffLeg1MC[x];
                        else if (lepFiredLeg1[y])
                            triggerWeight = lepTrigEffLeg1Data[y] / lepTrigEffLeg1MC[y];
                    }
                    else if (trig2l)
                    {
                        if ((lepFiredLeg1[x] || lepFiredLeg2[x]) && (lepFiredLeg1[y] || lepFiredLeg2[y]))
                        {
                            float trigEffData, trigEffMC;

                            trigEffData  = lepTrigEffLeg1Data[x] * lepTrigEffLeg2Data[y];
                            trigEffData += lepTrigEffLeg1Data[y] * lepTrigEffLeg2Data[x];
                            trigEffData -= lepTrigEffLeg1Data[x] * lepTrigEffLeg1Data[y];

                            trigEffMC  = lepTrigEffLeg1MC[x] * lepTrigEffLeg2MC[y];
                            trigEffMC += lepTrigEffLeg1MC[y] * lepTrigEffLeg2MC[x];
                            trigEffMC -= lepTrigEffLeg1MC[x] * lepTrigEffLeg1MC[y];
                            if (trigEffMC == 0)
                                trigEffMC = 1;

                            triggerWeight = trigEffData / trigEffMC;
                        }
                    }


                    // ID efficiency
                    idWeight *= lepIDSF[x] * lepIDSF[y] * lepIDSF[xx] * lepIDSF[yy];
                }
                float eventWeight = genWeight * PUWeight * triggerWeight * idWeight;



                //--- HISTOGRAMS ---//

                hAcceptedEvents->Fill(3, eventWeight);
                hTotalEvents->Fill(6);
                hTotalEvents->Fill(7, eventWeight);


                //--- CONTAINERS ---//

                evtNum = *evtNum_;                          nPV = *nPV_;        
                met = *met_;                                weight = eventWeight;
                z1p4 = lepP4[x] + lepP4[y];                 z2p4 = lepP4[xx] + lepP4[yy]; 
                zzp4 = z1p4 + z2p4;
                tie(l1p4, l1pdg, l1iso) = lepTuple[0];      tie(l2p4, l2pdg, l2iso) = lepTuple[1];
                tie(l3p4, l3pdg, l3iso) = lepTuple[2];      tie(l4p4, l4pdg, l4iso) = lepTuple[3];

                z1l1p4 = lepP4[x];  z1l2p4 = lepP4[y];      z2l1p4 = lepP4[xx]; z2l2p4 = lepP4[yy];
                tlp4 = l2p4 + l3p4 + l4p4;
                angle = z1l2p4.Angle(z2p4.Vect());

                tree->Fill();
            } 
        }
    }




    ///////////////////////
    //    TWO FLAVORS    //
    ///////////////////////


    else
    {


        //--- BRANCHES ---//

        TString Lep1 = muPair1 ? "Muon" : "Electron",   Lep2 = muPair2 ? "Muon" : "Electron";
        TString lep1 = muPair1 ? "muon" : "electron",   lep2 = muPair2 ? "muon" : "electron";
        Short_t pdg1 = muPair1 ? 13 : 11,               pdg2 = muPair2 ? 13 : 11;
        Float_t LEP1_MASS = muPair1 ? MU_MASS : ELE_MASS;
        Float_t LEP2_MASS = muPair2 ? MU_MASS : ELE_MASS;
        Float_t LEP1_PT_MIN = muPair1 ? MU_PT_MIN : ELE_PT_MIN;
        Float_t LEP2_PT_MIN = muPair2 ? MU_PT_MIN : ELE_PT_MIN;
        TString lep1IDSFStr = muPair1 ? muonIDSFStr : elecIDSFStr;
        TString lep2IDSFStr = muPair2 ? muonIDSFStr : elecIDSFStr;
        TString th2lep1Name = muPair1 ? "FINAL" : "EGamma_SF2D";
        TString th2lep2Name = muPair2 ? "FINAL" : "EGamma_SF2D";


        // Event
        TTreeReaderValue<UInt_t>    evtNum_(reader, "evtNumber.eventNumber");
        TTreeReaderValue<UShort_t>  nPV_(reader, "nPV");
        TTreeReaderValue<Float_t>   met_(reader, "met");

        TTreeReaderValue<Bool_t>    passMuonTrig_(reader, "evtMuonTriggered");
        TTreeReaderValue<Bool_t>    passElecTrig_(reader, "evtElectronTriggered");
        TTreeReaderValue<UShort_t>  nLep1s_(reader, "n"+Lep1+"s");
        TTreeReaderValue<UShort_t>  nLep2s_(reader, "n"+Lep2+"s");
        TTreeReaderValue<UShort_t>  nHZZLep1s_(reader, "nHZZ"+Lep1+"s");
        TTreeReaderValue<UShort_t>  nHZZLep2s_(reader, "nHZZ"+Lep2+"s");
        TTreeReaderValue<UShort_t>  nHZZAll_(reader, "nHZZLeptons");
        TTreeReaderValue<UShort_t>  nPartons_(reader, "nPartons");

        TTreeReaderValue<Float_t>   genWeight_(reader, genWeightStr);
        TTreeReaderValue<Float_t>   PUWeight_(reader, "PUWeight");


        // Leading flavor leptons
        TTreeReaderArray<TLorentzVector>    lep1P4_(reader, lep1+"P4");
        TTreeReaderValue<vector<Short_t>>   lep1Q_(reader, lep1+"Q");
        TTreeReaderValue<vector<Bool_t>>    lep1FiredLeg1_(reader, lep1+FiredLeg1Str);
        TTreeReaderValue<vector<Bool_t>>    lep1FiredLeg2_(reader, lep1+FiredLeg2Str);
        TTreeReaderValue<vector<Bool_t>>    lep1IsHZZ_(reader, lep1+"IsHZZ");

        TTreeReaderValue<vector<Float_t>>   lep1EnergySF_(reader, lep1+EnergySFStr);
        TTreeReaderValue<vector<Float_t>>   lep1IDSF_(reader, lep1IDSFStr);
        TTreeReaderValue<vector<Float_t>>   lep1TrigEffLeg1Data_(reader, lep1+TrigEffLeg1DataStr);
        TTreeReaderValue<vector<Float_t>>   lep1TrigEffLeg1MC_(reader, lep1+TrigEffLeg1MCStr);
        TTreeReaderValue<vector<Float_t>>   lep1TrigEffLeg2Data_(reader, lep1+TrigEffLeg2DataStr);
        TTreeReaderValue<vector<Float_t>>   lep1TrigEffLeg2MC_(reader, lep1+TrigEffLeg2MCStr);

        TTreeReaderValue<vector<Float_t>>   lep1Iso_(reader, lep1+"CombIso");


        // Subleading flavor leptons
        TTreeReaderArray<TLorentzVector>    lep2P4_(reader, lep2+"P4");
        TTreeReaderValue<vector<Short_t>>   lep2Q_(reader, lep2+"Q");
        TTreeReaderValue<vector<Bool_t>>    lep2FiredLeg1_(reader, lep2+FiredLeg1Str);
        TTreeReaderValue<vector<Bool_t>>    lep2FiredLeg2_(reader, lep2+FiredLeg2Str);
        TTreeReaderValue<vector<Bool_t>>    lep2IsHZZ_(reader, lep2+"IsHZZ");

        TTreeReaderValue<vector<Float_t>>   lep2EnergySF_(reader, lep2+EnergySFStr);
        TTreeReaderValue<vector<Float_t>>   lep2IDSF_(reader, lep2IDSFStr);

        TTreeReaderValue<vector<Float_t>>   lep2Iso_(reader, lep2+"CombIso");



        //--- HISTOGRAMS ---//

        // Fix for single-lepton trigger files
        TFile *file_id1 = TFile::Open("hzz_"+lep1+"_id_sf.root");
        TH2 *h_id1, *h_unc1;
        file_id1->GetObject(th2lep1Name, h_id1);
        h_id1->SetDirectory(0);
        file_id1->Close();

        TFile *file_id2 = TFile::Open("hzz_"+lep2+"_id_sf.root");
        TH2 *h_id2, *h_unc2;
        file_id2->GetObject(th2lep2Name, h_id2);
        h_id2->SetDirectory(0);
        file_id2->Close();
        

        // Lepton ID smear
        // FIXME
        if (systOn && smearID)
        {
            TFile *file_unc1 = TFile::Open("hzz_"+lep1+"_id_smear.root");
            file_unc1->GetObject("SMEAR", h_unc1);
            h_unc1->SetDirectory(0);
            file_unc1->Close();

            TFile *file_unc2 = TFile::Open("hzz_"+lep2+"_id_smear.root");
            file_unc2->GetObject("SMEAR", h_unc2);
            h_unc2->SetDirectory(0);
            file_unc2->Close();
        }




        //////////////////////
        //    EVENT LOOP    //
        //////////////////////


        unsigned event = 0;
        while (reader.Next() && event < 100)
        {
            //  event++;

            Bool_t   passTrigger;
            if (!isData)
            {
                if (muPair1)
                    passTrigger = *passMuonTrig_ && !*passElecTrig_;
                else
                    passTrigger = *passElecTrig_;
            }
            else
                passTrigger = kTRUE;

            UShort_t nLep1s = *nLep1s_,         nHZZLep1s = *nHZZLep1s_;
            UShort_t nLep2s = *nLep2s_,         nHZZLep2s = *nHZZLep2s_;
            UShort_t nHZZAll = *nHZZAll_,       nPartons = *nPartons_;
            Float_t  genWeight = *genWeight_,   PUWeight = *PUWeight_;

            vector<TLorentzVector> lep1P4, lep2P4;
            vector<Short_t> lep1Q, lep2Q;
            vector<Bool_t> lep1FiredLeg1, lep1FiredLeg2, lep1IsHZZ;
            vector<Bool_t> lep2FiredLeg1, lep2FiredLeg2, lep2IsHZZ;
            vector<Float_t> lep1EnergySF, lep1IDSF, lep2EnergySF, lep2IDSF;
            vector<Float_t> lep1TrigEffLeg1Data, lep1TrigEffLeg1MC;
            vector<Float_t> lep1TrigEffLeg2Data, lep1TrigEffLeg2MC;
            vector<Float_t> lep1Iso, lep2Iso;



            //--- PRESELECTION ---//

            // Remove parton events from inclusive DY sample
            if (isIncDY && nPartons != 0)
                continue;

            if (!passTrigger)
                continue;

            if (nHZZLep1s < 2 || nHZZLep2s < 2)
                continue;



            //--- LEPTONS ---//

            // Leading flavor
            for (const TLorentzVector& lep1P4__: lep1P4_)
                lep1P4.push_back(lep1P4__);
            for (unsigned i = 0; i < nLep1s; i++)
            {
                lep1Q.push_back((*lep1Q_)[i]);
                lep1FiredLeg1.push_back((*lep1FiredLeg1_)[i]);
                lep1FiredLeg2.push_back((*lep1FiredLeg2_)[i]);
                lep1IsHZZ.push_back((*lep1IsHZZ_)[i]);

                lep1EnergySF.push_back((*lep1EnergySF_)[i]);
                lep1IDSF.push_back((*lep1IDSF_)[i]);
                lep1TrigEffLeg1Data.push_back((*lep1TrigEffLeg1Data_)[i]);
                lep1TrigEffLeg1MC.push_back((*lep1TrigEffLeg1MC_)[i]);
                lep1TrigEffLeg2Data.push_back((*lep1TrigEffLeg2Data_)[i]);
                lep1TrigEffLeg2MC.push_back((*lep1TrigEffLeg2MC_)[i]);

                lep1Iso.push_back((*lep1Iso_)[i]);
            }

            // Subleading flavor
            for (const TLorentzVector& lep2P4__: lep2P4_)
                lep2P4.push_back(lep2P4__);
            for (unsigned i = 0; i < nLep2s; i++)
            {
                lep2Q.push_back((*lep2Q_)[i]);
                lep2FiredLeg1.push_back((*lep2FiredLeg1_)[i]);
                lep2FiredLeg2.push_back((*lep2FiredLeg2_)[i]);
                lep2IsHZZ.push_back((*lep2IsHZZ_)[i]);

                lep2EnergySF.push_back((*lep2EnergySF_)[i]);
                lep2IDSF.push_back((*lep2IDSF_)[i]);

                lep2Iso.push_back((*lep2Iso_)[i]);
            }



            //--- CORRECTIONS ---//

            // Leading flavor
            for (unsigned i = 0; i < nLep1s; i++)
            {
                // Energy correction
                if (systOn)
                {
                    if (smearPtMC && !isData)
                        lep1EnergySF[i] += PT_UNC * lep1EnergySF[i];
                    else if (smearPtData && isData)
                        lep1EnergySF[i] += PT_UNC * lep1EnergySF[i];
                }
                lep1P4[i].SetPtEtaPhiM(lep1P4[i].Pt() * lep1EnergySF[i], lep1P4[i].Eta(),
                                        lep1P4[i].Phi(), LEP1_MASS);


                // Remove leptons with Pt below threshold
                if (lep1IsHZZ[i] && lep1P4[i].Pt() < LEP1_PT_MIN)
                {
                    lep1IsHZZ[i] = kFALSE;      nHZZLep1s--;
                }


                // ID scale factor
                if (!isData && trig1l)
                    lep1IDSF[i] = GetBinContentPtEta(h_id1, lep1P4[i]);
                if (systOn)
                {
                    if (smearID && !isData)
                        lep1IDSF[i] += GetBinContentPtEta(h_unc1, lep1P4[i]);
                }
            }

            // Subleading flavor
            for (unsigned i = 0; i < nLep2s; i++)
            {
                // Energy correction
                if (systOn)
                {
                    if (smearPtMC && !isData)
                        lep2EnergySF[i] += PT_UNC * lep2EnergySF[i];
                    else if (smearPtData && isData)
                        lep2EnergySF[i] += PT_UNC * lep2EnergySF[i];
                }
                lep2P4[i].SetPtEtaPhiM(lep2P4[i].Pt() * lep2EnergySF[i], lep2P4[i].Eta(),
                                        lep2P4[i].Phi(), LEP2_MASS);


                // Remove leptons with Pt below threshold
                if (lep2IsHZZ[i] && lep2P4[i].Pt() < LEP2_PT_MIN)
                {
                    lep2IsHZZ[i] = kFALSE;      nHZZLep2s--;
                }


                // ID scale factor
                if (!isData && trig1l)
                    lep2IDSF[i] = GetBinContentPtEta(h_id2, lep2P4[i]);
                if (systOn)
                {
                    if (smearID && !isData)
                        lep2IDSF[i] += GetBinContentPtEta(h_unc2, lep2P4[i]);
                }
            }




            /////////////////////
            //    SELECTION    //
            /////////////////////


            //--- Z SELECTION ---//

            vector<pair<unsigned, unsigned>> z1Cand, z2Cand;

            // Z1 candidates
            for (unsigned j = 1; j < nLep1s; j++)
            {
                if (!lep1IsHZZ[j])
                    continue;

                for (unsigned i = 0; i < j; i++)
                {
                    if (!lep1IsHZZ[i])
                        continue;
                    if (lep1Q[i] == lep1Q[j])
                        continue;
                    //  if (trig2l && !lep1FiredLeg2[j])
                    //      continue;

                    Float_t mll = (lep1P4[i] + lep1P4[j]).M();
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

                    Float_t mll = (lep2P4[i] + lep2P4[j]).M();
                    if (mll > MZ_MIN && mll < MZ_MAX)
                        z2Cand.push_back(make_pair(i, j));
                }
            }
            if (z2Cand.size() < 1)
                continue;



            //--- ZZ SELECTION ---//

            vector<pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>>> zzCand;
            vector<Float_t> zzCand_m4l, zzCand_mZ1;
            for (unsigned i = 0; i < z1Cand.size(); i++)
            {
                for (unsigned j = 0; j < z2Cand.size(); j++)
                {
                    // Pair Pt, Q for easier manipulation
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
                    if (ZZ_p4.M() < M_MIN || ZZ_p4.M() > M_MAX)
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
                                foundGhost = kTRUE;     break;
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



            //--- BEST ZZ ---//

            unsigned zz = 0;
            Float_t massDiff = 1000.;
            for (unsigned m = 0; m < zzCand.size(); m++)
            {
                Float_t massDiff_ = fabs(zzCand_mZ1[m] - Z_MASS);
                if (massDiff_ < massDiff)
                {
                    massDiff = massDiff_;       zz = m;
                }
            }


            // Assemble indices of selected pairs
            pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>> ZZ = zzCand[zz];
            pair<unsigned, unsigned> Z1 = ZZ.first, Z2 = ZZ.second;
            unsigned x = Z1.first, y = Z1.second, xx = Z2.first, yy = Z2.second;


            // Sort lepton PDG ID and iso by Pt
            vector<tuple<TLorentzVector, Short_t, Float_t>> lepTuple;
            lepTuple.push_back(make_tuple(lep1P4[x], lep1Q[x]*pdg1, lep1Iso[x]/lep1P4[x].Pt()));
            lepTuple.push_back(make_tuple(lep1P4[y], lep1Q[y]*pdg1, lep1Iso[y]/lep1P4[y].Pt()));
            lepTuple.push_back(make_tuple(lep2P4[xx], lep2Q[xx]*pdg2, lep2Iso[xx]/lep2P4[xx].Pt()));
            lepTuple.push_back(make_tuple(lep2P4[yy], lep2Q[yy]*pdg2, lep2Iso[yy]/lep2P4[yy].Pt()));
            sort(lepTuple.begin(), lepTuple.end(), SortDecPt);



            //--- WEIGHTING ---//

            float triggerWeight = 1, idWeight = 1;

            if (!isData)
            {
                // Trigger efficiency
                // Uses only leading leptons, FIXME
                if (trig1l)
                {
                    if (lep1FiredLeg1[x] && lep1FiredLeg1[y])
                    {
                        if (lep1TrigEffLeg1MC[x] > 0 || lep1TrigEffLeg1MC[y] > 0)
                        {
                            triggerWeight *= 1. - (1. - lep1TrigEffLeg1Data[x]) * (1. - lep1TrigEffLeg1Data[y]);
                            triggerWeight /= 1. - (1. - lep1TrigEffLeg1MC[x]) * (1. - lep1TrigEffLeg1MC[y]);
                        }
                    }
                    else if (lep1FiredLeg1[x])
                        triggerWeight = lep1TrigEffLeg1Data[x] / lep1TrigEffLeg1MC[x];
                    else if (lep1FiredLeg1[y])
                        triggerWeight = lep1TrigEffLeg1Data[y] / lep1TrigEffLeg1MC[y];
                }
                else if (trig2l)
                {
                    if ((lep1FiredLeg1[x] || lep1FiredLeg2[x]) && (lep1FiredLeg1[y] || lep1FiredLeg2[y]))
                    {
                        float trigEffData, trigEffMC;

                        trigEffData  = lep1TrigEffLeg1Data[x] * lep1TrigEffLeg2Data[y];
                        trigEffData += lep1TrigEffLeg1Data[y] * lep1TrigEffLeg2Data[x];
                        trigEffData -= lep1TrigEffLeg1Data[x] * lep1TrigEffLeg1Data[y];

                        trigEffMC  = lep1TrigEffLeg1MC[x] * lep1TrigEffLeg2MC[y];
                        trigEffMC += lep1TrigEffLeg1MC[y] * lep1TrigEffLeg2MC[x];
                        trigEffMC -= lep1TrigEffLeg1MC[x] * lep1TrigEffLeg1MC[y];
                        if (trigEffMC == 0)
                            trigEffMC = 1;

                        triggerWeight = trigEffData / trigEffMC;
                    }
                }


                // ID efficiency
                idWeight *= lep1IDSF[x] * lep1IDSF[y] * lep2IDSF[xx] * lep2IDSF[yy];
            }
            float eventWeight = genWeight * PUWeight * triggerWeight * idWeight;



            //--- HISTOGRAMS ---//

            hAcceptedEvents->Fill(3, eventWeight);
            hTotalEvents->Fill(6);
            hTotalEvents->Fill(7, eventWeight);


            //--- CONTAINERS ---//

            evtNum = *evtNum_;                          nPV = *nPV_;        
            met = *met_;                                weight = eventWeight;
            z1p4 = lep1P4[x] + lep1P4[y];               z2p4 = lep2P4[xx] + lep2P4[yy]; 
            zzp4 = z1p4 + z2p4;
            tie(l1p4, l1pdg, l1iso) = lepTuple[0];      tie(l2p4, l2pdg, l2iso) = lepTuple[1];
            tie(l3p4, l3pdg, l3iso) = lepTuple[2];      tie(l4p4, l4pdg, l4iso) = lepTuple[3];

            z1l1p4 = lep1P4[x]; z1l2p4 = lep1P4[y];     z2l1p4 = lep2P4[xx];    z2l2p4 = lep2P4[yy];
            tlp4 = l2p4 + l3p4 + l4p4;
            angle = z1l2p4.Angle(z2p4.Vect());

            tree->Fill();
        }
    }



    //--- CLEANUP ---//

    // Write to file
    outFile->cd();
    tree->Write();
    hTotalEvents->Write();
    hAcceptedEvents->Write();
    outFile->Purge();
    outFile->Close();
    file->Close();



//  delete outFile;
//  delete tree;
//  delete hAcceptedEvents;
}




//////////////////////
// HELPER FUNCTIONS //
//////////////////////


bool SortDecPt(const tuple<TLorentzVector, Short_t, Float_t> &i_,
                const tuple<TLorentzVector, Short_t, Float_t> &j_)
{
    TLorentzVector i = get<0>(i_), j = get<0>(j_);
    return (i.Pt() > j.Pt());
}

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
