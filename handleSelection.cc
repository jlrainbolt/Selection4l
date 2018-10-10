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

bool            SortDecPt(  const tuple<TLorentzVector, Short_t, Float_t> &i_,
                            const tuple<TLorentzVector, Short_t, Float_t> &j_);

bool            SortDecP(   const tuple<TLorentzVector, Short_t, Float_t> &i_,
                            const tuple<TLorentzVector, Short_t, Float_t> &j_);

TLorentzVector  GetBoosted( const TLorentzVector &p4_,  const TVector3 &beta);

TLorentzVector  GetP4Sum(   const vector<tuple<TLorentzVector, Short_t, Float_t>> &leps);


/*
double GetBinContentPtEta(const TH2 *hist, const TLorentzVector &p4);
int GetXbin(const TH2 *hist, const double xval);
int GetYbin(const TH2 *hist, const double yval);
*/



void handleSelection(const TString suffix, const TString id, const TString systematics)
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



    //--- SELECTION ---//

    // Dataset
    int year = 2017;
    TString yearStr = TString::Format("%i", year);
    bool isData = suffix.Contains(yearStr);
    TString dir = yearStr; // + "_old";


    // Cuts & consts
    Float_t M_MIN = 80,         M_MAX = 100;
    Float_t MU_PT1_MIN = 20,    ELE_PT1_MIN = 25;
    Float_t MU_PT2_MIN = 10,    ELE_PT2_MIN = 15;
    Float_t MU_PT_MIN = 5,      ELE_PT_MIN = 7;
    Float_t MLL_MIN = 4,        MZ1_MIN = 12,   MZ_MIN = 4,     MZ_MAX = 120;
    Float_t DR_MIN = 0.02;
    Float_t Z_MASS = 91.2,      MU_MASS = 0.105658369,  ELE_MASS = 0.000511;



    //--- OUTPUT ---//

    // File
    TString output  = "trees_" + suffix + "_" + id + ".root";
    TFile *outFile  = new TFile(output, "RECREATE");


    // Selection
    const unsigned N = 6;
    unsigned                MM = 0, EE = 1, M4 = 2, ME = 3, EM = 4, E4 = 5; // Indices
    TString selection[N] = {"mumu", "ee",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned binIdx[N]   = {3,      4,      6,      7,      8,      9};


    // Trees
    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);


    // Initialize branch variables
    UInt_t          evtNum;
    UShort_t        nPV;
    Float_t         met, weight;
    TLorentzVector  z1p4, z2p4, zzp4, l1p4, l2p4, l3p4, l4p4;
    Short_t         z1pdg, z2pdg, l1pdg, l2pdg, l3pdg, l4pdg;
    Float_t         l1iso, l2iso, l3iso, l4iso;

    Short_t         z1l1pdg, z1l2pdg, z2l1pdg, z2l2pdg;
    Float_t         z1l1iso, z1l2iso, z2l1iso, z2l2iso;

    TLorentzVector  b_z1p4, b_z2p4, b_ttp4;
    TLorentzVector  b_l1p4, b_l2p4, b_l3p4, b_l4p4;
    Short_t         b_l1pdg, b_l2pdg, b_l3pdg, b_l4pdg;

    Float_t         b_theta, b_phi, b_z1alpha, b_z2alpha;
    Float_t         bb_z1theta, bb_z2theta;


    // Branches
    for (unsigned i = 0; i < N; i++)
    {
        // Event info
        tree[i]->Branch("evtNum",   &evtNum);       tree[i]->Branch("weight",   &weight);
        tree[i]->Branch("nPV",      &nPV);          tree[i]->Branch("met",      &met);


        // Pair/group momenta
        if (i > EE)
            tree[i]->Branch("zzp4",     &zzp4);

        tree[i]->Branch("z1p4",     &z1p4);

        if (i > EE)
            tree[i]->Branch("z2p4",     &z2p4);


        // Lepton momenta, id, iso
        tree[i]->Branch("l1p4",     &l1p4);         
        tree[i]->Branch("l1pdg",    &l1pdg);        tree[i]->Branch("l1iso",     &l1iso);

        tree[i]->Branch("l2p4",     &l2p4);
        tree[i]->Branch("l2pdg",    &l2pdg);        tree[i]->Branch("l2iso",     &l2iso);

        if (i > EE)
        {
        tree[i]->Branch("l3p4",     &l3p4);         
        tree[i]->Branch("l3pdg",    &l3pdg);        tree[i]->Branch("l3iso",     &l3iso);

        tree[i]->Branch("l4p4",     &l4p4);
        tree[i]->Branch("l4pdg",    &l4pdg);        tree[i]->Branch("l4iso",     &l4iso);


            // Boosted quantities for differential distributions
            tree[i]->Branch("b_z1p4",   &b_z1p4);       tree[i]->Branch("b_z2p4",   &b_z2p4);
            tree[i]->Branch("b_ttp4",   &b_ttp4);
            tree[i]->Branch("b_l1p4",   &b_l1p4);       tree[i]->Branch("b_l1pdg",  &b_l1pdg);
            tree[i]->Branch("b_l2p4",   &b_l2p4);       tree[i]->Branch("b_l2pdg",  &b_l2pdg);
            tree[i]->Branch("b_l3p4",   &b_l3p4);       tree[i]->Branch("b_l3pdg",  &b_l3pdg);
            tree[i]->Branch("b_l4p4",   &b_l4p4);       tree[i]->Branch("b_l4pdg",  &b_l4pdg);

            tree[i]->Branch("b_theta",  &b_theta);      tree[i]->Branch("b_phi",    &b_phi);
            tree[i]->Branch("b_z1alpha", &b_z1alpha);   tree[i]->Branch("b_z2alpha", &b_z2alpha);
            tree[i]->Branch("bb_z1theta", &bb_z1theta); tree[i]->Branch("bb_z2theta", &bb_z2theta);
        }
    }




    //////////////////////
    //    INPUT FILE    //
    //////////////////////


    // Path to directory
//  TString path = "root://cmsxrootd.fnal.gov//store/user/jrainbol/Trees/";
    TString path = "root://cmseos.fnal.gov//store/user/jrainbol/Trees/";


    // Input ROOT file
    TString input = path + dir + "/" + suffix + "/" + suffix + "_" + id + ".root";
    cout << input << endl;
    TFile *file = TFile::Open(input);
    TTreeReader reader("tree_" + suffix, file);



    //--- HISTOGRAMS ---//
    
    TH1D *hTotalEvents, *hPhaseSpaceEvents, *hFiducialEvents, *hSelectedEvents;

    file->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);
    hTotalEvents->Sumw2(kTRUE);

    file->GetObject("PhaseSpaceEvents_" + suffix, hPhaseSpaceEvents);
    if (!hPhaseSpaceEvents) 
        hPhaseSpaceEvents = new TH1D("PhaseSpaceEvents_" + suffix, "PhaseSpaceEvents", 10, 0.5, 10.5);
    hPhaseSpaceEvents->SetDirectory(outFile);
    hPhaseSpaceEvents->Sumw2(kTRUE);

    file->GetObject("FiducialEvents_" + suffix, hFiducialEvents);
    if (!hFiducialEvents) 
        hFiducialEvents = new TH1D("FiducialEvents_" + suffix, "FiducialEvents", 10, 0.5, 10.5);
    hFiducialEvents->SetDirectory(outFile);
    hFiducialEvents->Sumw2(kTRUE);

    hSelectedEvents = new TH1D("SelectedEvents_" + suffix, "SelectedEvents", 10, 0.5, 10.5);
    hSelectedEvents->SetDirectory(outFile);
    hSelectedEvents->Sumw2(kTRUE);




    /////////////////////////
    //    LOAD BRANCHES    //
    /////////////////////////


    // Event
    TTreeReaderValue<UInt_t>            evtNum_(reader,             "evtNumber.eventNumber");
    TTreeReaderValue<UShort_t>          nPV_(reader,                "nPV");
    TTreeReaderValue<Float_t>           met_(reader,                "met");

    TTreeReaderValue<Bool_t>            muonTrig_(reader,           "evtMuonTriggered");
    TTreeReaderValue<Bool_t>            elecTrig_(reader,           "evtElectronTriggered");
    TTreeReaderValue<UShort_t>          nMuons_(reader,             "nMuons");
    TTreeReaderValue<UShort_t>          nElecs_(reader,             "nElectrons");
    TTreeReaderValue<UShort_t>          nHZZMuons_(reader,          "nHZZMuons");
    TTreeReaderValue<UShort_t>          nHZZElecs_(reader,          "nHZZElectrons");
//  TTreeReaderValue<UShort_t>          nPartons_(reader,           "nPartons");

    TTreeReaderValue<Float_t>           genWeight_(reader,          "genWeight");
    TTreeReaderValue<Float_t>           PUWeight_(reader,           "PUWeight");


    // Muons
    TTreeReaderArray<TLorentzVector>    muonP4_(reader,             "muonP4");
    TTreeReaderValue<vector<Short_t>>   muonQ_(reader,              "muonQ");
    TTreeReaderValue<vector<Float_t>>   muonIso_(reader,            "muonCombIso");
    TTreeReaderValue<vector<Bool_t>>    muonFiredLeg1_(reader,      "muonFiredLeg1");
    TTreeReaderValue<vector<Bool_t>>    muonFiredLeg2_(reader,      "muonFiredLeg2");
    TTreeReaderValue<vector<Bool_t>>    muonIsHZZ_(reader,          "muonIsHZZ");
    TTreeReaderValue<vector<Bool_t>>    muonIsGhost_(reader,        "muonIsGhost");

    TTreeReaderValue<vector<Float_t>>   muonEnergySF_(reader,       "muonEnergySF");
    TTreeReaderValue<vector<Float_t>>   muonIDSF_(reader,           "muonHZZIDSF");
//  TTreeReaderValue<vector<Float_t>>   muonEffL1Data_(reader,      "muonTrigEffLeg1Data");
//  TTreeReaderValue<vector<Float_t>>   muonEffL1MC_(reader,        "muonTrigEffLeg1MC");
//  TTreeReaderValue<vector<Float_t>>   muonEffL2Data_(reader,      "muonTrigEffLeg2Data");
//  TTreeReaderValue<vector<Float_t>>   muonEffL2MC_(reader,        "muonTrigEffLeg2MC");


    // Electrons
    TTreeReaderArray<TLorentzVector>    elecP4_(reader,             "electronP4");
    TTreeReaderValue<vector<Short_t>>   elecQ_(reader,              "electronQ");
    TTreeReaderValue<vector<Float_t>>   elecIso_(reader,            "electronCombIso");
    TTreeReaderValue<vector<Bool_t>>    elecFiredLeg1_(reader,      "electronFiredLeg1");
    TTreeReaderValue<vector<Bool_t>>    elecFiredLeg2_(reader,      "electronFiredLeg2");
    TTreeReaderValue<vector<Bool_t>>    elecIsHZZ_(reader,          "electronIsHZZ");
    TTreeReaderValue<vector<Bool_t>>    elecPassMVA_(reader,        "electronPassNoIsoMVA");
    TTreeReaderValue<vector<Bool_t>>    elecIsGhost_(reader,        "electronIsGhost");

    TTreeReaderValue<vector<Float_t>>   elecEnergySF_(reader,       "electronEnergySF");
    TTreeReaderValue<vector<Float_t>>   elecIDSF_(reader,           "electronHZZIDSF");
//  TTreeReaderValue<vector<Float_t>>   elecEffL1Data_(reader,      "elecTrigEffLeg1Data");
//  TTreeReaderValue<vector<Float_t>>   elecEffL1MC_(reader,        "elecTrigEffLeg1MC");
//  TTreeReaderValue<vector<Float_t>>   elecEffL2Data_(reader,      "elecTrigEffLeg2Data");
//  TTreeReaderValue<vector<Float_t>>   elecEffL2MC_(reader,        "elecTrigEffLeg2MC");




    //////////////////////
    //    EVENT LOOP    //
    //////////////////////


    unsigned event = 0;
    while (reader.Next() && event < 100)
    {
//      event++;



        //--- INITIALIZE ---//

        Bool_t                  muonTrig    = *muonTrig_,           elecTrig    = *elecTrig_;
        UShort_t                nMuons      = *nMuons_,             nElecs      = *nElecs_;
        UShort_t                nHZZMuons   = *nHZZMuons_,          nHZZElecs   = *nHZZElecs_;
        Float_t                 genWeight   = *genWeight_,          PUWeight    = *PUWeight_;

        vector<TLorentzVector>  muonP4,                             elecP4;
        vector<Short_t>         muonQ,                              elecQ;
        vector<Float_t>         muonIso,                            elecIso;
        vector<Bool_t>          muonIsHZZ,                          elecIsHZZ,      elecPassMVA;
        vector<Bool_t>          muonIsGhost,                        elecIsGhost;
        vector<Bool_t>          muonFiredLeg1,  muonFiredLeg2,      elecFiredLeg1,  elecFiredLeg2;

        vector<Float_t>         muonEnergySF,                       elecEnergySF;
        vector<Float_t>         muonIDSF,                           elecIDSF;
//      vector<Float_t>         muonEffL1Data,  muonEffL2Data,      elecEffL1Data,  elecEffL2Data;
//      vector<Float_t>         muonEffL1MC,    muonEffL2MC,        elecEffL1MC,    elecEffL2MC;



        //--- PRESELECTION ---//

        // Make sure data events aren't triggered twice
        // (should really be done at BLT level...FIXME)
        if (isData && suffix.Contains("muon"))
            elecTrig = kFALSE;
        if (isData && suffix.Contains("electron"))
            muonTrig = kFALSE;

        if (!muonTrig && !elecTrig)
            continue;

        if (nHZZMuons < 2 && nHZZElecs < 2)
            continue;



        //--- FILL CONTAINERS ---//

        // Muons
        for (const TLorentzVector& muonP4__: muonP4_)
            muonP4.push_back(muonP4__);

        for (unsigned i = 0; i < nMuons; i++)
        {
            muonQ.push_back(            (*muonQ_)[i]);
            muonIso.push_back(          (*muonIso_)[i]);
            muonFiredLeg1.push_back(    (*muonFiredLeg1_)[i]);
            muonFiredLeg2.push_back(    (*muonFiredLeg2_)[i]);
            muonIsHZZ.push_back(        (*muonIsHZZ_)[i]);
            muonIsGhost.push_back(      (*muonIsGhost_)[i]);

            muonEnergySF.push_back(     (*muonEnergySF_)[i]);
            muonIDSF.push_back(         (*muonIDSF_)[i]);
//          muonEffL1Data.push_back(    (*muonEffL1Data_)[i]);
//          muonEffL1MC.push_back(      (*muonEffL1MC_)[i]);
//          muonEffL2Data.push_back(    (*muonEffL2Data_)[i]);
//          muonEffL2MC.push_back(      (*muonEffL2MC_)[i]);
        }


        // Electrons
        for (const TLorentzVector& elecP4__: elecP4_)
            elecP4.push_back(elecP4__);

        for (unsigned i = 0; i < nElecs; i++)
        {
            elecQ.push_back(            (*elecQ_)[i]);
            elecIso.push_back(          (*elecIso_)[i]);
            elecFiredLeg1.push_back(    (*elecFiredLeg1_)[i]);
            elecFiredLeg2.push_back(    (*elecFiredLeg2_)[i]);
            elecIsHZZ.push_back(        (*elecIsHZZ_)[i]);
            elecPassMVA.push_back(      (*elecPassMVA_)[i]);
            elecIsGhost.push_back(      (*elecIsGhost_)[i]);

            elecEnergySF.push_back(     (*elecEnergySF_)[i]);
            elecIDSF.push_back(         (*elecIDSF_)[i]);
//          elecEffL1Data.push_back(    (*elecEffL1Data_)[i]);
//          elecEffL1MC.push_back(      (*elecEffL1MC_)[i]);
//          elecEffL2Data.push_back(    (*elecEffL2Data_)[i]);
//          elecEffL2MC.push_back(      (*elecEffL2MC_)[i]);
        }



        //--- SYSTEMATICS ---//

        // FIXME
/*
        for (unsigned i = 0; i < nMuons; i++)
        {
            if (systOn)
            {
                if (smearPtMC && !isData)
                    lepEnergySF[i] += PT_UNC * lepEnergySF[i];
                else if (smearPtData && isData)
                    lepEnergySF[i] += PT_UNC * lepEnergySF[i];
            }
        }
*/


        //--- CORRECTIONS ---//

        // Muons
        for (unsigned i = 0; i < nMuons; i++)
        {
            // Apply energy correction
            float muonPt = muonP4[i].Pt() * muonEnergySF[i];
            muonP4[i].SetPtEtaPhiM(muonPt, muonP4[i].Eta(), muonP4[i].Phi(), MU_MASS);


            // Reevaluate Pt threshold
            if (muonP4[i].Pt() < MU_PT_MIN && muonIsHZZ[i])
            {
                muonIsHZZ[i] = kFALSE;
                nHZZMuons--;
            }


            // Remove "ghost" muons
            // (based on Pt...FIXME)
            if (muonIsGhost[i] && muonIsHZZ[i])
            {
                muonIsHZZ[i] = kFALSE;
                nHZZMuons--;
            }
        }


        // Electrons
        for (unsigned i = 0; i < nElecs; i++)
        {
            // Apply energy correction
            elecP4[i] *= elecEnergySF[i];


            // Reevaluate Pt threshold
            if (elecP4[i].Pt() < ELE_PT_MIN && elecIsHZZ[i])
            {
                elecIsHZZ[i] = kFALSE;
                nHZZElecs--;
            }


            // Remove electrons which do not pass MVA
            if (!elecPassMVA[i] && elecIsHZZ[i])
            {
                elecIsHZZ[i] = kFALSE;
                nHZZElecs--;
            }


            // Remove "ghost" electrons (cross cleaning)
            if (elecIsGhost[i] && elecIsHZZ[i])
            {
                elecIsHZZ[i] = kFALSE;
                nHZZElecs--;
            }
        }



        //--- CATEGORIZE ---//
/*
        // Gather HZZ leptons
        vector<unsigned> muonIdx, elecIdx;

        for (unsigned i = 0; i < nMuons; i++)
        {
            if (muonIsHZZ[i])
                muonIdx.push_back(i);
        }

        for (unsigned i = 0; i < nElecs; i++)
        {
            if (elecIsHZZ[i])
                elecIdx.push_back(i);
        }

        if ((muonIdx.size() != nHZZMuons) || (elecIdx.size() != nHZZElecs))
        {
            cout << "Wrong number of HZZ leptons :(" << endl;
            continue;
        }
*/

        // Assign category
        unsigned C;
        Bool_t sel2l = kFALSE, sel4l = kFALSE;
        Bool_t muLeads, muTrails;    // TRUE for muon, FALSE for electron

        // Prioritize electron-triggered events
        if (elecTrig)
        {
            muLeads = kFALSE;

            if      (nHZZElecs == 2 && nHZZMuons == 0)
            {
                C           = EE;
                sel2l       = kTRUE;
                sel4l       = kFALSE;
                muTrails    = kFALSE;   // to be safe?
            }
            else if (nHZZElecs == 2 && nHZZMuons == 2)
            {
                C           = EM;
                sel2l       = kFALSE;
                sel4l       = kTRUE;
                muTrails    = kTRUE;
            }
            else if (nHZZElecs == 4 && nHZZMuons == 0)
            {
                C           = E4;
                sel2l       = kFALSE;
                sel4l       = kTRUE;
                muTrails    = kFALSE;
            }
            else
                continue;
        }
        else if (muonTrig)
        {
            muLeads = kTRUE;

            if      (nHZZMuons == 2 && nHZZElecs == 0)
            {
                C           = MM;
                sel2l       = kTRUE;
                sel4l       = kFALSE;
                muTrails    = kTRUE;    // to be safe?
            }
            else if (nHZZMuons == 2 && nHZZElecs == 2)
            {
                C           = ME;
                sel2l       = kFALSE;
                sel4l       = kTRUE;
                muTrails    = kFALSE;
            }
            else if (nHZZMuons == 4 && nHZZElecs == 0)
            {
                C           = M4;
                sel2l       = kFALSE;
                sel4l       = kTRUE;
                muTrails    = kTRUE;
            }
            else
                continue;
        }
        else
            continue;


        // Choose cuts based on selection category
        Short_t                 lPDG            = muLeads   ? 13            : 11;                   
        Short_t                 tPDG            = muTrails  ? 13            : 11;

        Float_t                 PT1_MIN         = muLeads   ? MU_PT1_MIN    : ELE_PT1_MIN;
        Float_t                 PT2_MIN         = muLeads   ? MU_PT2_MIN    : ELE_PT2_MIN;


        // Fill NEW set of containers, wow, this code is disgusting

        // Leading pair flavor
        UShort_t                nLLeps          = muLeads   ? nMuons        : nElecs;
        UShort_t                nHZZLLeps       = muLeads   ? nHZZMuons     : nHZZElecs;
        vector<TLorentzVector>  lLepP4          = muLeads   ? muonP4        : elecP4;
        vector<Short_t>         lLepQ           = muLeads   ? muonQ         : elecQ;
        vector<Float_t>         lLepIso         = muLeads   ? muonIso       : elecIso;
        vector<Bool_t>          lLepIsHZZ       = muLeads   ? muonIsHZZ     : elecIsHZZ;
        vector<Bool_t>          lLepFiredLeg1   = muLeads   ? muonFiredLeg1 : elecFiredLeg1;
        vector<Bool_t>          lLepFiredLeg2   = muLeads   ? muonFiredLeg2 : elecFiredLeg2;
        vector<Float_t>         lLepIDSF        = muLeads   ? muonIDSF      : elecIDSF;
//      vector<Float_t>         lLepEffL1Data   = muLeads   ? muonEffL1Data : elecEffL1Data;
//      vector<Float_t>         lLepEffL2Data   = muLeads   ? muonEffL2Data : elecEffL2Data;
//      vector<Float_t>         lLepEffL1MC     = muLeads   ? muonEffL1MC   : elecEffL1MC;
//      vector<Float_t>         lLepEffL2MC     = muLeads   ? muonEffL2MC   : elecEffL2MC;

        // Subleading pair flavor
        UShort_t                nTLeps          = muTrails  ? nMuons        : nElecs;
        UShort_t                nHZZTLeps       = muTrails  ? nHZZMuons     : nHZZElecs;
        vector<TLorentzVector>  tLepP4          = muTrails  ? muonP4        : elecP4;
        vector<Short_t>         tLepQ           = muTrails  ? muonQ         : elecQ;
        vector<Float_t>         tLepIso         = muTrails  ? muonIso       : elecIso;
        vector<Bool_t>          tLepIsHZZ       = muTrails  ? muonIsHZZ     : elecIsHZZ;
        vector<Float_t>         tLepIDSF        = muTrails  ? muonIDSF      : elecIDSF;
        // Don't need subleading flavor trigger info?




        ////////////////////
        //    DILEPTON    //
        ////////////////////


        if (sel2l)
        {

            //--- LEADING LEPTON ---//

            // (for some reason I decided it needed to be the leading lepton?)
            // (I hope this isn't wrong...)
            if (!lLepIsHZZ[0])
                continue;
            if (lLepP4[0].Pt() < PT1_MIN)
                continue;
            if (!lLepFiredLeg1[0])
                continue;


            // Pair leading lepton
            unsigned x = 0, y = 0;
            bool foundPair = kFALSE;
            for (unsigned j = 1; j < nLLeps; j++)
            {
                if (!lLepIsHZZ[j])
                    continue;
                if (lLepP4[j].Pt() < PT2_MIN)
                    continue;
                if (!lLepFiredLeg2[j])
                    continue;

                // Charge requirement
                if (lLepQ[x] == lLepQ[j])
                    continue;

                // Mass cut
                Float_t mll = (lLepP4[x] + lLepP4[j]).M();
                if (mll > M_MIN && mll < M_MAX)
                {
                    y = j;
                    foundPair = kTRUE;
                    break;
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
                /*
                // Trigger efficiency
                //                  if ((lLepFiredLeg1[x] || lLepFiredLeg2[x]) && (lLepFiredLeg1[y] || lLepFiredLeg2[y]))
                if (kFALSE)
                {
                float trigEffData, trigEffMC;

                trigEffData  = lLepEffLeg1Data[x] * lLepEffLeg2Data[y];
                trigEffData += lLepEffLeg1Data[y] * lLepEffLeg2Data[x];
                trigEffData -= lLepEffLeg1Data[x] * lLepEffLeg1Data[y];

                trigEffMC  = lLepEffLeg1MC[x] * lLepEffLeg2MC[y];
                trigEffMC += lLepEffLeg1MC[y] * lLepEffLeg2MC[x];
                trigEffMC -= lLepEffLeg1MC[x] * lLepEffLeg1MC[y];
                if (trigEffMC == 0)
                trigEffMC = 1;

                triggerWeight = trigEffData / trigEffMC;
                }
                */
                // ID efficiency
                idWeight = lLepIDSF[x] * lLepIDSF[y];
            }
            weight = genWeight * PUWeight * triggerWeight * idWeight;



            //--- HISTOGRAMS ---//

            hSelectedEvents->Fill(1, weight);
            hSelectedEvents->Fill(binIdx[C], weight);

            hTotalEvents->Fill(6);
            hTotalEvents->Fill(7, weight);



            //--- CONTAINERS ---//

            evtNum  = *evtNum_;
            nPV     = *nPV_;                            met = *met_;
            l1p4    = lLepP4[x];                        l2p4    = lLepP4[y];
            l1pdg   = lLepQ[x] * lPDG;                  l2pdg   = lLepQ[y] * lPDG;
            l1iso   = lLepIso[x] / l1p4.Pt();           l2iso   = lLepIso[y] / l2p4.Pt();
            z1p4    = l1p4 + l2p4;

            tree[C]->Fill();
        }




        ////////////////////
        //    4-LEPTON    //
        ////////////////////


        else if (sel4l)
        {
            // Containers
            vector<tuple<TLorentzVector, Short_t, Float_t>> allLeps, z1leps, z2leps;




            /////////////////////////
            //    SINGLE FLAVOR    //
            /////////////////////////


            if (muLeads == muTrails)
            {

                //--- Z SELECTION ---//

                vector<pair<unsigned, unsigned>> zCand;
                for (unsigned j = 1; j < nLLeps; j++)
                {
                    if (!lLepIsHZZ[j])
                        continue;

                    for (unsigned i = 0; i < j; i++)
                    {
                        if (!lLepIsHZZ[i])
                            continue;
                        if (lLepQ[i] == lLepQ[j])
                            continue;

                        Float_t mll = (lLepP4[i] + lLepP4[j]).M();
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
                        if (
                                Zi.first == Zj.first    || Zi.first == Zj.second
                                ||  Zi.second == Zj.first   || Zi.second == Zj.second
                           )
                            continue;


                        // Choose Z1 closest to nominal Z mass
                        Float_t mll_i = (lLepP4[Zi.first] + lLepP4[Zi.second]).M();
                        Float_t mll_j = (lLepP4[Zj.first] + lLepP4[Zj.second]).M();
                        pair<unsigned, unsigned> Z1, Z2;

                        if (fabs(mll_i - Z_MASS) < fabs(mll_j - Z_MASS))
                        {
                            Z1 = Zi;
                            Z2 = Zj;
                        }
                        else
                        {
                            Z1 = Zj;
                            Z2 = Zi;
                        }


                        // Pair Pt, Q for easier manipulation
                        pair<TLorentzVector, TLorentzVector> Z1_p4s, Z2_p4s;
                        Z1_p4s = make_pair(lLepP4[Z1.first], lLepP4[Z1.second]);
                        Z2_p4s = make_pair(lLepP4[Z2.first], lLepP4[Z2.second]);

                        pair<Short_t, Short_t> Z1_qs, Z2_qs;
                        Z1_qs = make_pair(lLepQ[Z1.first], lLepQ[Z1.second]);
                        Z2_qs = make_pair(lLepQ[Z2.first], lLepQ[Z2.second]);

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
                        mll_x = (lLepP4[Zx.first] + lLepP4[Zx.second]).M();
                        mll_y = (lLepP4[Zy.first] + lLepP4[Zy.second]).M();

                        if (fabs(mll_x - Z_MASS) < fabs(mll_y - Z_MASS))
                        {
                            Za      = Zx;               Zb      = Zy;
                            mll_a   = mll_x;            mll_b   = mll_y;
                        }
                        else
                        {
                            Za      = Zy;               Zb      = Zx;
                            mll_a   = mll_y;            mll_b   = mll_x;
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
                            selLepP4.push_back(lLepP4[z[k]]);
                            selLepQ.push_back(lLepQ[z[k]]);
                        }

                        if (selLepP4[0].Pt() < PT1_MIN || selLepP4[1].Pt() < PT2_MIN)
                            continue;


                        // Trigger requirement
                        // (double check?)

                        bool firedLeg1 = kFALSE,    firedLeg2 = kFALSE;

                        for (unsigned k = 0; k < z.size(); k++)
                        {
                            if (lLepFiredLeg1[z[k]])
                                firedLeg1 = kTRUE;
                            if (lLepFiredLeg2[z[k]])
                                firedLeg2 = kTRUE;
                        }
                        if (!firedLeg1 || !firedLeg2)
                        {
                            cout << "Failed trigger requirement" << endl;
                            continue;
                        }


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

                // FIXME i don't think this actually does anything
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

                // Assemble indices of selected pairs
                pair<pair<unsigned, unsigned>, pair<unsigned, unsigned>> ZZ = zzCand[zz];
                pair<unsigned, unsigned> Z1 = ZZ.first, Z2 = ZZ.second;
                unsigned x = Z1.first, y = Z1.second, xx = Z2.first, yy = Z2.second;

                z1leps.push_back(make_tuple(lLepP4[x],  lLepQ[x]*lPDG,  lLepIso[x]/lLepP4[x].Pt()));
                z1leps.push_back(make_tuple(lLepP4[y],  lLepQ[y]*lPDG,  lLepIso[y]/lLepP4[y].Pt()));

                z2leps.push_back(make_tuple(lLepP4[xx], lLepQ[xx]*lPDG, lLepIso[xx]/lLepP4[xx].Pt()));
                z2leps.push_back(make_tuple(lLepP4[yy], lLepQ[yy]*lPDG, lLepIso[yy]/lLepP4[yy].Pt()));


                // Sort lepton PDG ID and iso by Pt
                allLeps.push_back(make_tuple(lLepP4[x],  lLepQ[x]*lPDG,  lLepIso[x]/lLepP4[x].Pt()));
                allLeps.push_back(make_tuple(lLepP4[y],  lLepQ[y]*lPDG,  lLepIso[y]/lLepP4[y].Pt()));
                allLeps.push_back(make_tuple(lLepP4[xx], lLepQ[xx]*lPDG, lLepIso[xx]/lLepP4[xx].Pt()));
                allLeps.push_back(make_tuple(lLepP4[yy], lLepQ[yy]*lPDG, lLepIso[yy]/lLepP4[yy].Pt()));
                sort(allLeps.begin(), allLeps.end(), SortDecPt);



                //--- WEIGHTING ---//

                float triggerWeight = 1, idWeight = 1;

                if (!isData)
                {
                    // Trigger efficiency
                    // Uses only leading leptons, FIXME
/*
//                  if ((lLepFiredLeg1[x] || lLepFiredLeg2[x]) && (lLepFiredLeg1[y] || lLepFiredLeg2[y]))
                    if (kFALSE)
                    {
                    float trigEffData, trigEffMC;

                    trigEffData  = lLepEffLeg1Data[x] * lLepEffLeg2Data[y];
                    trigEffData += lLepEffLeg1Data[y] * lLepEffLeg2Data[x];
                    trigEffData -= lLepEffLeg1Data[x] * lLepEffLeg1Data[y];

                    trigEffMC  = lLepEffLeg1MC[x] * lLepEffLeg2MC[y];
                    trigEffMC += lLepEffLeg1MC[y] * lLepEffLeg2MC[x];
                    trigEffMC -= lLepEffLeg1MC[x] * lLepEffLeg1MC[y];
                    if (trigEffMC == 0)
                    trigEffMC = 1;

                    triggerWeight = trigEffData / trigEffMC;
                    }
*/

                    // ID efficiency
                    idWeight *= lLepIDSF[x] * lLepIDSF[y] * lLepIDSF[xx] * lLepIDSF[yy];
                }
                weight = genWeight * PUWeight * triggerWeight * idWeight;
            }




            ////////////////////////
            //    MIXED FLAVOR    //
            ////////////////////////


            else if (muLeads != muTrails)
            {


                //--- Z SELECTION ---//

                vector<pair<unsigned, unsigned>> z1Cand, z2Cand;

                // Z1 candidates
                for (unsigned j = 1; j < nLLeps; j++)
                {
                    if (!lLepIsHZZ[j])
                        continue;
                    if (!lLepFiredLeg2[j])
                        continue;

                    for (unsigned i = 0; i < j; i++)
                    {
                        if (!lLepIsHZZ[i])
                            continue;
                        if (!lLepFiredLeg1[i])
                            continue;
                        if (lLepQ[i] == lLepQ[j])
                            continue;

                        Float_t mll = (lLepP4[i] + lLepP4[j]).M();
                        if (mll > MZ1_MIN && mll < MZ_MAX)
                            z1Cand.push_back(make_pair(i, j));
                    }
                }
                if (z1Cand.size() < 1)
                    continue;


                // Z2 candidates
                for (unsigned j = 1; j < nTLeps; j++)
                {
                    if (!tLepIsHZZ[j])
                        continue;

                    for (unsigned i = 0; i < j; i++)
                    {
                        if (!tLepIsHZZ[i])
                            continue;
                        if (tLepQ[i] == tLepQ[j])
                            continue;

                        Float_t mll = (tLepP4[i] + tLepP4[j]).M();
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
                        Z1_p4s = make_pair(lLepP4[Z1.first], lLepP4[Z1.second]);
                        Z2_p4s = make_pair(tLepP4[Z2.first], tLepP4[Z2.second]);

                        pair<Short_t, Short_t> Z1_qs, Z2_qs;
                        Z1_qs = make_pair(lLepQ[Z1.first], lLepQ[Z1.second]);
                        Z2_qs = make_pair(tLepQ[Z2.first], tLepQ[Z2.second]);

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

                z1leps.push_back(make_tuple(lLepP4[x],  lLepQ[x]*lPDG,  lLepIso[x]/lLepP4[x].Pt()));
                z1leps.push_back(make_tuple(lLepP4[y],  lLepQ[y]*lPDG,  lLepIso[y]/lLepP4[y].Pt()));

                z2leps.push_back(make_tuple(tLepP4[xx], tLepQ[xx]*lPDG, tLepIso[xx]/tLepP4[xx].Pt()));
                z2leps.push_back(make_tuple(tLepP4[yy], tLepQ[yy]*lPDG, tLepIso[yy]/tLepP4[yy].Pt()));


                // Sort lepton PDG ID and iso by Pt
                allLeps.push_back(make_tuple(lLepP4[x], lLepQ[x]*lPDG, lLepIso[x]/lLepP4[x].Pt()));
                allLeps.push_back(make_tuple(lLepP4[y], lLepQ[y]*lPDG, lLepIso[y]/lLepP4[y].Pt()));
                allLeps.push_back(make_tuple(tLepP4[xx], tLepQ[xx]*tPDG, tLepIso[xx]/tLepP4[xx].Pt()));
                allLeps.push_back(make_tuple(tLepP4[yy], tLepQ[yy]*tPDG, tLepIso[yy]/tLepP4[yy].Pt()));
                sort(allLeps.begin(), allLeps.end(), SortDecPt);



                //--- WEIGHTING ---//

                float triggerWeight = 1, idWeight = 1;

                if (!isData)
                {
                    // Trigger efficiency
                    // Uses only leading leptons, FIXME
                    // (oops, looks like there isn't even anything to fix)


                    // ID efficiency
                    idWeight *= lLepIDSF[x] * lLepIDSF[y] * tLepIDSF[xx] * tLepIDSF[yy];
                }
                weight = genWeight * PUWeight * triggerWeight * idWeight;
            }
            else
                continue;




            /////////////////////
            //    FILL TREE    //
            /////////////////////


            if (allLeps.size() != 4)
                continue;



            //--- HISTOGRAMS ---//

            hSelectedEvents->Fill(1, weight);
            hSelectedEvents->Fill(binIdx[C], weight);

            hTotalEvents->Fill(6);
            hTotalEvents->Fill(7, weight);



            //--- CONTAINERS ---//


            evtNum  = *evtNum_;
            nPV     = *nPV_;                            met = *met_;
            z1p4    = GetP4Sum(z1leps);                 z2p4 = GetP4Sum(z2leps); 
            zzp4    = z1p4 + z2p4;
            tie(l1p4, l1pdg, l1iso) = allLeps[0];       tie(l2p4, l2pdg, l2iso) = allLeps[1];
            tie(l3p4, l3pdg, l3iso) = allLeps[2];       tie(l4p4, l4pdg, l4iso) = allLeps[3];



            //--- BOOST ---//

            // Sort leptons by magnitude of momentum (P) in Z COM frame

            TVector3 zzb3 = zzp4.BoostVector();
            Float_t iso;    // this is really just a dummy variable

            vector<tuple<TLorentzVector, Short_t, Float_t>> b_allLeps, b_z1leps, b_z2leps;
            for (unsigned i = 0; i < allLeps.size(); i++)
            {
                TLorentzVector p4;
                Short_t pdg;
                tie(p4, pdg, iso) = allLeps[i];
                p4 = GetBoosted(p4, zzb3);
                b_allLeps.push_back(make_tuple(p4, pdg, iso));
            }
            sort(b_allLeps.begin(), b_allLeps.end(), SortDecP);

            for (unsigned i = 0; i < z1leps.size(); i++)
            {
                TLorentzVector p4;
                Short_t pdg;
                tie(p4, pdg, iso) = z1leps[i];
                p4 = GetBoosted(p4, zzb3);
                b_z1leps.push_back(make_tuple(p4, pdg, iso));
            }
            sort(b_z1leps.begin(), b_z1leps.end(), SortDecP);

            for (unsigned i = 0; i < z2leps.size(); i++)
            {
                TLorentzVector p4;
                Short_t pdg;
                tie(p4, pdg, iso) = z2leps[i];
                p4 = GetBoosted(p4, zzb3);
                b_z2leps.push_back(make_tuple(p4, pdg, iso));
            }
            sort(b_z2leps.begin(), b_z2leps.end(), SortDecP);

            tie(b_l1p4, b_l1pdg, iso) = b_allLeps[0];       tie(b_l2p4, b_l2pdg, iso) = b_allLeps[1];
            tie(b_l3p4, b_l3pdg, iso) = b_allLeps[2];       tie(b_l4p4, b_l4pdg, iso) = b_allLeps[3];

            TLorentzVector b_z1l1p4 = get<0>(b_z1leps[0]),  b_z1l2p4 = get<0>(b_z1leps[1]);
            TLorentzVector b_z2l1p4 = get<0>(b_z2leps[0]),  b_z2l2p4 = get<0>(b_z2leps[1]);


            // Pairs
            b_z1p4  = GetBoosted(z1p4, zzb3);               b_z2p4  = GetBoosted(z2p4, zzb3);

            // "Trailing trio"
            b_ttp4  = b_l2p4 + b_l3p4 + b_l4p4;


            // ANGLES
            // "phi":   angle between pairs
            TLorentzVector  z1lpp4  = get<1>(z1leps[0]) > 0 ? get<0>(z1leps[0]) : get<0>(z1leps[1]);
            TLorentzVector  z1lmp4  = get<1>(z1leps[0]) < 0 ? get<0>(z1leps[0]) : get<0>(z1leps[1]);
            TLorentzVector  z2lpp4  = get<1>(z2leps[0]) > 0 ? get<0>(z2leps[0]) : get<0>(z2leps[1]);
            TLorentzVector  z2lmp4  = get<1>(z2leps[0]) < 0 ? get<0>(z2leps[0]) : get<0>(z2leps[1]);
            TLorentzVector  b_z1lpp4 = GetBoosted(z1lpp4, zzb3);
            TLorentzVector  b_z1lmp4 = GetBoosted(z1lmp4, zzb3);
            TLorentzVector  b_z2lpp4 = GetBoosted(z2lpp4, zzb3);
            TLorentzVector  b_z2lmp4 = GetBoosted(z2lmp4, zzb3);

            // (normal to z1,z2 decay plane)
            TVector3        b_z1n   = b_z1lpp4.Vect().Cross(b_z1lmp4.Vect());
            TVector3        b_z2n   = b_z2lpp4.Vect().Cross(b_z2lmp4.Vect());
            b_phi   = b_z1n.Angle(b_z2n);

            // "theta": angle between trailing pair 1 lepton and low-mass pair
            b_theta = b_z1l2p4.Angle(b_z2p4.Vect());

            // "alpha": angle between paired leptons
            b_z1alpha = b_z1l1p4.Angle(b_z1l2p4.Vect());
            b_z2alpha = b_z2l1p4.Angle(b_z2l2p4.Vect());



            // OTHERS
            TVector3        z1b3    = z1p4.BoostVector(),       z2b3    = z2p4.BoostVector();

            TLorentzVector  bb_z1p4 = GetBoosted(z1p4, z2b3),   bb_z2p4 = GetBoosted(z2p4, z1b3);
            TLorentzVector  bb_z1lpp4 = GetBoosted(z1lpp4, z1b3);
            TLorentzVector  bb_z2lpp4 = GetBoosted(z2lpp4, z2b3);

            bb_z1theta  = bb_z1lpp4.Angle(bb_z2p4.Vect());
            bb_z2theta  = bb_z2lpp4.Angle(bb_z1p4.Vect());


            tree[C]->Fill();

        }   // END 4l selection

//      event++;

    }   // END event loop



    //--- CLEANUP ---//

    // Write to file
    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();
    hTotalEvents->Write();
    hSelectedEvents->Write();
    hFiducialEvents->Write();
    hPhaseSpaceEvents->Write();

    outFile->Purge();
    outFile->Close();
    file->Close();



//  delete outFile;
//  delete tree;
//  delete hSelectedEvents;
}




////////////////////////////
//    HELPER FUNCTIONS    //
////////////////////////////


bool SortDecPt( const tuple<TLorentzVector, Short_t, Float_t> &i_,
                const tuple<TLorentzVector, Short_t, Float_t> &j_)
{
    TLorentzVector i = get<0>(i_), j = get<0>(j_);
    return (i.Pt() > j.Pt());
}


bool SortDecP(  const tuple<TLorentzVector, Short_t, Float_t> &i_,
                const tuple<TLorentzVector, Short_t, Float_t> &j_)
{
    TLorentzVector i = get<0>(i_), j = get<0>(j_);
    return (i.P() > j.P());
}


TLorentzVector GetBoosted(const TLorentzVector &p4_, const TVector3 &beta)
{
    TLorentzVector p4(p4_);
    p4.Boost(-beta);
    return p4;
}


TLorentzVector GetP4Sum(const vector<tuple<TLorentzVector, Short_t, Float_t>> &leps)
{
    TLorentzVector p4sum;
    for (unsigned i = 0; i < leps.size(); i++)
        p4sum += get<0>(leps[i]);
    return p4sum;
}



/*
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
    else if (xval < hist->GetXaxis()->GetBinLowEdge(1))
        bin = 1;
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
    else if (yval < hist->GetYaxis()->GetBinLowEdge(1))
        bin = 1;
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
*/
