// STL
#include <vector>
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
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"
#include "SelectionTools.hh"

// Cuts
#include "Cuts2017.hh"
//#include "Cuts2016.hh"

using namespace std;



/* 
**  RecoSelection
**
**
**  This is "the" selection script.  Selects ll and 4l events based on cuts from CutsXXXX.hh.
**  It is reasonably self contained, but some algorithms are imported from SelectionTools.hh.
**
**  Reads from post-BLT analyzer (MultileptonAnalyzer) ntuples.  Output is split into ll (mumu, ee)
**  and 4l (4m, 2m2e, 2e2m, 4e) channels.  Only information relevant to each (ll, 4l) channel is 
**  included.  Copies over gen-level particle information for signal (zz_4l) events.
**
**  Integrates systematics analysis using histogeams in /data directory.
*/

void RecoSelection( const TString suffix,           const TString id,
                    const TString systematics = "", const TString idH = "0")
{


    //
    //  OPTIONS
    //

    bool debug = kFALSE;
    long selectEvents = 100;
    int  printEvery = 30000;

    // Systematics toggle
    const bool smearMuonID = systematics.EqualTo("muonID");
    const bool smearElecID = systematics.EqualTo("electronID");
    const bool smearElecReco = systematics.EqualTo("electronReco");
    const bool smearMuonPt = systematics.EqualTo("muonPt");
    const bool smearElecPt = systematics.EqualTo("electronPt");

    const bool smearOn = smearMuonID || smearElecID || smearElecReco || smearMuonPt || smearElecPt;



    //
    //  SAMPLE INFO
    //

    const bool isData   = suffix.Contains(YEAR_STR);
    const bool isSignal = suffix.EqualTo("zz_4l");

    const unsigned N = 8;   // Channel indices
    unsigned                   LL = 0, MM = 1, EE = 2, L4 = 3, M4 = 4, ME = 5, EM = 6, E4 = 7;
    TString selection[N]    = {"ll",   "mumu", "ee",   "4l",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = {2,      3,      4,      5,      6,      7,      8,      9};



    //
    //  OUTPUT FILE
    //

    TString prefix  = smearOn ? systematics + idH : "selected";
    TString outName = prefix + "_" + suffix + "_" + id + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);



    //
    //  OUTPUT BRANCHES
    //

    // Event info
    Int_t               runNum,     evtNum,     lumiSec;
    UShort_t            nPV;
    Float_t             weight,     genWeight,  qtWeight,   puWeight,   ecalWeight;
    Float_t             trigWeight, idWeight,   recoWeight;
    UInt_t              channel;

    // Pairs
    TLorentzVector      z1p4,       z2p4,       zzp4;
    Short_t             z1pdg,      z2pdg;

    // Leptons
    TLorentzVector      l1p4,       l2p4,       l3p4,       l4p4;
    Short_t             l1pdg,      l2pdg,      l3pdg,      l4pdg;
    Float_t             l1iso,      l2iso,      l3iso,      l4iso;
    UShort_t            l1z,        l2z,        l3z,        l4z;


    // Gen-level info (only written for signal)
    UShort_t            nHardProcMuons,             nHardProcElectrons,         nHardProcLeptons;
    UShort_t            nFinalStateMuons,           nFinalStateElectrons,       nFinalStateLeptons;

    TClonesArray        *hardProcMuonP4         =   new TClonesArray("TLorentzVector");
    TClonesArray        *hardProcElectronP4     =   new TClonesArray("TLorentzVector");
    vector<Short_t>     *hardProcMuonQ          = 0,        *hardProcElectronQ          = 0;
    vector<UShort_t>    *hardProcMuonZIndex     = 0,        *hardProcElectronZIndex     = 0;
    TLorentzVector      *hardProcLeptonsP4      = 0;

    TClonesArray        *finalStateMuonP4       =   new TClonesArray("TLorentzVector");
    TClonesArray        *finalStateElectronP4   =   new TClonesArray("TLorentzVector");
    vector<Short_t>     *finalStateMuonQ        = 0,        *finalStateElectronQ        = 0;
    vector<UShort_t>    *finalStateMuonZIndex   = 0,        *finalStateElectronZIndex   = 0;
    TLorentzVector      *finalStateLeptonsP4    = 0;

    TLorentzVector      u_l1p4,     u_l2p4,     u_l3p4,     u_l4p4;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("runNum",       &runNum);       tree[i]->Branch("evtNum",       &evtNum);               
        tree[i]->Branch("lumiSec",      &lumiSec);      tree[i]->Branch("nPV",          &nPV);
        tree[i]->Branch("weight",       &weight);       tree[i]->Branch("genWeight",    &genWeight);
        tree[i]->Branch("qtWeight",     &qtWeight);     tree[i]->Branch("puWeight",     &puWeight);
        tree[i]->Branch("ecalWeight",   &ecalWeight);   tree[i]->Branch("trigWeight",   &trigWeight);
        tree[i]->Branch("idWeight",     &idWeight);     tree[i]->Branch("recoWeight",   &recoWeight);
        tree[i]->Branch("channel",      &channel);

        if (i >= L4) {  tree[i]->Branch("zzp4", &zzp4);}
                        tree[i]->Branch("z1p4", &z1p4);     tree[i]->Branch("z1pdg",    &z1pdg);
        if (i >= L4) {  tree[i]->Branch("z2p4", &z2p4);     tree[i]->Branch("z2pdg",    &z2pdg);}

        tree[i]->Branch("l1p4",     &l1p4);                 tree[i]->Branch("l1pdg",    &l1pdg);
        tree[i]->Branch("l1iso",    &l1iso);   if (i >= L4) tree[i]->Branch("l1z",      &l1z);
        tree[i]->Branch("l2p4",     &l2p4);                 tree[i]->Branch("l2pdg",    &l2pdg);
        tree[i]->Branch("l2iso",    &l2iso);   if (i >= L4) tree[i]->Branch("l2z",      &l2z);

        if (i >= L4)
        {
            tree[i]->Branch("l3p4",     &l3p4);             tree[i]->Branch("l3pdg",    &l3pdg);
            tree[i]->Branch("l3iso",    &l3iso);            tree[i]->Branch("l3z",      &l3z);
            tree[i]->Branch("l4p4",     &l4p4);             tree[i]->Branch("l4pdg",    &l4pdg);
            tree[i]->Branch("l4iso",    &l4iso);            tree[i]->Branch("l4z",      &l4z);
        }


        if (isSignal && i >= L4)
        {
            tree[i]->Branch("nFinalStateMuons",         &nFinalStateMuons);
            tree[i]->Branch("nFinalStateElectrons",     &nFinalStateElectrons);
            tree[i]->Branch("nFinalStateLeptons",       &nFinalStateLeptons);
            tree[i]->Branch("nHardProcMuons",           &nHardProcMuons);
            tree[i]->Branch("nHardProcElectrons",       &nHardProcElectrons);
            tree[i]->Branch("nHardProcLeptons",         &nHardProcLeptons);

            tree[i]->Branch("hardProcMuonP4",           &hardProcMuonP4,            32000,      1);
            tree[i]->Branch("hardProcMuonQ",            &hardProcMuonQ);
            tree[i]->Branch("hardProcMuonZIndex",       &hardProcMuonZIndex);
            tree[i]->Branch("hardProcElectronP4",       &hardProcElectronP4,        32000,      1);
            tree[i]->Branch("hardProcElectronQ",        &hardProcElectronQ);
            tree[i]->Branch("hardProcElectronZIndex",   &hardProcElectronZIndex);
            tree[i]->Branch("hardProcLeptonsP4",        &hardProcLeptonsP4);

            tree[i]->Branch("finalStateMuonP4",         &finalStateMuonP4,          32000,      1);
            tree[i]->Branch("finalStateMuonQ",          &finalStateMuonQ);
            tree[i]->Branch("finalStateMuonZIndex",     &finalStateMuonZIndex);
            tree[i]->Branch("finalStateElectronP4",     &finalStateElectronP4,      32000,      1);
            tree[i]->Branch("finalStateElectronQ",      &finalStateElectronQ);
            tree[i]->Branch("finalStateElectronZIndex", &finalStateElectronZIndex);
            tree[i]->Branch("finalStateLeptonsP4",      &finalStateLeptonsP4);

            tree[i]->Branch("uncorr_l1p4",  &u_l1p4);   tree[i]->Branch("uncorr_l2p4",  &u_l2p4);
            tree[i]->Branch("uncorr_l3p4",  &u_l3p4);   tree[i]->Branch("uncorr_l4p4",  &u_l4p4);
        }
    }



    //
    //  INPUT FILE
    //

    TString inName  = suffix + "_" + id + ".root";
    TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/";
    TFile   *inFile = TFile::Open(inPath + inName);

    cout << endl << endl << "Opened " << inPath + inName << endl;



    //
    //  INPUT BRANCHES
    //

    TTreeReader reader("tree_" + suffix, inFile);

    TTreeReaderValue    <Int_t>             runNum_         (reader,    "runNumber");
    TTreeReaderValue    <Int_t>             evtNum_         (reader,    "evtNumber");
    TTreeReaderValue    <Int_t>             lumiSec_        (reader,    "lumiSection");
    TTreeReaderValue    <Float_t>           genWeight_      (reader,    "genWeight");
    TTreeReaderValue    <Float_t>           ecalWeight_     (reader,    "ECALWeight");
    TTreeReaderValue    <Float_t>           puWeight_       (reader,    "PUWeight");
    TTreeReaderValue    <UShort_t>          nPV_            (reader,    "nPV");

    TTreeReaderValue    <Bool_t>            muonTrig_       (reader,    "evtMuonTriggered");
    TTreeReaderValue    <Bool_t>            elecTrig_       (reader,    "evtElectronTriggered");
    TTreeReaderValue    <UShort_t>          nMuons_         (reader,    "nLooseMuons");
    TTreeReaderValue    <UShort_t>          nElecs_         (reader,    "nLooseElectrons");
    TTreeReaderValue    <UShort_t>          nTightMuons_    (reader,    "nTightMuons");
    TTreeReaderValue    <UShort_t>          nTightElecs_    (reader,    "nTightElectrons");

    TTreeReaderArray    <TLorentzVector>    muonP4_         (reader,    "muonP4");
    TTreeReaderArray    <TLorentzVector>    muonUncorrP4_   (reader,    "muonUncorrectedP4");
    TTreeReaderValue    <vector<Short_t>>   muonQ_          (reader,    "muonQ");
    TTreeReaderValue    <vector<Float_t>>   muonIso_        (reader,    "muonIsolation");
    TTreeReaderValue    <vector<Bool_t>>    muonIsTight_    (reader,    "muonIsTight");
    TTreeReaderValue    <vector<Float_t>>   muonIDSF_       (reader,    "muonIDSF");
    TTreeReaderValue    <vector<Bool_t>>    muonFiredLeg1_  (reader,    "muonFiredLeg1");
    TTreeReaderValue    <vector<Float_t>>   muonEffL1Data_  (reader,    "muonTrigEffLeg1Data");
    TTreeReaderValue    <vector<Float_t>>   muonEffL1MC_    (reader,    "muonTrigEffLeg1MC");
    TTreeReaderValue    <vector<Bool_t>>    muonFiredLeg2_  (reader,    "muonFiredLeg2");
    TTreeReaderValue    <vector<Float_t>>   muonEffL2Data_  (reader,    "muonTrigEffLeg2Data");
    TTreeReaderValue    <vector<Float_t>>   muonEffL2MC_    (reader,    "muonTrigEffLeg2MC");

    TTreeReaderArray    <TLorentzVector>    elecP4_         (reader,    "electronP4");
    TTreeReaderArray    <TLorentzVector>    elecUncorrP4_   (reader,    "electronUncorrectedP4");
    TTreeReaderValue    <vector<Short_t>>   elecQ_          (reader,    "electronQ");
    TTreeReaderValue    <vector<Float_t>>   elecIso_        (reader,    "electronIsolation");
    TTreeReaderValue    <vector<Bool_t>>    elecIsTight_    (reader,    "electronIsTight");
    TTreeReaderValue    <vector<Float_t>>   elecIDSF_       (reader,    "electronIDSF");
    TTreeReaderValue    <vector<Float_t>>   elecRecoSF_     (reader,    "electronRecoSF");
    TTreeReaderValue    <vector<Bool_t>>    elecFiredLeg1_  (reader,    "electronFiredLeg1");
    TTreeReaderValue    <vector<Float_t>>   elecEffL1Data_  (reader,    "electronTrigEffLeg1Data");
    TTreeReaderValue    <vector<Float_t>>   elecEffL1MC_    (reader,    "electronTrigEffLeg1MC");
    TTreeReaderValue    <vector<Bool_t>>    elecFiredLeg2_  (reader,    "electronFiredLeg2");
    TTreeReaderValue    <vector<Float_t>>   elecEffL2Data_  (reader,    "electronTrigEffLeg2Data");
    TTreeReaderValue    <vector<Float_t>>   elecEffL2MC_    (reader,    "electronTrigEffLeg2MC");

    cout << "Loaded branches" << endl;


    // Gen particles

    // Only read for signal events, so we have to do it the old-fashioned way

    TTree *inTree;

    if (isSignal && YEAR_STR.EqualTo("2017") && !smearOn)
    {
        inFile->GetObject("tree_" + suffix, inTree);

        inTree->SetBranchAddress(   "nFinalStateMuons",             &nFinalStateMuons);
        inTree->SetBranchAddress(   "nFinalStateElectrons",         &nFinalStateElectrons);
        inTree->SetBranchAddress(   "nFinalStateLeptons",           &nFinalStateLeptons);
        inTree->SetBranchAddress(   "nHardProcMuons",               &nHardProcMuons);
        inTree->SetBranchAddress(   "nHardProcElectrons",           &nHardProcElectrons);
        inTree->SetBranchAddress(   "nHardProcLeptons",             &nHardProcLeptons);

        inTree->GetBranch(          "hardProcMuonP4")               ->SetAutoDelete(kFALSE);
        inTree->SetBranchAddress(   "hardProcMuonP4",               &hardProcMuonP4);
        inTree->SetBranchAddress(   "hardProcMuonQ",                &hardProcMuonQ);
        inTree->SetBranchAddress(   "hardProcMuonZIndex",           &hardProcMuonZIndex);
        inTree->GetBranch(          "hardProcElectronP4")           ->SetAutoDelete(kFALSE);
        inTree->SetBranchAddress(   "hardProcElectronP4",           &hardProcElectronP4);
        inTree->SetBranchAddress(   "hardProcElectronQ",            &hardProcElectronQ);
        inTree->SetBranchAddress(   "hardProcElectronZIndex",       &hardProcElectronZIndex);
        inTree->SetBranchAddress(   "hardProcLeptonsP4",            &hardProcLeptonsP4);

        inTree->GetBranch(          "finalStateMuonP4")             ->SetAutoDelete(kFALSE);
        inTree->SetBranchAddress(   "finalStateMuonP4",             &finalStateMuonP4);
        inTree->SetBranchAddress(   "finalStateMuonQ",              &finalStateMuonQ);
        inTree->SetBranchAddress(   "finalStateMuonZIndex",         &finalStateMuonZIndex);
        inTree->GetBranch(          "finalStateElectronP4")         ->SetAutoDelete(kFALSE);
        inTree->SetBranchAddress(   "finalStateElectronP4",         &finalStateElectronP4);
        inTree->SetBranchAddress(   "finalStateElectronQ",          &finalStateElectronQ);
        inTree->SetBranchAddress(   "finalStateElectronZIndex",     &finalStateElectronZIndex);
        inTree->SetBranchAddress(   "finalStateLeptonsP4",          &finalStateLeptonsP4);

        cout << "Loaded gen branches" << endl;
    }



    //
    //  HISTOGRAMS
    //

    TH1D *hTotalEvents, *hSelectedEvents;
    
    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);    // we are currently in the domain of inFile
    hTotalEvents->Sumw2();

    hSelectedEvents = new TH1D("SelectedEvents_" + suffix, "SelectedEvents", 10, 0.5, 10.5);
    hSelectedEvents->SetDirectory(outFile);
    hSelectedEvents->Sumw2();

    // First bin of SelectedEvents is number of generated events
    hSelectedEvents->SetBinContent(1,
                        hTotalEvents->GetBinContent(1) - 2 * hTotalEvents->GetBinContent(10));
/*
    // Systematics
    TH2 *hSystematics;
    if (smearOn)
    {
        TString histName = "../data/" + systematics + "_smear_" + YEAR_STR + ".root";
        TFile *histFile = TFile::Open(histName);

        cout << "Opened " << histName << endl;

        histFile->GetObject("SMEAR" + idH, hSystematics);
        hSystematics->SetDirectory(outFile);
    }
*/





    ////
    ////
    ////    EVENT LOOP
    ////
    ////


    int nEvents = reader.GetEntries(kTRUE);

    cout << endl;
    cout << "Running over " << nEvents << " total events" << endl;
    if (debug)  cout << "Will select " << selectEvents << " events" << endl;
    cout << endl;

    long count = 0;     // track total number of selected events

    while (reader.Next())
    {
        
        //
        //  DEBUG PRINTOUT
        //

        const long currentEntry = reader.GetCurrentEntry();
        bool print = kFALSE;

        if (debug)
        {
            if (count >= selectEvents)
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

        // Quantities copied directly to output tree
        runNum      = *runNum_;         evtNum      = *evtNum_;         lumiSec     = *lumiSec_;
        genWeight   = *genWeight_;      qtWeight    = 1;                ecalWeight  = *ecalWeight_;
        puWeight    = *puWeight_;       nPV         = *nPV_;

        // Quantities used in analysis, but not written out
        bool        muonTrig    = *muonTrig_,       elecTrig    = *elecTrig_;
        unsigned    nMuons      = *nMuons_,         nElecs      = *nElecs_;
        unsigned    nTightMuons = *nTightMuons_,    nTightElecs = *nTightElecs_;



        //
        //  PRESELECTION
        //

        // Make sure data events aren't triggered twice (hmm...should be done at BLT level?)

        if (isData && suffix.Contains("muon"))
            elecTrig = kFALSE;
        if (isData && suffix.Contains("electron"))
            muonTrig = kFALSE;


        // Sanity checks (these are both implemented at BLT level)

        if (!muonTrig && !elecTrig)         // event does not pass trigger
            continue;

        if (nTightMuons < 2 && nTightElecs < 2) // not enough HZZ-identified leptons
            continue;


        hTotalEvents->Fill(6);

        if (print)
            cout << "Passed preselection" << endl;




        //
        //  SYSTEMATICS
        //
/*
        // Get smudge factors from histograms and apply them before creating objects
        if (smearOn)
        {
            if      (smearMuonID)
            {
                for (unsigned i = 0; i < nMuons; i++)
                {
                    TLorentzVector p4 = muonP4_.At(i);
                    int bin = hSystematics->FindBin(p4.Eta(), p4.Pt());
                    (*muonIDSF_)[i] += hSystematics->GetBinContent(bin);
                }
            }
            else if (smearElecID || smearElecReco) 
            {
                for (unsigned i = 0; i < nElecs; i++)
                {
                    TLorentzVector p4 = elecP4_.At(i);
                    int bin = hSystematics->FindBin(p4.Eta(), p4.Pt());
                    (*elecIDSF_)[i] += hSystematics->GetBinContent(bin);
                }
            }
            else if (smearMuonPt)
            {
                inTree->GetEntry(currentEntry);

                for (unsigned i = 0; i < nMuons; i++)
                {
                    TLorentzVector p4 = muonP4_.At(i);
                    int bin = hSystematics->FindBin(p4.Eta(), p4.Pt());
                    (*muonEnergySF_)[i] *= 1 + hSystematics->GetBinContent(bin);
                }
            }
            else if (smearElecPt)
            {
                for (unsigned i = 0; i < nElecs; i++)
                    (*elecEnergySF_)[i] *= 1 + ELEC_PT_SHIFT;
            }
        }
*/




        ////
        ////
        ////    LEPTONS
        ////
        ////


        // Muons

        vector<Lepton> muons;

        for (unsigned i = 0; i < nMuons; i++)
        {
            if (!(*muonIsTight_)[i])                // failed ID *after* Rochester correction
                continue;


            // Build and store lepton object

            Lepton  muon;

            muon.p4         = muonP4_.At(i);
            muon.u_p4       = muonUncorrP4_.At(i);
            muon.q          = (*muonQ_)[i];
            muon.pdg        = -13 * muon.q;
            muon.iso        = (*muonIso_)[i] / muon.p4.Pt();    // actually rel iso
            muon.id_sf      = make_pair((*muonIDSF_)[i],        1);
            muon.fired      = make_pair((*muonFiredLeg1_)[i],   (*muonFiredLeg2_)[i]);
            muon.te_data    = make_pair((*muonEffL1Data_)[i],   (*muonEffL2Data_)[i]);
            muon.te_mc      = make_pair((*muonEffL1MC_)[i],     (*muonEffL2MC_)[i]);

            muons.push_back(muon);
        }


        // Electrons

        vector<Lepton> elecs;

        for (unsigned i = 0; i < nElecs; i++)
        {
            if (!(*elecIsTight_)[i])                // failed ID *after* energy correction
                continue;


            // Build and store lepton object

            Lepton  elec;

            elec.p4         = elecP4_.At(i);
            elec.u_p4       = elecUncorrP4_.At(i);
            elec.q          = (*elecQ_)[i];
            elec.pdg        = -11 * elec.q;
            elec.iso        = (*elecIso_)[i] / elec.p4.Pt();    // actually rel iso
            elec.id_sf      = make_pair((*elecIDSF_)[i],        (*elecRecoSF_)[i]);
            elec.fired      = make_pair((*elecFiredLeg1_)[i],   (*elecFiredLeg2_)[i]);
            elec.te_data    = make_pair((*elecEffL1Data_)[i],   (*elecEffL2Data_)[i]);
            elec.te_mc      = make_pair((*elecEffL1MC_)[i],     (*elecEffL2MC_)[i]);

            elecs.push_back(elec);
        }



        //
        //  LEPTON COUNT
        //

        sort(muons.begin(), muons.end(), DecreasingPt);
        nTightMuons = muons.size();
        sort(elecs.begin(), elecs.end(), DecreasingPt);
        nTightElecs = elecs.size();
        vector<Lepton> tmp_leps = muons;
        tmp_leps.insert(tmp_leps.end(), elecs.begin(), elecs.end());
        unsigned nTightLeps = tmp_leps.size();

        if (nTightLeps != 2 && nTightLeps != 4) // wrong number of leptons
            continue;
        

        hTotalEvents->Fill(7);

        if (print)
            cout << "Passed lepton count requirement (" << nTightLeps << " leps)" << endl;



        //
        //  MASS WINDOW
        //

        TLorentzVector tmp_p4 = TotalP4(tmp_leps);

        if ((tmp_p4.M() < M_MIN) || (tmp_p4.M() > M_MAX))   // total mass below/above Z window
            continue;

        if (print)
            cout << "Passed total mass cut (" << tmp_p4.M() << " GeV)" << endl;



        //
        //  DIVERGENCE CUT
        //
        //  &
        //
        //  GHOST REMOVAL
        //

        bool failedDivCut = kFALSE, ghostBusted = kFALSE;

        // "exitloops" at the end of the following three outer loops

        // Muons
        for (unsigned j = 1; j < nTightMuons; j++)
        {
            for (unsigned i = 0; i < j; i++)
            {
                TLorentzVector dilep_p4 = muons[i].p4 + muons[j].p4;

                if ((muons[i].q != muons[j].q) && (dilep_p4.M() < MLL_MIN)) // failed divergence
                {
                    failedDivCut = kTRUE;
                    goto exitloops;
                }
                if (muons[i].p4.DeltaR(muons[j].p4) < SF_DR_MIN)    // failed same-flavor DeltaR
                {
                    ghostBusted = kTRUE;
                    goto exitloops;
                }
            }
        }


        // Electrons
        for (unsigned j = 1; j < nTightElecs; j++)
        {
            for (unsigned i = 0; i < j; i++)
            {
                TLorentzVector dilep_p4 = elecs[i].p4 + elecs[j].p4;

                if ((elecs[i].q != elecs[j].q) && (dilep_p4.M() < MLL_MIN)) // failed divergence
                {
                    failedDivCut = kTRUE;
                    goto exitloops;
                }
                if (elecs[i].p4.DeltaR(elecs[j].p4) < SF_DR_MIN)    // failed same-flavor DeltaR
                {
                    ghostBusted = kTRUE;
                    goto exitloops;
                }
            }
        }


        // Apply (larger) minimum DeltaR cut to opposite-flavor leptons
        for (unsigned i = 0; i < nTightElecs; i++)
        {
            for (unsigned j = 0; j < nTightMuons; j++)
            {
                if (elecs[i].p4.DeltaR(muons[j].p4) < OF_DR_MIN)    // failed opp-flavor DeltaR
                {
                    ghostBusted = kTRUE;
                    goto exitloops;
                }
            }
        }

        exitloops:


        if (failedDivCut)
            continue;

        if (print)
            cout << "Passed divergence cut" << endl;


        if (ghostBusted)
            continue;

        if (print)
            cout << "Passed DeltaR requirements" << endl;


        hTotalEvents->Fill(8);



        //
        //  CATEGORIZE EVENT
        //

        unsigned C;     // index for filling trees by channel
        LeptonPair z1, z2;

        if (elecTrig)                       // higher-priority electron trigger (reevaluate?)
        {
            if      (nTightElecs == 2 && nTightMuons == 0)      // ee
            {
                C = EE;

                z1.SetMembers(elecs[0], elecs[1]);
            }
            else if (nTightElecs == 2 && nTightMuons == 2)      // 2e2m
            {
                C = EM;

                // Choose z1 as triggered pair (reevaluate?)
                z1.SetMembers(elecs[0], elecs[1]); 
                z2.SetMembers(muons[0], muons[1]);
            }
            else if (nTightElecs == 4 && nTightMuons == 0)      // 4e
            {
                C = E4;

                // Choose (opposite-sign) pairs such that their mass difference is maximized
                MakePairsMaxDiff(elecs, &z1, &z2);
            } 
            else
                continue;
        }

        else if (muonTrig)                  // lower-priority muon trigger
        {
            if      (nTightMuons == 2 && nTightElecs == 0)      // mumu
            {
                C = MM;

                z1.SetMembers(muons[0], muons[1]);
            }
            else if (nTightMuons == 2 && nTightElecs == 2)      // 2m2e
            {
                C = ME;

                // Choose z1 as triggered pair (reevaluate?)
                z1.SetMembers(muons[0], muons[1]);
                z2.SetMembers(elecs[0], elecs[1]);
            }
            else if (nTightMuons == 4 && nTightElecs == 0)      // 4m
            {
                C = M4;

                // Choose (opposite-sign) pairs such that their mass difference is maximized
                MakePairsMaxDiff(muons, &z1, &z2);
            }
            else
                continue;
        }


        if (print)
            cout << "Passed pair requirement for " << selection[C] << endl;



        //
        //  SET CUTS
        //

        // false => elecPairLeads
        const bool  muonPairLeads   = ((C == MM) || (C == ME) || (C == M4));

        const float LEP_PT1_MIN     = muonPairLeads ? MUON_PT1_MIN  : ELEC_PT1_MIN;
        const float LEP_PT2_MIN     = muonPairLeads ? MUON_PT2_MIN  : ELEC_PT2_MIN;


        const bool  isDilepton      = ((C == MM) || (C == EE));
        const bool  isFourLepton    = ((C == M4) || (C == ME) || (C == EM) || (C == E4));
        const bool  isMixedFlavor   = ((C == ME) || (C == EM));






        ////
        ////
        ////    DILEPTON SELECTION
        ////
        ////

        if (isDilepton)
        {

            //
            //  CHARGE REQUIREMENT
            //

            if ((z1.Plus().q != 1) || (z1.Minus().q != -1)) // lepton charges are not +1 and -1
                continue;

            if (print)
                cout << "Passed charge requirement" << endl;



            //
            //  PT REQUIREMENTS
            //

            if (z1.First().p4.Pt() < LEP_PT1_MIN)           // leading lepton Pt too low
                continue;

            if (z1.Second().p4.Pt() < LEP_PT2_MIN)          // trailing lepton Pt too low
                continue;

            if (print)
                cout << "Passed Pt requirements" << endl;



            //
            //  EVENT WEIGHT
            //

            idWeight    = z1.First().id_sf.first * z1.Second().id_sf.first;
            recoWeight  = z1.First().id_sf.second * z1.Second().id_sf.second;
            trigWeight  = 1;    // FIXME


            hTotalEvents->Fill(9);

            if (print)
                cout << "PASSED DILEPTON SELECTION" << endl;

        } // END dilepton case






        ////
        ////
        ////    FOUR-LEPTON SELECTION
        ////
        ////

        else if (isFourLepton)
        {

            //
            //  CHARGE REQUIREMENT
            //

            if ((z1.Plus().q != 1) || (z1.Minus().q != -1)) // lepton charges are not +1 and -1 
                continue;

            if ((z2.Plus().q != 1) || (z2.Minus().q != -1)) // (both extraneous for same-flavor)
                continue;

            if (print)
                cout << "Passed charge requirements" << endl;



            //
            //  PT REQUIREMENTS
            //

            if (isMixedFlavor)  // z1 leptons are triggered
            {
                if (z1.First().p4.Pt() < LEP_PT1_MIN)
                    continue;

                if (z1.Second().p4.Pt() < LEP_PT2_MIN)
                    continue;
            }
            else                // triggered leptons could be in either pair
            {
                vector<Lepton> all_leps = z1.GetMembers();
                vector<Lepton> z2_leps = z2.GetMembers();
                all_leps.insert(all_leps.end(), z2_leps.begin(), z2_leps.end());
                sort(all_leps.begin(), all_leps.end(), DecreasingPt);

                if (all_leps[0].p4.Pt() < LEP_PT1_MIN)  // no lepton passes Pt1 threshold
                    continue;

                if (all_leps[1].p4.Pt() < LEP_PT2_MIN)  // no lepton passes Pt2 threshold
                    continue;
            }

            if (print)
                cout << "Passed Pt requirement" << endl;



            //
            //  Z1, Z2 MASS REQUIREMENTS
            //

            if (z1.p4.M() < z2.p4.M())                          // mixed-flavor pairs are swapped
                swap(z1, z2);

            if ((z1.p4.M() < Z1_M_MIN) || (z1.p4.M() > Z_M_MAX))// z1 failed pair mass requirements
                continue;                                       // (z2's mass is bound by z1)

            z1.SetMothers(1);       z2.SetMothers(2);           // bookkeeping

            if (print)
                cout << "Passed z1, z2 mass requirements" << endl;



            //
            //  EVENT WEIGHT
            //

            idWeight    = z1.First().id_sf.first * z1.Second().id_sf.first;
            idWeight   *= z2.First().id_sf.first * z2.Second().id_sf.first;
            recoWeight  = z1.First().id_sf.second * z1.Second().id_sf.second;
            recoWeight *= z2.First().id_sf.second * z2.Second().id_sf.second;
            trigWeight  = 1;    // FIXME


            hTotalEvents->Fill(9);

            if (print)
                cout << "PASSED FOUR-LEPTON SELECTION" << endl;

        } // END four-lepton case






        ////
        ////
        ////    FILL
        ////
        ////


        // Event has been selected!


        // Assemble leptons
        vector<Lepton> leps = z1.GetMembers();

        if (isFourLepton)
        {
            vector<Lepton> z2_leps = z2.GetMembers();
            leps.insert(leps.end(), z2_leps.begin(), z2_leps.end());
        }

        sort(leps.begin(), leps.end(), DecreasingPt);


        // Set branch variables
        weight  = genWeight * qtWeight * ecalWeight * puWeight * trigWeight * idWeight * recoWeight;
        channel = chanIdx[C];

        z1p4    = z1.p4;                z1pdg   = z1.pdg;
        l1p4    = leps[0].p4;           l1pdg   = leps[0].pdg;          l1iso   = leps[0].iso;
        l2p4    = leps[1].p4;           l2pdg   = leps[1].pdg;          l2iso   = leps[1].iso;

        if (isFourLepton)
        { 
            z2p4    = z2.p4;                z2pdg   = z2.pdg;
            zzp4    = z1p4 + z2p4;
            l1z     = leps[0].mother;       l2z     = leps[1].mother;
            l3p4    = leps[2].p4;           l3pdg   = leps[2].pdg;
            l3iso   = leps[2].iso;          l3z     = leps[2].mother;
            l4p4    = leps[3].p4;           l4pdg   = leps[3].pdg;
            l4iso   = leps[3].iso;          l4z     = leps[3].mother;
        }

        // Get gen particle info
        if (isFourLepton && isSignal && YEAR_STR.EqualTo("2017") && !smearOn)
        {
            u_l1p4  = leps[0].u_p4;         u_l2p4  = leps[1].u_p4;
            u_l3p4  = leps[2].u_p4;         u_l4p4  = leps[3].u_p4;

            hardProcMuonP4->Delete();       finalStateMuonP4->Delete();
            hardProcElectronP4->Delete();   finalStateElectronP4->Delete();

            inTree->GetEntry(currentEntry);

            if (print)
                cout << nHardProcLeptons  << " hard process leptons" << endl;
        }

 
        // Histograms
        unsigned D = (C > L4) ? L4 : LL;    // Index for ll, 4l channel

        hSelectedEvents->Fill(chanIdx[C], weight);
        hSelectedEvents->Fill(chanIdx[D], weight);


        // Trees
        tree[C]->Fill();

        if (C > L4)                         // don't write ll tree because it's HUGE and useless
            tree[D]->Fill();


        count++;

    } // END event loop



    //
    //  PRINT RESULTS
    //

    int nSelected = tree[LL]->GetEntries() + tree[L4]->GetEntries(); 

    cout << endl << endl;
    cout << "Done!" << endl;
    cout << endl << endl;
    cout << "Yields (unweighted):" << endl << endl;
    for (unsigned i = 0; i < N; i++)
    {
        if (i != LL && i != L4)
            cout << "\t\t";
        if (i == L4)
            cout << endl;

        cout << selection[i] << ":\t";
        cout << tree[i]->GetEntries() << endl;
    }
    cout << endl << endl;



    //
    //  WRITE FILE
    //

    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();

    hTotalEvents->Write();
    hSelectedEvents->Write();
/*
    if (smearOn)
        hSystematics->Write();
*/
    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
