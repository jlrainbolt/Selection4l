// STL
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"
#include "SelectionTools.hh"

// Cuts
//#include "Cuts2018.hh"
//#include "Cuts2017.hh"
#include "Cuts2016.hh"
//#include "Cuts2012.hh"

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
**  Integrates systematics analysis using extra branches.
*/

void RecoSelection( const TString suffix,           const TString id,
                    const TString systematics = "", const TString var = "")
{


    //
    //  OPTIONS
    //

    bool debug = kFALSE;
    long selectEvents = 10000;
    int  printEvery = 30000;

    // Systematics toggle
    const bool smearMuonPt = systematics.EqualTo("muonPt");
    const bool smearElecPt = systematics.EqualTo("electronPt");
    const bool allowUntriggered = systematics.EqualTo("triggerEff");

    const bool systOn = allowUntriggered || smearMuonPt || smearElecPt;

    if ((smearMuonPt || smearElecPt) && !var.EqualTo("Up") && !var.EqualTo("Down"))
    {
        cout << "Incorrect value for \"var\"" << endl;
        return;
    }



    //
    //  SAMPLE INFO
    //

    const bool isData       = suffix.Contains(YEAR_STR);

    const bool isSignal     = suffix.EqualTo("zz_4l");
    const bool isDrellYan   = suffix.EqualTo("zjets_m-50");

    const unsigned N = 6;   // Channel indices
    unsigned                   MM = 0, EE = 1, L4 = 2, M4 = 3, ME = 4, E4 = 5;
    TString selection[N]    = {"mumu", "ee",   "4l",   "4m",   "2m2e", "4e"};
    unsigned chanIdx[N]     = {3,      4,      5,      6,      7,      9};

    const TString UncorrP4 = YEAR_STR.EqualTo("2012") ? "UncorrP4" : "UncorrectedP4";





    ////
    ////
    ////    OUTPUT
    ////
    ////


    //
    //  FILE
    //

    TString prefix  = systOn ? systematics + var : "selected";
    TString outName = prefix + "_" + suffix + "_" + id + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);



    //
    //  BRANCHES
    //

    // Event info
    Int_t               runNum,     evtNum,     lumiSec;
    UShort_t            nPV;
    Float_t             weight,     genWeight,      qtWeight,       trigWeight,     idWeight;
    Float_t             ecalWeight, ecalWeightUp,   ecalWeightDown;
    Float_t             puWeight,   puWeightUp,     puWeightDown,   nPU;
    UShort_t            channel;
    Bool_t              hasTauDecay;
    Bool_t              muonTrig,   siMuTrig,   diMuTrig,   elecTrig,   siElTrig,   diElTrig;

    // Pairs
    TLorentzVector      z1p4,       z2p4,       zzp4;
    UShort_t            z1pdg,      z2pdg;

    // Leptons
    TLorentzVector      l1p4,       l2p4,       l3p4,       l4p4;
    Short_t             l1pdg,      l2pdg,      l3pdg,      l4pdg;
    Float_t             l1iso,      l2iso,      l3iso,      l4iso;
    UShort_t            l1z,        l2z,        l3z,        l4z;


    // Gen-level info (only written for signal)
    UShort_t            nDressedMuons,          nDressedElectrons,      nDressedLeptons;
    TClonesArray        *dressedMuonP4          = new TClonesArray("TLorentzVector");
    TClonesArray        *dressedElectronP4      = new TClonesArray("TLorentzVector");
    TLorentzVector      *dressedLeptonsP4 = 0;

    vector<Short_t>     *dressedMuonQ = 0,          *dressedElectronQ = 0;
    vector<UShort_t>    *dressedMuonZIndex = 0,     *dressedElectronZIndex = 0;

    TLorentzVector      u_l1p4,     u_l2p4,     u_l3p4,     u_l4p4;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("runNum",       &runNum);       tree[i]->Branch("evtNum",       &evtNum);               
        tree[i]->Branch("lumiSec",      &lumiSec);      tree[i]->Branch("nPV",          &nPV);
        tree[i]->Branch("weight",       &weight);       tree[i]->Branch("genWeight",    &genWeight);
        tree[i]->Branch("qtWeight",     &qtWeight);     tree[i]->Branch("puWeight",     &puWeight);
        tree[i]->Branch("ecalWeight",   &ecalWeight);   tree[i]->Branch("trigWeight",   &trigWeight);
        tree[i]->Branch("idWeight",     &idWeight);
        tree[i]->Branch("channel",      &channel);      tree[i]->Branch("hasTauDecay",  &hasTauDecay);

        if (isData || isSignal || isDrellYan)
        {
            tree[i]->Branch("muonTrig",     &muonTrig); tree[i]->Branch("doubleMuTrig", &diMuTrig);
            tree[i]->Branch("singleMuTrig", &siMuTrig); tree[i]->Branch("electronTrig", &elecTrig);
            tree[i]->Branch("doubleElTrig", &diElTrig); tree[i]->Branch("singleElTrig", &siElTrig);
        }

        if (isSignal || isDrellYan)
        {
            tree[i]->Branch("nPU",              &nPU);

            if      (YEAR_STR.EqualTo("2012"))
            {
                tree[i]->Branch("puWeightUp",       &puWeightUp);
                tree[i]->Branch("puWeightDown",     &puWeightDown);
            }
            else if (YEAR_STR.EqualTo("2016") || YEAR_STR.EqualTo("2017"))
            {
                tree[i]->Branch("ecalWeightUp",     &ecalWeightUp);
                tree[i]->Branch("ecalWeightDown",   &ecalWeightDown);
            }
        }

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
            tree[i]->Branch("nDressedMuons",            &nDressedMuons);
            tree[i]->Branch("nDressedElectrons",        &nDressedElectrons);
            tree[i]->Branch("nDressedLeptons",          &nDressedLeptons);

            tree[i]->Branch("dressedMuonP4",            &dressedMuonP4,         32000,      1);
            tree[i]->Branch("dressedMuonQ",             &dressedMuonQ);
            tree[i]->Branch("dressedMuonZIndex",        &dressedMuonZIndex);
            tree[i]->Branch("dressedElectronP4",        &dressedElectronP4,     32000,      1);
            tree[i]->Branch("dressedElectronQ",         &dressedElectronQ);
            tree[i]->Branch("dressedElectronZIndex",    &dressedElectronZIndex);
            tree[i]->Branch("dressedLeptonsP4",         &dressedLeptonsP4);

            tree[i]->Branch("uncorr_l1p4",  &u_l1p4);   tree[i]->Branch("uncorr_l2p4",  &u_l2p4);
            tree[i]->Branch("uncorr_l3p4",  &u_l3p4);   tree[i]->Branch("uncorr_l4p4",  &u_l4p4);
        }
    }





    ////
    ////
    ////    INPUT
    ////
    ////


    //
    //  FILE
    //

    TString inName  = suffix + "_" + id + ".root";
    TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "_new/";
    TFile   *inFile = TFile::Open(inPath + inName);

    cout << endl << endl << "Opened " << inPath + inName << endl;



    //
    //  BRANCHES
    //

    TTreeReader reader("tree_" + suffix, inFile);

    TTreeReaderValue    <Int_t>             runNum_         (reader,    "runNumber");
    TTreeReaderValue    <Int_t>             evtNum_         (reader,    "evtNumber");
    TTreeReaderValue    <Int_t>             lumiSec_        (reader,    "lumiSection");
    TTreeReaderValue    <Float_t>           genWeight_      (reader,    "genWeight");
    TTreeReaderValue    <Float_t>           ecalWeight_     (reader,    "ECALWeight");
    TTreeReaderValue    <Float_t>           puWeight_       (reader,    "PUWeight");
    TTreeReaderValue    <UShort_t>          nPV_            (reader,    "nPV");
    TTreeReaderValue    <Bool_t>            hasTauDecay_    (reader,    "hasTauDecay");

    TTreeReaderValue    <Bool_t>            muonTrig_       (reader,    "evtMuonTriggered");
    TTreeReaderValue    <Bool_t>            diMuTrig_       (reader,    "evtDoubleMuTriggered");
    TTreeReaderValue    <Bool_t>            siMuTrig_       (reader,    "evtSingleMuTriggered");
    TTreeReaderValue    <Bool_t>            elecTrig_       (reader,    "evtElectronTriggered");
    TTreeReaderValue    <Bool_t>            diElTrig_       (reader,    "evtDoubleElTriggered");
    TTreeReaderValue    <Bool_t>            siElTrig_       (reader,    "evtSingleElTriggered");

    TTreeReaderValue    <UShort_t>          nMuons_         (reader,    "nLooseMuons");
    TTreeReaderValue    <UShort_t>          nElecs_         (reader,    "nLooseElectrons");
    TTreeReaderValue    <UShort_t>          nTightMuons_    (reader,    "nTightMuons");
    TTreeReaderValue    <UShort_t>          nTightElecs_    (reader,    "nTightElectrons");

    TTreeReaderArray    <TLorentzVector>    muonP4_         (reader,    "muonP4");
    TTreeReaderArray    <TLorentzVector>    muonUncorrP4_   (reader,    "muon" + UncorrP4);
    TTreeReaderValue    <vector<Float_t>>   muonEnergySF_   (reader,    "muonEnergySF" + var);
    TTreeReaderValue    <vector<Short_t>>   muonQ_          (reader,    "muonQ");
    TTreeReaderValue    <vector<Float_t>>   muonIso_        (reader,    "muonIsolation");
    TTreeReaderValue    <vector<Bool_t>>    muonIsTight_    (reader,    "muonIsTight");
    TTreeReaderValue    <vector<Float_t>>   muonIDSF_       (reader,    "muonIDSF");
    TTreeReaderValue    <vector<Bool_t>>    muonFiredLeg1_  (reader,    "muonFiredLeg1");
    TTreeReaderValue    <vector<Bool_t>>    muonFiredLeg2_  (reader,    "muonFiredLeg2");
    TTreeReaderValue    <vector<Bool_t>>    muonFiredSing_  (reader,    "muonFiredSingle");

    TTreeReaderArray    <TLorentzVector>    elecP4_         (reader,    "electronP4");
    TTreeReaderArray    <TLorentzVector>    elecUncorrP4_   (reader,    "electron" + UncorrP4);
    TTreeReaderValue    <vector<Float_t>>   elecEnergySF_   (reader,    "electronEnergySF" + var);
    TTreeReaderValue    <vector<Short_t>>   elecQ_          (reader,    "electronQ");
    TTreeReaderValue    <vector<Float_t>>   elecIso_        (reader,    "electronIsolation");
    TTreeReaderValue    <vector<Bool_t>>    elecIsTight_    (reader,    "electronIsTight");
    TTreeReaderValue    <vector<Float_t>>   elecIDSF_       (reader,    "electronIDSF");
    TTreeReaderValue    <vector<Float_t>>   elecRecoSF_     (reader,    "electronRecoSF");
    TTreeReaderValue    <vector<Bool_t>>    elecFiredLeg1_  (reader,    "electronFiredLeg1");
    TTreeReaderValue    <vector<Bool_t>>    elecFiredLeg2_  (reader,    "electronFiredLeg2");
    TTreeReaderValue    <vector<Bool_t>>    elecFiredSing_  (reader,    "electronFiredSingle");

    cout << "Loaded branches" << endl;


    // Gen particles

    // Only read for signal events, so we have to do it the old-fashioned way

    TTree *inTree;

    if (isSignal || isDrellYan)
    {
        inFile->GetObject("tree_" + suffix, inTree);

        inTree->SetBranchAddress(   "nPU",  &nPU);

        if      (YEAR_STR.EqualTo("2012"))
        {
            inTree->SetBranchAddress(   "PUWeightUp",       &puWeightUp);
            inTree->SetBranchAddress(   "PUWeightDown",     &puWeightDown);
        }
        else if (YEAR_STR.EqualTo("2016") || YEAR_STR.EqualTo("2017"))
        {
            inTree->SetBranchAddress(   "ECALWeightUp",     &ecalWeightUp);
            inTree->SetBranchAddress(   "ECALWeightDown",   &ecalWeightDown);
        }
    }

    if (isSignal && !systOn)
    {
        inTree->SetBranchAddress(   "nDressedMuons",            &nDressedMuons);
        inTree->SetBranchAddress(   "nDressedElectrons",        &nDressedElectrons);
        inTree->SetBranchAddress(   "nDressedLeptons",          &nDressedLeptons);
        inTree->SetBranchAddress(   "dressedLeptonsP4",         &dressedLeptonsP4);

        inTree->GetBranch(          "dressedMuonP4")    ->SetAutoDelete(kFALSE);
        inTree->SetBranchAddress(   "dressedMuonP4",            &dressedMuonP4);
        inTree->SetBranchAddress(   "dressedMuonQ",             &dressedMuonQ);
        inTree->SetBranchAddress(   "dressedMuonZIndex",        &dressedMuonZIndex);

        inTree->GetBranch(          "dressedElectronP4")->SetAutoDelete(kFALSE);
        inTree->SetBranchAddress(   "dressedElectronP4",        &dressedElectronP4);
        inTree->SetBranchAddress(   "dressedElectronQ",         &dressedElectronQ);
        inTree->SetBranchAddress(   "dressedElectronZIndex",    &dressedElectronZIndex);

        cout << "Loaded gen branches" << endl;
    }



    //
    //  DATA
    //

    // Dilepton Qt reweighting
    TString graphName = "../data/qt_weights_" + YEAR_STR + ".root";

    TGraphAsymmErrors *qtGraph[2];

    TFile *graphFile = TFile::Open(graphName);
    graphFile->GetObject(selection[EE] + "_weight", qtGraph[0]);    // ee: muonPairLeads = 0
    graphFile->GetObject(selection[MM] + "_weight", qtGraph[1]);    // mumu: muonPairLeads = 1

    graphFile->Close();


    // Trigger scale factors (only single muon for now)
    TString histName[2];
    TH2 *sfHist[2];

    if (allowUntriggered)
    {
        for (unsigned i = 0; i < N_TRIG; i++)
        {
            histName[i] = "../data/muon_trig_sf_" + YEAR_STR + "_" + TRIG_PD[i] + ".root";

            TFile *histFile = TFile::Open(histName[i]);
            histFile->GetObject(TRIG_NAME + "_PtEtaBins/abseta_pt_ratio", sfHist[i]);
            sfHist[i]->SetDirectory(0);

            histFile->Close();
        }
    }

    // RNG to choose file
    TRandom3 rng;



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

        // Copied directly to output tree
        runNum      = *runNum_;         evtNum      = *evtNum_;         lumiSec     = *lumiSec_;
        genWeight   = *genWeight_;      ecalWeight  = *ecalWeight_;     puWeight    = *puWeight_;
        nPV         = *nPV_;            hasTauDecay = *hasTauDecay_;
        muonTrig    = *muonTrig_;       diMuTrig    = *diMuTrig_;       siMuTrig    = *siMuTrig_;
        elecTrig    = *elecTrig_;       diElTrig    = *diElTrig_;       siElTrig    = *siElTrig_;

        // Needing initialization
        trigWeight  = 1;                qtWeight    = 1;                idWeight    = 1;

        // Used in analysis, but not written out
        unsigned    nMuons      = *nMuons_,             nElecs      = *nElecs_;
        unsigned    nTightMuons = *nTightMuons_,        nTightElecs = *nTightElecs_;



        //
        //  PRESELECTION
        //

        if (!muonTrig && !elecTrig && !allowUntriggered)    // event does not pass trigger
            continue;

        if (nTightMuons < 2 && nTightElecs < 2)             // not enough tight leptons
            continue;


        hTotalEvents->Fill(6);

        if (print)
            cout << "Passed preselection" << endl;

        // Throw a random number to decide which trigger SF file to pull from
        const unsigned TI = rng.Rndm() < TRIG_FRAC ? 0 : 1;
        TH2 *hist = sfHist[TI];






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
            muon.di_hlt     = make_pair((*muonFiredLeg1_)[i],   (*muonFiredLeg2_)[i]);
            muon.si_hlt     = (*muonFiredSing_)[i];


            // Correct energy (if necessary)
            if (smearMuonPt)
            {
                float corrPt = muon.u_p4.Pt() * (*muonEnergySF_)[i];
                muon.p4.SetPtEtaPhiM(corrPt, muon.u_p4.Eta(), muon.u_p4.Phi(), muon.u_p4.M());
            }
            
            // Get trigger scale factor
            if (!isData and !YEAR_STR.EqualTo("2012") and allowUntriggered)
            {
                int bin = hist->FindBin(fabs(muon.p4.Eta()), muon.p4.Pt());
                if (hist->IsBinUnderflow(bin) || hist->IsBinOverflow(bin))
                    muon.trig_sf = 1;
                else
                    muon.trig_sf = hist->GetBinContent(bin);
            }
            else
                muon.trig_sf = 1;

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
            elec.di_hlt     = make_pair((*elecFiredLeg1_)[i],   (*elecFiredLeg2_)[i]);
            elec.si_hlt     = (*elecFiredSing_)[i];


            // Correct energy (if necessary)
            if (smearElecPt)
                elec.p4 = elec.u_p4 * (*elecEnergySF_)[i];

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
        //  TRIGGER MATCHING
        //

        if (muonTrig)
        {
            bool matchedLeg1 = kFALSE, matchedLeg2 = kFALSE, matchedSingle = kFALSE;

            for (unsigned i = 0; i < nTightMuons; i++)
            {
                if (muons[i].di_hlt.first && (muons[i].p4.Pt() > MUON_LEG1_PT))
                    matchedLeg1 = kTRUE;
                if (muons[i].di_hlt.second && (muons[i].p4.Pt() > MUON_LEG2_PT))
                    matchedLeg2 = kTRUE;
                if (muons[i].si_hlt && (muons[i].p4.Pt() > MUON_SINGLE_PT))
                    matchedSingle = kTRUE;
            }

            diMuTrig = diMuTrig && matchedLeg1 && matchedLeg2;
            siMuTrig = siMuTrig && matchedSingle;
        }

        if (elecTrig)
        {
            bool matchedLeg1 = kFALSE, matchedLeg2 = kFALSE, matchedSingle = kFALSE;

            for (unsigned i = 0; i < nTightElecs; i++)
            {
                if (elecs[i].di_hlt.first && (elecs[i].p4.Pt() > ELEC_LEG1_PT))
                    matchedLeg1 = kTRUE;
                if (elecs[i].di_hlt.second && (elecs[i].p4.Pt() > ELEC_LEG2_PT))
                    matchedLeg2 = kTRUE;
                if (elecs[i].si_hlt && (elecs[i].p4.Pt() > ELEC_SINGLE_PT))
                    matchedSingle = kTRUE;
            }

            diElTrig = diElTrig && matchedLeg1 && matchedLeg2;
            siElTrig = siElTrig && matchedSingle;
        }

        if (allowUntriggered)
        {
            bool matchedLeg1 = kFALSE, matchedLeg2 = kFALSE, matchedSingle = kFALSE;
            for (unsigned i = 0; i < nTightMuons; i++)
            {
                if (muons[i].p4.Pt() > MUON_LEG1_PT)
                    matchedLeg1 = kTRUE;
                if (muons[i].p4.Pt() > MUON_LEG2_PT)
                    matchedLeg2 = kTRUE;
                if (muons[i].p4.Pt() > MUON_SINGLE_PT)
                    matchedSingle = kTRUE;
            }
            bool matchedMuTrig = (matchedLeg1 && matchedLeg2) || matchedSingle;

            matchedLeg1 = kFALSE; matchedLeg2 = kFALSE; matchedSingle = kFALSE;
            for (unsigned i = 0; i < nTightElecs; i++)
            {
                if (elecs[i].p4.Pt() > ELEC_LEG1_PT)
                    matchedLeg1 = kTRUE;
                if (elecs[i].p4.Pt() > ELEC_LEG2_PT)
                    matchedLeg2 = kTRUE;
                if (elecs[i].p4.Pt() > ELEC_SINGLE_PT)
                    matchedSingle = kTRUE;
            }
            bool matchedElTrig = (matchedLeg1 && matchedLeg2) || matchedSingle;

            if (!matchedMuTrig && !matchedElTrig)
                continue;
        }

        muonTrig = diMuTrig || siMuTrig;
        elecTrig = diElTrig || siElTrig;

        if (!muonTrig && !elecTrig && !allowUntriggered)
            continue;



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
                TLorentzVector dilep_p4 = elecs[i].p4 + muons[j].p4;

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

        if      (nTightMuons == 2 && nTightElecs == 0)      // mumu
        {
            C = MM;
            z1.SetMembers(muons[0], muons[1]);
        }
        else if (nTightMuons == 0 && nTightElecs == 2)      // ee
        {
            C = EE;
            z1.SetMembers(elecs[0], elecs[1]);
        }
        else if (nTightMuons == 4 && nTightElecs == 0)      // 4m
        {
            C = M4;
            MakePairsMaxDiff(muons, &z1, &z2);
        }
        else if (nTightMuons == 2 && nTightElecs == 2)      // 2m2e
        {
            C = ME;
            z1.SetMembers(muons[0], muons[1]);
            z2.SetMembers(elecs[0], elecs[1]);

            if (z1.p4.M() < z2.p4.M())
                swap(z1, z2);
        }
        else if (nTightMuons == 0 && nTightElecs == 4)      // 4e
        {
            C = E4;
            MakePairsMaxDiff(elecs, &z1, &z2);
        }
        else
            continue;


        if (print)
            cout << "Passed pair requirement for " << selection[C] << endl;

        const bool  muonPairLeads   = z1.pdg == 13;
        const bool  isDilepton      = (C == MM) || (C == EE);
        const bool  isFourLepton    = (C == M4) || (C == ME) || (C == E4);

        const float LEP_PT1_MIN     = FID_PT1_MIN;
        const float LEP_PT2_MIN     = FID_PT2_MIN;





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

            unsigned Q  = muonPairLeads;
            if (isDrellYan)
                qtWeight = qtGraph[Q]->Eval(z1.p4.Pt());

            idWeight    *= z1.First().id_sf.first * z1.Second().id_sf.first;
            idWeight    *= z1.First().id_sf.second * z1.Second().id_sf.second;


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
            //  Z1, Z2 MASS REQUIREMENTS
            //

            if ((z1.p4.M() < Z1_M_MIN) || (z1.p4.M() > Z_M_MAX))// z1 failed pair mass requirements
                continue;                                       // (z2's mass is bound by z1)

            z1.SetMothers(1);       z2.SetMothers(2);           // bookkeeping

            if (print)
                cout << "Passed z1, z2 mass requirements" << endl;



            //
            //  PT REQUIREMENTS
            //

            vector<Lepton> all_leps = z1.GetMembers();
            vector<Lepton> z2_leps = z2.GetMembers();
            all_leps.insert(all_leps.end(), z2_leps.begin(), z2_leps.end());
            sort(all_leps.begin(), all_leps.end(), DecreasingPt);

            if (all_leps[0].p4.Pt() < LEP_PT1_MIN)  // no lepton passes Pt1 threshold
                continue;

            if (all_leps[1].p4.Pt() < LEP_PT2_MIN)  // no lepton passes Pt2 threshold
                continue;

            if (print)
                cout << "Passed Pt requirement" << endl;



            //
            //  EVENT WEIGHT
            //

            TLorentzVector zzp4 = z1.p4 + z2.p4;
            unsigned Q  = muonPairLeads;
            if (isDrellYan)
                qtWeight = qtGraph[Q]->Eval(zzp4.Pt());

            idWeight    *= z1.First().id_sf.first * z1.Second().id_sf.first;
            idWeight    *= z2.First().id_sf.first * z2.Second().id_sf.first;
            idWeight    *= z1.First().id_sf.second * z1.Second().id_sf.second;
            idWeight    *= z2.First().id_sf.second * z2.Second().id_sf.second;


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

        // Adjust trigger weights
        if (!isData)
        {
            if      (muonTrig)
            {
                ecalWeight = 1;
                trigWeight = 1;
                for (unsigned i = 0; i < nTightMuons; i++)
                    trigWeight *= muons[i].trig_sf;
            }
            else if (elecTrig)
            {
                // check for events with more than one electron with eta > 2.1
                unsigned n_pref = 0;
                for (unsigned i = 0; i < nTightElecs; i++)
                {
                    if (elecs[i].p4.Eta() > ELEC_ETA_PREF)
                        n_pref++;
                }
                if (n_pref > 1)
                    ecalWeight = 1;

                trigWeight = ELEC_TRIG_SF;
            }
            else // untriggered
                ecalWeight = 1;
        }
        weight = genWeight * puWeight * ecalWeight * trigWeight * qtWeight * idWeight;


        // Assemble leptons
        vector<Lepton> leps = z1.GetMembers();

        if (isFourLepton)
        {
            vector<Lepton> z2_leps = z2.GetMembers();
            leps.insert(leps.end(), z2_leps.begin(), z2_leps.end());
        }

        sort(leps.begin(), leps.end(), DecreasingPt);


        // Set branch variables
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
        if (isFourLepton && isSignal && !systOn)
        {
            u_l1p4  = leps[0].u_p4;         u_l2p4  = leps[1].u_p4;
            u_l3p4  = leps[2].u_p4;         u_l4p4  = leps[3].u_p4;

            dressedMuonP4->Delete();        dressedElectronP4->Delete();

            inTree->GetEntry(currentEntry);

            if (print)
                cout << nDressedLeptons  << " dressed leptons" << endl;
        }

 
        // Fill trees and histograms

        tree[C]->Fill();
        hSelectedEvents->Fill(chanIdx[C], weight);

        if (C > L4)
        {
            tree[L4]->Fill();
            hSelectedEvents->Fill(chanIdx[L4], weight);
        }

        hSelectedEvents->Fill(1, weight);


        count++;

    } // END event loop



    //
    //  PRINT RESULTS
    //

    cout << endl << endl;
    cout << "Done!" << endl;
    cout << endl << endl;
    cout << "Yields (unweighted):" << endl << endl;
    for (unsigned i = 0; i < N; i++)
    {
        if (i == L4)
            cout << endl;
        else
            cout << "\t\t";

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

    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
