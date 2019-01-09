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
#include "TVector3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

// Cuts
#include "Cuts2017.hh"

using namespace std;



/*
**  GenRescaling
**
**
**  All "selected" events are within the phase space region.  Events which also fall into the
**  fiducial region are flagged (isFiducial) based on cuts from header file (CutsXXXX.hh).  The
**  input argument determines whether all or only fiducial events are written out.
**
**  Reads from post-BLT analyzer (PhaseSpaceAnalyzer) ntuples and write out gen-level files that
**  are identical to the RescaledAnalysis output format ("rescaled_").
**
**  Currently only implemented for Drell-Yan (zjets_m-50) events...and only uses hard leptons.
*/

void GenRescaling(const bool fidOnly = kFALSE)
{

    //
    //    OPTIONS
    //

    // fidOnly = kTRUE  => only events which fall into the fiducial region are written out
    //         = kFALSE => all input (phase space region) events are written out


    bool debug = kFALSE;
    int selectEvents = 1000;
    int printEvery = 30000;



    //
    //  SAMPLE INFO
    //

    const unsigned N = 3;
    unsigned                   LL = 0,  MM = 1, EE = 2;     // Indices
    TString selection[N]    = {"ll",    "mumu", "ee"};
    unsigned chanIdx[N]     = {2,       3,      4};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "rescaled";
    TString suffix  = fidOnly ? "fiducial" : "phase_space"; 
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);



    //
    //  OUTPUT BRANCHES
    //

    // Event info
    Int_t               runNum,     evtNum,     lumiSec;
    Float_t             met,        weight,     rescale;
    UInt_t              channel,    nPV;
    Bool_t              isFiducial;


    // Lab frame objects
    TLorentzVector      z1p4,       l1p4,       l2p4;
    Short_t             z1pdg,      l1pdg,      l2pdg;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("runNum",   &runNum);               tree[i]->Branch("evtNum",   &evtNum);
        tree[i]->Branch("lumiSec",  &lumiSec);              tree[i]->Branch("nPV",      &nPV);
        tree[i]->Branch("met",      &met);                  tree[i]->Branch("channel",  &channel);
        tree[i]->Branch("weight",   &weight);               tree[i]->Branch("rescale",  &rescale);
        tree[i]->Branch("isFiducial", &isFiducial);

        tree[i]->Branch("z1p4",     &z1p4);                 tree[i]->Branch("z1pdg",    &z1pdg);
        tree[i]->Branch("l1p4",     &l1p4);                 tree[i]->Branch("l1pdg",    &l1pdg);
        tree[i]->Branch("l2p4",     &l2p4);                 tree[i]->Branch("l2pdg",    &l2pdg);
    }



    //
    //  INPUT FILE
    //

    TString inDir   = "gen_zjets_m-50";
    TString inName  = "genHardProc_zjets_m-50.root";
    TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/" + inDir + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;



    //
    //  WEIGHT UTILS
    //

    TH1 *hRescale[N];

    TString sfName = "../data/z1pt_rescale_" + YEAR_STR + ".root";
    TFile*  sfFile = new TFile(sfName, "OPEN");
    for (unsigned i = 0; i < N; i++)
    {
        sfFile->GetObject(selection[i] + "/hist_z1pt_" + selection[i], hRescale[i]);
        hRescale[i]->SetDirectory(0);
    }



    //
    //  INPUT BRANCHES
    //

    TTreeReader reader("tree_zjets_m-50", inFile);

    TTreeReaderValue    <Int_t>                 runNum_     (reader,    "runNumber");
    TTreeReaderValue    <Int_t>                 evtNum_     (reader,    "evtNumber.eventNumber");
    TTreeReaderValue    <Int_t>                 lumiSec_    (reader,    "lumiSection");
    TTreeReaderValue    <Float_t>               genWeight_  (reader,    "genWeight");
    TTreeReaderValue    <UShort_t>              channel_    (reader,    "decayChannel");
    TTreeReaderValue    <UShort_t>              nMuons_     (reader,    "nHardProcMuons");
    TTreeReaderValue    <UShort_t>              nElecs_     (reader,    "nHardProcElectrons");
    TTreeReaderArray    <TLorentzVector>        muonP4_     (reader,    "hardProcMuonP4");
    TTreeReaderValue    <vector<Short_t>>       muonQ_      (reader,    "hardProcMuonQ");
    TTreeReaderValue    <vector<UShort_t>>      muonZ_      (reader,    "hardProcMuonZIndex");
    TTreeReaderArray    <TLorentzVector>        elecP4_     (reader,    "hardProcElectronP4");
    TTreeReaderValue    <vector<Short_t>>       elecQ_      (reader,    "hardProcElectronQ");
    TTreeReaderValue    <vector<UShort_t>>      elecZ_      (reader,    "hardProcElectronZIndex");
    TTreeReaderValue    <TLorentzVector>        lepsP4_     (reader,    "hardProcLeptonsP4");

    cout << "Loaded branches" << endl;



    //
    //  HISTOGRAMS
    //

    TH1D *hPhaseSpaceEvents, *hFiducialEvents;

    hPhaseSpaceEvents = new TH1D("PhaseSpaceEvents_" + suffix, "PhaseSpaceEvents", 10, 0.5, 10.5);
    hPhaseSpaceEvents->SetDirectory(outFile);
    hPhaseSpaceEvents->Sumw2();

    hFiducialEvents = new TH1D("FiducialEvents_" + suffix, "FiducialEvents", 10, 0.5, 10.5);
    hFiducialEvents->SetDirectory(outFile);
    hFiducialEvents->Sumw2();






    ////
    ////
    ////    EVENT LOOP
    ////
    ////


    int nEvents = reader.GetEntries(kTRUE);

    cout << endl;
    cout << "Running over " << nEvents << " total phase space events" << endl;
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
        runNum  = *runNum_;         evtNum  = *evtNum_;         lumiSec     = *lumiSec_;
        weight  = *genWeight_;      channel = *channel_;

        // Quantities used in analysis, but not written out
        unsigned    nMuons  = *nMuons_,         nElecs  = *nElecs_;



        //
        //  LEPTONS
        //

        vector<Lepton> muons, elecs;

        // Muons
        for (unsigned i = 0; i < nMuons; i++)
        {
            Lepton      muon;

            muon.p4     = muonP4_.At(i);
            muon.q      = (*muonQ_)[i];
            muon.pdg    = -13 * muon.q;

            muons.push_back(muon);
        }
        sort(muons.begin(), muons.end(), DecreasingPt);


        // Electrons
        for (unsigned i = 0; i < nElecs; i++)
        {
            Lepton      elec;

            elec.p4     = elecP4_.At(i);
            elec.q      = (*elecQ_)[i];
            elec.pdg    = -11 * elec.q;

            elecs.push_back(elec);
        }
        sort(elecs.begin(), elecs.end(), DecreasingPt);

        if ((muons.size() != 2) && (elecs.size() != 2))
        {
            cout << "Event does not have two leptons" << endl;
            continue;
        }



        //
        //  CATEGORIZE EVENT
        //

        unsigned C;     // index for filling trees by channel
        LeptonPair z1;
        

        if (print)
            cout << channel << endl;


        if      (channel == 3)                      // mumu
        {
            C = MM;
            z1.SetMembers(muons[0], muons[1]);
        }
        else if (channel == 4)                      // ee
        {
            C = EE;
            z1.SetMembers(elecs[0], elecs[1]);
        }
        else
        {
            cout << "Invalid decay channel" << endl;
            continue;
        }
/*
        if (print)
        {
            cout << "Lab Pt:\t";
            for (unsigned i = 0; i < leps.size(); i++)
                cout << leps[i].p4.Pt() << "\t";
            cout << endl;
        }
*/


        //
        //  RESCALE
        //

        rescale = hRescale[C]->GetBinContent(hRescale[C]->FindBin(z1.p4.Pt()));

        hPhaseSpaceEvents->Fill(chanIdx[C], weight * rescale);
        hPhaseSpaceEvents->Fill(chanIdx[LL], weight * rescale);



        //
        //  FIDUCIAL CHECK
        //

        isFiducial = kTRUE;     // innocent until proven guilty?

        // Pt requirement
        if (z1.First().p4.Pt() < FID_PT1_MIN)
            isFiducial = kFALSE;

        if (z1.Second().p4.Pt() < FID_PT2_MIN)
            isFiducial = kFALSE;


        // Eta requirement
        if (fabs(z1.First().p4.Eta()) > FID_ETA_MAX)
            isFiducial = kFALSE;

        if (fabs(z1.Second().p4.Eta()) > FID_ETA_MAX)
            isFiducial = kFALSE;


        hFiducialEvents->Fill(1, weight * rescale);

        if (isFiducial)
        {
            hFiducialEvents->Fill(chanIdx[C], weight * rescale);
            hFiducialEvents->Fill(chanIdx[LL], weight * rescale);
        }

        if (fidOnly && !isFiducial)
            continue;



        //
        //  FILL TREE
        //

        // Sort all leptons by P in Z CM frame
        z1p4    = z1.p4;                z1pdg   = z1.pdg;
        l1p4    = z1.First().p4;        l1pdg   = z1.First().pdg;
        l2p4    = z1.Second().p4;       l2pdg   = z1.Second().pdg;


        tree[C]->Fill();
        tree[LL]->Fill();

        count++;

    } // END event loop



    //
    //  PRINT RESULTS
    //

    int nFiducial = tree[LL]->GetEntries("isFiducial");

    cout << endl << endl;
    cout << "Done!" << endl;
    cout << "Found " << nFiducial << "/" << reader.GetCurrentEntry() << " fiducial events" << endl;

    if (!fidOnly)
    {
        cout << endl << endl;
        cout << "Yields for phase space region:" << endl;
        for (unsigned i = 0; i < N; i++)
            cout << "\t\t" << selection[i] << ":\t" << tree[i]->GetEntries() << endl;
    }
    cout << endl << endl;

    cout << "Yields for fiducial region:" << endl;
    for (unsigned i = 0; i < N; i++)
            cout << "\t\t" << selection[i] << ":\t" << tree[i]->GetEntries("isFiducial") << endl;
    cout << endl << endl;



    //
    //  WRITE FILE
    //

    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();

    hPhaseSpaceEvents->Write();
    hFiducialEvents->Write();

    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
