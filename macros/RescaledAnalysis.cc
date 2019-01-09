// STL
#include <vector>
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"

// Cuts
#include "Cuts2017.hh"

using namespace std;



/*
**  RescaledAnalysis
**
**  Reads dilepton events from a "Selected" level tree.  Calculates event rescale factors using
**  dilepton transverse momentum ratio plot.
*/

void RescaledAnalysis(const TString suffix, bool doWeights = kTRUE)
{

    //
    //  SAMPLE INFO
    //

    const bool isData = suffix.Contains(YEAR_STR);

    const unsigned N = 3;
    unsigned                   LL = 0,  MM = 1, EE = 2;     // Indices
    TString selection[N]    = {"ll",    "mumu", "ee"};
    unsigned chanIdx[N]     = {2,       3,      4};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "rescaled";
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


    // Lab frame objects
    TLorentzVector      z1p4,       l1p4,       l2p4;
    Short_t             z1pdg,      l1pdg,      l2pdg;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("runNum",   &runNum);               tree[i]->Branch("evtNum",   &evtNum);
        tree[i]->Branch("lumiSec",  &lumiSec);              tree[i]->Branch("nPV",      &nPV);
        tree[i]->Branch("met",      &met);                  tree[i]->Branch("channel",  &channel);
        tree[i]->Branch("weight",   &weight);
        
        if (doWeights)  tree[i]->Branch("rescale",  &rescale);

        tree[i]->Branch("z1p4",     &z1p4);                 tree[i]->Branch("z1pdg",    &z1pdg);
        tree[i]->Branch("l1p4",     &l1p4);                 tree[i]->Branch("l1pdg",    &l1pdg);
        tree[i]->Branch("l2p4",     &l2p4);                 tree[i]->Branch("l2pdg",    &l2pdg);
    }



    //
    //  INPUT FILE
    //

    TString inName  = "selected_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl;



    //
    //  HISTOGRAMS
    //

    TH1D *hTotalEvents, *hSelectedEvents;

    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);
    hTotalEvents->Sumw2();

    hSelectedEvents = new TH1D("SelectedEvents_" + suffix, "SelectedEvents", 10, 0.5, 10.5);
    hSelectedEvents->SetDirectory(outFile);
    hSelectedEvents->Sumw2();

    // First bin of SelectedEvents is number of generated events
    hSelectedEvents->SetBinContent(1,
                        hTotalEvents->GetBinContent(1) - 2 * hTotalEvents->GetBinContent(10));



    //
    //  WEIGHT UTILS
    //

    TH1 *hRescale[N];

    if (doWeights)
    {
        TString sfName = "../data/z1pt_rescale_" + YEAR_STR + ".root";  
        TFile*  sfFile = new TFile(sfName, "OPEN");
        for (unsigned i = 0; i < N; i++)
        {
            sfFile->GetObject(selection[i] + "/hist_z1pt_" + selection[i], hRescale[i]);
            hRescale[i]->SetDirectory(0);
        }
    }



    //
    //  INPUT BRANCHES
    //

    for (unsigned i = 1; i < N; i++)
    {
        TTreeReader reader(selection[i] + "_" + suffix, inFile);

        TTreeReaderValue    <Int_t>                 runNum_     (reader,    "runNum");
        TTreeReaderValue    <Int_t>                 evtNum_     (reader,    "evtNum");
        TTreeReaderValue    <Int_t>                 lumiSec_    (reader,    "lumiSec");
        TTreeReaderValue    <UShort_t>              nPV_        (reader,    "nPV");
        TTreeReaderValue    <Float_t>               met_        (reader,    "met");
        TTreeReaderValue    <Float_t>               weight_     (reader,    "weight");
        TTreeReaderValue    <UInt_t>                channel_    (reader,    "channel");
        TTreeReaderValue    <TLorentzVector>        z1p4_       (reader,    "z1p4");
        TTreeReaderValue    <Short_t>               z1pdg_      (reader,    "z1pdg");
        TTreeReaderValue    <TLorentzVector>        l1p4_       (reader,    "l1p4");
        TTreeReaderValue    <Short_t>               l1pdg_      (reader,    "l1pdg");
        TTreeReaderValue    <TLorentzVector>        l2p4_       (reader,    "l2p4");
        TTreeReaderValue    <Short_t>               l2pdg_      (reader,    "l2pdg");






        ////
        ////
        ////    EVENT LOOP
        ////
        ////


        int nEvents = reader.GetEntries(kTRUE);

        cout << endl;
        cout << "Running over " << nEvents << " selected " + selection[i] + " events..." << flush;

        while (reader.Next())
        {

            //
            //  EVENT INFO
            //                

            // Quantities copied directly to output tree
            runNum  = *runNum_;         evtNum  = *evtNum_;         lumiSec = *lumiSec_;
            weight  = *weight_;         nPV     = *nPV_;            met     = *met_;
            channel = *channel_;
            z1p4    = *z1p4_;           z1pdg   = *z1pdg_;
            l1p4    = *l1p4_;           l1pdg   = *l1pdg_;
            l2p4    = *l2p4_;           l2pdg   = *l2pdg_;



            //
            //  RESCALE
            //

            if (doWeights)
            {
                if (!isData)
                    rescale = hRescale[i]->GetBinContent(hRescale[i]->FindBin(z1p4.Pt()));
                else
                    rescale = 1;
            }



            //
            //  FILL
            //

            hSelectedEvents->Fill(chanIdx[i], weight * rescale);
            hSelectedEvents->Fill(chanIdx[LL], weight * rescale);

            tree[i]->Fill();
            tree[LL]->Fill();

        } // END event loop

        cout << "done!" << endl;

    } // END tree loop


    cout << endl << endl;
    cout << "Processed all trees in " << inName << endl;



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
