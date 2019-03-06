// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TError.h"

// Custom
#include "Cuts2016.hh"

using namespace std;


/*
**  DrawRunRange
*/ 

void DrawRunRange(const TString suffix)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 3;
    unsigned                   LL = 0,  MM = 1, EE = 2;     // Indices
    TString selection[N]    = {"ll",    "mumu", "ee"};
    unsigned chanIdx[N]     = {2,       3,      4};
    TString lepChan[N]      = {_l,      _mu,    _e};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "rr";
    TString outName = prefix + "_" + YEAR_STR + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  INPUT FILES
    //

    TString inName  = "selected_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;






    ////
    ////
    ////    CHANNEL LOOP
    ////
    ////


    float run[14] = {   272006, 275377,     // Period B
                        275656, 276284,     // Period C
                        276314, 276812,     // Period D
                        276830, 277421,     // Period E
                        277771, 278809,     // Period F
                        278819, 280396,     // Period G
                        281612, 284045  };  // Period H

    TString weight = suffix.EqualTo("zjets_m-50") ? "weight/trigWeight" : "weight/trigWeight/qtWeight";


    for (unsigned i = 1; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        TTree *tree;
        inFile->GetObject(selection[i] + "_" + suffix, tree);

        cout << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;



        //
        //  DRAW
        //

        cout << "\t" << "Drawing histograms..." << flush;

        // Create and draw histogram
        TH1D *h = new TH1D("runNum_" + suffix, "", 13, run);
//      TH1D *h = new TH1D("runNum", "", 284045 - 272007, 272007, 284045);
        tree->Draw("runNum >>+ runNum_" + suffix, weight);

        h->GetXaxis()->SetTitle("Run");
        h->Sumw2(kTRUE);
        h->Write();

        cout << "done!" << endl;
    }

    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
