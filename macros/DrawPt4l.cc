// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"

// Custom
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;


/*
**  DrawPt4l
**
**  Draws lepton momenta by flavor for a "selected_" sample
*/ 

void DrawPt4l(const TString suffix, const TString year)
{
    if (!year.EqualTo(YEAR_STR))
    {
        cout << "Wrong year in header file!" << endl;
        return;
    }

    //
    //  SAMPLE INFO
    //

    const unsigned N = 2;
    unsigned                    M = 0,  E = 1;     // Indices
    TString selection[N]    = {"mu",    "e"};
    TString lepChan[N]      = {"#mu",  "e"};
    int     lepPDG[N]       = { 13,     11  };

    int lColor[4]           = {lBlue, lOrange, lPurple, lGreen};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "pt4l";
    TString outName = prefix + "_" + year + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");


    //
    //  INPUT FILE
    //

//  TString inName  = "selected_" + suffix + ".root";
    TString inName  = "boosted_" + suffix + ".root";
//  TString inPath  = EOS_PATH + "/Selected/" + year + "_new/" + inName;
    TString inPath  = EOS_PATH + "/Boosted/" + year + "_new/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;






    ////
    ////
    ////    CHANNEL LOOP
    ////
    ////


    TTree *tree;
    inFile->GetObject("4l_" + suffix, tree);

    cout << "4l tree has " << tree->GetEntries() << " events." << flush;

    TH1D *h[N][4];
    TCanvas *c[N];

    for (unsigned i = 0; i < N; i++)
    {
        cout << "\t" << "Drawing " << selection[i] << " histograms..." << flush;

        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

//      TString evtWeight = "(weight/trigWeight/qtWeight) * ";
        TString evtWeight = "weight * isFiducial * ";

        for (unsigned j = 0; j < 4; j++)
        {
            TString index = TString::Format("%i", j+1);
            TString name = "l" + index + "pt_" + selection[i];
            TString pdg = TString::Format("%i", lepPDG[i]);

            TString weight = evtWeight + "(abs(l" + index + "pdg) == " + pdg + ")";
            cout << weight << endl;

            h[i][j] = new TH1D(name, "", 60, 0, 60);
            h[i][j]->Sumw2(kTRUE);
            h[i][j]->SetStats(0);
            h[i][j]->SetMinimum(0);

            tree->Draw("l" + index + "p4.Pt()>>+" + name, weight);

            Facelift(h[i][j]);
            h[i][j]->SetXTitle("p_{T}^{" + lepChan[i] + "} (GeV)");
            h[i][j]->SetYTitle("Events/GeV");
            h[i][j]->SetLineColor(j+6);
            h[i][j]->SetLineWidth(2);

            h[i][j]->Write();
        }

        c[i] = new TCanvas("c_lpt_" + selection[i], "", lCanvasSize, lCanvasSize);
        Facelift(c[i]);
        c[i]->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c[i]->SetMargin(lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);

        c[i]->cd();

        h[i][3]->Draw("HIST");

        for (unsigned j = 0; j < 4; j++)
            h[i][j]->Draw("HIST SAME");

        c[i]->Write();
        cout << "done!" << endl;
    }

    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
