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
#include "TLine.h"

// Custom
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;


/*
**  DrawMeanPileup
**
**  Draws nPV distributions and calculates mean for Z->mumu samples
*/ 

void DrawMeanPileup(TString selection = "mumu")
{

    //
    //  INPUT FILE
    //

    TString suffix = MU_SUFF;
    TString inName = "selected_" + suffix + ".root";
    TString inPath = EOS_PATH + "/Selected/" + YEAR_STR + "_update/" + inName;
//  TString inPath = EOS_PATH + "/Selected/" + YEAR_STR + "_new/" + inName;
    TFile *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;

    TTree *tree;



    //
    //  OUTPUT FILE
    //

    TString prefix  = "mean_pileup";
    TString outName = prefix + "_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  DRAW
    //

    inFile->GetObject(selection + "_" + suffix, tree);

    cout << "\t" << "Drawing histograms..." << flush;
    cout << selection + " tree has " << tree->GetEntries() << " events." << flush;

    // Get parameters
    TString hname = "nPV",  quantity = "nPV";
    int     bins = 34;
    float   xmin = 2,       xmax = 70;
    TString weight = "weight";

    // Create and draw histogram
    TH1D *h = new TH1D(YEAR_STR, "", bins, xmin, xmax);
    h->SetStatOverflows(TH1::kIgnore);
    tree->Draw(quantity + ">>+" + YEAR_STR, weight);

    TString xlabel = "n_{PV}";
    h->GetXaxis()->SetTitle(xlabel);
    h->GetYaxis()->SetTitle("Events");
    h->Sumw2(kTRUE);

    TCanvas *canvas = new TCanvas(YEAR_STR + "_" + selection + "_" + hname, "", 100, 100);
    canvas->cd();
    Facelift(canvas);
    canvas->SetMargin(1.2*lCanvasMargin, 0.9*lCanvasMargin, 0.9*lCanvasMargin, 0.6*lCanvasMargin);

    h->SetMaximum(1.25 * h->GetMaximum());
    h->SetMarkerStyle(kFullCircle);
    h->SetMarkerSize(2);
    h->SetMarkerColor(kBlack);
    h->SetLineWidth(2);
    Facelift(h);
    h->SetStats(1);

//  if (i == L4)
//      h->Fit("gaus", "", "SAME", -3, 3);
//  else
//      h->Fit("gaus", "", "SAME", -4, 4);

    h->Draw("E P");
//  h->GetListOfFunctions()->FindObject("gaus")->Draw("SAME");
    h->Write();

    canvas->Update();
    TPaveStats *stats = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
    stats->SetOptStat(2211);
    stats->SetTextFont(lHelveticaMediumR);
    stats->SetTextSize(lSmall);
    stats->SetX1NDC(0.5); stats->SetY1NDC(0.8);

    canvas->SaveAs(".pdf");
    canvas->Write();

    cout << "done!" << endl;



    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms, canvases to " << outName << endl << endl << endl;
}
