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
//#include "TError.h"

// Custom
#include "Cuts2017.hh"

using namespace std;


/*
**  DrawResolution
**
**  Draws resolution of distributions for a "matched_" sample
*/ 

void DrawAccEff(const bool scale = kFALSE)
{

    //
    //  OPTIONS
    //

//  gErrorIgnoreLevel = kError;
    const double PI = TMath::Pi();



    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                   L4 = 0,  M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = {"4l",    "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = {5,       6,      7,      8,      9};
    TString lepChan[N]      = {_l,      _mu,    _l,     _l,     _e};
    int lColor[N]           = {lPurple, lBlue,  lPurple, lPurple, lRed};



    //
    //  INPUT FILE
    //

    TString inPath  = "acc_x_eff.root";//"unscaled4l_zz_4l.root";// "unscaled4l_phase_space.root";
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;



    //
    //  OUTPUT FILE
    //

    TString outName = "example.root";
    TFile *outFile  = new TFile(outName, "UPDATE");
    outFile->cd();


    // draw

    TH1D *h;
    TString suffix = "4l";  //"zz_4l";   // "phase_space";
    inFile->GetObject("4l/b_l1p_" + suffix, h);


    TCanvas *canvas = new TCanvas(suffix, "", lCanvasSize, lCanvasSize);
    canvas->cd();
    Facelift(canvas);
    canvas->SetCanvasSize(lCanvasSize, 0.625*lCanvasSize);
    canvas->SetMargin(lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);
    h->SetTitle("");
    h->SetLineWidth(2);     h->SetLineColor(lPurple);
//  h->SetLineColor(lRed);    h->SetFillColor(lRed);
//  h->SetLineColor(lBlue);    h->SetFillColor(lBlue);
    Facelift(h);
//  h->SetStats(kTRUE);
    h->Draw("E");    //h->Draw("HIST");
/*
    canvas->Update();
    TPaveStats *stats = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
    stats->SetOptStat(1000001110);
    stats->SetTextFont(lHelveticaMediumR);
    stats->SetTextSize(lSmall);
    stats->SetX1NDC(1.5 * lCanvasMargin); stats->SetX2NDC(0.4);
    stats->SetY1NDC(0.6); stats->SetY2NDC(1 - 1.5 * lCanvasMargin);
*/
    gPad->Modified();
    gPad->Update();

    canvas->Write();



    outFile->Close();
    inFile->Close();
}
