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
//#include "Cuts2018.hh"
//#include "Cuts2017.hh"
#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;


/*
**  DrawMassResolution
**
**  Draws mass resolution of distributions for a "matched_" sample
*/ 

void DrawMassResolution()
{

    //
    //  INPUT FILES
    //

    const unsigned M = 2;
    unsigned                   DY = 0,        ZZ = 1;
    TString suffix[2]       = {"zjets_m-50",  "zz_4l"};
    TString inName[2],  inPath[2];
    TFile *inFile[2];

    for (unsigned i = 0; i < M; i++)
    {
        inName[i]   = "matched_" + suffix[i] + ".root";
        inPath[i]   = EOS_PATH + "/Boosted/" + YEAR_STR + "_update/" + inName[i];
        inFile[i]   = TFile::Open(inPath[i]);

        cout << endl << endl << "Opened " << inPath[i] << endl << endl;
    }

    TTree *tree;




    //
    //  SAMPLE INFO
    //

    const unsigned N = 6;
    unsigned                   MM = 0,  EE = 1,  L4 = 2,  M4 = 3,  ME = 4,  E4 = 5;     // Indices
    TString selection[N]    = {"mumu",  "ee",    "4l",    "4m",    "2m2e",  "4e",   };



    //
    //  OUTPUT FILE
    //

    TString prefix  = "mass_resolution";
    TString outName = prefix + "_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  DRAW
    //


    for (unsigned i = 0; i < N; i++)
    {
        bool isDY = i < 2;
        unsigned I = isDY ? DY : ZZ;

        inFile[I]->GetObject(selection[i] + "_" + suffix[I], tree);

        cout << "\t" << "Drawing histograms..." << flush;
        cout << selection[i] + " tree has " << tree->GetEntries() << " events." << flush;

        // Get parameters
        TString hname = isDY ? "z1m" : "zzm";
        TString quantity = isDY ? "z1p4.M()" : "zzp4.M()";
        int     bins = 100;
        float   xmin = -20,     xmax = 20,      width = 1;
        float   ngen = isDY ? NGEN_ZJETS : NGEN_ZZ_4L;
        float   xsec = isDY ? XSEC_ZJETS : XSEC_ZZ_4L;

        TString weight = "weight";

        // Add subtraction to quantity
        quantity = "gen_" + quantity + " - " + quantity;

        // Create and draw histogram
        TH1D *h = new TH1D(hname + "_" + selection[i], "", bins, xmin, xmax);
        tree->Draw(quantity + ">>+" + hname + "_" + selection[i], weight);

        TString xlabel;
        if (i == MM)
            xlabel = "\\Delta m_{\\mu^{+}\\mu^{-}}\\mbox{ (GeV)}";
        else if (i == EE)
            xlabel = "\\Delta m_{\\mbox{e}^{+}\\mbox{e}^{-}}\\mbox{ (GeV)}";
        else if (i == L4)
            xlabel = "\\Delta m_{4l}\\mbox{ (GeV)}";
        else if (i == M4)
            xlabel = "\\Delta m_{4\\mu}\\mbox{ (GeV)}";
        else if (i == ME)
            xlabel = "\\Delta m_{2\\mu2\\mbox{e}}\\mbox{ (GeV)}";
        else if (i == E4)
            xlabel = "\\Delta m_{4\\mbox{e}}\\mbox{ (GeV)}";

        h->GetXaxis()->SetTitle(xlabel);
        h->GetYaxis()->SetTitle("Events");
        h->Sumw2(kFALSE);
        h->Scale(INT_LUMI * 1000. * xsec / ngen);

        TCanvas *canvas = new TCanvas(YEAR_STR + "_" + selection[i] + "_" + hname + "_resolution", "", 100, 100);
        canvas->cd();
        Facelift(canvas);
        canvas->SetMargin(1.2*lCanvasMargin, 0.9*lCanvasMargin, 0.9*lCanvasMargin, 0.6*lCanvasMargin);

        h->SetMaximum(1.25 * h->GetMaximum());
        h->SetFillColor(kGray);
        h->SetLineColor(kGray);
        h->SetLineWidth(4);
        Facelift(h);
        h->SetStats(1);

        if (i == L4)
            h->Fit("gaus", "", "SAME", -3, 3);
        else
            h->Fit("gaus", "", "SAME", -4, 4);

        h->Draw("HIST");
        h->GetListOfFunctions()->FindObject("gaus")->Draw("SAME");
        h->Write();

        canvas->Update();
        TPaveStats *stats = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
        stats->SetOptFit(1);
        stats->SetOptStat(0);
//      stats->SetOptStat(10);
        stats->SetTextFont(lHelveticaMediumR);
        stats->SetTextSize(lSmall);
        stats->SetX1NDC(0.5); stats->SetY1NDC(0.8);

//      float x = 0.5 * width;
//      float y1 = 0, y2 = gPad->GetUymax();
//      TLine *line[2] = {new TLine(-x, y1, -x, y2), new TLine(x, y1, x, y2)};
//      for (unsigned l = 0; l < 2; l++)
//      {
//          line[l]->SetLineColor(kBlack);
//          line[l]->SetLineStyle(kDashed);
//          line[l]->SetLineWidth(2);
//          line[l]->Draw();
//      }

        canvas->SaveAs(".pdf");
        canvas->Write();

        cout << "done!" << endl;
    }



    outFile->Close();
    inFile[DY]->Close();
    inFile[ZZ]->Close();

    cout << endl << "Wrote histograms, canvases to " << outName << endl << endl << endl;
}
