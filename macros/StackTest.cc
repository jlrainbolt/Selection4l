// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TColor.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TMathText.h"

// Cuts
//#include "Cuts2017.hh"
#include "Cuts2016.hh"

using namespace std;


/*
 **  StackTest
 **
 **  Scales and stacks all "unscaled2l_" distributions for Bacon sanity check
 */

void StackTest()
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 3;
    unsigned                    LL = 0, MM = 1, EE = 2;     // Indices
    TString selection[N]    = { "ll",   "mumu", "ee"};
    unsigned chanIdx[N]     = { 2,      3,      4};



    //
    //  DATA
    //

    TString prefix  = "unscaled2l";
    cout << endl << endl;

    // Muon file
    TString muName = prefix + "_" + MU_SUFF + ".root";
    TFile *muFile = TFile::Open(muName);
    cout << "Opened " << muName << endl;

    // Electron file
    TString elName = prefix + "_" + EL_SUFF + ".root";
    TFile *elFile = TFile::Open(elName);
    cout << "Opened " << elName << endl;


    // Get histogram names
    vector<TString> hname;
    TDirectory *keyDir = muFile->GetDirectory("/" + selection[MM], kTRUE, "GetDirectory");
    TKey *histKey;
    TIter next(keyDir->GetListOfKeys());
    while ((histKey = (TKey*) next()))
    {   
        TString hname_ = histKey->GetName(); 
        hname_.Resize(hname_.Length() - (1 + MU_SUFF.Length()));    // truncate before suffix
        hname.push_back(hname_);
    }

    // Now get the histograms
    const unsigned H = hname.size();
    TH1* data[N][H];

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        TH1 *hist;
        for (unsigned i = 1; i < N; i++)    // channel loop
        {
            if      (i == MM)
                muFile->GetObject(selection[i] + "/" + hname[h] + "_" + MU_SUFF, hist);
            else if (i == EE)
                elFile->GetObject(selection[i] + "/" + hname[h] + "_" + EL_SUFF, hist);

            hist->SetDirectory(0);

            hist->Rebin(2);

            data[i][h] = hist;
        }

        // Placeholder for 2l histogram
        data[LL][h] = (TH1*) hist->Clone();
        data[LL][h]->Reset();
    }
    muFile->Close();
    elFile->Close();

    cout << "Got data histograms" << endl << endl;



    //
    //  MONTE CARLO
    //

    TH1 *mc[N][H];
    TString inName = prefix + "_DYJetsToLL.root";
    TFile *inFile = TFile::Open(inName);

    cout << "Opened " << inName << endl;

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        TH1 *hist;
        for (unsigned i = 1; i < N; i++)    // channel loop
        {
            inFile->GetObject(selection[i] + "/" + hname[h] + "_DYJetsToLL", hist);
            hist->SetDirectory(0);

            hist->Rebin(2);
            hist->Scale(data[i][h]->Integral() / hist->Integral());

            mc[i][h] = hist;
        }

        // Placeholder for 2l histogram
        mc[LL][h] = (TH1*) hist->Clone();
        mc[LL][h]->Reset();
    }
    inFile->Close();
    cout << "Got MC histograms" << endl << endl;



    //
    //  ADD CHANNELS
    //

    cout << endl << "Combining..." << flush;
    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        // Data
        TH1 *data_ = (TH1*) data[1][h]->Clone();

        for (unsigned i = 2; i < N; i++)
            data_->Add(data[i][h]);

        data[LL][h] = data_;


        // Monte Carlo
        TH1 *mc_ = (TH1*) mc[1][h]->Clone();

        for (unsigned i = 2; i < N; i++)
            mc_->Add(mc[i][h]);

        mc[LL][h] = mc_;
    }



    //
    //  MAKE STACKS
    //

    THStack *stack[N][H];

    for (unsigned i = 0; i < N; i++)    // channel loop
    {
        for (unsigned h = 0; h < H; h++)    // distribution loop
        {
            data[i][h]->SetTitle("");
            data[i][h]->SetStats(0);
            data[i][h]->SetMarkerColor(kBlack);
            data[i][h]->SetMarkerStyle(kFullCircle);
            data[i][h]->SetMarkerSize(2);
            data[i][h]->SetLineWidth(2);
            data[i][h]->SetLineColor(kBlack);


            // Create and fill stack object
            stack[i][h] = new THStack(hname[h] + "_" + selection[i], "");

            mc[i][h]->SetTitle("");
            mc[i][h]->SetStats(0);
            mc[i][h]->SetFillColor(COLOR[DY]);
            mc[i][h]->SetLineColor(COLOR[DY]);
            mc[i][h]->SetLineWidth(0);
            stack[i][h]->Add(mc[i][h]);
        }
    }



    //
    //  OUTPUT FILE
    //

    TString tag     = "test";
    TString outName = "stacks2l_" + tag + "_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  DRAW CANVASES
    //

    TCanvas *canvas[N][H];
    TRatioPlot *ratio[N][H];

    for (unsigned i = 0; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        for (unsigned h = 0; h < H; h++)
        {
            canvas[i][h] = new TCanvas(hname[h] + "_" + selection[i],"", lCanvasSize, lCanvasSize);

            Facelift(canvas[i][h]);
            canvas[i][h]->cd();

            ratio[i][h] = new TRatioPlot(data[i][h], mc[i][h], "divsym");
            ratio[i][h]->SetH1DrawOpt("E");
            ratio[i][h]->SetH2DrawOpt("E");
            ratio[i][h]->SetSeparationMargin(0.0);
            ratio[i][h]->Draw();

            TPad *upper = ratio[i][h]->GetUpperPad(), *lower = ratio[i][h]->GetLowerPad();
            upper->cd();

            stack[i][h]->Draw("HIST SAME");
            stack[i][h]->GetXaxis()->SetTitle(data[i][h]->GetXaxis()->GetTitle());
            Facelift(stack[i][h]);
            data[i][h]->Draw("E SAME");

            float maximum = mc[i][h]->GetMaximum();
            if (data[i][h]->GetMaximum() > maximum)
                maximum = data[i][h]->GetMaximum();
            data[i][h]->SetMaximum(1.1 * maximum);
  
            Facelift(ratio[i][h]->GetLowerRefXaxis());
            Facelift(ratio[i][h]->GetLowerRefYaxis());
            ratio[i][h]->GetLowerRefGraph()->SetMinimum(0.8);
            ratio[i][h]->GetLowerRefGraph()->SetMaximum(1.2);
            lower->SetBottomMargin(3 * lCanvasMargin);
            lower->Modified();

            canvas[i][h]->Write();
        }
    }
    cout << "done!" << endl << endl;
    outFile->Close();

    cout << "Wrote canvases to " << outName << endl << endl << endl;
}
