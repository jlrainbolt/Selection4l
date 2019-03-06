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
 **  StackDists2l
 **
 **  Scales and stacks all "unscaled2l_" distributions
 */

void StackRuns()
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

    TString prefix  = "rr_" + YEAR_STR;
    cout << endl << endl;

    // Muon file
    TString muName = prefix + "_" + MU_SUFF + ".root";
    TFile *muFile = TFile::Open(muName);
    cout << "Opened " << muName << endl;

    // Electron file
    TString elName = prefix + "_" + EL_SUFF + ".root";
    TFile *elFile = TFile::Open(elName);
    cout << "Opened " << elName << endl;


    // Now get the histograms
    TH1* data[N];

    for (unsigned i = 1; i < N; i++)    // channel loop
    {
        if      (i == MM)
            muFile->GetObject(selection[i] + "/runNum_" + MU_SUFF, data[i]);
        else if (i == EE)
            elFile->GetObject(selection[i] + "/runNum_" + EL_SUFF, data[i]);

        data[i]->SetDirectory(0);
    }
    muFile->Close();
    elFile->Close();

    cout << "Got data histograms" << endl << endl;



    //
    //  MONTE CARLO
    //

    float LUMI[14] = {0, 5.75, 0, 2.573, 0, 4.242, 0, 4.025, 0, 3.105, 0, 7.576, 0, 8.651};

    TH1 *mc[N][N_MC];
    for (unsigned j = 0; j < N_MC; j++) // sample loop
    {   
        TString inName = prefix + "_" + MC_SUFF[j] + ".root";
        TFile *inFile = TFile::Open(inName);

        cout << "Opened " << inName << endl;

        for (unsigned i = 1; i < N; i++)    // channel loop
        {
            inFile->GetObject(selection[i] + "/runNum_" + MC_SUFF[j], mc[i][j]);
            mc[i][j]->SetDirectory(0);

            for (unsigned k = 1; k <= mc[i][j]->GetNbinsX(); k++)
            {
                float sf = LUMI[k] * 1000 * XSEC[j] / NGEN[j];
                mc[i][j]->SetBinContent(k, mc[i][j]->GetBinContent(0) * sf);
            }
        }

        inFile->Close();
    }
    cout << "Got MC histograms" << endl << endl;



    //
    //  ADD CHANNELS
    //

    cout << endl << "Combining..." << flush;

    // Data
    TH1 *data_ = (TH1*) data[1]->Clone();

    for (unsigned i = 2; i < N; i++)
        data_->Add(data[i]);

    data[LL] = data_;


    // Monte Carlo
    for (unsigned j = 0; j < N_MC; j++)
    {
        TH1 *mc_ = (TH1*) mc[1][j]->Clone();

        for (unsigned i = 2; i < N; i++)
            mc_->Add(mc[i][j]);

        mc[LL][j] = mc_;
    }



    //
    //  ADD SAMPLES
    //

    TH1 *total[N];
    for (unsigned i = 0; i < N; i++)    // channel loop
    {
        TH1 *total_ = (TH1*) mc[i][0]->Clone();

        for (unsigned j = 1; j < N_MC; j++)     // sample loop
            total_->Add(mc[i][j]);

        total[i] = total_;
        total[i]->Sumw2();
        total[i]->SetLineColor(0);
    }



    //
    //  MAKE STACKS
    //

    THStack *stack[N];

    for (unsigned i = 0; i < N; i++)    // channel loop
    {
        data[i]->SetTitle("");
        data[i]->SetStats(0);
        data[i]->SetMarkerColor(kBlack);
        data[i]->SetMarkerStyle(kFullCircle);
        data[i]->SetMarkerSize(2);
        data[i]->SetLineWidth(2);
        data[i]->SetLineColor(kBlack);


        // Create and fill stack object
        stack[i] = new THStack("runNum_" + selection[i], "");

        mc[i][ZZ]->SetTitle("");
        mc[i][ZZ]->SetStats(0);
        mc[i][ZZ]->SetFillColor(COLOR[ZZ]);
        mc[i][ZZ]->SetLineColor(COLOR[ZZ]);
        stack[i]->Add(mc[i][ZZ]);

        for (unsigned j_ = N_MC; j_ > DY+1; j_--) // sample loop
        {
            unsigned j = j_ - 1;

            mc[i][j]->SetTitle("");
            mc[i][j]->SetStats(0);
            mc[i][j]->SetFillColor(COLOR[j]);
            mc[i][j]->SetLineColor(COLOR[j]);

            stack[i]->Add(mc[i][j]);
        }

        for (unsigned j = DY; j <= DY+1; j++) // sample loop
        {
            mc[i][j]->SetTitle("");
            mc[i][j]->SetStats(0);
            mc[i][j]->SetFillColor(COLOR[j]);
            mc[i][j]->SetLineColor(COLOR[j]);

            stack[i]->Add(mc[i][j]);
        }
    }



    //
    //  LEGEND
    //

    float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
    float RightPosition = 1,        RightMargin = -lLegendMargin;
    float TopPosition = 1,          TopMargin = -lLegendMargin;
    float BottomPosition = TopPosition - 0.085 * 6;//0.065 * 6.;
    float BottomMargin = 2. * lCanvasMargin - lLegendMargin;
    TLegend *legend = new TLegend(LeftPosition + LeftMargin, BottomPosition - TopMargin,
                                    TopPosition + TopMargin, TopPosition + TopMargin);

    TString dentry = _sp+"\\mbox{Data}";
    TString lentry[5] = {_sp+_Z+_to+_ll, _sp+_H, _sp+_ttbar, _sp+_V+_V, _sp+_ZZ+_to+_4l};
    Int_t lfill[5] = {lYellow, lPurple, lGreen, lOrange, lLightBlue};

    TH1D* dummy = new TH1D(dentry, "", 1, 0, 1);
    dummy->SetMarkerColor(kBlack);
    dummy->SetMarkerStyle(kFullCircle);
    dummy->SetMarkerSize(2);
    dummy->SetLineWidth(2);
    dummy->SetLineColor(kBlack);
    legend->AddEntry(dummy, dentry, "LP");

    for (unsigned h = 0; h < 5; h++)
    {
        TH1D* hist = new TH1D(lentry[h], "", 1, 0, 1);
        hist->SetFillColor(lfill[h]);
        hist->SetLineColor(lfill[h]);
        legend->AddEntry(hist, lentry[h], "F");
    }
    Facelift(legend);



    //
    //  OUTPUT FILE
    //

    TString outName = "stacksRun_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  DRAW CANVASES
    //

    TCanvas *canvas[N];
    TRatioPlot *ratio[N];

    for (unsigned i = 0; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        canvas[i] = new TCanvas("runNum_" + selection[i],"", lCanvasSize, lCanvasSize);

        Facelift(canvas[i]);
        canvas[i]->cd();

        data[i]->SetMinimum(0);
        total[i]->SetMinimum(0);

        data[i]->GetXaxis()->SetTitle(mc[i][DY]->GetXaxis()->GetTitle()); //FIXME

        ratio[i] = new TRatioPlot(data[i], total[i], "divsym");
        ratio[i]->SetH1DrawOpt("E");
        ratio[i]->SetH2DrawOpt("E");
        ratio[i]->SetSeparationMargin(0.0);
        ratio[i]->Draw();

        TPad *upper = ratio[i]->GetUpperPad(), *lower = ratio[i]->GetLowerPad();
        upper->cd();
        upper->SetLeftMargin(1.5 * lCanvasMargin);

        stack[i]->Draw("HIST SAME");
        stack[i]->GetYaxis()->SetTitle(mc[i][DY]->GetYaxis()->GetTitle());
        Facelift(stack[i]);
        data[i]->Draw("E SAME");

        float maximum = total[i]->GetMaximum();
        if (data[i]->GetMaximum() > maximum)
            maximum = data[i]->GetMaximum();
        data[i]->SetMaximum(1.1 * maximum);
        upper->Modified();

        ratio[i]->GetLowerRefYaxis()->SetTitle("Data/MC");
        Facelift(ratio[i]->GetLowerRefXaxis());
        Facelift(ratio[i]->GetLowerRefYaxis());
        ratio[i]->GetLowerRefYaxis()->SetTitleOffset(lTitleOffsetY);
        ratio[i]->GetLowerRefGraph()->SetMinimum(0.8);
        ratio[i]->GetLowerRefGraph()->SetMaximum(1.2);
        lower->SetBottomMargin(3 * lCanvasMargin);
        lower->SetLeftMargin(1.5 * lCanvasMargin);
        lower->Modified();

        legend->Draw();

        canvas[i]->Write();
    }
    cout << "done!" << endl << endl;
    outFile->Close();

    cout << "Wrote canvases to " << outName << endl << endl << endl;
}
