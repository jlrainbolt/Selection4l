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

void StackDists2l(const TString tag, bool useLog = kFALSE)
{
    const unsigned N_MC_ = 2;

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

    TString prefix  = "2l_" + tag;
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

    TH1 *mc[N][H][N_MC_];
    for (unsigned j = 0; j < N_MC_; j++) // sample loop
    {   
        TString inName = prefix + "_" + MC_SUFF[j] + ".root";
        TFile *inFile = TFile::Open(inName);

        cout << "Opened " << inName << endl;

        for (unsigned h = 0; h < H; h++)    // distribution loop
        {
            TH1 *hist;
            for (unsigned i = 1; i < N; i++)    // channel loop
            {
                inFile->GetObject(selection[i] + "/" + hname[h] + "_" + MC_SUFF[j], hist);
                hist->SetDirectory(0);

                float LUMI;
                if      (i == MM)
                    LUMI = MUON_TRIG_LUMI;
                else if (i == EE)
                    LUMI = ELEC_TRIG_LUMI * ELEC_TRIG_SF;
                float sf = LUMI * 1000 * XSEC[j] / NGEN[j];

                hist->Scale(sf);

                mc[i][h][j] = hist;
            }

            // Placeholder for 4l histogram
            mc[LL][h][j] = (TH1*) hist->Clone();
            mc[LL][h][j]->Reset();
        }
        inFile->Close();
    }
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
        for (unsigned j = 0; j < N_MC_; j++)
        {
            TH1 *mc_ = (TH1*) mc[1][h][j]->Clone();

            for (unsigned i = 2; i < N; i++)
                mc_->Add(mc[i][h][j]);

            mc[LL][h][j] = mc_;
        }
    }



    //
    //  ADD SAMPLES
    //

    TH1 *total[N][H];
    for (unsigned i = 0; i < N; i++)    // channel loop
    {
        for (unsigned h = 0; h < H; h++)    // distribution loop
        {
            TH1 *total_ = (TH1*) mc[i][h][0]->Clone();

            for (unsigned j = 1; j < N_MC_; j++)     // sample loop
                total_->Add(mc[i][h][j]);

            total[i][h] = total_;
            total[i][h]->Sumw2();
            total[i][h]->SetLineColor(0);
        }
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

            mc[i][h][ZZ]->SetTitle("");
            mc[i][h][ZZ]->SetStats(0);
            mc[i][h][ZZ]->SetFillColor(COLOR[ZZ]);
            mc[i][h][ZZ]->SetLineColor(COLOR[ZZ]);
            stack[i][h]->Add(mc[i][h][ZZ]);

            for (unsigned j_ = N_MC_; j_ > DY+1; j_--) // sample loop
            {
                unsigned j = j_ - 1;

                mc[i][h][j]->SetTitle("");
                mc[i][h][j]->SetStats(0);
                mc[i][h][j]->SetFillColor(COLOR[j]);
                mc[i][h][j]->SetLineColor(COLOR[j]);

                stack[i][h]->Add(mc[i][h][j]);
            }

            for (unsigned j = DY; j <= DY+1; j++) // sample loop
            {
                mc[i][h][j]->SetTitle("");
                mc[i][h][j]->SetStats(0);
                mc[i][h][j]->SetFillColor(COLOR[j]);
                mc[i][h][j]->SetLineColor(COLOR[j]);

                stack[i][h]->Add(mc[i][h][j]);
            }
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

            if (useLog)
            {
                data[i][h]->SetMinimum(1);
                total[i][h]->SetMinimum(1);
            }
            else
            {
                data[i][h]->SetMinimum(0);
                total[i][h]->SetMinimum(0);
            }

            data[i][h]->GetXaxis()->SetTitle(mc[i][h][DY]->GetXaxis()->GetTitle()); //FIXME

            ratio[i][h] = new TRatioPlot(data[i][h], total[i][h], "divsym");
            ratio[i][h]->SetH1DrawOpt("E");
            ratio[i][h]->SetH2DrawOpt("E");
            ratio[i][h]->SetSeparationMargin(0.0);
            ratio[i][h]->Draw();

            TPad *upper = ratio[i][h]->GetUpperPad(), *lower = ratio[i][h]->GetLowerPad();
            upper->cd();
            upper->SetLeftMargin(1.5 * lCanvasMargin);

            stack[i][h]->Draw("HIST SAME");
            stack[i][h]->GetYaxis()->SetTitle(mc[i][h][DY]->GetYaxis()->GetTitle());
            Facelift(stack[i][h]);
            data[i][h]->Draw("E SAME");

            float maximum = total[i][h]->GetMaximum();
            if (data[i][h]->GetMaximum() > maximum)
                maximum = data[i][h]->GetMaximum();
            data[i][h]->SetMaximum(1.1 * maximum);
            upper->Modified();
  
            ratio[i][h]->GetLowerRefYaxis()->SetTitle("Data/MC");
            Facelift(ratio[i][h]->GetLowerRefXaxis());
            Facelift(ratio[i][h]->GetLowerRefYaxis());
            ratio[i][h]->GetLowerRefYaxis()->SetTitleOffset(lTitleOffsetY);
            ratio[i][h]->GetLowerRefGraph()->SetMinimum(0.8);
            ratio[i][h]->GetLowerRefGraph()->SetMaximum(1.2);
            lower->SetBottomMargin(3 * lCanvasMargin);
            lower->SetLeftMargin(1.5 * lCanvasMargin);
            lower->Modified();

            legend->Draw();

            if (useLog)
            {
                stack[i][h]->SetMinimum(1);
                upper->SetLogy();
            }

            canvas[i][h]->Write();
        }
    }
    cout << "done!" << endl << endl;
    outFile->Close();

    cout << "Wrote canvases to " << outName << endl << endl << endl;
}
