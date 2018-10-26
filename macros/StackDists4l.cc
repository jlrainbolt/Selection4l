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
#include "Cuts2017.hh"

using namespace std;


/*
 **  StackDists4l
 **
 **  Scales and stacks all "unscaled4l_" distributions
 */

void StackDists4l()
{

    //
    //  OPTIONS
    //

//  gErrorIgnoreLevel = kError;
    const int nRebin = 5;



    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                    L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = { "4l",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = { 5,      6,      7,      8,      9};



    //
    //  DATA
    //

    TString prefix  = "unscaled4l";
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
    TDirectory *keyDir = muFile->GetDirectory("/" + selection[M4], kTRUE, "GetDirectory");
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
            if      (i == M4 || i == ME)
                muFile->GetObject(selection[i] + "/" + hname[h] + "_" + MU_SUFF, hist);
            else if (i == E4 || i == EM)
                elFile->GetObject(selection[i] + "/" + hname[h] + "_" + EL_SUFF, hist);

            hist->SetDirectory(0);

            data[i][h] = hist;
        }

        // Placeholder for 4l histogram
        data[L4][h] = (TH1*) hist->Clone();
        data[L4][h]->Reset();
    }
    muFile->Close();
    elFile->Close();

    cout << "Got data histograms" << endl << endl;



    //
    //  MONTE CARLO
    //

    TH1 *mc[N][H][N_MC];
    for (unsigned j = 0; j < N_MC; j++) // sample loop
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
                if      (i == M4 || i == ME)
                    LUMI = MUON_TRIG_LUMI;
                else if (i == E4 || i == EM)
                    LUMI = ELEC_TRIG_LUMI;
                float sf = LUMI * 1000 * XSEC[j] / NGEN[j];

                hist->Scale(sf);

                mc[i][h][j] = hist;
            }

            // Placeholder for 4l histogram
            mc[L4][h][j] = (TH1*) hist->Clone();
            mc[L4][h][j]->Reset();
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

        data[L4][h] = data_;


        // Monte Carlo
        for (unsigned j = 0; j < N_MC; j++)
        {
            TH1 *mc_ = (TH1*) mc[1][h][j]->Clone();

            for (unsigned i = 2; i < N; i++)
                mc_->Add(mc[i][h][j]);

            mc[L4][h][j] = mc_;
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

            for (unsigned j = 1; j < N_MC; j++)     // sample loop
                total_->Add(mc[i][h][j]);

            total[i][h] = total_;
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

            for (unsigned j = 0; j < N_MC; j++) // sample loop
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
    TString lentry[5] = {_sp+_ZZ+_to+_4l, _sp+_Z+_to+_ll, _sp+_H, _sp+_ttbar, _sp+_V+_V};
    Int_t lfill[5] = {lLightBlue, lYellow, lPurple, lGreen, lRed};

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

    TString outName = "stacks4l_raw_2017.root";
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

            ratio[i][h] = new TRatioPlot(data[i][h], total[i][h], "divsym");
//          ratio[i][h]->SetGraphDrawOpt("Z");
            ratio[i][h]->Draw();

            TPad *upper = ratio[i][h]->GetUpperPad(), *lower = ratio[i][h]->GetLowerPad();
            upper->cd();

            data[i][h]->Draw("E");

            // Stacks must be drawn in order for their axes to exist
            if (data[i][h]->GetMaximum() > stack[i][h]->GetMaximum())
            {
                data[i][h]->Draw("E");
                Facelift(data[i][h]);
                upper->Modified();

                stack[i][h]->Draw("HIST SAME");
                Facelift(stack[i][h]);
                upper->Modified();

                data[i][h]->Draw("E SAME");
            }
            else
            { 
                stack[i][h]->Draw("HIST");
                Facelift(stack[i][h]);
                upper->Modified();

                data[i][h]->Draw("E SAME");
                Facelift(data[i][h]);
                upper->Modified();
            }

            upper->RedrawAxis();

            Facelift(ratio[i][h]->GetLowerRefXaxis());
            Facelift(ratio[i][h]->GetLowerRefYaxis());
            lower->Modified();

//          canvas[i][h]->cd();
            legend->Draw();
            gPad->Modified();

            canvas[i][h]->Write();
        }
    }
    cout << "done!" << endl << endl;
    outFile->Close();

    cout << "Wrote canvases to " << outName << endl << endl << endl;
}
