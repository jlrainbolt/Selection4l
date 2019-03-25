// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TMathText.h"

// Custom
//#include "Cuts2016.hh"
#include "Cuts2017.hh"

using namespace std;


/*
**  StackBackground
**
**  Draws ratio plots from background studies
*/ 

void StackBackground(const bool signalOnly = kFALSE)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 7;   // Channel indices
    unsigned                    L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4, M3 = 5, E3 = 6;
    TString selection[N]    = { "4l",   "4m",   "2m2e", "2e2m", "4e",   "3m1e", "1m3e"  };
    unsigned chanIdx[N]     = { 10,     4,      6,      7,      9,      5,      8       };



    //
    //  INPUT FILES
    //

    TString prefix = signalOnly ? "bkg" : "bkg_all";

    TString muName = prefix + "_" + YEAR_STR + "_" + MU_SUFF + ".root";
    TFile   *muFile = TFile::Open(muName);
    cout << endl << endl << "Opened " << muName << endl;

    // Get directory keys
    vector<TString> hnames;
    TDirectory *keyDir = muFile->GetDirectory("/" + selection[M4], kTRUE, "GetDirectory");
    TKey *histKey;
    TIter next(keyDir->GetListOfKeys());
    while ((histKey = (TKey*) next()))
    {
        TString hname = histKey->GetName();
        hname.Resize(hname.Length() - (1 + MU_SUFF.Length()));  // truncate before suffix
        hnames.push_back(hname);
    }
    cout << "Got directory keys from " << muName << endl;

    const unsigned H = hnames.size();

    TString elName = prefix + "_" + YEAR_STR + "_" + EL_SUFF + ".root";
    TFile   *elFile = TFile::Open(elName);
    cout << endl << "Opened " << elName << endl << endl;



    //
    //  GET HISTOGRAMS
    //

    // Data histograms
    TH1 *data[N][H];

    for (unsigned i = 1; i < N; i++)
    {
        for (unsigned h = 0; h < H; h++)
        {
            if      ((i == M4) || (i == ME) || (i == M3))
                muFile->GetObject(selection[i] + "/" + hnames[h] + "_" + MU_SUFF, data[i][h]);
            else if ((i == E4) || (i == EM) || (i == E3))
                elFile->GetObject(selection[i] + "/" + hnames[h] + "_" + EL_SUFF, data[i][h]);

            data[i][h]->SetDirectory(0);
        }
    }
    cout << "Got data histograms" << endl;

    muFile->Close();
    cout << "Closed " << muName << endl;
    elFile->Close();
    cout << "Closed " << elName << endl;


    // MC histograms
    TH1 *mc[N][H][N_MC];


    for (unsigned j = 0; j < N_MC; j++)
    {
        TString mcName  = prefix + "_" + YEAR_STR + "_" + MC_SUFF[j] + ".root";
        TFile   *mcFile = TFile::Open(mcName);
        cout << "Opened " << mcName << endl;

        for (unsigned i = 1; i < N; i++)
        {
            float LUMI;

            if      ((i == M4) || (i == ME) || (i == M3))
                LUMI = MUON_TRIG_LUMI;
            else if ((i == E4) || (i == EM) || (i == E3))
                LUMI = ELEC_TRIG_LUMI * ELEC_TRIG_SF;
            
            float sf = LUMI * 1000 * XSEC[j] / NGEN[j];
            if (i == 1)
                cout << "Scale factor: " << sf << endl;

            for (unsigned h = 0; h < H; h++)
            {
                mcFile->GetObject(selection[i] + "/" + hnames[h] + "_" + MC_SUFF[j], mc[i][h][j]);
                mc[i][h][j]->SetDirectory(0);
                mc[i][h][j]->Scale(sf);
            }
        }
        mcFile->Close();
        cout << "Closed " << mcName << endl << endl;
    }



    //
    //  ADD CHANNELS
    //

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        data[ME][h]->Add(data[EM][h]);
        data[L4][h] = (TH1*) data[EM][h]->Clone();

        for (unsigned j = 0; j < N_MC; j++)
        {
            mc[ME][h][j]->Add(mc[EM][h][j]);
            mc[L4][h][j] = (TH1*) mc[EM][h][j]->Clone();
        }


        for (unsigned i = 1; i < N; i++)
        {
            if ((i != EM) && (i != ME))
            {
                data[L4][h]->Add(data[i][h]);
        
                for (unsigned j = 0; j < N_MC; j++)
                    mc[L4][h][j]->Add(mc[i][h][j]);
            }
        }
    }

    TH1 *total[N][H];
    for (unsigned i = 0; i < N; i++)    // channel loop
    {
        if (i == EM)    // skip "2e2m" since we added those (gross)
            continue;

        for (unsigned h = 0; h < H; h++)    // distribution loop
        {
            total[i][h] = (TH1*) mc[i][h][0]->Clone();
            total[i][h]->Sumw2();
            total[i][h]->SetLineColor(0);

            for (unsigned j = 1; j < N_MC; j++)     // sample loop
                total[i][h]->Add(mc[i][h][j]);
        }
    }



    //
    //  MAKE STACKS
    //

    THStack *stack[N][H];

    for (unsigned i = 0; i < N; i++)    // channel loop
    {
        if (i == EM)    // skip "2e2m" since we added those (gross)
            continue;

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
            stack[i][h] = new THStack(hnames[h] + "_" + selection[i], "");

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

    TLegend *legend = new TLegend(0.05, 0.55, 0.3, 0.95);

    TString lentry[6] = {_sp+_ZZ+_to+_4l, _sp+_Z+_to+_ll, _sp+_ttbar, _sp+_V+_V, _sp+_V+_V+_V, _sp+_H};
    Int_t lfill[6] = {lLightBlue, lYellow, lGreen, lOrange, lBlue, lPurple};

    legend->AddEntry(data[0][0], "Data", "LP");
    for (unsigned h = 0; h < 6; h++)
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

    TString outName = prefix + "_overlays_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    // Draw
    TCanvas *canvas[N][H];
    for (unsigned i = 0; i < N; i++)
    {
        if (signalOnly && ((i == E3) || (i == M3)))
            continue;

        if (i == EM)
            continue;

        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        for (unsigned h = 0; h < H; h++)
        {
            canvas[i][h] = new TCanvas(hnames[h] +"_"+ selection[i], "", lCanvasSize, lCanvasSize);

            Facelift(canvas[i][h]);
            canvas[i][h]->cd();

            data[i][h]->SetMinimum(0);
            stack[i][h]->SetMinimum(0);

            if (stack[i][h]->GetMaximum() > data[i][h]->GetMaximum())
                data[i][h]->SetMaximum(1.1 * stack[i][h]->GetMaximum());

            TRatioPlot *ratio = new TRatioPlot(data[i][h], total[i][h], "divsym");
            ratio->SetH1DrawOpt("E");
            ratio->SetH2DrawOpt("E");
            ratio->SetSeparationMargin(0.01);
            ratio->Draw();

            TPad *upper = ratio->GetUpperPad(), *lower = ratio->GetLowerPad();
            upper->cd();

            stack[i][h]->Draw("HIST SAME");
            Facelift(stack[i][h]);
            data[i][h]->Draw("E SAME");

            Facelift(ratio->GetUpperRefYaxis());
            ratio->GetUpperRefYaxis()->SetTitle("");
            ratio->GetUpperRefYaxis()->SetTitleOffset(lTitleOffsetY);
            upper->SetLeftMargin(1.2 * lCanvasMargin);
            upper->Modified();

            Facelift(ratio->GetLowerRefXaxis());
            Facelift(ratio->GetLowerRefYaxis());
            ratio->GetLowerRefGraph()->SetMinimum(0);
            ratio->GetLowerRefGraph()->SetMaximum(5);
            ratio->GetLowerRefYaxis()->SetTitle("Data/MC");
            lower->SetBottomMargin(3 * lCanvasMargin);
            lower->SetLeftMargin(1.2 * lCanvasMargin);
            lower->Modified();

            canvas[i][h]->Update();
            legend->Draw();

            canvas[i][h]->Write();
        }
    }

    outFile->Close();
    cout << "Wrote output to " << outName << endl << endl;
}
