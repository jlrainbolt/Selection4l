#include <vector>
#include <sstream>
#include <tuple>

#include "TString.h"
#include "TFile.h"
#include "TParameter.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMathText.h"


using namespace std;



void overlayHists(const TString suffix)
{
    bool drawNorm = kTRUE;
    bool drawLog  = kFALSE;


    // Selection
    const unsigned N = 5;
    unsigned                M4 = 0, ME = 1, EM = 2, E4 = 3, L4 = 4;
    TString selection[N] = {"4m",   "2m2e", "2e2m", "4e",   "4l"};


    // Levels
    const unsigned M = 4;
    unsigned                PS = 0,        FR = 1,        GS = 2,           RS = 3;            
    TString inPrefix[M] = { "PhaseSpace",  "Fiducial",    "GenSelected",    "RecoSelected"  };
    int color[M]        = { 46,            30,            9,                lLightBlue      };
    int fill[M]         = { 0,             0,             0,                3944            };
    TString entry[M]    = { "Phase Sp.",   "Fiducial",    "Selected",       "Reco Sel."     };



    //--- HIST KEYS ---//
    
    // Get histogram names from one of the files
    TFile *inFile = TFile::Open("hists_" + inPrefix[0] + "_" + suffix + ".root");

    vector<TString> hname_4l;
    TKey *histKey;

    TDirectory *dir_4l = inFile->GetDirectory("/"+selection[L4], kTRUE, "GetDirectory");
    TIter next_4l(dir_4l->GetListOfKeys());

    while ((histKey = (TKey*)next_4l()))
        hname_4l.push_back(histKey->GetName());

    inFile->Close();



    //--- STORE HISTS ---//
    
    vector<TH1*> hists[M][N];

    for (unsigned h = 0; h < M; h++)
    {
        inFile = TFile::Open("hists_" + inPrefix[h] + "_" + suffix + ".root");

        for (unsigned i = 0; i < N; i++)
        {
            vector<TString> hname = hname_4l;

            for (unsigned j = 0; j < hname.size(); j++)
            {
                TH1 *hist;

                inFile->GetObject(selection[i] + "/" + hname[j], hist);
                hist->SetDirectory(0);
                hist->SetStats(0);
                hist->SetLineWidth(2);
                hist->SetFillStyle(fill[h]);
                hists[h][i].push_back(hist);
            }
        }
        inFile->Close();
    }



    //--- FILL STACKS ---//

    vector<THStack*> stack[N];

    for (unsigned i = 0; i < N; i++)
    {
        vector<TString> hname = hname_4l;
        for (unsigned j = 0; j < hname.size(); j++)
        {
            stack[i].push_back(new THStack(hname[j], hname[j]));

            for (unsigned h = 0; h < M; h++)
            {
                TH1 *hist = hists[h][i][j];

                if (drawNorm)
                    hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));

                stack[i][j]->Add(hist);
            }
        }
    }



    //--- LEGEND ---//

    float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
    float RightPosition = 1,        RightMargin = -lLegendMargin;
    float TopPosition = 1,          TopMargin = -lLegendMargin;
    float BottomPosition = TopPosition - 0.065 * 4.;
    float BottomMargin = 2. * lCanvasMargin - lLegendMargin;
    TLegend *legend = new TLegend(LeftPosition + LeftMargin, BottomPosition - TopMargin,
                                    TopPosition + TopMargin, TopPosition + TopMargin);

    for (unsigned h = 0; h < M; h++)
    {
        TH1D* hist = new TH1D(inPrefix[h], "", 1, 0, 1);
        hist->SetFillColor(color[h]);
        hist->SetFillStyle(fill[h]);
        hist->SetLineColor(color[h]);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, entry[h], "LE");
    }
    Facelift(legend);



    //--- CANVASES ---//

    TString outPrefix = "overlays";
    if (drawNorm)
        outPrefix += "_norm";
    if (drawLog)
        outPrefix += "_log";
    TFile *outFile = new TFile(outPrefix + "_" + suffix + ".root", "RECREATE");
    vector<TCanvas*> canvas[N];

    for (unsigned i = 0; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        vector<TString> hname = hname_4l;

        for (unsigned j = 0; j < hname.size(); j++)
        {
            canvas[i].push_back(new TCanvas(hname[j] + "_" + selection[i], "", lCanvasSize, lCanvasSize));

            Facelift(canvas[i][j]);
            canvas[i][j]->cd();
            stack[i][j]->Draw("E NOSTACK");
            Facelift(stack[i][j]);
/*
            for (unsigned h = 0; h < M; h++)
            {
                TString option;
                if      (h == PS)
                    option = "HIST";
                else if (h == RS)
                    option = "E SAME";
                else
                    option = "HIST SAME";

                if (drawNorm)
                    hists[h][i][j]->DrawNormalized(option);
                else
                    hists[h][i][j]->Draw(option);

                Facelift(hists[h][i][j]);
            }
*/
            if (drawLog)
            {
                stack[i][j]->SetMinimum(1);
                gPad->SetLogy();
            }
            gPad->Modified();
            gPad->RedrawAxis();

            legend->Draw();
            gPad->Modified();

            canvas[i][j]->Write();
        }
    }
    outFile->Close();
}
