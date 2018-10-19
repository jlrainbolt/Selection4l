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



void stackHists()
{
    // Selection
    const unsigned N = 8;
    unsigned                MM = 0, EE = 1, LL = 2, M4 = 3, ME = 4, EM = 5, E4 = 6, L4 = 7;
    TString selection[N] = {"mumu", "ee",   "ll",   "4m",   "2m2e", "2e2m", "4e",   "4l"};


    // Samples
    const unsigned M = 9;
    TString suffix[M]   = { "zjets_m-50", "ttbar", "ww_2l2nu", "wz_2l2q", "wz_3lnu", "zz_2l2q",
                            "zz_4l", "ggH_zz_4l", "vbfH_zz_4l"};



    //--- HIST KEYS ---//
    
    // Get histogram names from data file
    TFile *dataFile = TFile::Open("hists_data.root");

    vector<TString> hname_2l, hname_4l;
    TKey *histKey;

    TDirectory *dir_2l = dataFile->GetDirectory("/"+selection[LL], kTRUE, "GetDirectory");
    TIter next_2l(dir_2l->GetListOfKeys());

    while ((histKey = (TKey*)next_2l()))
        hname_2l.push_back(histKey->GetName());


    TDirectory *dir_4l = dataFile->GetDirectory("/"+selection[L4], kTRUE, "GetDirectory");
    TIter next_4l(dir_4l->GetListOfKeys());

    while ((histKey = (TKey*)next_4l()))
        hname_4l.push_back(histKey->GetName());



    //--- STORE HISTS ---//
    

    vector<TH1*> dataHist[N];
    vector<THStack*> stack[N];

    for (unsigned i = 0; i < N; i++)
    {
        vector<TString> hname = i < M4 ? hname_2l : hname_4l;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1 *hist;
            dataFile->GetObject(selection[i] + "/" + hname[h], hist);
            hist->SetDirectory(0);

            dataHist[i].push_back(hist);
            stack[i].push_back(new THStack(hname[h], ""));
        }
    }
    dataFile->Close();



    //--- FILL STACKS ---//

    for (unsigned j_ = 0; j_ < M; j_++)
    {
        unsigned j = M - 1 - j_;
        TFile *mcFile = TFile::Open("hists_" + suffix[j] + ".root");

        for (unsigned i = 0; i < N; i++)
        {
            vector<TString> hname = i < M4 ? hname_2l : hname_4l;

            for (unsigned h = 0; h < hname.size(); h++)
            {
                TH1 *hist;
                mcFile->GetObject(selection[i] + "/" + hname[h], hist);
                hist->SetDirectory(0);
                stack[i][h]->Add(hist);
            }
        }
        mcFile->Close();
    }



    //--- LEGEND ---//

    float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
    float RightPosition = 1,        RightMargin = -lLegendMargin;
    float TopPosition = 1,          TopMargin = -lLegendMargin;
    float BottomPosition = TopPosition - 0.065 * 6.;
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



    //--- CANVASES ---//

    TFile *outFile = new TFile("stacks.root", "RECREATE");
    vector<TCanvas*> canvas[N];

    for (unsigned i = 0; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        vector<TString> hname = i < M4 ? hname_2l : hname_4l;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            canvas[i].push_back(new TCanvas(hname[h], "", lCanvasSize, lCanvasSize));

            Facelift(canvas[i][h]);
            canvas[i][h]->cd();
            dataHist[i][h]->Draw("E");

            // Stacks must be drawn in order for their axes to exist (Rene Brun, WHY??)
            if (dataHist[i][h]->GetMaximum() > stack[i][h]->GetMaximum())
            {
                dataHist[i][h]->Draw("E");
                Facelift(dataHist[i][h]);
                gPad->Modified();

                stack[i][h]->Draw("HIST SAME");
                Facelift(stack[i][h]);
                gPad->Modified();

                dataHist[i][h]->Draw("E SAME");
            }
            else
            { 
                stack[i][h]->Draw("HIST");
                Facelift(stack[i][h]);
                gPad->Modified();

                dataHist[i][h]->Draw("E SAME");
                Facelift(dataHist[i][h]);
                gPad->Modified();
            }

            gPad->RedrawAxis();

            legend->Draw();
            gPad->Modified();

            canvas[i][h]->Write();
        }
    }
    outFile->Close();
}
