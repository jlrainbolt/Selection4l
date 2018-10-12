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



void subtractHists(const TString suffix)
{
    // Selection
    const unsigned N = 5;
    unsigned                M4 = 0, ME = 1, EM = 2, E4 = 3, L4 = 4;
    TString selection[N] = {"4m",   "2m2e", "2e2m", "4e",   "4l"};



    //--- HIST KEYS ---//
    
    // Get histogram names from gen file
    TFile *genFile = TFile::Open("hists_gen_" + suffix + ".root");

    vector<TString> hname_4l;
    TKey *histKey;

    TDirectory *dir_4l = genFile->GetDirectory("/"+selection[L4], kTRUE, "GetDirectory");
    TIter next_4l(dir_4l->GetListOfKeys());

    while ((histKey = (TKey*)next_4l()))
        hname_4l.push_back(histKey->GetName());



    //--- STORE HISTS ---//
    

    TFile *recoFile = TFile::Open("hists_" + suffix + ".root");

    vector<TH1*> genHist[N], recoHist[N];

    for (unsigned i = 0; i < N; i++)
    {
        vector<TString> hname = hname_4l;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1 *_genHist, *_recoHist;

            genFile->GetObject(selection[i] + "/" + hname[h], _genHist);
            _genHist->SetDirectory(0);
            _genHist->SetLineWidth(2);
            _genHist->SetFillStyle(0);
            genHist[i].push_back(_genHist);

            recoFile->GetObject(selection[i] + "/" + hname[h], _recoHist);
            _recoHist->SetDirectory(0);
            _recoHist->SetLineWidth(2);
            recoHist[i].push_back(_recoHist);
        }
    }
    genFile->Close();
    recoFile->Close();



    //--- SUBTRACT ---//

    vector<TH1*> diffHist[N];

    for (unsigned i = 0; i < N; i++)
    {
        vector<TString> hname = hname_4l;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1 *_diffHist = (TH1*) genHist[i][h]->Clone();
            _diffHist->Add(recoHist[i][h], -1);

            diffHist[i].push_back(_diffHist);
        }
    }


/*
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
*/


    //--- WRITE ---//

    TFile *diffFile = new TFile("diff_" + suffix + ".root", "RECREATE");
    for (unsigned i = 0; i < N; i++)
    {
        diffFile->mkdir(selection[i]);
        diffFile->cd(selection[i]);

        vector<TString> hname = hname_4l;

        for (unsigned h = 0; h < hname.size(); h++)
            diffHist[i][h]->Write();
    }
    diffFile->Close();



    //--- CANVASES ---//

    TFile *compFile = new TFile("comp_" + suffix + ".root", "RECREATE");
    vector<TCanvas*> canvas[N];

    for (unsigned i = 0; i < N; i++)
    {
        compFile->mkdir(selection[i]);
        compFile->cd(selection[i]);

        vector<TString> hname = hname_4l;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            canvas[i].push_back(new TCanvas(hname[h], "", lCanvasSize, lCanvasSize));

            Facelift(canvas[i][h]);
            canvas[i][h]->cd();

            // Stacks must be drawn in order for their axes to exist (Rene Brun, WHY??)
            if (genHist[i][h]->GetMaximum() > recoHist[i][h]->GetMaximum())
            {
                genHist[i][h]->Draw("HIST");
                Facelift(genHist[i][h]);
                gPad->Modified();

                recoHist[i][h]->Draw("E SAME");
                Facelift(recoHist[i][h]);
                gPad->Modified();
            }
            else
            { 
                recoHist[i][h]->Draw("E");
                Facelift(recoHist[i][h]);
                gPad->Modified();

                genHist[i][h]->Draw("HIST SAME");
                Facelift(genHist[i][h]);
                gPad->Modified();
            }

            gPad->RedrawAxis();

//          legend->Draw();
//          gPad->Modified();

            canvas[i][h]->Write();
        }
    }
    compFile->Close();
}
