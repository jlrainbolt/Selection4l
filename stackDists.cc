#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"

using namespace std;

void ApplyAlpha(float rgb[][3], const int, const float);
void FaceliftTH1(TH1*);
void FaceliftLegend(TLegend*);


void stackDists(const TString filePath, const TString fileSuffix)
{


    //--- CHANNEL OPTIONS ---//

    // Channel names and scale factors
    const unsigned N = 4;       TString prefix[N];          Double_t scale[N];  
    unsigned M4 = 0;            prefix[M4] = "4m";          scale[M4] = 15.8508;
    unsigned ME = 1;            prefix[ME] = "2m2e";        scale[ME] = 34.0485;
    unsigned EM = 2;            prefix[EM] = "2e2m";        scale[EM] = 30.3738;
    unsigned E4 = 3;            prefix[E4] = "4e";          scale[E4] = 55.9601;
    Double_t totalScale = 5.34227e-10;


    // Channel colors
    Float_t parula[N][3] = {{0.2422, 0.1504, 0.6603},   // Indigo
                            {0.1380, 0.6276, 0.8973},   // Teal
                            {0.5616, 0.7942, 0.3045},   // Orange
                            {0.9956, 0.7434, 0.2371}};  // Yellow
    ApplyAlpha(parula, N, 0.75);
    Int_t color[N];   TColor *color_[N];
    for (unsigned i = 0; i < N; i++)
    {
        color[i] = TColor::GetFreeColorIndex();
        color_[i] = new TColor(color[i], parula[i][0], parula[i][1], parula[i][2]);
    }

    // Channel legend labels
    TString sp = "\\,", mu = "\\mu", el = "\\mbox{e}", plus = "^{+}", minus = "^{-}";
    TString label[N];
    label[M4] = sp + mu+plus + mu+minus + mu+plus + mu+minus;
    label[ME] = sp + mu+plus + mu+minus + el+plus + el+minus;
    label[EM] = sp + el+plus + el+minus + mu+plus + mu+minus;
    label[E4] = sp + el+plus + el+minus + el+plus + el+minus;



    //--- DISTRIBUTION OPTIONS ---//

    vector<TString> dname = {"l1pt_d", "z2m_d", "m3l_d", "angle_d"};

    // Create objects
    vector<THStack*> stack;
    vector<TCanvas*> canvas;    int cdim = 800;
    vector<TH1*> sumHist, hists[N];
    for (unsigned j = 0; j < dname.size(); j++)
    {
        canvas.push_back(new TCanvas(dname[j], "", cdim, cdim));
        stack.push_back(new THStack(dname[j], ""));
    }

    // Axis titles
    TString BF = "\\mathcal{B}(\\mbox{Z} \\to 4\\ell)";
    TString pT1 = "p_{\\mbox{T} 1}", Z2 = "\\mbox{Z}_{2}";
    vector<TString> dtitle = {pT1, "m_{"+Z2+"}", "m_{3\\ell}", "\\theta_{"+pT1+","+Z2+"}"};
    vector<TString> dunit = {"10 GeV", "5 GeV", "10 GeV", "0.4 rad"};



    //--- LOAD HISTOGRAMS ---//

    for (unsigned i = 0; i < N; i++)
    {
        TString rootFile = filePath + prefix[i] + "_" + fileSuffix + ".root";
        TFile *file = TFile::Open(rootFile, "READ");
        TDirectory *dir = file->GetDirectory("/Calculation/Signal - Background", kTRUE, "GetDirectory");
        Double_t weight = scale[i] * totalScale;

        for (unsigned j = 0; j < dname.size(); j++)
        {
            TH1F *hist;
            dir->GetObject(dname[j], hist);
            hist->SetDirectory(0);
            hist->SetTitle("");
            hist->Sumw2(kTRUE);
            hist->Scale(weight); 
            hist->SetLineColorAlpha(color[i], 0.75);       //  hist->SetLineWidth(2);
            hist->SetFillColorAlpha(color[i], 0.75);
            stack[j]->Add(hist);
            hists[i].push_back(hist);

            if (sumHist.size() <= j)
            {
                sumHist.push_back((TH1*) hist->Clone());
                sumHist.back()->SetDirectory(0);
            }
            else
                sumHist[j]->Add(hist);
        }

        file->Close();
    }



    //--- LEGEND ---//

    float c_mrg = .12, l_mrg = 0.04;
    TLegend *legend = new TLegend(.5 +2.*c_mrg-l_mrg, .5 + 2.*c_mrg-l_mrg, 1.-l_mrg, 1.-l_mrg);
    for (unsigned i_ = dname.size(); i_ > 0; i_--)
    {
        unsigned i = i_ - 1;
        legend->AddEntry(hists[i][0], label[i], "F");
    }
    FaceliftLegend(legend);



    //--- DRAW ---//

    for (unsigned j = 0; j < dname.size(); j++)
    {
        canvas[j]->cd();
        canvas[j]->SetMargin(c_mrg, c_mrg, c_mrg, c_mrg);

        sumHist[j]->SetYTitle(BF + "\\mbox{ / " + dunit[j] + "}");
        sumHist[j]->SetXTitle(dtitle[j]);
        sumHist[j]->SetLineColor(kBlack);
        sumHist[j]->SetLineWidth(2);
        FaceliftTH1(sumHist[j]);

        sumHist[j]->Draw("E");
        stack[j]->Draw("HIST SAME");
        sumHist[j]->Draw("E SAME");

        legend->Draw();

        gPad->RedrawAxis();
        gPad->Update();
    }



    //--- WRITE OUTPUT ---//

    TFile *outFile = new TFile("dists_" + fileSuffix + ".root", "RECREATE");

    TDirectory *c_dir = outFile->mkdir("Canvases");
    c_dir->cd();
    for (unsigned j = 0; j < dname.size(); j++)
    {
        canvas[j]->Write();
        canvas[j]->SaveAs("../2016_12a/" + dname[j] + ".png");
    }

    TDirectory *s_dir = outFile->mkdir("Stacks");
    s_dir->cd();
    for (unsigned j = 0; j < dname.size(); j++)
        stack[j]->Write();

    TDirectory *h_dir = outFile->mkdir("Sum Histograms");
    h_dir->cd();
    for (unsigned j = 0; j < dname.size(); j++)
        sumHist[j]->Write();

    TDirectory *o_dir = outFile->mkdir("Histograms");
    TDirectory *o_subdir[N];
    for (unsigned i = 0; i < N; i++)
    {
        o_subdir[i] = o_dir->mkdir(prefix[i]);
        o_subdir[i]->cd();
        for (unsigned j = 0; j < dname.size(); j++)
            hists[i][j]->Write();
    }

    outFile->Close();

}



void ApplyAlpha(float rgb[][3], const int size, const float alpha)
{
    for (unsigned i = 0; i < size; i++)
    {
        for (unsigned j = 0; j < 3; j++)
            rgb[i][j] += (1 - rgb[i][j]) * (1 - alpha);
    }
}

void FaceliftTH1(TH1 *hist)
{
    int helv = 43;      // font code 3 is expressed in pixels
    // Point font sizes
    float axlabel = 24, axtitle = 36, axoff = 1.1;

    hist->SetStats(0);
    hist->SetMinimum(0);

    hist->GetXaxis()->SetLabelFont(helv);
    hist->GetXaxis()->SetLabelSize(axlabel);
    hist->GetXaxis()->SetTitleFont(helv);
    hist->GetXaxis()->SetTitleSize(axtitle);
//  hist->GetXaxis()->SetTitleOffset(axoff);

    hist->GetYaxis()->SetLabelFont(helv);
    hist->GetYaxis()->SetLabelSize(axlabel);
    hist->GetYaxis()->SetTitleFont(helv);
    hist->GetYaxis()->SetTitleSize(axtitle);
    hist->GetYaxis()->SetTitleOffset(axoff);
}

void FaceliftLegend(TLegend *legend)
{
    int helv = 43;      // font code 3 is expressed in pixels
    // Point font sizes
    float textsize = 36;

//  legend->SetFillStyle(0);
//  legend->SetBorderSize(0);

    legend->SetTextFont(helv);
    legend->SetTextSize(textsize);

//  legend->SetEntrySeparation();
}
