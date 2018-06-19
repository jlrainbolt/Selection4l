#include <vector>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TParameter.h"


using namespace std;

void stackResults(const TString inFile, const TString selection)
{   
    // Choose electron or muon
    Bool_t selMuMu;
    if (selection == "mumu")
        selMuMu = kTRUE;
    else if (selection == "ee")
        selMuMu = kFALSE;

    TFile *file = TFile::Open(inFile, "UPDATE");

    // Get lists of samples and their contents
    TDirectory *dir = file->GetDirectory("/Histograms", kTRUE, "GetDirectory");
    TIter next(dir->GetListOfKeys());
    TKey *subdirKey;
   
    vector<TDirectory*> dataSubdir, mcSubdir;
    vector<TString> dataSuffix, mcSuffix;
    vector<Int_t> dataColor, mcColor;
    Float_t lumi;   vector<Float_t> xsec;

    while ((subdirKey = (TKey*)next()))
    {
        TString suffix = subdirKey->GetName();
        TDirectory *subdir;             subdir = dir->GetDirectory(suffix);
        TParameter<Bool_t> *isData;     subdir->GetObject("isData", isData);
        TParameter<Int_t> *color;       subdir->GetObject("color", color);
        TParameter<Float_t> *factor;
        if (isData->GetVal())
        {
            dataSubdir.push_back(subdir);
            dataSuffix.push_back(suffix);       dataColor.push_back(color->GetVal());
            subdir->GetObject("lumi", factor);  lumi = factor->GetVal();
        }
        else
        {
            mcSubdir.push_back(subdir);
            mcSuffix.push_back(suffix);         mcColor.push_back(color->GetVal());
            subdir->GetObject("xsec", factor);  xsec.push_back(factor->GetVal());
        }
    }


    // Create stacks, canvases, legend
    vector<TString> hname = {"DileptonMass", "DileptonPt", "Lepton1Pt", "Lepton2Pt", "Lepton1Eta", 
        "Lepton2Eta", "nPV", "nExtraLeps", "xLepPt", "xLepEta", "xLepQ", "xLepIso", "xLepTrkIso",
        "xLepD0", "xLepDz", "xLepIsTrig", "TotalEvents"};
    if (selMuMu)
    {
        vector<TString> hname2 = {"xLepIsGLB", "xLepNMatchStn", "xLepNPixHits","xLepNTkLayers",
            "xLepNValidHits", "xLepMuNChi2"};
        hname.insert(hname.end(), hname2.begin(), hname2.end());
    }
    vector<THStack*> dataStack, mcStack;
    vector<TCanvas*> canvas;
    for (unsigned h = 0; h < hname.size(); h++)
    {
        canvas.push_back(new TCanvas(hname[h], hname[h], 640, 480));
        dataStack.push_back(new THStack(hname[h], hname[h]));
        mcStack.push_back(new THStack(hname[h], hname[h]));
    }
    TLegend *legend = new TLegend(.75, .55, .95, .95);


    // Fill stacks and legend
    for (int i = 0; i < dataSubdir.size(); i++)
    {
        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1* hist;
            dataSubdir[i]->GetObject(hname[h] + "_" + dataSuffix[i], hist);
            hist->SetLineColor(dataColor[i]);   hist->SetLineWidth(2);
            hist->Sumw2();          dataStack[h]->Add(hist);

            if (i == 0)
                dataStack[h]->SetTitle(hist->GetTitle());
            if (h == 0)
                legend->AddEntry(hist, dataSuffix[i], "E");
        }
    }

    for (int j_ = mcSubdir.size(); j_ > 0; j_--)
    {
        unsigned j = j_ - 1;
        TH1* hTotalEvents;
        mcSubdir[j]->GetObject("TotalEvents_" + mcSuffix[j], hTotalEvents);
        Float_t ngen = hTotalEvents->GetBinContent(1) - 2*hTotalEvents->GetBinContent(10);
        Float_t weight = lumi * 1000. * xsec[j] / ngen;
        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1* hist;
            mcSubdir[j]->GetObject(hname[h] + "_" + mcSuffix[j], hist);
            hist->SetFillColor(mcColor[j]); hist->SetLineColor(mcColor[j]);
            hist->Scale(weight);    mcStack[h]->Add(hist);

            if (j == 0)
                mcStack[h]->SetTitle(hist->GetTitle());
            if (h == 0)
                legend->AddEntry(hist, mcSuffix[j], "F");
        }
    }


    // Draw on canvases
    for (unsigned h = 0; h < hname.size(); h++)
    {
        canvas[h]->cd();
        if (dataStack[h]->GetMaximum() > mcStack[h]->GetMaximum())
        {
            dataStack[h]->Draw("E");
            mcStack[h]->Draw("HIST SAME");
        }
        else
            mcStack[h]->Draw("HIST");
        dataStack[h]->Draw("E SAME");
        legend->Draw();
    }

    // Write to file
//  TFile *outFile = new TFile(output, "RECREATE");
    TDirectory *c_dir = file->mkdir("Canvases");
    c_dir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        canvas[h]->Write();
    TDirectory *s_dir = file->mkdir("Stacks");
    TDirectory *s_mcSubdir = s_dir->mkdir("Monte Carlo");
    s_mcSubdir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        mcStack[h]->Write();
    TDirectory *s_dataSubdir = s_dir->mkdir("Data");
    s_dataSubdir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        dataStack[h]->Write();
    file->Close();
// DELETE YOUR SHIT
    
}
