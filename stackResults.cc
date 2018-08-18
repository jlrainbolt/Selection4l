#include <vector>
#include <sstream>

#include "TString.h"
#include "TFile.h"
#include "TMacro.h"
#include "TH1.h"
#include "TColor.h"
#include "TParameter.h"


using namespace std;

void ApplyAlpha(float rgb[][3], const int, const float);



void stackResults(const TString inFile, const TString selection)
{
    // Choose selection type
    Bool_t sel2l = kFALSE, sel4l = kFALSE;
    TString sigSuffix, selName;
    if (selection == "mumu" || selection == "2m" || selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;      selName = "DY";     sigSuffix = "dy_m-50";
    }
    else if (selection == "2m2e" || selection == "2mu2e" || selection == "2e2m" || selection== "2e2mu"
             || selection == "4mu" || selection == "4m" || selection == "4e")
    {
        sel4l = kTRUE;      selName = "4l";     sigSuffix = "zz_4l";
    }
    TFile *file = TFile::Open(inFile, "UPDATE");


    // Create colors
    const unsigned L = 8;
    Float_t lines[L][3] = { {0, 0, 0},                  // Black    0
                            {0, 0.4470, 0.7410},        // Blue     1
                            {0.8500, 0.3250, 0.0980},   // "Red"    2
                            {0.9290, 0.6940, 0.1250},   // Yellow   3
                            {0.4940, 0.1840, 0.5560},   // Purple   4
                            {0.4660, 0.6740, 0.1880},   // Green    5
                            {0.3010, 0.7450, 0.9330},   // Lt blue  6
                            {0.6350, 0.0780, 0.1840}};  // Dk red   7
    Float_t alpha = 0.75;   ApplyAlpha(lines, L, alpha);
    TColor *lines_[L];
    for (unsigned i = 0; i < L; i++)
    {
        lines_[i] = new TColor(1179 + (Int_t) i, lines[i][0], lines[i][1], lines[i][2]);
    }


    // Get lists of samples and their contents
    TDirectory *dir = file->GetDirectory("/Histograms", kTRUE, "GetDirectory");
    TIter next(dir->GetListOfKeys());
    TKey *subdirKey;
   
    vector<TDirectory*> dataSubdir, mcSubdir;
    vector<TString> dataSuffix, mcSuffix;
    vector<Int_t> dataColor, mcColor;
    vector<Bool_t> mcSignal;
    Float_t lumi;   vector<Float_t> xsec;

    while ((subdirKey = (TKey*)next()))
    {
        TString suffix = subdirKey->GetName();
        TDirectory *subdir;             subdir = dir->GetDirectory(suffix);
        TParameter<Bool_t> *isData;     subdir->GetObject("isData", isData);
        TParameter<Int_t> *color;       subdir->GetObject("color", color);
        TParameter<Bool_t> *isSignal;   subdir->GetObject("isSignal", isSignal);
        TParameter<Float_t> *factor;

        if (isData->GetVal())
        {
            dataSubdir.push_back(subdir);
            dataSuffix.push_back(suffix);
            dataColor.push_back(color->GetVal());
            subdir->GetObject("lumi", factor);  lumi = factor->GetVal();
        }
        else
        {
            mcSubdir.push_back(subdir);
            mcSuffix.push_back(suffix);
            mcColor.push_back(color->GetVal());
            subdir->GetObject("xsec", factor);  xsec.push_back(factor->GetVal());
            mcSignal.push_back(isSignal->GetVal());
        }
    }


    // Get accepted events hist before scaling
    unsigned si;    TH1* hAcceptedEvents;
    for (unsigned j = 0; j < mcSubdir.size(); j++)
    {
        if (mcSuffix[j] == sigSuffix)
        {
            mcSubdir[j]->GetObject("AcceptedEvents_" + mcSuffix[j], hAcceptedEvents);
            hAcceptedEvents->SetName("AcceptedEvents");
            break;
        }
    }


    // Create stacks, canvases, legend
    vector<TString> hname;
    if (sel2l)
        hname = {"nPV", "met", "z1m", "z1pt",
                    "l1pt", "l1eta", "l1iso", "l1pdg", "l2pt", "l2eta", "l2iso", "l2pdg",
                    "TotalEvents", "AcceptedEvents"};
    else if (sel4l)
        hname = {"nPV", "met", "zzm", "zzpt", "z1m", "z1pt", "z2m", "z2pt",
                    "l1pt", "l1eta", "l1iso", "l1pdg", "l2pt", "l2eta", "l2iso", "l2pdg",
                    "l3pt", "l3eta", "l3iso", "l3pdg", "l4pt", "l4eta", "l4iso", "l4pdg",
                    "TotalEvents", "AcceptedEvents",
                    "l1pt_d", "z2m_d", "m3l_d", "angle_d"};

    vector<THStack*> dataStack, mcStack;
    vector<TH1*> dataSum, mcSum, sigSum, bgSum;
    vector<TCanvas*> canvas;
    for (unsigned h = 0; h < hname.size(); h++)
    {
        canvas.push_back(new TCanvas(hname[h], hname[h], 640, 480));
        dataStack.push_back(new THStack(hname[h], hname[h]));
        mcStack.push_back(new THStack(hname[h], hname[h]));
    }
    TLegend *legend = new TLegend(.75, .25, .95, .95);


    // Macro to write out yields
    TMacro log("printYields", "printYields");
    stringstream buff;
    buff << "void printYields()" << endl;
    buff << "{" << endl;
    buff << "    cout << endl << endl;" << endl;
    buff << "    cout << \" SAMPLE YIELDS\" << endl;" << endl;
    buff << "    cout << \"------------------------------------------------------\" << endl;"<<endl;


    // Fill stacks
    for (int i = 0; i < dataSubdir.size(); i++)
    {
        TH1* hTotalEvents;
        dataSubdir[i]->GetObject("TotalEvents_" + dataSuffix[i], hTotalEvents);

        buff << "    cout << \"" << setw(16) << dataSuffix[i] << "\\t";
        buff << setw(10) << hTotalEvents->GetBinContent(7) << "\\t+- ";
        buff << setw(10) << hTotalEvents->GetBinError(7) << "\" << endl;" << endl;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1* hist;
            dataSubdir[i]->GetObject(hname[h] + "_" + dataSuffix[i], hist);
            hist->Sumw2(kTRUE);
            hist->SetLineColor(dataColor[i]);   hist->SetLineWidth(2);
            dataStack[h]->Add(hist);

            if (i == 0)
                dataStack[h]->SetTitle(hist->GetTitle());
            if (h == 0)
                legend->AddEntry(hist, dataSuffix[i], "E");

            if (dataSum.size() <= h)
                dataSum.push_back((TH1*) hist->Clone(hname[h]));
            else
                dataSum[h]->Add(hist);
        }
    }
    buff << "    cout << endl;" << endl;

    for (unsigned j_ = 0; j_ < mcSubdir.size(); j_++)
    {
        unsigned j = j_;
        if (sel2l)
            j = mcSubdir.size() - 1 - j_;

        TH1* hTotalEvents;
        mcSubdir[j]->GetObject("TotalEvents_" + mcSuffix[j], hTotalEvents);
        Float_t ngen = hTotalEvents->GetBinContent(1) - 2 * hTotalEvents->GetBinContent(10);
        Float_t weight = lumi * 1000. * xsec[j] / ngen;

        buff << "    cout << \"" << setw(16) << mcSuffix[j] << "\\t";
        buff << setw(10) << weight * hTotalEvents->GetBinContent(7) << "\\t+- ";
        buff << setw(10) << weight * hTotalEvents->GetBinError(7) << "\" << endl;" << endl;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1* hist;
            mcSubdir[j]->GetObject(hname[h] + "_" + mcSuffix[j], hist);
            hist->SetFillColor(mcColor[j]); hist->SetLineColor(mcColor[j]);
            hist->Sumw2(kTRUE);
            hist->Scale(weight);    mcStack[h]->Add(hist);

            if (j == 0)
                mcStack[h]->SetTitle(hist->GetTitle());

            if (mcSum.size() <= h)
                mcSum.push_back((TH1*) hist->Clone(hname[h]));
            else
                mcSum[h]->Add(hist);

            if (mcSignal[j])
            {
                if (sigSum.size() <= h)
                    sigSum.push_back((TH1*) hist->Clone(hname[h]));
                else
                    sigSum[h]->Add(hist);
            }
            else
            {
                if (bgSum.size() <= h)
                    bgSum.push_back((TH1*) hist->Clone(hname[h]));
                else
                    bgSum[h]->Add(hist);
            }
        }
    }
    buff << "    cout << \"------------------------------------------------------\" << endl;"<<endl;
    buff << "    cout << endl << endl;" << endl;


    // Fill legend
    for (unsigned j_ = 0; j_ < mcSubdir.size(); j_++)
    {
        unsigned j = j_;
        if (sel2l)
            j = mcSubdir.size() - 1 - j_;

        TH1* hist;
        mcSubdir[j]->GetObject(hname[0] + "_" + mcSuffix[j], hist);
        legend->AddEntry(hist, mcSuffix[j], "F");
    }


    // Draw on canvases
    for (unsigned h = 0; h < hname.size(); h++)
    {
        dataStack[h]->SetMinimum(0.001);
        mcStack[h]->SetMinimum(0.001);
/*
        if (sel4l && h == hname.size() - 1)
        {
            dataStack[h]->SetMaximum(1000);
            mcStack[h]->SetMaximum(1000);
        }
*/
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
   

    // Create signal - bg histograms 
    vector<TH1*> diffHist;
    for (unsigned h = 0; h < hname.size(); h++)
    {
        diffHist.push_back((TH1*) dataSum[h]->Clone());
        diffHist.back()->Add(bgSum[h], -1);
    }


    // Get signal/background MC yields
    unsigned ti = hname.size() - 2;
    TH1* hDiff = (TH1*) dataSum[ti]->Clone();   hDiff->Add(bgSum[ti], -1);
    buff << "    cout << \" TOTALS\" << endl;" << endl;
    buff << "    cout << \"------------------------------------------------------\" << endl;"<<endl;
    buff << "    cout << \"" << setw(16) << "Data" << "\\t";
    buff << setw(10) << dataSum[ti]->GetBinContent(7) << "\\t+- ";
    buff << setw(10) << dataSum[ti]->GetBinError(7) << "\" << endl;" << endl;
    buff << "    cout << \"" << setw(16) << "MC" << "\\t";
    buff << setw(10) << mcSum[ti]->GetBinContent(7) << "\\t+- ";
    buff << setw(10) << mcSum[ti]->GetBinError(7) << "\" << endl << endl;" << endl;
    buff << "    cout << \"" << setw(16) << "Signal MC" << "\\t";
    buff << setw(10) << sigSum[ti]->GetBinContent(7) << "\\t+- ";
    buff << setw(10) << sigSum[ti]->GetBinError(7) << "\" << endl;" << endl;
    buff << "    cout << \"" << setw(16) << "Background MC" << "\\t";
    buff << setw(10) << bgSum[ti]->GetBinContent(7) << "\\t+- ";
    buff << setw(10) << bgSum[ti]->GetBinError(7) << "\" << endl << endl;" << endl;
    buff << "    cout << \"" << setw(16) << "Data - Bgnd" << "\\t";
    buff << setw(10) << hDiff->GetBinContent(7) << "\\t+- ";
    buff << setw(10) << hDiff->GetBinError(7) << "\" << endl;" << endl;
    buff << "    cout << \"------------------------------------------------------\" << endl;"<<endl;
    buff << "    cout << endl << endl;" << endl;


    // Get acceptance/efficiency factors
    Double_t nSamp = hAcceptedEvents->GetBinContent(1);
    Double_t nSpc = hAcceptedEvents->GetBinContent(2);
    Double_t nSel = hAcceptedEvents->GetBinContent(3);
    Double_t accEff = nSel / nSpc, eUnc = accEff * sqrt(1./nSel + 1./nSpc);
    buff << "    cout << \" ACCEPTANCE & EFFICIENCY (UNSCALED)\" << endl;" << endl;
    buff << "    cout << \"------------------------------------------------------\" << endl;"<<endl;
    buff << "    cout << \"" << setw(16) << "Total " + selName << "\\t";
    buff << setw(10) << nSamp << "\" << endl;" << endl;
    buff << "    cout << \"" << setw(16) << "Phase space " + selection << "\\t";
    buff << setw(10) << nSpc << "\\t+- " << setw(10) << sqrt(nSpc) << "\" << endl;" << endl;
    buff << "    cout << \"" << setw(16) << "Selected " + selection << "\\t";
    buff << setw(10) << nSel << "\\t+- " << setw(10) << sqrt(nSel) << "\" << endl << endl;" << endl;
    buff << "    cout << \"" << setw(16) << selection + " acc * eff" << "\\t";
    buff << setw(10) << nSel / nSpc << "\\t+- " << setw(10) << eUnc << "\" << endl;" << endl;
    buff << "    cout << \"------------------------------------------------------\" << endl;"<<endl;



    // Write to file

    // Sum histograms
    TDirectory *x_dir = file->mkdir("Calculation");
    x_dir->cd();
    hAcceptedEvents->Write();

    TDirectory *x_diffSubdir = x_dir->mkdir("Signal - Background");
    x_diffSubdir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        diffHist[h]->Write();

    TDirectory *x_dataSubdir = x_dir->mkdir("All Data");
    x_dataSubdir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        dataSum[h]->Write();

    TDirectory *x_mcSubdir = x_dir->mkdir("All MC");
    x_mcSubdir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        mcSum[h]->Write();

    TDirectory *x_sigSubdir = x_dir->mkdir("Signal MC");
    x_sigSubdir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        sigSum[h]->Write();

    TDirectory *x_bgSubdir = x_dir->mkdir("Background MC");
    x_bgSubdir->cd();
    for (unsigned h = 0; h < hname.size(); h++)
        bgSum[h]->Write();

    // Canvases and stacks
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

    // Yield macro
    buff << "    cout << endl;" << endl;
    buff << "}" << endl;
    log.AddLine(buff.str().c_str());
    file->cd();
    log.Write();

    file->Close();
// DELETE YOUR SHIT
}




void ApplyAlpha(float rgb[][3], const int size, const float alpha)
{
    for (unsigned i = 0; i < size; i++)
    {
        for (unsigned j = 0; j < 3; j++)
            rgb[i][j] += (1 - rgb[i][j]) * (1 - alpha);
    }
}
