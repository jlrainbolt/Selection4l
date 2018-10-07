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



void stackResults(const TString inFile, const TString selection)
{
    // Choose selection type
    Bool_t sel2l = kFALSE, sel4l = kFALSE, muPair1, muPair2;
    TString sigSuffix, selName;
    if (selection == "mumu" || selection == "2m" || selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;      selName = "DY";     sigSuffix = "zjets_m-50";
    }
    else if (selection == "2m2e" || selection == "2mu2e" || selection == "2e2m" || selection== "2e2mu"
             || selection == "4mu" || selection == "4m" || selection == "4e" || selection == "4l")
    {
        sel4l = kTRUE;      selName = "4l";     sigSuffix = "zz_4l";
    }
    TFile *file = TFile::Open(inFile, "UPDATE");


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

    if (selection == "4l")
    {
        for (unsigned i = 0; i < dataColor.size(); i++)
            dataColor[i] /= 2;
        for (unsigned i = 0; i < mcColor.size(); i++)
            mcColor[i] /= 4;
        for (unsigned i = 0; i < xsec.size(); i++)
            xsec[i] /= 4;
        lumi = 41.529 * 2 + 36.735 * 2;    // hTotalEvents got hadded over the four channels...
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
    vector<TString> hname, hxtitle; //, hytitle;
    if (sel2l)
    {
        hname = {"nPV", "met", "z1m", "z1pt",
                    "l1pt", "l1eta", "l1iso", "l1pdg", "l2pt", "l2eta", "l2iso", "l2pdg",
                    "TotalEvents", "AcceptedEvents"};

        hxtitle = {"n_{\\mbox{PV}}", _met, _m(_ll), _pT(_ll),
                    _pT(1), _eta(1), _iso+"_{1}", "id 1",
                    _pT(2), _eta(2), _iso+"_{2}", "id 2",
                    "", ""};
/*
        hytitle = {_EventsPer(_unit), _EventsPer(.5, _GeV),
                    _EventsPer(.2, _GeV), _EventsPer(_GeV),
                    _EventsPer(_GeV), _EventsPer(.05, _rad), "", "",
                    _EventsPer(_GeV), _EventsPer(.05, _rad), "", "",
                    "", ""};
*/
    }
    else if (sel4l)
    {
/*
        hname = {"zzm", "zzpt", "l1pt_d", "z2m_d", "TotalEvents", "AcceptedEvents"};
        hxtitle = {_m(_4l), _pT(_4l), _pT(1), _m(_Z2), "", ""};
*/
        hname = {"nPV", "met",
                    "zzm", "zzpt", "z1m", "z1pt", "z2m", "z2pt",
                    "l1pt", "l1eta", "l1iso", "l1pdg", "l2pt", "l2eta", "l2iso", "l2pdg",
                    "l3pt", "l3eta", "l3iso", "l3pdg", "l4pt", "l4eta", "l4iso", "l4pdg",
                    // "l1pt_d", "z2m_d", "m3l_d", "angle_d",
                    "TotalEvents", "AcceptedEvents"};
        hxtitle = {"N_{\\mbox{PV}}", _met, 
                    _m(_4l), _pT(_4l), _m(_Z1), _pT(_Z1), _m(_Z2), _pT(_Z2), 
                    _pT(1), _eta(1), _iso+"_{1}", "id 1",
                    _pT(2), _eta(2), _iso+"_{2}", "id 2",
                    _pT(3), _eta(3), _iso+"_{3}", "id 3",
                    _pT(4), _eta(4), _iso+"_{4}", "id 4",
                    // _pT(1), _m(_Z2), _m("234"), "\\theta_{"+_p4(2)+_and+_p4(_Z2)+"}",
                    "", ""};
/*
        hytitle = {_EventsPer(_unit), _EventsPer(.5, _GeV),
                    _EventsPer(_GeV), _EventsPer(10, _GeV),
                    _EventsPer(5, _GeV), _EventsPer(10, _GeV),
                    _EventsPer(2.5, _GeV), _EventsPer(5, _GeV),
                    _EventsPer(10, _GeV), _EventsPer(.2, _rad), "", "",
                    _EventsPer(5, _GeV), _EventsPer(.2, _rad), "", "",
                    _EventsPer(5, _GeV), _EventsPer(.2, _rad), "", "",
                    _EventsPer(5, _GeV), _EventsPer(.2, _rad), "", "",
                    "", "",
                    _EventsPer(10, _GeV), _EventsPer(5, _GeV), _EventsPer(10, _GeV),
                    _EventsPer(.2, _rad)};
*/
    }

    vector<THStack*> dataStack, mcStack;
    vector<TH1*> dataSum, mcSum, sigSum, bgSum;
    vector<TCanvas*> canvas;
    for (unsigned h = 0; h < hname.size(); h++)
    {
        canvas.push_back(new TCanvas(hname[h], "", lCanvasSize, lCanvasSize));
        dataStack.push_back(new THStack(hname[h], ""));
        mcStack.push_back(new THStack(hname[h], ""));
    }



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
            //  hist->Sumw2(kTRUE);
            hist->SetMarkerColor(dataColor[i]);
            hist->SetMarkerStyle(kFullCircle);  hist->SetMarkerSize(2); 
            hist->SetLineColor(dataColor[i]);   hist->SetLineWidth(2);

            if (dataSum.size() <= h)
                dataSum.push_back((TH1*) hist->Clone(hname[h]));
            else
                dataSum[h]->Add(hist);

            if (selection != "4l")
                dataStack[h]->Add(hist);
        }
    }

    if (selection == "4l")
    {
        for (unsigned h = 0; h < hname.size(); h++)
            dataStack[h]->Add(dataSum[h]);
    }

    buff << "    cout << endl;" << endl;

    for (unsigned j_ = 0; j_ < mcSubdir.size(); j_++)
    {
        unsigned j = mcSubdir.size() - 1 - j_;

        TH1* hTotalEvents;
        mcSubdir[j]->GetObject("TotalEvents_" + mcSuffix[j], hTotalEvents);
        Float_t ngen = hTotalEvents->GetBinContent(1) - 2 * hTotalEvents->GetBinContent(10);
        Float_t weight = lumi * 1000. * xsec[j] / ngen;

        buff << "    cout << \"" << setw(16) << mcSuffix[j] << "\\t";
        buff << setw(10) << weight * hTotalEvents->GetBinContent(7) << "\\t+- ";
        buff << setw(10) << weight * hTotalEvents->GetBinError(7) << "\" << endl;" << endl;

        cout << ngen << endl;
        cout << lumi << endl;
        cout << xsec[j] << endl;
        cout << weight << endl;

        for (unsigned h = 0; h < hname.size(); h++)
        {
            TH1* hist;
            mcSubdir[j]->GetObject(hname[h] + "_" + mcSuffix[j], hist);
            hist->SetFillColor(mcColor[j]); hist->SetLineColor(mcColor[j]);
            //  hist->Sumw2(kTRUE);
            hist->Scale(weight);    mcStack[h]->Add(hist);

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


    // Fill legend with dummy histograms
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


    // Draw on canvases
    for (unsigned h = 0; h < hname.size(); h++)
    {
        Facelift(canvas[h]);
        canvas[h]->cd();
/*
        // ROOT sucks so you gotta draw the legend before changing its coordinates 
        legend->Draw();
        Facelift(legend);
        gPad->Modified();
*/
        dataStack[h]->Draw("E");

        // Also stacks must be drawn in order for their axes to exist (Rene Brun, WHY??)
        if (dataStack[h]->GetMaximum() > mcStack[h]->GetMaximum())
        {
            dataStack[h]->Draw("E");
            Facelift(dataStack[h]);
            dataStack[h]->GetXaxis()->SetTitle(hxtitle[h]);
//          dataStack[h]->GetYaxis()->SetTitle(hytitle[h]);
            gPad->Modified();
        
            mcStack[h]->Draw("HIST SAME");
            Facelift(mcStack[h]);
            gPad->Modified();

            dataStack[h]->Draw("E SAME");
        }
        else
        { 
            mcStack[h]->Draw("HIST");
            Facelift(mcStack[h]);
            mcStack[h]->GetXaxis()->SetTitle(hxtitle[h]);
//          mcStack[h]->GetYaxis()->SetTitle(hytitle[h]);
            gPad->Modified();

            dataStack[h]->Draw("E SAME");
            Facelift(dataStack[h]);
            gPad->Modified();
        }

        gPad->RedrawAxis();

        legend->Draw();
        gPad->Modified();
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
    {
        canvas[h]->Write();
//      canvas[h]->SaveAs("../2016_12a/" + selection + "_" + hname[h] + ".png");
    }

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
