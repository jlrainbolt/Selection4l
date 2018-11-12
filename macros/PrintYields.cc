// STL
#include <iostream>
#include <fstream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"

// Cuts
#include "Cuts2017.hh"

using namespace std;


/*
**  PrintYields
**
**  Prints event yields for each channel.
*/

void PrintYields(bool useDY = kTRUE)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 8;   // Channel indices
    unsigned                LL = 0, MM = 1, EE = 2, L4 = 3, M4 = 4, ME = 5, EM = 6, E4 = 7;
    TString selection[N] = {"ll",   "mumu", "ee",   "4l",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N] = { 2,      3,      4,      5,      6,      7,      8,      9};






    ////
    ////
    ////    FILL HISTOGRAMS
    ////
    ////


    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/";
    TString prefix = "selected";



    //
    //  DATA
    //

    // Muon file
    TString muName = inPath + prefix + "_" + MU_SUFF + ".root";
    TFile *muFile = TFile::Open(muName);
    cout << "Opened " << muName << endl;

    TH1 *muHist;
    muFile->GetObject("SelectedEvents_" + MU_SUFF, muHist);
    muHist->SetDirectory(0);    
    muHist->Sumw2();
    muFile->Close();


    // Electron file
    TString elName = inPath + prefix + "_" + EL_SUFF + ".root";
    TFile *elFile = TFile::Open(elName);
    cout << "Opened " << elName << endl;

    TH1 *elHist;
    elFile->GetObject("SelectedEvents_" + EL_SUFF, elHist);
    elHist->SetDirectory(0);    
    elHist->Sumw2();
    elFile->Close();


    TH1 *dataHist = (TH1*) muHist->Clone("SelectedEvents_data_" + YEAR_STR);
    dataHist->Add(elHist);



    //
    //  MONTE CARLO
    //

    TH1 *mcHist[N_MC];
    for (unsigned j = 0; j < N_MC; j++) // sample loop
    {
        TString inName = inPath + prefix + "_" + MC_SUFF[j] + ".root";
        TFile *inFile = TFile::Open(inName);

        cout << "Opened " << inName << endl;

        TH1 *hist_;
        inFile->GetObject("SelectedEvents_" + MC_SUFF[j], hist_);
        hist_->SetDirectory(0);
        hist_->Sumw2();

        mcHist[j] = hist_;

        inFile->Close();
    }






    ////
    ////
    ////    CALCULATIONS
    ////
    ////


    //
    //  SCALING
    //

    // Make a histogram of the appropriate lumi for each channel
    TH1D *lumiHist = new TH1D("Lumi", "Lumi", 10, 0.5, 10.5);

    lumiHist->SetBinContent(chanIdx[MM], MUON_TRIG_LUMI);
    lumiHist->SetBinContent(chanIdx[M4], MUON_TRIG_LUMI);
    lumiHist->SetBinContent(chanIdx[ME], MUON_TRIG_LUMI);

    lumiHist->SetBinContent(chanIdx[EE], ELEC_TRIG_LUMI * ELEC_TRIG_SF);
    lumiHist->SetBinContent(chanIdx[EM], ELEC_TRIG_LUMI * ELEC_TRIG_SF);
    lumiHist->SetBinContent(chanIdx[E4], ELEC_TRIG_LUMI * ELEC_TRIG_SF);

    
    // Scale by appropriate lumi, etc.
    for (unsigned j = 0; j < N_MC; j++) // sample loop
    {
        float sf = 1000 * XSEC[j] / NGEN[j];
        mcHist[j]->Multiply(lumiHist);
        mcHist[j]->Scale(sf);
    }



    //
    //  SUMS
    //

    TH1 *mcTotal;
    mcTotal = (TH1*) mcHist[0]->Clone("ExpectedEvents");
    for (unsigned j = 1; j < N_MC; j++)
        mcTotal->Add(mcHist[j]);

    TH1 *sigBackground = (TH1*) mcHist[2]->Clone("BackgroundEvents");
    if (useDY)  // otherwise, omit DY background from signal calculation
        sigBackground->Add(mcHist[DY]);

    TH1 *dyBackground = (TH1*) mcHist[ZZ]->Clone("DYBackgroundEvents");
    dyBackground->Add(mcHist[2]);
    for (unsigned j = 3; j < N_MC; j++)
    {
        sigBackground->Add(mcHist[j]);
        dyBackground->Add(mcHist[j]);
    }
    for (unsigned i = MM; i < L4; i++)
        sigBackground->SetBinContent(chanIdx[i], 0);
    for (unsigned i = M4; i < N; i++)
        dyBackground->SetBinContent(chanIdx[i], 0);

    TH1 *obsMinusBg = (TH1*) dataHist->Clone("ObsMinusBgEvents");
    obsMinusBg->Add(sigBackground, -1);
    obsMinusBg->Add(dyBackground, -1);




    //
    //  LATEX
    //

    TString selLaTeX[N] = {"\\lplm", "\\mumu", "\\ee", "4\\ell", "4\\mu", "2\\mu2\\el","","4\\el"};
    for (unsigned i = MM; i < N; i++)
    {
        if (i == L4)
            continue;
        if (i == EM)
            continue;

        unsigned SIG = (i < L4) ? DY : ZZ;


        int nObserved = dataHist->GetBinContent(chanIdx[i]);
        float nExpected = mcTotal->GetBinContent(chanIdx[i]);
        float uExpected = mcTotal->GetBinError(chanIdx[i]);
        float nSignal = mcHist[SIG]->GetBinContent(chanIdx[i]);
        float uSignal = mcHist[SIG]->GetBinError(chanIdx[i]);
        float nBackground = (i > L4) ? sigBackground->GetBinContent(chanIdx[i])
                                     : dyBackground->GetBinContent(chanIdx[i]);
        float uBackground = (i > L4) ? sigBackground->GetBinError(chanIdx[i])
                                     : dyBackground->GetBinError(chanIdx[i]);

        if (i == ME)
        {
            nObserved += dataHist->GetBinContent(chanIdx[EM]);
            nExpected += mcTotal->GetBinContent(chanIdx[EM]);
            uExpected = sqrt(pow(uExpected, 2) + pow(mcTotal->GetBinError(chanIdx[EM]), 2));
            nSignal += mcHist[ZZ]->GetBinContent(chanIdx[EM]);
            uSignal = sqrt(pow(uSignal, 2) + pow(mcHist[ZZ]->GetBinError(chanIdx[EM]), 2));
            nBackground += sigBackground->GetBinContent(chanIdx[EM]);
            uBackground = sqrt(pow(uBackground,2)+pow(sigBackground->GetBinError(chanIdx[EM]),2));
        }
        if (i < L4)
        {
            nObserved = round(nObserved);
            nExpected = round(nExpected);       uExpected = round(uExpected);
            nSignal = round(nSignal);           uSignal = round(uSignal);
            nBackground = round(nBackground);   uBackground = round(uBackground);
        }

        TString texName = "2017_" + selection[i] + ".tex";

        ofstream texFile;
        texFile.open(texName);
        if (i > L4)
        {
            texFile.precision(2);
            texFile << fixed;
            texFile << "\\begin{tabular}{l l l S[table-format=3.2] l S[table-format=2.2]}" << endl; 
        }
        else
        {
            texFile.precision(8);
            texFile << "\\sisetup{group-minimum-digits=4}" << endl << endl;
            texFile << "\\begin{tabular}{l l l S[table-format=8.0] l S[table-format=5.0]}" << endl; 
        }

        texFile << "\\toprule" << endl << "\t\\multicolumn{3}{l}{Type} & ";
        texFile << "\\multicolumn{3}{l}{$N_{" << selLaTeX[i] << "}$ (events)} \\\\" << endl;
        texFile << "\\midrule" << endl;
        texFile << "\t\\multicolumn{3}{l}{Observed} & " << nObserved << " \\\\" << endl;
        texFile << "\\addlinespace\\addlinespace" << endl;
        texFile << "\t\\multicolumn{3}{l}{Expected} & ";
        texFile << nExpected << " & $\\pm$ & " << uExpected << " \\\\" << endl;
        texFile << "\\addlinespace" << endl << "\t& \\multicolumn{2}{l}{Signal} & ";
        texFile << nSignal << " & $\\pm$ & " << uSignal << " \\\\" << endl;
        texFile << "\\addlinespace" << endl << "\t& \\multicolumn{2}{l}{Background} & ";
        texFile << nBackground << " & $\\pm$ & " << uBackground << " \\\\" << endl;
        texFile << "\\addlinespace" << endl;

        for (unsigned j = ZZ; j < N_MC; j++)
        {
            if (j == SIG)
                continue;

            float nSample = mcHist[j]->GetBinContent(chanIdx[i]);
            float uSample = mcHist[j]->GetBinError(chanIdx[i]);

            if (i < L4)
            {
                nSample = round(nSample);       uSample = round(uSample);
            }

            texFile << "\t\t& & $" << MC_TEX[j] << "$ & ";
            texFile << nSample << " & $\\pm$ & " << uSample << " \\\\" << endl;
        }
        texFile << "\\bottomrule" << endl << "\\end{tabular}" << endl;

        texFile.close();
        cout << "Wrote LaTeX table to " << texName << endl;
    }
}
