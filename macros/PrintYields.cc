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

    // Containers
    double nObserved[N] = { 0,      0,      0,      0,      0,      0,      0,      0};
    double nExpected[N] = { 0,      0,      0,      0,      0,      0,      0,      0};
    double nSignal[N] = {   0,      0,      0,      0,      0,      0,      0,      0};
    double nBackground[N]={ 0,      0,      0,      0,      0,      0,      0,      0};
    double nSample[j][N];






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

/*
 *
 *
 *              UGH
 *
 */

    for (unsigned i = MM; i < N; i++)
    {
        if      (i > L4)
        {
            if (i == EM)
            {
                nObserved[ME] += dataHist->GetBinContent(chanIdx[i]);
                nExpected[ME] += mcTotal->GetBinContent(chanIdx[i]);
                nSignal[ME] += mcHist[ZZ]->GetBinContent(chanIdx[i]);
                nBackground[ME] += dyBackground->GetBinContent(chanIdx[i]);
                for (unsigned j = 0; j < mcHist
            }
            else
        }
        else if (i > L4)
        {
            nObserved[i] = dataHist->GetBinContent(chanIdx[i]);
            nObserved[L4] += nObserved[i];

            nExpected[i] = mcTotal->GetBinContent(chanIdx[i]);
            nExpected[L4] += nExpected[i];

            nSignal[i] = mcHist[ZZ]->GetBinContent(chanIdx[i]);
            nSignal[L4] += nSignal[i];

            nBackground[i] = sigBackground->GetBinContent(chanIdx[i]);
            nBackground[L4] += nBackground[i];
        }

        if (i < L4)
        {
            nObserved[i] = dataHist->GetBinContent(chanIdx[i]);
            nObserved[LL] += nObserved[i];

            nExpected[i] = mcTotal->GetBinContent(chanIdx[i]);
            nExpected[LL] += nExpected[i];

            nSignal[i] = mcHist[DY]->GetBinContent(chanIdx[i]);
            nSignal[LL] += nSignal[i];

            nBackground[i] = dyBackground->GetBinContent(chanIdx[i]);
            nBackground[LL] += nBackground[i];
        }
    }



    //
    //  LATEX
    //

    TString selLaTeX[N] = {"\\lplm", "\\mumu", "\\ee", "4\\ell", "4\\mu", "2\\mu2\\el","","4\\el"};
    for (unsigned i = M4; i < N; i++)
    {
        if (i == EM)
            continue;


        int nObserved = dataHist->GetBinContent(chanIdx[i]);
        float nExpected = mcTotal->GetBinContent(chanIdx[i]);
        float uExpected = mcTotal->GetBinError(chanIdx[i]);
        float nSignal = mcHist[ZZ]->GetBinContent(chanIdx[i]);
        float uSignal = mcHist[ZZ]->GetBinError(chanIdx[i]);
        float nBackground = sigBackground->GetBinContent(chanIdx[i]);
        float uBackground = sigBackground->GetBinError(chanIdx[i]);

        TString texName = "2017_yields_" + selection[i] + ".tex";

        ofstream texFile;
        texFile.open(texName);
        texFile.precision(2);
        texFile << fixed;

        texFile << "\\begin{tabular}{l l l S l S}" << endl << "\\toprule" << endl;
        texFile << "\t\\multicolumn{3}{l}{Type} & ";
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

        for (unsigned j = DY; j < N_MC; j++)
        {
            float nSample = mcHist[j]->GetBinContent(chanIdx[i]);
            float uSample = mcHist[j]->GetBinError(chanIdx[i]);

            texFile << "\t\t& & $" << MC_TEX[j] << "$ & ";
            texFile << nSample << " & $\\pm$ & " << uSample << " \\\\" << endl;
        }
        texFile << "\\bottomrule" << endl << "\\end{tabular}" << endl;

        texFile.close();
        cout << endl << endl << "Wrote LaTeX table to " << texName << endl;
    }
}
