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
**  CalculateBF
**
**  Calculates Z -> 4l branching fractions using yields from SelectedEvents histograms
*/

void CalculateBF(bool useDY = kFALSE)
{
    // Constants
    Double_t Zto2l = 0.03366, f_nr = 0.96;



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
    double nBackground[N]= {0,      0,      0,      0,      0,      0,      0,      0};
    double nObsMinusBg[N]= {0,      0,      0,      0,      0,      0,      0,      0};





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



    //
    //  GEN LEVEL
    //

    inPath = EOS_PATH + "/BLT/" + YEAR_STR + "/";
    prefix = "genHardProc";

    // Make phase space histogram by hand
    TH1D *psHist = new TH1D("PhaseSpace", "PhaseSpace", 10, 0.5, 10.5);

    // Start with Drell-Yan
    TString dyName = inPath + "gen_" + MC_SUFF[DY] + "/" + prefix + "_" + MC_SUFF[DY] + ".root";
    TFile *dyFile = TFile::Open(dyName);

    cout << "Opened " << dyName << endl;

    TH1 *dyHist;
    dyFile->GetObject("PhaseSpaceEvents_" + MC_SUFF[DY], dyHist);

    psHist->SetBinContent(chanIdx[MM], dyHist->GetBinContent(chanIdx[MM]));
    psHist->SetBinContent(chanIdx[EE], dyHist->GetBinContent(chanIdx[EE]));

    dyFile->Close();


    // Now with ZZ sample
    TString zzName = inPath + "gen_" + MC_SUFF[ZZ] + "/" + prefix + "_" + MC_SUFF[ZZ] + ".root";
    TFile *zzFile = TFile::Open(zzName);

    cout << "Opened " << zzName << endl;

    TH1 *zzHist;
    zzFile->GetObject("PhaseSpaceEvents_" + MC_SUFF[ZZ], zzHist);

    psHist->SetBinContent(chanIdx[M4], zzHist->GetBinContent(chanIdx[M4]));
    psHist->SetBinContent(chanIdx[E4], zzHist->GetBinContent(chanIdx[E4]));

    // Combine 2m2e channels
    float mixedPS = zzHist->GetBinContent(chanIdx[ME]);
    mixedPS += zzHist->GetBinContent(chanIdx[EM]);
    psHist->SetBinContent(chanIdx[ME], mixedPS);
    psHist->SetBinContent(chanIdx[EM], mixedPS);

    zzFile->Close();






    ////
    ////
    ////    CALCULATIONS
    ////
    ////

    
    //
    //  ACCEPTANCE * EFFICIENCY
    //

    // Make selected events histogram by hand
    TH1D *totalAE = new TH1D("AccXEff", "AccXEff", 10, 0.5, 10.5);

    totalAE->SetBinContent(chanIdx[MM], mcHist[DY]->GetBinContent(chanIdx[MM]));
    totalAE->SetBinContent(chanIdx[EE], mcHist[DY]->GetBinContent(chanIdx[EE]));
    totalAE->SetBinContent(chanIdx[M4], mcHist[ZZ]->GetBinContent(chanIdx[M4]));
    totalAE->SetBinContent(chanIdx[E4], mcHist[ZZ]->GetBinContent(chanIdx[E4]));
    float mixedSel = mcHist[ZZ]->GetBinContent(chanIdx[ME]);
    mixedSel += mcHist[ZZ]->GetBinContent(chanIdx[EM]);
    totalAE->SetBinContent(chanIdx[ME], mixedSel);
    totalAE->SetBinContent(chanIdx[EM], mixedSel);

    totalAE->Divide(psHist);



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
    //  PRINT
    //

    for (unsigned i = MM; i < N; i++)
    {
        if      (i < L4)
        {
            nObserved[i] = dataHist->GetBinContent(chanIdx[i]);
            nObserved[LL] += nObserved[i];

            nExpected[i] = mcTotal->GetBinContent(chanIdx[i]);
            nExpected[LL] += nExpected[i];

            nSignal[i] = mcHist[DY]->GetBinContent(chanIdx[i]);
            nSignal[LL] += nSignal[i];

            nBackground[i] = dyBackground->GetBinContent(chanIdx[i]);
            nBackground[LL] += nBackground[i];

            nObsMinusBg[i] = obsMinusBg->GetBinContent(chanIdx[i]);
            nObsMinusBg[LL] += nObsMinusBg[i];
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

            nObsMinusBg[i] = obsMinusBg->GetBinContent(chanIdx[i]);
            nObsMinusBg[L4] += nObsMinusBg[i];
        }
    }

    cout << endl << endl;
    cout << "BEFORE CORRECTIONS" << endl;
    cout << "\t\t" << "Observed" << "\t" << "Expected" << "\t" << "Signal  " << "\t";
    cout << "Background" << "\t" << "Obs - Bkg" << endl;
    for (unsigned i = 0; i < N; i++)
    {
        cout << selection[i] << "\t\t" << setw(8) << nObserved[i] << "\t";
        cout << setw(8) << nExpected[i] << "\t" << setw(8) << nSignal[i] << "\t";
        cout << setw(8) << nBackground[i] << "\t" << setw(8) << nObsMinusBg[i] << endl;

        if (i == EE)
            cout << endl;
    }






    ////
    ////
    ////    LUMI CORRECTION
    ////
    ////

    TH1D *lumiCorr = new TH1D("LumiCorr", "LumiCorr", 10, 0.5, 10.5);

    lumiCorr->SetBinContent(chanIdx[MM], ELEC_TRIG_LUMI / MUON_TRIG_LUMI);
    lumiCorr->SetBinContent(chanIdx[M4], ELEC_TRIG_LUMI / MUON_TRIG_LUMI);
    lumiCorr->SetBinContent(chanIdx[ME], ELEC_TRIG_LUMI / MUON_TRIG_LUMI);

    lumiCorr->SetBinContent(chanIdx[EE], 1);
    lumiCorr->SetBinContent(chanIdx[EM], 1);
    lumiCorr->SetBinContent(chanIdx[E4], 1);

    
    // Apply

    for (unsigned j = 0; j < N_MC; j++)
        mcHist[j]->Multiply(lumiCorr);

    mcTotal->Multiply(lumiCorr);
    dataHist->Multiply(lumiCorr);
    obsMinusBg->Multiply(lumiCorr);



    //
    //  PRINT
    //

    nObserved[LL] = 0;  nExpected[LL] = 0;  nSignal[LL] = 0;    nBackground[LL] = 0;
    nObsMinusBg[LL] = 0;
    nObserved[L4] = 0;  nExpected[L4] = 0;  nSignal[L4] = 0;    nBackground[L4] = 0;
    nObsMinusBg[L4] = 0;

    for (unsigned i = MM; i < N; i++)
    {
        if      (i < L4)
        {
            nObserved[i] = dataHist->GetBinContent(chanIdx[i]);
            nObserved[LL] += nObserved[i];

            nExpected[i] = mcTotal->GetBinContent(chanIdx[i]);
            nExpected[LL] += nExpected[i];

            nSignal[i] = mcHist[DY]->GetBinContent(chanIdx[i]);
            nSignal[LL] += nSignal[i];

            nBackground[i] = dyBackground->GetBinContent(chanIdx[i]);
            nBackground[LL] += nBackground[i];

            nObsMinusBg[i] = obsMinusBg->GetBinContent(chanIdx[i]);
            nObsMinusBg[LL] += nObsMinusBg[i];
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

            nObsMinusBg[i] = obsMinusBg->GetBinContent(chanIdx[i]);
            nObsMinusBg[L4] += nObsMinusBg[i];
        }
    }

    cout << endl << endl;
    cout << "CORRECTED FOR LUMI" << endl;
    cout << "\t\t" << "Observed" << "\t" << "Expected" << "\t" << "Signal  " << "\t";
    cout << "Background" << "\t" << "Obs - Bkg" << endl;
    for (unsigned i = 0; i < N; i++)
    {
        cout << selection[i] << "\t\t" << setw(8) << nObserved[i] << "\t";
        cout << setw(8) << nExpected[i] << "\t" << setw(8) << nSignal[i] << "\t";
        cout << setw(8) << nBackground[i] << "\t" << setw(8) << nObsMinusBg[i] << endl;

        if (i == EE)
            cout << endl;
    }






    ////
    ////
    ////    ACCEPTANCE * EFFICIENCY
    ////
    ////


    for (unsigned j = 0; j < N_MC; j++)
        mcHist[j]->Divide(totalAE);

    mcTotal->Divide(totalAE);
    dataHist->Divide(totalAE);
    obsMinusBg->Divide(totalAE);



    //
    //  PRINT
    //

    nObserved[LL] = 0;  nExpected[LL] = 0;  nSignal[LL] = 0;    nBackground[LL] = 0;
    nObsMinusBg[LL] = 0;
    nObserved[L4] = 0;  nExpected[L4] = 0;  nSignal[L4] = 0;    nBackground[L4] = 0;
    nObsMinusBg[L4] = 0;

    for (unsigned i = MM; i < N; i++)
    {
        if      (i < L4)
        {
            nObserved[i] = dataHist->GetBinContent(chanIdx[i]);
            nObserved[LL] += nObserved[i];

            nExpected[i] = mcTotal->GetBinContent(chanIdx[i]);
            nExpected[LL] += nExpected[i];

            nSignal[i] = mcHist[DY]->GetBinContent(chanIdx[i]);
            nSignal[LL] += nSignal[i];

            nBackground[i] = dyBackground->GetBinContent(chanIdx[i]);
            nBackground[LL] += nBackground[i];

            nObsMinusBg[i] = obsMinusBg->GetBinContent(chanIdx[i]);
            nObsMinusBg[LL] += nObsMinusBg[i];
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

            nObsMinusBg[i] = obsMinusBg->GetBinContent(chanIdx[i]);
            nObsMinusBg[L4] += nObsMinusBg[i];
        }
    }
    nObserved[ME] += nObserved[EM];
    nExpected[ME] += nExpected[EM];
    nSignal[ME] += nSignal[EM];
    nBackground[ME] += nBackground[EM];
    nObsMinusBg[ME] += nObsMinusBg[EM];

    cout << endl << endl;
    cout << "CORRECTED FOR ACC * EFF, LUMI" << endl;
    cout << "\t\t" << "Observed" << "\t" << "Expected" << "\t" << "Signal  " << "\t";
    cout << "Background" << "\t" << "Obs - Bkg" << endl;
    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)
            continue;

        cout << selection[i] << "\t\t" << setw(8) << nObserved[i] << "\t";
        cout << setw(8) << nExpected[i] << "\t" << setw(8) << nSignal[i] << "\t";
        cout << setw(8) << nBackground[i] << "\t" << setw(8) << nObsMinusBg[i] << endl;

        if (i == EE)
            cout << endl;
    }






    ////
    ////
    ////    BRANCHING FRACTIONS
    ////
    ////


    double fracUnc[N];

    for (unsigned i = 0; i < N; i++)
    {
        if ((i == LL) || (i == L4) || (i == EM))
            continue;

        fracUnc[i] = dataHist->GetBinError(chanIdx[i]) / dataHist->GetBinContent(chanIdx[i]);
    }
    fracUnc[LL] = sqrt(dataHist->GetBinError(chanIdx[MM]) * dataHist->GetBinError(chanIdx[MM])
                    + dataHist->GetBinError(chanIdx[EE]) * dataHist->GetBinError(chanIdx[EE]))
                    / nObsMinusBg[LL];
    fracUnc[L4] = 0;
    for (unsigned i = M4; i < N; i++)
    {
        if (i == EM)
            continue;
        fracUnc[L4] += dataHist->GetBinError(chanIdx[i]) * dataHist->GetBinError(chanIdx[i]);
    }
    fracUnc[L4] = sqrt(fracUnc[L4]) / nObsMinusBg[L4];


    double branchingFraction[N], bfFracUnc[N], bfUnc[N];
    for (unsigned i = L4; i < N; i++)
    {
        if (i == EM)
            continue;

        branchingFraction[i] = 2 * Zto2l * f_nr * nObsMinusBg[i] / nObsMinusBg[LL];
        bfFracUnc[i] = sqrt(fracUnc[i] * fracUnc[i] + fracUnc[LL] * fracUnc[LL]);
        bfUnc[i] = bfFracUnc[i] * branchingFraction[i];
    }

    cout << endl << endl;
    cout << "BRANCHING FRACTIONS" << endl;
    cout << "\t\t" << "Value" << "\t\t\t\t\t" << "Uncertainty" << endl;
    for (unsigned i = L4; i < N; i++)
    {
        if (i == EM)
            continue;

        cout << selection[i] << "\t\t" << setw(6) << branchingFraction[i];
        cout << " +- " << "\t" << bfUnc[i] << "\t\t";
        cout << bfFracUnc[i] << endl;
    }



    //
    //  LATEX
    //

    TString selLaTeX[N] = {"", "", "", "4\\ell", "4\\mu", "2\\mu2\\el", "", "4\\el"};
    TString texName = "2017_prelim.tex";

    ofstream texFile;
    texFile.open(texName);
    texFile.precision(3);

    texFile << "\\begin{tabular}{l l l l S}" << endl << "\\toprule" << endl;
    texFile << "\tChannel & \\multicolumn{3}{l}{$\\BF$ ($\\times 10^{-6}$)} & ";
    texFile << "\\multicolumn{1}{l}{Unc.~(\\%)} \\\\" << endl << "\\midrule" << endl;

//  texFile << fixed;
    for (unsigned i = L4; i < N; i++)
    {
        if (i == EM)
            continue;

        float val = branchingFraction[i] * 1000000;
        float unc = bfUnc[i] * 1000000;
        float frac = bfFracUnc[i] * 100;

        texFile << "\t$" << selLaTeX[i] << "$ & ";
        texFile << setw(2) << val << " & $\\pm$ & " << unc << " & " << frac;
        texFile << " \\\\";

        if (i == L4)
            texFile << " \\addlinespace";

        texFile << endl;
    }
    texFile << "\\bottomrule" << endl << "\\end{tabular}" << endl;

    texFile.close();
    cout << endl << endl << "Wrote LaTeX table to " << texName << endl;
}
