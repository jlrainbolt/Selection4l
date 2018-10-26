// STL
#include <iostream>

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

void CalculateBF()
{
    // Constants
    Double_t Zto2l = 0.03366, f_nr = 0.96;

    //
    //  SAMPLE INFO
    //

    // or whatever the heck this is

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
    ////    GET HISTOGRAMS
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

    TH1 *psHist[2];
    for (unsigned j = ZZ; j <= DY; j++)    // loop over only signal & Drell-Yan
    {
        TString inName = inPath + "gen_" + MC_SUFF[j] + "/" + prefix + "_" + MC_SUFF[j] + ".root";
        TFile *inFile = TFile::Open(inName);

        cout << "Opened " << inName << endl;

        TH1 *hist_;
        inFile->GetObject("PhaseSpaceEvents_" + MC_SUFF[j], hist_);
        hist_->SetDirectory(0);
        hist_->Sumw2();

        psHist[j] = hist_;

        inFile->Close();
    }






    ////
    ////
    ////    CALCULATIONS
    ////
    ////

    
    //
    //  ACCEPTANCE * EFFICIENCY
    //

    TH1 *aeHist[2];
    for (unsigned j = ZZ; j <= DY; j++)    // loop over only signal & Drell-Yan
    {
        // Divide selected by phase space yields
        aeHist[j] = (TH1*) mcHist[j]->Clone("AcceptanceXEfficiency");
        aeHist[j]->Divide(psHist[j]);
    }
    TH1 *totalAE = (TH1*) aeHist[ZZ]->Clone();
    totalAE->Add(aeHist[DY]);



    //
    //  SCALING
    //

    // Make a histogram of the appropriate lumi for each channel...
    TH1D *lumiHist = new TH1D("Lumi", "Lumi", 10, 0.5, 10.5);

    lumiHist->SetBinContent(chanIdx[MM], MUON_TRIG_LUMI);
    lumiHist->SetBinContent(chanIdx[M4], MUON_TRIG_LUMI);
    lumiHist->SetBinContent(chanIdx[ME], MUON_TRIG_LUMI);

    lumiHist->SetBinContent(chanIdx[EE], ELEC_TRIG_LUMI);
    lumiHist->SetBinContent(chanIdx[EM], ELEC_TRIG_LUMI);
    lumiHist->SetBinContent(chanIdx[E4], ELEC_TRIG_LUMI);

    
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

    TH1 *sigBackground = (TH1*) mcHist[DY]->Clone("BackgroundEvents");
    TH1 *dyBackground = (TH1*) mcHist[ZZ]->Clone("DYBackgroundEvents");
    for (unsigned j = 2; j < N_MC; j++)
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

    cout << endl << endl;
    cout << "CORRECTED FOR ACC & EFF" << endl;
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
    ////    BRANCHING FRACTIONS
    ////
    ////


    double fracUnc[N];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == LL)
            continue;
        if (i == L4)
            continue;

        fracUnc[i] = obsMinusBg->GetBinError(chanIdx[i]) / obsMinusBg->GetBinContent(chanIdx[i]);
    }
    fracUnc[LL] = sqrt(obsMinusBg->GetBinError(chanIdx[MM]) * obsMinusBg->GetBinError(chanIdx[MM])
                    + obsMinusBg->GetBinError(chanIdx[EE]) * obsMinusBg->GetBinError(chanIdx[EE]))
                    / nObsMinusBg[LL];
    fracUnc[L4] = 0;
    for (unsigned i = M4; i < N; i++)
        fracUnc[L4] += obsMinusBg->GetBinError(chanIdx[i]) * obsMinusBg->GetBinError(chanIdx[i]);
    fracUnc[L4] = sqrt(fracUnc[L4]) / nObsMinusBg[L4];


    double branchingFraction[N], bfFracUnc[N], bfUnc[N];
/*
    branchingFraction[L4] = 2. * Zto2l * nObsMinusBg[L4] / nObsMinusBg[LL];
    bfFracUnc[L4] = sqrt(fracUnc[L4] * fracUnc[L4] + fracUnc[LL] * fracUnc[LL]);
    bfUnc[L4] = bfFracUnc[L4] * branchingFraction[L4];
*/
    for (unsigned i = M4; i <= ME; i++)
    {
        branchingFraction[i] = Zto2l * nObsMinusBg[i] / nObsMinusBg[MM];
        bfFracUnc[i] = sqrt(fracUnc[i] * fracUnc[i] + fracUnc[MM] * fracUnc[MM]);
        bfUnc[i] = bfFracUnc[i] * branchingFraction[i];
    }
    for (unsigned i = EM; i <= E4; i++)
    {
        branchingFraction[i] = Zto2l * nObsMinusBg[i] / nObsMinusBg[EE];
        bfFracUnc[i] = sqrt(fracUnc[i] * fracUnc[i] + fracUnc[EE] * fracUnc[EE]);
        bfUnc[i] = bfFracUnc[i] * branchingFraction[i];
    }

    branchingFraction[L4] = 0;
    bfUnc[L4] = 0;
    for (unsigned i = M4; i < N; i++)
    {
        branchingFraction[L4] += branchingFraction[i];
        bfUnc[L4] += bfUnc[i] * bfUnc[i];
    }
    bfUnc[L4] = sqrt(bfUnc[L4]);
    bfFracUnc[L4] = bfUnc[L4] / branchingFraction[L4];

    cout << endl << endl;
    cout << "BRANCHING FRACTIONS" << endl;
    cout << "\t\t" << "Value" << "\t\t\t\t\t" << "Uncertainty" << endl;
    for (unsigned i = L4; i < N; i++)
    {
        cout << selection[i] << "\t\t" << setw(6) << branchingFraction[i];
        cout << " +- " << "\t" << bfUnc[i] << "\t\t";
        cout << bfFracUnc[i] << endl;
    }
}
