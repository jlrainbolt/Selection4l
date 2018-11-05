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
**  CalculateAsymmetry
**
**  Calculates triple product asymmetry observable using sin(phi) from unweighted4l histograms
*/

void CalculateAsymmetry(bool useDY = kFALSE)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;   // Channel indices
    unsigned                L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;
    TString selection[N] = {"4l",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N] = { 5,      6,      7,      8,      9};


    TH1 *hObserved[N], *hExpected[N], *hBackground[N], *hObsMinusBg[N];






    ////
    ////
    ////    FILL HISTOGRAMS
    ////
    ////


//  TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/";
    TString inPath = "";
    TString prefix = "unscaled4l";



    //
    //  DATA
    //

    // Muon file
    TString muName = inPath + prefix + "_" + MU_SUFF + ".root";
    TFile *muFile = TFile::Open(muName);
    cout << "Opened " << muName << endl;

    // Electron file
    TString elName = inPath + prefix + "_" + EL_SUFF + ".root";
    TFile *elFile = TFile::Open(elName);
    cout << "Opened " << elName << endl;

    for (unsigned i = 1; i < N; i++)
    {
        TH1 *hist_;
        if      ((i == M4) || (i == ME))
            muFile->GetObject(selection[i] + "/sin_phi_2_" + MU_SUFF, hist_);
        else if ((i == E4) || (i == EM))
            elFile->GetObject(selection[i] + "/sin_phi_2_" + EL_SUFF, hist_);

        hist_->SetDirectory(0);
        hist_->Sumw2();
        hist_->SetName("sin_phi_2_data");

        hObserved[i] = hist_;
    }
    muFile->Close();
    elFile->Close();



    //
    //  MONTE CARLO
    //

    TH1 *mcHist[N][N_MC];
    for (unsigned j = 0; j < N_MC; j++) // sample loop
    {
        TString inName = inPath + prefix + "_" + MC_SUFF[j] + ".root";
        TFile *inFile = TFile::Open(inName);

        cout << "Opened " << inName << endl;

        for (unsigned i = 1; i < N; i++)
        {
            TH1 *hist_;
            inFile->GetObject(selection[i] + "/sin_phi_2_" + MC_SUFF[j], hist_);

            hist_->SetDirectory(0);
            hist_->Sumw2();

            // Scale by appropriate lumi, etc.
            float LUMI;
            if      ((i == M4) || (i == ME))
                LUMI = MUON_TRIG_LUMI;
            else if ((i == E4) || (i == EM))
                LUMI = ELEC_TRIG_LUMI;
            float sf = LUMI * 1000 * XSEC[j] / NGEN[j];
            hist_->Scale(sf);

            mcHist[i][j] = hist_;
        }
        inFile->Close();
    }


/*
    //
    //  ACCEPTANCE * EFFICIENCY
    //

    TH1 *aeHist[N];
    TString aeName = inPath + "acc_x_eff.root";
    TFile *aeFile = TFile::Open(aeName);

    cout << "Opened " << aeName << endl;

    for (unsigned i = 1; i < N; i++)
    {
        TH1 *hist_;
        aeFile->GetObject(selection[i] + "/sin_phi_2_" + selection[i], hist_);

        hist_->SetDirectory(0);
        hist_->Sumw2();
        hist_->SetName("sin_phi_2_axe");

        aeHist[i] = hist_;
    }
    aeFile->Close();
*/





    ////
    ////
    ////    CALCULATIONS
    ////
    ////

 
    //
    //  SUMS
    //

    // 4l
    hObserved[0] = (TH1*) hObserved[1]->Clone();
    for (unsigned i = 2; i < N; i++)
        hObserved[0]->Add(hObserved[i]);
    for (unsigned j = 0; j < N_MC; j++)
    {
        mcHist[0][j] = (TH1*) mcHist[1][j]->Clone();
        for (unsigned i = 2; i < N; i++)
            mcHist[0][j]->Add(mcHist[i][j]);
    }

    // 2m2e
    hObserved[ME]->Add(hObserved[EM]);
    for (unsigned j = 0; j < N_MC; j++)
        mcHist[ME][j]->Add(mcHist[EM][j]);

    // Expected
    for (unsigned i = 0; i < N; i++)
    {
        hExpected[i] = (TH1*) mcHist[i][0]->Clone("sin_phi_2_exp");
        for (unsigned j = 1; j < N_MC; j++)
            hExpected[i]->Add(mcHist[i][j]);
    }

    // Signal, background
    for (unsigned i = 0; i < N; i++)
    {
        hBackground[i] = (TH1*) mcHist[i][1]->Clone("sin_phi_2_bkg");
        for (unsigned j = 2; j < N_MC; j++)
            hBackground[i]->Add(mcHist[i][j]);

        // Observed minus background
        hObsMinusBg[i] = (TH1*) hObserved[i]->Clone("sin_phi_2_omb");
        hObsMinusBg[i]->Add(hBackground[i], -1);
    }



    //
    //  PRINT
    //

    cout << endl << endl;
    cout << "BEFORE CORRECTIONS" << endl;
    for (unsigned b = 1; b <= 2; b++)
    {
        if (b == 1)
            cout << "NEGATIVE" << endl;
        else
            cout << "POSITIVE" << endl;
        cout << "\t\t" << "Observed" << "\t" << "Expected" << "\t" << "Signal  " << "\t";
        cout << "Background" << "\t" << "Obs - Bkg" << endl;
        for (unsigned i = 0; i < N; i++)
        {
            if (i == EM)
                continue;

            cout << selection[i] << "\t\t";
            cout << setw(8) << hObserved[i]->GetBinContent(b) << "\t";
            cout << setw(8) << hExpected[i]->GetBinContent(b) << "\t";
            cout << setw(8) << mcHist[i][0]->GetBinContent(b) << "\t";
            cout << setw(8) << hBackground[i]->GetBinContent(b) << "\t";
            cout << setw(8) << hObsMinusBg[i]->GetBinContent(b) << endl;
        }
    }




/*

    ////
    ////
    ////    ACCEPTANCE * EFFICIENCY
    ////
    ////


    // Scale by appropriate acc * eff histogram
    for (unsigned i = 1; i < N; i++)
    {
        for (unsigned j = 0; j < N_MC; j++)
            mcHist[i][j]->Divide(aeHist[i]);
        hObserved[i]->Divide(aeHist[i]);
    }

    // Repeat previous calculations...
    hObserved[0] = (TH1*) hObserved[1]->Clone();
    for (unsigned i = 2; i < N; i++)
        hObserved[0]->Add(hObserved[i]);
    for (unsigned j = 0; j < N_MC; j++)
    {
        mcHist[0][j] = (TH1*) mcHist[1][j]->Clone();
        for (unsigned i = 2; i < N; i++)
            mcHist[0][j]->Add(mcHist[i][j]);
    }

    for (unsigned i = 0; i < N; i++)
    {
        hExpected[i] = (TH1*) mcHist[i][0]->Clone("sin_phi_2_exp");
        for (unsigned j = 1; j < N_MC; j++)
            hExpected[i]->Add(mcHist[i][j]);

        // Choose whether DY is counted among "background" MC samples 
        unsigned idx1, idx2;
        if (useDY)
        {
            idx1 = 1;
            idx2 = 2;
        }
        else
        {
            idx1 = 2;
            idx2 = 3;
        }
        hBackground[i] = (TH1*) mcHist[i][idx1]->Clone("sin_phi_2_bkg");
        for (unsigned j = idx2; j < N_MC; j++)
            hBackground[i]->Add(mcHist[i][j]);

        // Observed minus background
        hObsMinusBg[i] = (TH1*) hObserved[i]->Clone("sin_phi_2_omb");
        hObsMinusBg[i]->Add(hBackground[i], -1);
    }



    //
    //  PRINT
    //

    cout << endl << endl;
    cout << "CORRECTED FOR ACC & EFF" << endl;
    for (unsigned b = 1; b <= 2; b++)
    {
        if (b == 1)
            cout << "NEGATIVE" << endl;
        else
            cout << "POSITIVE" << endl;
        cout << "\t\t" << "Observed" << "\t" << "Expected" << "\t" << "Signal  " << "\t";
        cout << "Background" << "\t" << "Obs - Bkg" << endl;
        for (unsigned i = 0; i < N; i++)
        {
            cout << selection[i] << "\t\t";
            cout << setw(8) << hObserved[i]->GetBinContent(b) << "\t";
            cout << setw(8) << hExpected[i]->GetBinContent(b) << "\t";
            cout << setw(8) << mcHist[i][0]->GetBinContent(b) << "\t";
            cout << setw(8) << hBackground[i]->GetBinContent(b) << "\t";
            cout << setw(8) << hObsMinusBg[i]->GetBinContent(b) << endl;
        }
    }
*/





    ////
    ////
    ////    ASYMMETRIES
    ////
    ////


    float asymmetry[N], unc[N], fracUnc[N];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)
            continue;

        float neg = hObsMinusBg[i]->GetBinContent(1), negUnc = hObsMinusBg[i]->GetBinError(1);
        float pos = hObsMinusBg[i]->GetBinContent(2), posUnc = hObsMinusBg[i]->GetBinError(2);

        float sum = pos + neg,      diff = pos - neg;
        asymmetry[i] = diff / sum;

        unc[i] = sqrt(4*pos*neg / pow(sum, 3));
        fracUnc[i] = unc[i] / asymmetry[i];
    }

    cout << endl << endl;
    cout << "ASYMMETRIES" << endl;
    cout << "\t\t" << "Value" << "\t\t\t\t\t" << "Fractional Uncertainty" << endl;
    for (unsigned i = L4; i < N; i++)
    {
        if (i == EM)
            continue;

        cout << selection[i] << "\t\t" << setw(6) << asymmetry[i];
        cout << " +- " << "\t" << unc[i] << "\t\t";
        cout << fracUnc[i] << endl;
    }
}
