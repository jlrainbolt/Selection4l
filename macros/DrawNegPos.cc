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
**  DrawPosNeg
**
**  Draws correctly-weighted positive and negative sin(phi) histograms
*/

void DrawNumDenom(const TString suffix, bool useScaling = kTRUE)
{

    //
    //  SAMPLE INFO
    //

    const bool isData   = suffix.Contains(YEAR_STR);
    if (isData)
        useScaling = kFALSE;

    const unsigned N = 4;   // Channel indices
    unsigned                M4 = 2, ME = 3, EM = 4, E4 = 5;
    TString selection[N] = {"4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N] = { 6,      7,      8,      9};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "channels";
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    TH1D *hNumerator = new TH1D("Numerator_" + suffix, "Numerator", 5, 4.5, 9.5);
    hNumerator->Sumw2();
    TH1D *hDenominator = new TH1D("Denominator_" + suffix, "Denominator", 5, 4.5, 9.5);
    hDenominator->Sumw2();



    //
    //  SCALING
    //

    float xsec = 1, ngen = 1;

    if (useScaling) // find MC sample cross section
    {
        for (unsigned i = 0; i < N_MC; i++)
        {
            if (suffix.EqualTo(MC_SUFF[i]))
            {
                xsec = XSEC[i];
                ngen = NGEN[i];
                break;
            }
        }
    }



    //
    //  INPUT FILE
    //

    TString inName  = "selected_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;






    ////
    ////
    ////    CHANNEL LOOP
    ////
    ////


    for (unsigned i = 0; i < N; i++)
    {
        TTree *tree;
        inFile->GetObject(selection[i] + "_" + suffix, tree);

        cout << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;

        // Weight string
        float sf = lumi[i] * 1000 * xsec / ngen;
        TString scale = "";
        if (useScaling)
            scale.Form("weight*%f", sf);

        // Bin string
        TString bin;



        //
        //  DRAW
        //

        cout << "  Drawing histograms..." << endl;

        if      (i == MM)   // fill muon-triggered denominator channels
        {
            for (unsigned ii = M4; ii <= ME; i++)
            {
                bin.Form("%i", chanIdx[ii]);
                tree->Draw(bin + ">>+Denominator_" + suffix, scale);
            }
            tree->Draw("5>>+Denominator_" + suffix, scale);
        }
        else if (i == EE)   // fill electron-triggered denominator channels
        {
            for (unsigned ii = EM; ii <= E4; i++)
            {
                bin.Form("%i", chanIdx[ii]);
                tree->Draw(bin + ">>+Denominator_" + suffix, scale);
            }
            tree->Draw("5>>+Denominator_" + suffix, scale);
        }
        else    // fill appropriate numerator channel
        {
            bin.Form("%i", chanIdx[i]);
            tree->Draw(bin + ">>+Numerator_" + suffix, scale);
            tree->Draw("5>>+Numerator_" + suffix, scale);
        }
    }
    inFile->Close();



    //
    //  WRITE OUTPUT
    //

    outFile->cd();
    hNumerator->Write();
    hDenominator->Write();
    outFile->Close();

    cout << "Wrote histograms to " << outName << endl << endl << endl;
}
