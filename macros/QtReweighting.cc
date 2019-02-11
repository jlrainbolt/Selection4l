// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

// Custom
#include "Cuts2017.hh"
//#include "Cuts2016.hh"

using namespace std;


/*
**  QtReweighting
**
**  Creates graphs to use for dimuon and dielectron pt reweighting
*/ 

void QtReweighting()
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 2;
    unsigned                    MM = 0,     EE = 1;     // Indices
    TString selection[N]    = { "mumu",     "ee"    };
    unsigned chanIdx[N]     = { 3,          4       };
    TString lepChan[N]      = { _mu,        _e      };



    //
    //  INPUT FILES
    //

    TString prefix  = "selected";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/";
    cout << endl << endl;

    // Muon data
    TString muName  = prefix + "_" + MU_SUFF + ".root";
    TFile *muFile   = TFile::Open(inPath + muName);
    cout << "Opened " << inPath + muName << endl;

    // Electron data
    TString elName  = prefix + "_" + EL_SUFF + ".root";
    TFile *elFile   = TFile::Open(inPath + elName);
    cout << "Opened " << inPath + elName << endl;

    // Drell-Yan MC
    TString dyName  = prefix + "_" + MC_SUFF[DY] + ".root";
    TFile   *dyFile = TFile::Open(inPath + dyName);
    cout << "Opened " << inPath << dyName << endl << endl;



    //
    //  OUTPUT FILE
    //

    TString outName = "qt_weights_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");
    outFile->cd();



    //
    //  HISTOGRAMS
    //

    // Binning

    const unsigned  P = 10;
    float   binwidth[P] = { 1,  2,  5,  10, 20, 50, 100,200,400,1000};
    int     nbins[P]    = { 10, 5,  4,  2,  2,  2,  2,  1,  1,  1   };

    const unsigned  M = 29;
    double xval = 0, xbins[M+1];
    unsigned k = 0;             // counter for bin index
    for (unsigned i = 0; i < P; i++)    // loop over widths
    {
        for (unsigned j = 0; j < nbins[i]; j++)
        {
            xbins[k] = xval;
            xval += binwidth[i];
            k++;
        }
    }
    const unsigned B = xbins[M] - xbins[0];


    // Draw

    TH1D *dataGeV[N], *mcGeV[N];

    for (unsigned i = 0; i < N; i++)
    {
        // 1-GeV binning
        dataGeV[i]  = new TH1D(selection[i] + "_data_GeV",  "", B, xbins[0], xbins[M]);
        dataGeV[i]->Sumw2();
        mcGeV[i]    = new TH1D(selection[i] + "_mc_GeV",    "", B, xbins[0], xbins[M]);
        mcGeV[i]->Sumw2();

        // Data
        TTree *tree;
        if      (i == MM)
            muFile->GetObject(selection[i] + "_" + MU_SUFF, tree);
        else if (i == EE)
            elFile->GetObject(selection[i] + "_" + EL_SUFF, tree);

        cout << "Data " << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;
        cout << "\t" << "Drawing histograms..." << flush;
        
        tree->Draw("z1p4.Pt()>>+" + selection[i] + "_data_GeV", "");
        cout << "done!" << endl;


        // MC
        dyFile->GetObject(selection[i] + "_" + MC_SUFF[DY], tree);

        cout << "MC " << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;
        cout << "\t" << "Drawing histograms..." << flush;

        tree->Draw("z1p4.Pt()>>+" + selection[i] + "_mc_GeV", "weight/qtWeight");
        cout << "done!" << endl;
    }
    muFile->Close();
    elFile->Close();
    dyFile->Close();


    // Create ratio histogram

    TH1 *data[N], *mc[N];
    TH1D *ratio[N];

    for (unsigned i = 0; i < N; i++)
    {
        data[i] = dataGeV[i]->Rebin(M, selection[i] + "_data",  xbins);
        mc[i]   = mcGeV[i]->Rebin(  M, selection[i] + "_mc",    xbins);

        float dataScale = 1. / data[i]->Integral();
        data[i]->Scale(dataScale);
        dataGeV[i]->Scale(dataScale);

        float mcScale = 1. / mc[i]->Integral();
        mc[i]->Scale(mcScale);
        mcGeV[i]->Scale(mcScale);

        ratio[i] = new TH1D(selection[i] + "_ratio", "", M, xbins);
        ratio[i]->Divide(data[i], mc[i], 1, 1);
    }



    //
    //  CREATE GRAPHS
    //

    // Get x values: corrected bin centers

    float x[N][M];

    for (unsigned i = 0; i < N; i++)
    {
//      x[N][0] = 0;

        for (unsigned j = 0; j < M; j++)
        {
            float binContent = data[i]->GetBinContent(j+1) + mc[i]->GetBinContent(j+1);
            int binWidth = ratio[i]->GetBinWidth(j+1); 

            int k = mcGeV[i]->FindBin(mc[i]->GetBinLowEdge(j+1) + 0.5);

            // Find 1 GeV bin where center lies
            for (unsigned l = 0; l <= binWidth; l++)
            {
                float tempContent = dataGeV[i]->Integral(k, k+l) + mcGeV[i]->Integral(k, k+l);

                if (tempContent >= 0.5 * binContent)
                {
                    x[i][j] = dataGeV[i]->GetBinCenter(k+l);
                    break;
                }
            }
        }
    }


    // Get other values from histogram

    float y[N][M], exl[N][M], exh[N][M], ey[N][M];

    for (unsigned i = 0; i < N; i++)
    {
        for (unsigned j = 0; j <= M; j++)
        {
            y[i][j] = ratio[i]->GetBinContent(j+1);
            ey[i][j] = ratio[i]->GetBinError(j+1);

            exl[i][j] = x[i][j] - xbins[j];
            exh[i][j] = xbins[j+1] - x[i][j];
        }
    }


    // Finally, create graph objects

    TGraphAsymmErrors *graph[N];

    for (unsigned i = 0; i < N; i++)
    {
        graph[i] = new TGraphAsymmErrors(M, x[i], y[i], exl[i], exh[i], ey[i], ey[i]);
        graph[i]->SetName(selection[i] + "_weight");
    }





    //
    //  WRITE OUTPUT
    //

    for (unsigned i = 0; i < N; i++)
    {
        data[i]->Write();
        mc[i]->Write();
        ratio[i]->Write();
        graph[i]->Write();
    }
    outFile->Close();
    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
