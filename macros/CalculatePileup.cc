// STL
#include <vector>
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// Cuts
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2012.hh"

using namespace std;



void CalculatePileup()
{


    //
    //  CONTAINERS
    //

    const unsigned N = 2;    // Number of l+l- trees
    TString sel_2l[N] = {"mumu", "ee"};

    float   sel_4l_gen = 0, sel_4l_nom = 0, sel_4l_up = 0,  sel_4l_dn = 0;
    float   sel_2l_gen = 0, sel_2l_nom = 0, sel_2l_up = 0,  sel_2l_dn = 0;



    //
    //  HISTOGRAMS
    //

    TH1 *h_nom, *h_up, *h_dn;

    TString inName  = "pu_weights_" + YEAR_STR + ".root";
    TString inPath  = BLT_PATH + "/BLTAnalysis/data/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    inFile->GetObject("weights", h_nom);
    inFile->GetObject("weights_varUp", h_up);
    inFile->GetObject("weights_varDn", h_dn);
    h_nom->SetDirectory(0);
    h_up->SetDirectory(0);
    h_dn->SetDirectory(0);

    inFile->Close();
    cout << "Closed " << inPath << endl;



    //
    //  FOUR-LEPTON
    //

    inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/selected_zz_4l.root";
    TFile   *zzFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    TTreeReader reader("4l_zz_4l", zzFile);

    TTreeReaderValue    <Float_t>               genWeight_      (reader,    "genWeight");
    TTreeReaderValue    <Float_t>               nPU_            (reader,    "nPU");

    cout << "Loaded branches" << endl;


    // Event loop

    int nEvents = reader.GetEntries(kTRUE);

    cout << "Running over " << nEvents << " total 4l events" << endl;
    cout << endl;

    while (reader.Next())
    {
        float nPU = *nPU_,  genWeight = *genWeight_;

        float nomWeight = genWeight * h_nom->GetBinContent(h_nom->FindBin(nPU));
        float upWeight = genWeight * h_up->GetBinContent(h_up->FindBin(nPU));
        float dnWeight = genWeight * h_dn->GetBinContent(h_dn->FindBin(nPU));

//      if (fabs(dnWeight) < fabs(upWeight))
//          swap(upWeight, dnWeight);

        sel_4l_gen += genWeight;    sel_4l_nom += nomWeight;
        sel_4l_up += upWeight;      sel_4l_dn += dnWeight;
    }

    zzFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  DILEPTON
    //

    // Open one file

    inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/selected_zjets_m-50.root";
    TFile   *dyFile = TFile::Open(inPath);

    float mod = 100;

    if (YEAR_STR == "2018")
        mod = 50;

    cout << "Opened " << inPath << endl;

    for (unsigned i = 0; i < N; i++)
    {
        reader.SetTree(sel_2l[i] + "_zjets_m-50", dyFile);
        reader.Restart();

        cout << "Loaded branches" << endl;


        // Event loop

        int nEvents = reader.GetEntries(kTRUE) / mod;

        cout << endl;
        cout << "Running over " << nEvents << " total " << sel_2l[i] << " events" << endl;
        cout << endl;

        while (reader.Next() && (reader.GetCurrentEntry() < nEvents))
        {
            float nPU = *nPU_,  genWeight = *genWeight_;

            float nomWeight = genWeight * h_nom->GetBinContent(h_nom->FindBin(nPU));
            float upWeight = genWeight * h_up->GetBinContent(h_up->FindBin(nPU));
            float dnWeight = genWeight * h_dn->GetBinContent(h_dn->FindBin(nPU));

            sel_2l_gen += genWeight;    sel_2l_nom += nomWeight;
            sel_2l_up += upWeight;      sel_2l_dn += dnWeight;
        }
    }
    dyFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  PRINTOUT
    //

    cout << "year = " << YEAR_STR << endl << endl;
    cout << "%\t" << "Efficiency" << "\t\t\t" << "4l" << "\t\t\t" << "2l" << endl << endl;

    cout << setprecision(16);

    cout << "eff_nom = [\t" << sel_4l_gen / sel_4l_nom << "\t" << sel_2l_gen / sel_2l_nom << endl;
    cout << "\t" << "];" << endl << endl;

    cout << "eff_up = [\t" << sel_4l_gen / sel_4l_up << "\t" << sel_2l_gen / sel_2l_up << endl;
    cout << "\t" << "];" << endl << endl;

    cout << "eff_dn = [\t" << sel_4l_gen / sel_4l_dn << "\t" << sel_2l_gen / sel_2l_dn << endl;
    cout << "\t" << "];" << endl << endl;
}
