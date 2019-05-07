// STL
#include <vector>
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// Cuts
#include "Cuts2017.hh"
//#include "Cuts2012.hh"

using namespace std;



void CalculateECAL()
{


    //
    //  CONTAINERS
    //

    const unsigned N = 2;    // Number of l+l- trees
    TString sel_2l[N] = {"mumu", "ee"};

    float   sel_4l_gen = 0, sel_4l_nom = 0, sel_4l_up = 0,  sel_4l_dn = 0;
    float   sel_2l_gen = 0, sel_2l_nom = 0, sel_2l_up = 0,  sel_2l_dn = 0;



    //
    //  FOUR-LEPTON
    //

    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/selected_zz_4l.root";
    TFile   *zzFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    TTreeReader reader("4l_zz_4l", zzFile);

    TTreeReaderValue    <Float_t>   genWeight_      (reader,    "genWeight");
    TTreeReaderValue    <Float_t>   ecalWeight_     (reader,    "ecalWeight");
    TTreeReaderValue    <Float_t>   ecalWeightUp_   (reader,    "ecalWeightUp");
    TTreeReaderValue    <Float_t>   ecalWeightDown_ (reader,    "ecalWeightDown");

    cout << "Loaded branches" << endl;


    // Event loop

    int nEvents = reader.GetEntries(kTRUE);

    cout << "Running over " << nEvents << " total 4l events" << endl;
    cout << endl;

    while (reader.Next())
    {
        float genWeight = *genWeight_;

        float nomWeight = genWeight * *ecalWeight_;
        float upWeight = genWeight * *ecalWeightUp_;
        float dnWeight = genWeight * *ecalWeightDown_;

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
            float genWeight = *genWeight_;

            float nomWeight = genWeight * *ecalWeight_;
            float upWeight = genWeight * *ecalWeightUp_;
            float dnWeight = genWeight * *ecalWeightDown_;

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
