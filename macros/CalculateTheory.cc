// STL
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

// Cuts
#include "Cuts2017.hh"
//#include "Cuts2012.hh"

using namespace std;



void CalculateTheory()
{


    //
    //  CONTAINERS
    //

    const unsigned  N_QCD = 9,  N_PDF = 100;
    float           fid_4l_nom = 0, gen_4l_nom = 0, fid_2l_nom = 0, gen_2l_nom = 0;
    float           fid_4l[N_QCD],  gen_4l[N_QCD],  fid_2l[N_QCD],  gen_2l[N_QCD];

    for (unsigned i = 0; i < N_QCD; i++)
    {
        fid_4l[i] = 0;
        gen_4l[i] = 0;
        fid_2l[i] = 0;
        gen_2l[i] = 0;
    }
    



    //
    //  FOUR-LEPTON
    //

//  TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/gen_zz_4l_0.root";
    TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/gen_zjets_m-50_0.root";
    TFile   *inFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;



    //
    //  INPUT BRANCHES
    //

//  TTreeReader reader("tree_zz_4l", inFile);
    TTreeReader reader("tree_zjets_m-50", inFile);

    TTreeReaderValue    <Bool_t>                isFiducial_ (reader,    "isFiducial");
    TTreeReaderValue    <Float_t>               genWeight_  (reader,    "genWeight");
//  TTreeReaderValue    <vector<UShort_t>>      pdfID_      (reader,    "pdfID");
    TTreeReaderValue    <vector<Float_t>>       pdfWeight_  (reader,    "pdfWeight");
//  TTreeReaderValue    <vector<UShort_t>>      qcdID_      (reader,    "qcdID");
    TTreeReaderValue    <vector<Float_t>>       qcdWeight_  (reader,    "qcdWeight");

    cout << "Loaded branches" << endl;





    ////
    ////
    ////    EVENT LOOP
    ////
    ////


    int nEvents = reader.GetEntries(kTRUE);

    cout << endl;
    cout << "Running over " << nEvents << " total phase space events" << endl;
    cout << endl;

    while (reader.Next())
    {
        if (reader.GetCurrentEntry() % 10000 == 0)
        {
            cout << "Processed " << reader.GetCurrentEntry() << " of " << nEvents << " events";
            cout << endl;
        }
        
        gen_4l_nom += *genWeight_;

        if (*isFiducial_)
            fid_4l_nom += *genWeight_;

        for (unsigned i = 0; i < N_QCD; i++)
        {
            gen_4l[i] += (*qcdWeight_)[i];

            if (*isFiducial_)
                fid_4l[i] += (*qcdWeight_)[i];
        }
    }
    inFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  ACCEPTANCE
    //

    float   acc_4l_nom = fid_4l_nom / gen_4l_nom,   acc_2l_nom = 0; //fid_2l_nom / gen_2l_nom;
    float   acc_4l[N_QCD],  acc_2l[N_QCD];

    for (unsigned i = 0; i < N_QCD; i++)
    {
        acc_4l[i] = fid_4l[i] / gen_4l[i];
//      acc_2l[i] = fid_2l[i] / gen_2l[i];
    }



    cout << "Accepted " << fid_4l_nom << " of " << gen_4l_nom << " nominal events (";
    cout << acc_4l_nom * 100. << "%)" << endl << endl;

    for (unsigned i = 0; i < N_QCD; i++)
    {
        cout << "Accepted " << fid_4l[i] << " of " << gen_4l[i] << " events (";
        cout << acc_4l[i] * 100. << "%)" << endl;
    }
    cout << endl << endl;
}
