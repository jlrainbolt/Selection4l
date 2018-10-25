// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
//#include "TError.h"

// Custom
#include "Cuts2017.hh"

using namespace std;


/*
**  DrawResolution
**
**  Draws resolution of distributions for a "matched_" sample
*/ 

void DrawResolution(const bool scale = kTRUE, const bool rebin = kFALSE)
{

    //
    //  OPTIONS
    //

//  gErrorIgnoreLevel = kError;
    const int nRebin = 5;



    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                   L4 = 0,  M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = {"4l",    "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = {5,       6,      7,      8,      9};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "resolution";
    TString suffix1 = scale ? "scaled" : "unscaled";
    TString suffix2 = rebin ? "" : "autobin";
    TString outName = prefix + "_" + suffix1 + "_" + suffix2 + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, Int_t, Double_t, Double_t>> v;

    //                      name            quantity            bins    xmin        xmax
    // Lab frame kinematics
    v.push_back(make_tuple( "zzm",          "zzp4.M()",         100,    80,         100));
    v.push_back(make_tuple( "zzpt",         "zzp4.Pt()",        100,    0,          100));

    v.push_back(make_tuple( "z1m",          "z1p4.M()",         100,    0,          100));
    v.push_back(make_tuple( "z1pt",         "z1p4.Pt()",        100,    0,          120));

    v.push_back(make_tuple( "z2m",          "z2p4.M()",         105,    0,          35));
    v.push_back(make_tuple( "z2pt",         "z2p4.Pt()",        100,    0,          75));

    v.push_back(make_tuple( "l1pt",         "l1p4.Pt()",        100,    0,          100));
    v.push_back(make_tuple( "l1eta",        "l1p4.Eta()",       100,    -2.5,       2.5));

    v.push_back(make_tuple( "l2pt",         "l2p4.Pt()",        120,    0,          60));
    v.push_back(make_tuple( "l2eta",        "l2p4.Eta()",       100,    -2.5,       2.5));

    v.push_back(make_tuple( "l3pt",         "l3p4.Pt()",        105,    0,          35));
    v.push_back(make_tuple( "l3eta",        "l3p4.Eta()",       100,    -2.5,       2.5));

    v.push_back(make_tuple( "l4pt",         "l4p4.Pt()",        100,    0,          25));
    v.push_back(make_tuple( "l4eta",        "l4p4.Eta()",       100,    -2.5,       2.5));


    // Z rest frame kinematics
    v.push_back(make_tuple( "b_z1p",        "b_z1p4.P()",       100,    0,          50));
    v.push_back(make_tuple( "b_ttm",        "b_ttp4.M()",       100,    0,          65));

    v.push_back(make_tuple( "b_l1p",        "b_l1v3.Mag()",     100,    25,         50));

    v.push_back(make_tuple( "b_l2p",        "b_l2v3.Mag()",     100,    15,         50));

    v.push_back(make_tuple( "b_l3p",        "b_l3v3.Mag()",     100,    0,          30));

    v.push_back(make_tuple( "b_l4p",        "b_l4v3.Mag()",     100,    0,          20));


    // Observables
    v.push_back(make_tuple( "psi",              "psi",              100,    -5000,  5000));
    v.push_back(make_tuple( "sin_phi",          "sin_phi",          100,    -1,     1));
    v.push_back(make_tuple( "theta_z1",         "theta_z1",         100,    0,      TMath::Pi()));
    v.push_back(make_tuple( "theta_z2",         "theta_z2",         100,    0,      TMath::Pi()));
    v.push_back(make_tuple( "angle_z1leps",     "angle_z1leps",     100,    0,      TMath::Pi()));
    v.push_back(make_tuple( "angle_z2leps",     "angle_z2leps",     100,    0,      TMath::Pi()));
    v.push_back(make_tuple( "angle_z1l2_z2",    "angle_z1l2_z2",    100,    0,      TMath::Pi()));



    //
    //  INPUT FILE
    //

    TString suffix  = "zz_4l";
    TString inName  = "matched_" + suffix + ".root";
    TString inPath  = HOME_PATH + "/Boosted/" + YEAR_STR + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;






    ////
    ////
    ////    CHANNEL LOOP
    ////
    ////

//  outFile->mkdir(selection[L4]);  // put 4l directory at the top

    for (unsigned i = 1; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        TTree *tree;
        inFile->GetObject(selection[i] + "_" + suffix, tree);

        cout << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;

        // Overall scale factor
        float LUMI; 
        if      (i == M4 || i == ME)
            LUMI = MUON_TRIG_LUMI;
        else if (i == E4 || i == EM)
            LUMI = ELEC_TRIG_LUMI;

        float sf = LUMI * 1000 * XSEC_ZZ_4L / NGEN_ZZ_4L;



        //
        //  DRAW
        //

        cout << "\t" << "Drawing histograms..." << flush;

        for (unsigned j = 0; j < v.size(); j++)
        {
            // Get parameters from tuple and form string
            TString hname,  quantity;
            int     bins;
            float   xmin,   xmax;
            TString weight = "weight";
            tie(hname, quantity, bins, xmin, xmax) = v[j];
            
            // Add subtraction to quantity
            quantity = "gen_" + quantity + " - " + quantity;

            TString varexp;
//          varexp.Form(quantity + ">>" + hname + "(%i,%g,%g)", bins, xmin, xmax);
            varexp.Form(quantity + ">>" + hname);


            // Draw and retrieve tree
            tree->Draw(varexp, weight);

            TH1F* hist;
            gDirectory->GetObject(hname, hist);

            if (hist)
                hist->SetName(hname + "_" + selection[i]);
            else    // the tree was empty and we need to make a dummy histogram
                hist = new TH1F(hname + "_" + suffix, quantity+" {"+weight+"}", bins, xmin, xmax);


            // Scaling, etc.
            if (scale)
                hist->Scale(sf);

            if (rebin)
                hist->Rebin(nRebin);


//          hist->Sumw2(kTRUE);
            hist->Write();
        }
        cout << "done!" << endl;
    }



    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
