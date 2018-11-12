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

void DrawResolution(const bool scale = kFALSE)
{

    //
    //  OPTIONS
    //

//  gErrorIgnoreLevel = kError;
    const double PI = TMath::Pi();



    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                   L4 = 0,  M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = {"4l",    "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = {5,       6,      7,      8,      9};
    TString lepChan[N]      = {_l,      _mu,    _l,     _l,     _e};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "resolution";
    TString suffix2 = scale ? "scaled" : "unscaled";
    TString outName = prefix + "_" + suffix2 + "_100bins.root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, Int_t, Double_t, Double_t>> v = {

        //          name            quantity            axis label          bins    xmin    xmax
        // Z rest frame kinematics
        make_tuple( "b_ttm",        "b_ttp4.M()",       _m_(_l_("2,3,4")),  100,    -10,    10),
        make_tuple( "b_l1p",        "b_l1v3.Mag()",     _p_(_l_(1)),        100,    -10,    10),
        make_tuple( "b_l2p",        "b_l2v3.Mag()",     _p_(_l_(2)),        100,    -10,    10),
        make_tuple( "b_l3p",        "b_l3v3.Mag()",     _p_(_l_(3)),        100,    -5,     5),
        make_tuple( "b_l4p",        "b_l4v3.Mag()",     _p_(_l_(4)),        100,    -5,     5),


        // Observables
        make_tuple( "psi",              "psi",              _psi,           100,    -1000,  1000),
        make_tuple( "sin_phi",          "sin_phi",          _sinphi,        100,    -0.2,   0.2),
        make_tuple( "theta_z1",         "theta_z1",         _theta_(_Z1),   100,    -PI/6., PI/6.),
        make_tuple( "theta_z2",         "theta_z2",         _theta_(_Z2),   100,    -PI/6., PI/6.),
        make_tuple( "angle_z1leps",     "angle_z1leps",     _alpha_(_Z1),   100,    -PI/8., PI/8.),
        make_tuple( "angle_z2leps",     "angle_z2leps",     _alpha_(_Z2),   100,    -PI/8., PI/8.),
        make_tuple( "angle_z1l2_z2",    "angle_z1l2_z2",    _beta,          100,    -PI/6., PI/6.)
    };



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
            TString hname,  quantity,   xlabel;
            int     bins;
            float   xmin,   xmax;
            TString weight = "weight";
            tie(hname, quantity, xlabel, bins, xmin, xmax) = v[j];
            
            // Add subtraction to quantity
            quantity = "gen_" + quantity + " - " + quantity;

            // Create and draw histogram
            TH1D *h = new TH1D(hname + "_" + suffix, quantity+" {"+weight+"}", bins, xmin, xmax);
            tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

            xlabel.ReplaceAll(_l, lepChan[i]);
            xlabel = "\\Delta " + xlabel;
            h->GetXaxis()->SetTitle(xlabel);
            h->Sumw2(kTRUE);
            h->Write();
        }

        cout << "done!" << endl;
    }



    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
