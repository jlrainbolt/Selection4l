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
#include "TError.h"

// Custom
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;


/*
**  DrawDists2l
**
**  Draws UNSCALED distributions for a "selected_" sample
*/ 

void DrawDists2l(const TString suffix, const TString year)
{
    if (!year.EqualTo(YEAR_STR))
    {
        cout << "Wrong year in header file!" << endl;
        return;
    }

    //
    //  SAMPLE INFO
    //

    const unsigned N = 2;
    unsigned                   MM = 0, EE = 1;     // Indices
    TString selection[N]    = {"mumu", "ee"};
    unsigned chanIdx[N]     = {3,      4};
    TString lepChan[N]      = {_mu,    _e};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "2l";
    TString outName = prefix + "_" + YEAR_STR + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, TString, int, float, float>> v = {

        //          name        quantity        x label         unit        bins    xmin    xmax
        make_tuple( "nPV",      "nPV",          _nPV,           _unit,      60,     0,      60),
                                                                    
        // Lab frame kinematics                                     
        make_tuple( "z1m",      "z1p4.M()",     _m_(_ll),       _GeV,       40,     80,     100),
        make_tuple( "z1pt",     "z1p4.Pt()",    _pT_(_ll),      _GeV,       40,     0,      80),
        make_tuple( "z1y",      "z1p4.Rapidity()",_y_(_ll),     _units,     40,     -2.5,   2.5),
                                                                
        make_tuple( "l1pt",     "l1p4.Pt()",    _pT_(_l_(1)),   _GeV,       40,     20,     100),
        make_tuple( "l1eta",    "l1p4.Eta()",   _eta_(_l_(1)),  _units,     40,     -2.5,   2.5),
                                                                
        make_tuple( "l2pt",     "l2p4.Pt()",    _pT_(_l_(2)),   _GeV,       40,     10,     50),
        make_tuple( "l2eta",    "l2p4.Eta()",   _eta_(_l_(2)),  _units,     40,     -2.5,   2.5),

        make_tuple( "dphi",     "fabs(l1p4.Phi()-l2p4.Phi())/3.141492654",
                                                _dphi_(_ll),    _pirad,     40,     -0,     2)
    };



    //
    //  INPUT FILE
    //

    TString inName  = "selected_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "_v1/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;


    // Get all "Events" histograms
    vector<TH1*> hEvents;
    TKey *key;
    TIter next(inFile->GetListOfKeys());
    while ((key = (TKey*) next()))
    {
        TString hname = key->GetName();
        if (hname.Contains("Events"))
        {
            TH1 *hist;
            inFile->GetObject(hname, hist);
            hist->SetDirectory(0);
            hEvents.push_back(hist);

            cout << "Got histogram " << hname << endl;
        }
    }






    ////
    ////
    ////    CHANNEL LOOP
    ////
    ////


    for (unsigned i = 0; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        TTree *tree;
        inFile->GetObject(selection[i] + "_" + suffix, tree);

        cout << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;



        //
        //  DRAW
        //

        cout << "\t" << "Drawing histograms..." << flush;

        for (unsigned j = 0; j < v.size(); j++)
        {
            // Get parameters from tuple
            TString hname,  quantity,   xlabel, unit;
            int     bins;
            float   xmin,   xmax;
            TString weight = "weight";

            tie(hname, quantity, xlabel, unit, bins, xmin, xmax) = v[j];


            // Create and draw histogram
            TH1D *h = new TH1D(hname + "_" + suffix, quantity+" {"+weight+"}", bins, xmin, xmax);
            h->Sumw2(kTRUE);
            h->SetBinErrorOption(TH1::kPoisson);

            tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

            xlabel.ReplaceAll(_l, lepChan[i]);

            h->GetXaxis()->SetTitle(xlabel);
            h->GetYaxis()->SetTitle(unit);
            h->Write();
        }

        cout << "done!" << endl;

    }

    outFile->cd();
    for (unsigned k = 0; k < hEvents.size(); k++)
        hEvents[k]->Write();
    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
