// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

// Custom
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;


/*
**  DrawDists4l
**
**  Draws UNSCALED distributions for a "boosted_" sample
*/ 

void DrawDists4lCut(const TString suffix, const TString year, const TString cut = "",
        const bool isBkg = kFALSE)
{
    if (!year.EqualTo(YEAR_STR))
    {
        cout << "Wrong year in header file!" << endl;
        return;
    }

    //
    //  SAMPLE INFO
    //

    const unsigned N = 4;
    unsigned                   L4 = 0,  M4 = 1, ME = 2, E4 = 3;     // Indices
    TString selection[N]    = {"4l",    "4m",   "2m2e", "4e"};
    unsigned chanIdx[N]     = {5,       6,      7,      9};
    TString lepChan[N]      = {_l,      _mu,    _l,     _e};



    //
    //  OUTPUT FILE
    //

    TString prefix  = isBkg ? "bkg_" : "";
    TString outName = prefix + cut + "_" + year + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, TString, int, float, float>> v = {

        //          name        quantity         axis label     unit        bins    xmin    xmax
        make_tuple( "nPV",      "nPV",           _nPV,          _unit,      20,     0,      60),

        // Lab frame kinematics
        make_tuple( "zzm",      "zzp4.M()",      _m_(_4l),      _GeV,       20,     80,     100),
        make_tuple( "zzpt",     "zzp4.Pt()",     _pT_(_4l),     _GeV,       20,     0,      100),

        make_tuple( "z1m",      "z1p4.M()",      _m_(_Z1),      _GeV,       20,     12,     92),
        make_tuple( "z1pt",     "z1p4.Pt()",     _pT_(_Z1),     _GeV,       20,     0,      100),

        make_tuple( "z2m",      "z2p4.M()",      _m_(_Z2),      _GeV,       20,     4,      44),
        make_tuple( "z2pt",     "z2p4.Pt()",     _pT_(_Z2),     _GeV,       20,     0,      60),

        make_tuple( "l1pt",     "l1p4.Pt()",     _pT_(_l_(1)),  _GeV,       26,     20,     72),
        make_tuple( "l1eta",    "l1p4.Eta()",    _eta_(_l_(1)), _units,     20,     -2.5,   2.5),

        make_tuple( "l2pt",     "l2p4.Pt()",     _pT_(_l_(2)),  _GeV,       20,     10,     50),
        make_tuple( "l2eta",    "l2p4.Eta()",    _eta_(_l_(2)), _units,     20,     -2.5,   2.5),

        make_tuple( "l3pt",     "l3p4.Pt()",     _pT_(_l_(3)),  _GeV,       24,     5,      29),
        make_tuple( "l3eta",    "l3p4.Eta()",    _eta_(_l_(3)), _units,     20,     -2.5,   2.5),

        make_tuple( "l4pt",     "l4p4.Pt()",     _pT_(_l_(4)),  _GeV,       20,     5,      25),
        make_tuple( "l4eta",    "l4p4.Eta()",    _eta_(_l_(4)), _units,     20,     -2.5,   2.5)
    };



    //
    //  INPUT FILE
    //

    TString tag     = isBkg ? "background" : "selected";
    TString inName  = tag + "_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + year + "_trig/" + inName;
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
            TString weight = "weight * " + cut;

            tie(hname, quantity, xlabel, unit, bins, xmin, xmax) = v[j];

            if (selection[i].EqualTo("4m") && quantity.Contains("Eta"))
            {
                bins = 24;
                xmin = -2.4;
                xmax = 2.4;
            }


            // Create and draw histogram
            TH1D *h = new TH1D(hname + "_" + suffix, "", bins, xmin, xmax);
            h->Sumw2(kTRUE);
            h->SetBinErrorOption(TH1::kPoisson);

            tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

            xlabel.ReplaceAll(_l, lepChan[i]);

            h->GetXaxis()->SetTitle(xlabel);
            h->GetYaxis()->SetTitle(unit);
            h->SetStats(0);
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
