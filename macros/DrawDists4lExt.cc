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

void DrawDists4lExt(const TString suffix, const TString year, const bool isBkg = kFALSE)
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
    TString outName = prefix + "ext_" + year + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, TString, int, float, float>> v = {

        //          name        quantity         axis label     unit        bins    xmin    xmax

        // Lab frame kinematics
        make_tuple( "zzm_60_120",   "zzp4.M()",  _m_(_4l),      _GeV,       60,     60,     120),
        make_tuple( "zzm_70_110",   "zzp4.M()",  _m_(_4l),      _GeV,       40,     70,     110),
        make_tuple( "zzm_80_100",   "zzp4.M()",  _m_(_4l),      _GeV,       20,     80,     100),
        make_tuple( "zzm_75_105",   "zzp4.M()",  _m_(_4l),      _GeV,       30,     75,     105),
        make_tuple( "zzm_50_150",   "zzp4.M()",  _m_(_4l),      _GeV,       100,    50,     150),
        make_tuple( "zzm_50_250",   "zzp4.M()",  _m_(_4l),      _GeV,       100,    50,     250),
        make_tuple( "zzm_50_350",   "zzp4.M()",  _m_(_4l),      _GeV,       100,    50,     350),
        make_tuple( "zzm_50_550",   "zzp4.M()",  _m_(_4l),      _GeV,       100,    50,     550)
    };



    //
    //  INPUT FILE
    //

    TString tag     = isBkg ? "background" : "selected";
    TString inName  = tag + "_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Extended/" + year + "_new/" + inName;
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
