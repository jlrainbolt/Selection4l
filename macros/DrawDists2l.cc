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
#include "Cuts2017.hh"
//#include "Cuts2016.hh"

using namespace std;


/*
**  DrawDists2l
**
**  Draws UNSCALED distributions for a "selected_" sample
*/ 

void DrawDists2l(const TString suffix)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 3;
    unsigned                   LL = 0,  MM = 1, EE = 2;     // Indices
    TString selection[N]    = {"ll",    "mumu", "ee"};
    unsigned chanIdx[N]     = {2,       3,      4};
    TString lepChan[N]      = {_l,      _mu,    _e};



    //
    //  OUTPUT FILE
    //

//  TString prefix  = "unscaled2l";
    TString prefix  = "rescaled2l";
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, Int_t, Double_t, Double_t>> v = {

        //          name        quantity            axis label          bins    xmin    xmax
        make_tuple( "nPV",      "nPV",              _nPV,               60,     0,      60),
        make_tuple( "met",      "met",              _met,               40,     0,      100),
                                                                        
        // Lab frame kinematics                                         
        make_tuple( "z1m",      "z1p4.M()",         _m_(_ll),           40,     80,     100),
        make_tuple( "z1pt",     "z1p4.Pt()",        _pT_(_ll),          40,     0,      80),
        make_tuple( "z1y",      "z1p4.Rapidity()",  "y_{"+_ll+"}",      40,     -2.5,   2.5),
        make_tuple( "z1pdg",    "z1pdg",            "",                 3,      10.5,   13.5),
                                                                    
        make_tuple( "l1pt",     "l1p4.Pt()",        _pT_(_l_(1)),       40,     20,     100),
        make_tuple( "l1eta",    "l1p4.Eta()",       _eta_(_l_(1)),      40,     -2.5,   2.5),
        make_tuple( "l1pdg",    "l1pdg",            "",                 27,     -13.5,  13.5),
                                                                    
        make_tuple( "l2pt",     "l2p4.Pt()",        _pT_(_l_(2)),       40,     10,     50),
        make_tuple( "l2eta",    "l2p4.Eta()",       _eta_(_l_(2)),      40,     -2.5,   2.5),
        make_tuple( "l2pdg",    "l2pdg",            "",                 27,     -13.5,  13.5),

        make_tuple( "dphi",     "fabs(l1p4.Phi()-l2p4.Phi())/3.141492654",
                                "|\\Delta\\phi_{"+_ll+"}|/\\pi",        40,     -0,     2)
    };



    //
    //  INPUT FILE
    //

//  TString inName  = "selected_" + suffix + ".root";
//  TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/" + inName;
    TString inName  = "rescaled_" + suffix + ".root";
    TString inPath  = inName;
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


    for (unsigned i = 1; i < N; i++)
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
            TString hname,  quantity,   xlabel;
            int     bins;
            float   xmin,   xmax;
//          TString weight = "weight";
            TString weight = "weight*rescale";
            tie(hname, quantity, xlabel, bins, xmin, xmax) = v[j];


            // Create and draw histogram
            TH1D *h = new TH1D(hname + "_" + suffix, quantity+" {"+weight+"}", bins, xmin, xmax);
            tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

            xlabel.ReplaceAll(_l, lepChan[i]);
            h->GetXaxis()->SetTitle(xlabel);
            h->Sumw2(kTRUE);
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
