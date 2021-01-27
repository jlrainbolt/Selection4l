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

void DrawDiffDists(const TString suffix, const TString year, const bool isBkg = kFALSE)
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
    vector<TString> weight  = {"weight", "wtMuonIDUp", "wtMuonIDDn", "wtElecIDUp", "wtElecIDDn",
                                "wtElecRecoUp", "wtElecRecoDn"};



    //
    //  OUTPUT FILE
    //

    TString prefix  = isBkg ? "bkg_" : "";
    TString outName = "ddr_" + prefix + year + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, TString, int, float, float>> v = {

        //          name        quantity         axis label     unit        bins    xmin    xmax
        make_tuple( "b_z1m",    "b_z1p4.M()",    _m_(_Z1),      _GeV,       10,     12,     92),
        make_tuple( "b_z2m",    "b_z2p4.M()",    _m_(_Z2),      _GeV,       11,     4,      37),
        make_tuple( "b_ttm",    "b_ttp4.M()",    _m_(_l_("2,3,4")), _GeV,   11,     5,      60),

        make_tuple( "b_l1p",    "b_l1v3.Mag()",  _p_(_l_(1)),   _GeV,       10,     25,     50),

        make_tuple( "phi",      "phi/3.141592654",  _phi,           _pirad, 20,     -1,     1),
        make_tuple( "cos_phi",      "cos_phi",      _cosphi,        _units, 20,     -1,     1),
        make_tuple( "sin_phi",      "sin_phi",      _sinphi,        _units, 20,     -1,     1),
        make_tuple( "cos_theta_z1", "cos_theta_z1", _costheta_(_Z1),_units, 10,     -1,     1),
        make_tuple( "cos_theta_z2", "cos_theta_z2", _costheta_(_Z2),_units, 10,     -1,     1),
        make_tuple( "angle_z1leps", "angle_z1leps/3.141592654",
                                                    _alpha_(_Z1),   _pirad, 10,     0,      1),
        make_tuple( "angle_z2leps", "angle_z2leps/3.141592654",
                                                    _alpha_(_Z2),   _pirad, 10,     0,      1),
        make_tuple( "angle_z1l2_z2","angle_z1l2_z2/3.141592654",
                                                    _beta,          _pirad, 10,     0,      1)
    };



    //
    //  INPUT FILE
    //

    TString tag     = isBkg ? "boosted_bkg" : "boosted";
    TString inName  = tag + "_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Boosted/" + year + "_wt2/" + inName;
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


    for (unsigned w = 0; w < weight.size(); w++)
    {
        outFile->mkdir(weight[w]);

            for (unsigned i = 0; i < N; i++)
            {
                TString thisDir = weight[w] + "/" + selection[i];
                outFile->mkdir(thisDir);
                outFile->cd(thisDir);

                TTree *tree;
                inFile->GetObject(selection[i] + "_" + suffix, tree);

                cout << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;



                //
                //  DRAW
                //

                cout << "\t" << "Drawing '" << weight[w] << "' histograms..." << flush;

                for (unsigned j = 0; j < v.size(); j++)
                {
                    // Get parameters from tuple
                    TString hname,  quantity,   xlabel, unit;
                    int     bins;
                    float   xmin,   xmax;

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

                    tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight[w]);

                    xlabel.ReplaceAll(_l, lepChan[i]);

                    h->GetXaxis()->SetTitle(xlabel);
                    h->GetYaxis()->SetTitle(unit);
                    h->SetStats(0);
                    h->Write();
                }

                cout << "done!" << endl;

            }
    }

    outFile->cd();
    for (unsigned k = 0; k < hEvents.size(); k++)
        hEvents[k]->Write();
    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
