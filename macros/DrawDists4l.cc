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
//#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
#include "Cuts2012.hh"

using namespace std;


/*
**  DrawDists4l
**
**  Draws UNSCALED distributions for a "boosted_" sample
*/ 

void DrawDists4l(const TString suffix, const TString year)
{
    if (!year.EqualTo(YEAR_STR))
    {
        cout << "Wrong year in header file!" << endl;
        return;
    }

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

    TString prefix  = "4l";
    TString outName = prefix + "_" + year + "_" + suffix + ".root";
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

        make_tuple( "z1m",      "z1p4.M()",      _m_(_Z1),      _GeV,       16,     12,     92),
        make_tuple( "z1pt",     "z1p4.Pt()",     _pT_(_Z1),     _GeV,       16,     0,      100),

        make_tuple( "z2m",      "z2p4.M()",      _m_(_Z2),      _GeV,       10,     4,      36),
        make_tuple( "z2pt",     "z2p4.Pt()",     _pT_(_Z2),     _GeV,       10,     0,      60),

        make_tuple( "l1pt",     "l1p4.Pt()",     _pT_(_l_(1)),  _GeV,       20,     0,      100),
        make_tuple( "l1eta",    "l1p4.Eta()",    _eta_(_l_(1)), _units,     10,     -2.5,   2.5),

        make_tuple( "l2pt",     "l2p4.Pt()",     _pT_(_l_(2)),  _GeV,       20,     0,      100),
        make_tuple( "l2eta",    "l2p4.Eta()",    _eta_(_l_(2)), _units,     10,     -2.5,   2.5),

        make_tuple( "l3pt",     "l3p4.Pt()",     _pT_(_l_(3)),  _GeV,       20,     0,      100),
        make_tuple( "l3eta",    "l3p4.Eta()",    _eta_(_l_(3)), _units,     10,     -2.5,   2.5),

        make_tuple( "l4pt",     "l4p4.Pt()",     _pT_(_l_(4)),  _GeV,       20,     0,      100),
        make_tuple( "l4eta",    "l4p4.Eta()",    _eta_(_l_(4)), _units,     10,     -2.5,   2.5),
 
        // Z rest frame kinematics                                                
        make_tuple( "b_ttm",    "b_ttp4.M()",    _m_(_l_("2,3,4")), _GeV,   11,     5,      60),
                                
        make_tuple( "b_l1p",    "b_l1v3.Mag()",  _p_(_l_(1)),   _GeV,       10,     25,     50),
        make_tuple( "b_l2p",    "b_l2v3.Mag()",  _p_(_l_(2)),   _GeV,       10,     15,     45),
        make_tuple( "b_l3p",    "b_l3v3.Mag()",  _p_(_l_(3)),   _GeV,       10,     4,      24),
        make_tuple( "b_l4p",    "b_l4v3.Mag()",  _p_(_l_(4)),   _GeV,       10,     0,      20),
                        
                        
        // Observables      
        make_tuple( "psi",          "psi",          _psi,           "",     20,     -5000,  5000),
        make_tuple( "phi",      "phi/3.141592654",  _phi,           _pirad, 20,     -1,     1),
        make_tuple( "cos_phi",      "cos_phi",      _cosphi,        _units, 20,     -1,     1),
        make_tuple( "sin_phi",      "sin_phi",      _sinphi,        _units, 20,     -1,     1),
        make_tuple( "sin_phi_2",    "sin_phi",      _sinphi,        _units, 2,      -1,     1),
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

//  TString inName  = "boosted_bkg_" + suffix + ".root";
    TString inName  = "boosted_" + suffix + ".root";
    TString inPath  = HOME_PATH + "/Boosted/" + year + "/" + inName;
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
            TString weight = "weight/trigWeight/qtWeight";

            tie(hname, quantity, xlabel, unit, bins, xmin, xmax) = v[j];

            if      (suffix.EqualTo("zjets_m-50"))
                weight = "weight/trigWeight";
            else if (suffix.EqualTo("phase_space"))
            {
                weight = "weight";
//              bins = 2 * bins;
            }


            // Create and draw histogram
            TH1D *h = new TH1D(hname + "_" + suffix, "", bins, xmin, xmax);
            tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

            xlabel.ReplaceAll(_l, lepChan[i]);

            h->GetXaxis()->SetTitle(xlabel);
            h->GetYaxis()->SetTitle(unit);
            h->Sumw2(kTRUE);
            h->SetStats(0);
            h->Write();

//          if (suffix == "phase_space")
//          {
//              h->SetFillColor(lLightBlue);
//              h->SetLineColor(lLightBlue);

//              TCanvas *c = new TCanvas("c_" + hname + "_" + suffix, "", 800, 800);
//              c->cd();
//              h->Draw("HIST");
//              c->Write();
//          }
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
