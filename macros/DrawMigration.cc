// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

// Custom
#include "Cuts2017.hh"

using namespace std;


/*
**  DrawMigration
**
**  Draws migration matrices and gen/reco distributions for a "matched_" sample
*/ 

void DrawMigration()
{

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

    TString prefix  = "migration";
    TString outName = prefix + "_zz_4l.root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, Int_t, Double_t, Double_t>> v = {

        //          name            quantity        axis label              bins    xmin    xmax
        // Z rest frame kinematics
//      make_tuple( "b_ttm",        "b_ttp4.M()",    _m_(_l_("2,3,4")),     15,     0,      60),
        make_tuple( "b_l1p",        "b_l1v3.Mag()",  _p_(_l_(1)),           10,     25,     50),


        // Observables
//      make_tuple( "psi",          "psi",              _psi,               20,     -5000,  5000),
//      make_tuple( "sin_phi",      "sin_phi",          _sinphi,            20,     -1,     1),
//      make_tuple( "cos_theta_z1", "cos_theta_z1",     _costheta_(_Z1),    10,     -1,     1),
//      make_tuple( "cos_theta_z2", "cos_theta_z2",     _costheta_(_Z2),    10,     -1,     1),
//      make_tuple( "cos_zeta_z1",  "cos_zeta_z1",      _coszeta_(_Z1),     10,     -1,     1),
//      make_tuple( "cos_zeta_z2",  "cos_zeta_z2",      _coszeta_(_Z2),     10,     -1,     1),
//      make_tuple( "angle_z1leps",
//                          "angle_z1leps/3.14",        _alpha_(_Z1),       12,     0,      1),
//      make_tuple( "angle_z2leps",
//                          "angle_z2leps/3.14",        _alpha_(_Z2),       10,     0,      1),
//      make_tuple( "angle_z1l2_z2",
//                          "angle_z1l2_z2/3.14",       _beta,              15,     0,      1)
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
            // Get parameters from tuple and form string
            TString hname,  quantity,   xlabel;
            int     bins;
            float   xmin,   xmax;
            TString weight = "weight";
            tie(hname, quantity, xlabel, bins, xmin, xmax) = v[j];

            xlabel.ReplaceAll(_l, lepChan[i]);
            

            // Create and draw histograms
            float width = (xmax - xmin) / bins;

            TH1D *h_reco = new TH1D(hname + "_reco", quantity, bins + 3, xmin - width, xmax + width);
            tree->Draw(quantity + ">>+" + hname + "_reco", weight);
            h_reco->GetXaxis()->SetTitle(xlabel);
            h_reco->Sumw2(kTRUE);
            outFile->cd(selection[i]);
            h_reco->Write();

            TH1D *h_gen = new TH1D(hname + "_gen", "gen_" + quantity, bins, xmin, xmax);
            tree->Draw("gen_" + quantity + ">>+" + hname + "_gen", weight);
//          h_gen->ClearUnderflowAndOverflow();
            h_gen->GetXaxis()->SetTitle(xlabel);
            h_gen->Sumw2(kTRUE);
            outFile->cd(selection[i]);
            h_gen->Write();


            TH2D *h_2d = new TH2D(hname + "_2d", quantity + ":gen_" + quantity,
                    bins + 3, xmin - width, xmax + width, bins, xmin, xmax);
            tree->Draw("gen_"+quantity + ":" + quantity + ">>+" + hname + "_2d", weight);
/*
            for (unsigned b = 0; b <= bins + 1; b++)
            {
                h_2d->SetBinContent(0, b, 0);
                h_2d->SetBinContent(bins + 1, b, 0);
            }
*/
            h_2d->GetXaxis()->SetTitle("reco");
            h_2d->GetYaxis()->SetTitle("gen");
            h_2d->Sumw2(kTRUE);
            outFile->cd(selection[i]);
            h_2d->Write();
        }

        cout << "done!" << endl;
    }



    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
