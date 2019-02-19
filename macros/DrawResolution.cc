// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"

// Custom
#include "Cuts2016.hh"
//#include "Cuts2017.hh"

using namespace std;


/*
**  DrawResolution
**
**  Draws resolution of distributions for a "matched_" sample
*/ 

void DrawResolution()
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                   L4 = 0,  M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = {"4l",    "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = {5,       6,      7,      8,      9};
    TString lepChan[N]      = {_l,      _mu,    _l,     _l,     _e};
    int lColor[N]           = {lPurple, lBlue,  lPurple, lPurple, lRed};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "resolution";
    TString outName = prefix + "_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, int, float, float, float>> v = {

        //          name            quantity            axis label          bins xmin   xmax  width
        // Z rest frame kinematics
        make_tuple( "b_ttm",        "b_ttp4.M()",       _m_(_l_("2,3,4")),  100, -10,   10,     5),
        make_tuple( "b_l1p",        "b_l1v3.Mag()",     _p_(_l_(1)),        100, -10,   10,   2.5),
//      make_tuple( "b_l2p",        "b_l2v3.Mag()",     _p_(_l_(2)),        100, -10,   10,     2),
//      make_tuple( "b_l3p",        "b_l3v3.Mag()",     _p_(_l_(3)),        100, -5,    5,      2),
//      make_tuple( "b_l4p",        "b_l4v3.Mag()",     _p_(_l_(4)),        100, -5,    5,    1.5),


        // Observables
//      make_tuple( "psi",          "psi",              _psi,               100, -1000, 1000, 500),
        make_tuple( "sin_phi",      "sin_phi",          _sinphi,            100, -0.25, 0.25, 0.1),
        make_tuple( "sin_phi_2",    "sin_phi",          _sinphi,            100, -0.25, 0.25,   1),
        make_tuple( "cos_theta_z1", "cos_theta_z1",     _costheta_(_Z1),    100, -0.25, 0.25, 0.2),
        make_tuple( "cos_theta_z2", "cos_theta_z2",     _costheta_(_Z2),    100, -0.25, 0.25, 0.2),
        make_tuple( "cos_zeta_z1",  "cos_zeta_z1",      _coszeta_(_Z1),     100, -0.25, 0.25, 0.2),
        make_tuple( "cos_zeta_z2",  "cos_zeta_z2",      _coszeta_(_Z2),     100, -0.25, 0.25, 0.2),
        make_tuple( "angle_z1leps", "angle_z1leps/3.14",_alpha_(_Z1),       100, -0.25, 0.25, 0.1),
        make_tuple( "angle_z2leps", "angle_z2leps/3.14",_alpha_(_Z2),       100, -0.25, 0.25, 0.1),
        make_tuple( "angle_z1l2_z2","angle_z1l2_z2/3.14",   _beta,          100, -0.25, 0.25, 0.1)
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
            // Get parameters from tuple and form string
            TString hname,  quantity,   xlabel;
            int     bins;
            float   xmin,   xmax,       width;
            TString weight = "genWeight";
            tie(hname, quantity, xlabel, bins, xmin, xmax, width) = v[j];
            
            // Add subtraction to quantity
            quantity = "gen_" + quantity + " - " + quantity;

            // Create and draw histogram
            TH1D *h = new TH1D(hname + "_" + suffix, "", bins, xmin, xmax);
            tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

            xlabel.ReplaceAll(_l, lepChan[i]);
            xlabel = "\\Delta " + xlabel;
            h->GetXaxis()->SetTitle(xlabel);
            h->Sumw2(kTRUE);
            h->Write();

            TCanvas *canvas = new TCanvas(hname + "_" + selection[i] + "_canvas", "", 100, 100);
            canvas->cd();
            Facelift(canvas);
            canvas->SetCanvasSize(lCanvasSize, 0.625*lCanvasSize);
            canvas->SetMargin(lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);

            h->SetFillColor(lColor[i]);
            h->SetLineColor(lColor[i]);
            Facelift(h);
            h->SetStats(kTRUE);
            h->Draw("HIST");

            canvas->Update();
            TPaveStats *stats = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
            stats->SetOptStat(1000111110);
            stats->SetTextFont(lHelveticaMediumR);
            stats->SetTextSize(lSmall);
            stats->SetX1NDC(0.7); stats->SetY1NDC(0.5);

            float x = 0.5 * width;
            float y1 = gPad->GetUymin(), y2 = gPad->GetUymax();
            TLine *line[2] = {new TLine(-x, y1, -x, y2), new TLine(x, y1, x, y2)};
            for (unsigned l = 0; l < 2; l++)
            {
                line[l]->SetLineColor(kBlack);
                line[l]->SetLineStyle(kDashed);
                line[l]->SetLineWidth(2);
                line[l]->Draw();
            }

            gPad->Modified();
            gPad->Update();

            canvas->Write();
        }

        cout << "done!" << endl;
    }



    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms, canvases to " << outName << endl << endl << endl;
}
