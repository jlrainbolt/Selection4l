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
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;


/*
**  DrawResolution
**
**  Draws resolution of distributions for a "matched_" sample
*/ 

void DrawResolution(TString chan = "4l")
{

    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, int, float, float, float>> v = {

        //          name            quantity            axis label          bins xmin   xmax  width
        make_tuple( "zzm",          "zzp4.M()",         _m_(_4l),           100, -20,   20,    1),
        make_tuple( "b_z1m",        "b_z1p4.M()",       _m_(_Z1),           100, -10,   10,    8),
        make_tuple( "b_z2m",        "b_z2p4.M()",       _m_(_Z2),           100, -10,   10,    4),
        make_tuple( "b_ttm",        "b_ttp4.M()",       _m_(_l_("2,3,4")),  100, -10,   10,    5),
        make_tuple( "b_l1p",        "b_l1v3.Mag()",     _p_(_l_(1)),        100, -10,   10,  2.5),
        make_tuple( "angle_z1leps", "angle_z1leps/3.14",_alpha_(_Z1),       100, -0.25, 0.25, 0.1),
        make_tuple( "angle_z2leps", "angle_z2leps/3.14",_alpha_(_Z2),       100, -0.25, 0.25, 0.1),
        make_tuple( "angle_z1l2_z2","angle_z1l2_z2/3.14",   _beta,          100, -0.25, 0.25, 0.1),
        make_tuple( "cos_theta_z1", "cos_theta_z1",     _costheta_(_Z1),    100, -0.25, 0.25, 0.2),
        make_tuple( "cos_theta_z2", "cos_theta_z2",     _costheta_(_Z2),    100, -0.25, 0.25, 0.2),
        make_tuple( "sin_phi",      "sin_phi",          _sinphi,            100, -0.25, 0.25, 0.2)
    };



    //
    //  INPUT FILE
    //

    TString suffix  = "zz_4l";
    TString inName  = "matched_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Boosted/" + YEAR_STR + "_new/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;

    TTree *tree;
    inFile->GetObject(chan + "_" + suffix, tree);

    cout << chan + " tree has " << tree->GetEntries() << " events." << flush;



    //
    //  OUTPUT FILE
    //

    TString prefix  = "resolution";
    TString outName = prefix + "_" + YEAR_STR + "_" + chan + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



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
        TString weight = "weight";
        tie(hname, quantity, xlabel, bins, xmin, xmax, width) = v[j];
        if (chan.EqualTo("4e"))
            width *= 2;

        // Add subtraction to quantity
        quantity = "gen_" + quantity + " - " + quantity;

        // Create and draw histogram
        TH1D *h = new TH1D(hname + "_" + suffix, "", bins, xmin, xmax);
        tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

        xlabel.ReplaceAll("\\ ", " ");
        xlabel.ReplaceAll("\\ell", "l");
        xlabel.ReplaceAll("\\cos", "cos ");
        xlabel.ReplaceAll("\\mbox{", "");
        xlabel.ReplaceAll("})", ")");
        xlabel.ReplaceAll("}_", "_");
        xlabel.ReplaceAll("  ", " ");
        xlabel = "#Delta " + xlabel;
        h->GetXaxis()->SetTitle(xlabel);
        h->Sumw2(kFALSE);
        h->Scale(INT_LUMI * 1000. * XSEC_ZZ_4L / NGEN_ZZ_4L);

        TCanvas *canvas = new TCanvas(YEAR_STR + "_" + hname + "_resolution", "", 100, 100);
        canvas->cd();
        Facelift(canvas);
//      canvas->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
//      canvas->SetMargin(lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);

        h->SetFillColor(kMagenta);
        h->SetLineColor(kMagenta);
        Facelift(h);
        h->SetStats(kTRUE);
        h->Draw("HIST");
        h->SetStats(0);
        h->Write();

        canvas->Update();
//      TPaveStats *stats = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
//      stats->SetOptStat(1000111110);
//      stats->SetTextFont(lHelveticaMediumR);
//      stats->SetTextSize(lSmall);
//      stats->SetX1NDC(0.7); stats->SetY1NDC(0.5);

        float x = 0.5 * width;
        float y1 = 0, y2 = gPad->GetUymax();
        TLine *line[2] = {new TLine(-x, y1, -x, y2), new TLine(x, y1, x, y2)};
        for (unsigned l = 0; l < 2; l++)
        {
            line[l]->SetLineColor(kBlack);
            line[l]->SetLineStyle(kDashed);
            line[l]->SetLineWidth(2);
            line[l]->Draw();
        }

//      gPad->Modified();
//      gPad->Update();

//      canvas->SaveAs(".pdf");
        canvas->Write();
    }

    cout << "done!" << endl;



    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms, canvases to " << outName << endl << endl << endl;
}
