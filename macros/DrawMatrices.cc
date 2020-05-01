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
#include "TCanvas.h"
#include "TExec.h"
//#include "TMathText.h"

// Custom
#include "Cuts2018.hh"

using namespace std;


/*
**  DrawMatrices
**
**  Draws matrices from unfolding process
*/ 

void DrawMatrices()
{

    //
    //  INPUT FILE
    //

    TString inName  = "unfolding_comb.root";
    TFile   *inFile = TFile::Open(inName);

    cout << endl << endl << "Opened " << inName << endl << endl;

    // Get directory keys
    vector<TString> hnames;
    TKey *dirKey;
    TIter next(inFile->GetListOfKeys());
    while ((dirKey = (TKey*) next()))
    {
        TString hname = dirKey->GetName();
        hnames.push_back(hname);
    }
    cout << "Got directory keys from " << inName << endl;

    const unsigned H = hnames.size();



    //
    //  GET HISTOGRAMS
    //

    TH2 *resp[H], *cov[H], *var[H];
    TH1 *data[H];

    for (unsigned h = 0; h < H; h++)
    {
        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_response", resp[h]);
        resp[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_covariance", cov[h]);
        cov[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_data_cov", var[h]);
        var[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_data", data[h]);
        data[h]->SetDirectory(0);
    }
    cout << "Got histograms from " << inName << endl;

    inFile->Close();



    //
    //  OUTPUT FILE
    //

    TString outName = "matrices_comb.root";
    TFile *outFile  = new TFile(outName, "RECREATE");


    // Draw
    for (unsigned h = 0; h < H; h++)
    {
/*
        TString title = data[h]->GetXaxis()->GetTitle();
        TString width; width.Form("%g", data[h]->GetBinWidth(1));
        if (title.Contains("("))
        {
            title.ReplaceAll("(", "(" + width);
            title.ReplaceAll("})", " bins})");
            title.ReplaceAll("mbox{", "mbox{ ");
        }
        else
            title.Append("\\ (" + width + "\\mbox{ unit bins})"); 
*/


        // Response matrix
        TString title = "Response matrix";
        TCanvas *c_resp = new TCanvas("comb_" + hnames[h] + "_response", "", 
                lCanvasSize, lCanvasSize);
        c_resp->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_resp->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_resp->cd();
        c_resp->AddExec("t_resp", "gStyle->SetPaintTextFormat(\".4f\")");
        c_resp->AddExec("p_resp", "gStyle->SetPalette(kBird)");
        resp[h]->SetTitle(title);
        resp[h]->SetStats(0);
        resp[h]->SetMarkerSize(2);
        resp[h]->SetZTitle("P(reco in bin i | gen in bin j)");
        resp[h]->SetTitleSize(0.05, "xyz");
        resp[h]->SetTickLength(0, "xy");
        resp[h]->SetNdivisions(resp[h]->GetNbinsX(), "x");
        resp[h]->SetNdivisions(resp[h]->GetNbinsY(), "y");
        resp[h]->GetYaxis()->SetTitleOffset(0.6);
        resp[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            resp[h]->SetMarkerSize(1.25);

        resp[h]->SetMinimum(0);
        resp[h]->SetMaximum(1);
        resp[h]->DrawClone("COLZ");
        resp[h]->GetZaxis()->SetRangeUser(0, 0.5);
        resp[h]->SetMarkerColor(kWhite);
        resp[h]->DrawClone("TEXT SAME");
        resp[h]->GetZaxis()->SetRangeUser(0.5, 1);
        resp[h]->SetMarkerColor(kBlack);
        resp[h]->DrawClone("TEXT SAME");

        c_resp->AutoExec();
        c_resp->SaveAs(".pdf");


        // Covariance matrix
        title = "Unfolded covariance matrix";
        TCanvas *c_cov = new TCanvas("comb_" + hnames[h] + "_covariance", "", 
                lCanvasSize, lCanvasSize);
        c_cov->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_cov->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_cov->cd();
        c_cov->AddExec("t_cov", "gStyle->SetPaintTextFormat(\".2f\")");
        c_cov->AddExec("p_cov", "gStyle->SetPalette(kBird)");
        cov[h]->SetTitle(title);
        cov[h]->SetStats(0);
        cov[h]->SetMarkerSize(2);
        cov[h]->SetZTitle("\\sigma^{2}_{ij} (Events/bin)^{2}");
        cov[h]->SetTitleSize(0.05, "xyz");
        cov[h]->SetTickLength(0, "xy");
        cov[h]->SetNdivisions(resp[h]->GetNbinsX(), "x");
        cov[h]->SetNdivisions(resp[h]->GetNbinsY(), "y");
        cov[h]->GetYaxis()->SetTitleOffset(0.6);
        cov[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            cov[h]->SetMarkerSize(1.25);

        float cov_median = 0.5 * (cov[h]->GetMaximum() - cov[h]->GetMinimum()) + cov[h]->GetMinimum();
        cov[h]->DrawClone("COLZ");
        cov[h]->GetZaxis()->SetRangeUser(cov[h]->GetMinimum(), cov_median);
        cov[h]->SetMarkerColor(kWhite);
        cov[h]->DrawClone("TEXT SAME");
        cov[h]->GetZaxis()->SetRangeUser(cov_median, cov[h]->GetMaximum());
        cov[h]->SetMarkerColor(kBlack);
        cov[h]->DrawClone("TEXT SAME");

        c_cov->AutoExec();
        c_cov->SaveAs(".pdf");


        // Data "covariance" (variance) matrix
        title = "Data covariance matrix";
        TCanvas *c_var = new TCanvas("comb_" + hnames[h] + "_data_cov", "", 
                lCanvasSize, lCanvasSize);
        c_var->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_var->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_var->cd();
        c_var->AddExec("t_var", "gStyle->SetPaintTextFormat(\".2f\")");
        c_var->AddExec("p_var", "gStyle->SetPalette(kBird)");
        var[h]->SetTitle(title);
        var[h]->SetStats(0);
        var[h]->SetMarkerSize(2);
        var[h]->SetZTitle("\\sigma^{2}_{ij} (Events/bin)^{2}");
        var[h]->SetTitleSize(0.05, "xyz");
        var[h]->SetTickLength(0, "xy");
        var[h]->SetNdivisions(resp[h]->GetNbinsX(), "x");
        var[h]->SetNdivisions(resp[h]->GetNbinsY(), "y");
        var[h]->GetYaxis()->SetTitleOffset(0.6);
        var[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            var[h]->SetMarkerSize(1.25);

        float var_median = 0.5 * (var[h]->GetMaximum() - var[h]->GetMinimum()) + var[h]->GetMinimum();
        var[h]->DrawClone("COLZ");
        var[h]->GetZaxis()->SetRangeUser(var[h]->GetMinimum(), var_median);
        var[h]->SetMarkerColor(kWhite);
        var[h]->DrawClone("TEXT SAME");
        var[h]->GetZaxis()->SetRangeUser(var_median, var[h]->GetMaximum());
        var[h]->SetMarkerColor(kBlack);
        var[h]->DrawClone("TEXT SAME");

        c_var->AutoExec();
        c_var->SaveAs(".pdf");


        // Write to file
        outFile->mkdir(hnames[h]);
        outFile->cd(hnames[h]);
        c_resp->Write();
        c_cov->Write();
        c_var->Write();
    }

    outFile->Close();
}
