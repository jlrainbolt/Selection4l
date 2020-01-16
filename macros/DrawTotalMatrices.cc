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
**  DrawTotalMatrices
**
**  Draws systematic and total matrices from unfolding process
*/ 

void DrawTotalMatrices()
{

    //
    //  INPUT FILE
    //

    TString inName  = "ddr_uncertainty.root";
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

    TH2 *syst[H], *tot[H];

    for (unsigned h = 0; h < H; h++)
    {
        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_syst", syst[h]);
        syst[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_tot", tot[h]);
        tot[h]->SetDirectory(0);
    }
    cout << "Got histograms from " << inName << endl;

    inFile->Close();



    //
    //  OUTPUT FILE
    //

    TString outName = "matrices_total.root";
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


        // Systematic covariance matrix
        TString title = "Systematic covariance matrix";
        TCanvas *c_syst = new TCanvas("comb_" + hnames[h] + "_syst_cov", "", 
                lCanvasSize, lCanvasSize);
        c_syst->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_syst->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_syst->cd();
        c_syst->AddExec("t_syst", "gStyle->SetPaintTextFormat(\".2f\")");
        c_syst->AddExec("p_syst", "gStyle->SetPalette(kBird)");
        syst[h]->SetTitle(title);
        syst[h]->SetStats(0);
        syst[h]->SetMarkerSize(2);
        syst[h]->SetZTitle("\\sigma^{2}_{ij} (Events/bin)^{2}");
        syst[h]->SetTitleSize(0.05, "xyz");
        syst[h]->SetTickLength(0, "xy");
        syst[h]->SetNdivisions(syst[h]->GetNbinsX(), "x");
        syst[h]->SetNdivisions(syst[h]->GetNbinsY(), "y");
        syst[h]->GetYaxis()->SetTitleOffset(0.6);
        syst[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            syst[h]->SetMarkerSize(1.25);

        float syst_median = 0.5 * (syst[h]->GetMaximum() - syst[h]->GetMinimum());
        syst[h]->DrawClone("COLZ");
        syst[h]->GetZaxis()->SetRangeUser(syst[h]->GetMinimum(), syst_median);
        syst[h]->SetMarkerColor(kWhite);
        syst[h]->DrawClone("TEXT SAME");
        syst[h]->GetZaxis()->SetRangeUser(syst_median, syst[h]->GetMaximum());
        syst[h]->SetMarkerColor(kBlack);
        syst[h]->DrawClone("TEXT SAME");

        c_syst->AutoExec();
        c_syst->SaveAs(".pdf");


        // Total covariance matrix
        title = "Total covariance matrix";
        TCanvas *c_tot = new TCanvas("comb_" + hnames[h] + "_tot_cov", "", 
                lCanvasSize, lCanvasSize);
        c_tot->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_tot->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_tot->cd();
        c_tot->AddExec("t_tot", "gStyle->SetPaintTextFormat(\".2f\")");
        c_tot->AddExec("p_tot", "gStyle->SetPalette(kBird)");
        tot[h]->SetTitle(title);
        tot[h]->SetStats(0);
        tot[h]->SetMarkerSize(2);
        tot[h]->SetZTitle("\\sigma^{2}_{ij} (Events/bin)^{2}");
        tot[h]->SetTitleSize(0.05, "xyz");
        tot[h]->SetTickLength(0, "xy");
        tot[h]->SetNdivisions(tot[h]->GetNbinsX(), "x");
        tot[h]->SetNdivisions(tot[h]->GetNbinsY(), "y");
        tot[h]->GetYaxis()->SetTitleOffset(0.6);
        tot[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            tot[h]->SetMarkerSize(1.25);

        float tot_median = 0.5 * (tot[h]->GetMaximum() - tot[h]->GetMinimum());
        tot[h]->DrawClone("COLZ");
        tot[h]->GetZaxis()->SetRangeUser(tot[h]->GetMinimum(), tot_median);
        tot[h]->SetMarkerColor(kWhite);
        tot[h]->DrawClone("TEXT SAME");
        tot[h]->GetZaxis()->SetRangeUser(tot_median, tot[h]->GetMaximum());
        tot[h]->SetMarkerColor(kBlack);
        tot[h]->DrawClone("TEXT SAME");

        c_tot->AutoExec();
        c_tot->SaveAs(".pdf");


        // Write to file
        outFile->mkdir(hnames[h]);
        outFile->cd(hnames[h]);
        c_tot->Write();
        c_syst->Write();
    }

    outFile->Close();
}
