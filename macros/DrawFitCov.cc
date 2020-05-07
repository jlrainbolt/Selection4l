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
**  DrawFitCov
**
**  Draws covariance matrices from simultaneous fir
*/ 

void DrawFitCov()
{

    //
    //  INPUT FILE
    //

    TString inName  = "fit_covariance.root";
    TFile   *inFile = TFile::Open(inName);

    cout << endl << endl << "Opened " << inName << endl << endl;



    //
    //  GET HISTOGRAMS
    //

    TH2D *h_cov_total, *h_cov_syst, *h_cov_stat;

    inFile->GetObject("cov_syst", h_cov_syst);
    h_cov_syst->SetDirectory(0);

    inFile->GetObject("cov_stat", h_cov_stat);
    h_cov_stat->SetDirectory(0);

    inFile->GetObject("cov_total", h_cov_total);
    h_cov_total->SetDirectory(0);

    cout << "Got histograms from " << inName << endl;

    inFile->Close();



    //
    //  OUTPUT FILE
    //

    TString outName = "fit_matrices.root";
    TFile *outFile  = new TFile(outName, "RECREATE");


    // Draw

    // Systematic covariance matrix
    TString title = "Systematic covariance matrix";
    TCanvas *c_syst = new TCanvas("c_cov_syst", "", lCanvasSize, lCanvasSize);
    c_syst->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
    c_syst->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

    c_syst->cd();
    c_syst->AddExec("t_syst", "gStyle->SetPaintTextFormat(\".1e\")");
    c_syst->AddExec("p_syst", "gStyle->SetPalette(kBird)");
    h_cov_syst->SetTitle(title);
    h_cov_syst->SetStats(0);
    h_cov_syst->SetMarkerSize(1.5);
    h_cov_syst->SetZTitle("\\sigma^{2}_{ij}");
    h_cov_syst->SetTitleSize(0.05, "xyz");
    h_cov_syst->SetTickLength(0, "xy");
    h_cov_syst->SetNdivisions(h_cov_syst->GetNbinsX(), "x");
    h_cov_syst->SetNdivisions(h_cov_syst->GetNbinsY(), "y");
    h_cov_syst->GetYaxis()->SetTitleOffset(0.6);
    h_cov_syst->GetZaxis()->SetTitleOffset(0.8);

    float syst_median = 0.5 * (h_cov_syst->GetMaximum() - h_cov_syst->GetMinimum()) + h_cov_syst->GetMinimum();
    h_cov_syst->DrawClone("COLZ");
    h_cov_syst->GetZaxis()->SetRangeUser(h_cov_syst->GetMinimum(), syst_median);
    h_cov_syst->SetMarkerColor(kWhite);
    h_cov_syst->DrawClone("TEXT SAME");
    h_cov_syst->GetZaxis()->SetRangeUser(syst_median, h_cov_syst->GetMaximum());
    h_cov_syst->SetMarkerColor(kBlack);
    h_cov_syst->DrawClone("TEXT SAME");

    c_syst->AutoExec();
    c_syst->SaveAs(".pdf");


    // Statistical covariance matrix
    title = "Statistical covariance matrix";
    TCanvas *c_stat = new TCanvas("c_cov_stat", "", lCanvasSize, lCanvasSize);
    c_stat->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
    c_stat->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

    c_stat->cd();
    c_stat->AddExec("t_stat", "gStyle->SetPaintTextFormat(\".1e\")");
    c_stat->AddExec("p_stat", "gStyle->SetPalette(kBird)");
    h_cov_stat->SetTitle(title);
    h_cov_stat->SetStats(0);
    h_cov_stat->SetMarkerSize(1.5);
    h_cov_stat->SetZTitle("\\sigma^{2}_{ij}");
    h_cov_stat->SetTitleSize(0.05, "xyz");
    h_cov_stat->SetTickLength(0, "xy");
    h_cov_stat->SetNdivisions(h_cov_syst->GetNbinsX(), "x");
    h_cov_stat->SetNdivisions(h_cov_syst->GetNbinsY(), "y");
    h_cov_stat->GetYaxis()->SetTitleOffset(0.6);
    h_cov_stat->GetZaxis()->SetTitleOffset(0.8);

    float stat_median = 0.5 * (h_cov_stat->GetMaximum() - h_cov_stat->GetMinimum()) + h_cov_stat->GetMinimum();
    h_cov_stat->DrawClone("COLZ");
    h_cov_stat->GetZaxis()->SetRangeUser(h_cov_stat->GetMinimum(), stat_median);
    h_cov_stat->SetMarkerColor(kWhite);
    h_cov_stat->DrawClone("TEXT SAME");
    h_cov_stat->GetZaxis()->SetRangeUser(stat_median, h_cov_stat->GetMaximum());
    h_cov_stat->SetMarkerColor(kBlack);
    h_cov_stat->DrawClone("TEXT SAME");

    c_stat->AutoExec();
    c_stat->SaveAs(".pdf");


    // Total covariance matrix
    title = "Total covariance matrix";
    TCanvas *c_tot = new TCanvas("c_cov_total", "", lCanvasSize, lCanvasSize);
    c_tot->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
    c_tot->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

    c_tot->cd();
    c_tot->AddExec("t_tot", "gStyle->SetPaintTextFormat(\".1e\")");
    c_tot->AddExec("p_tot", "gStyle->SetPalette(kBird)");
    h_cov_total->SetTitle(title);
    h_cov_total->SetStats(0);
    h_cov_total->SetMarkerSize(1.5);
    h_cov_total->SetZTitle("\\sigma^{2}_{ij}");
    h_cov_total->SetTitleSize(0.05, "xyz");
    h_cov_total->SetTickLength(0, "xy");
    h_cov_total->SetNdivisions(h_cov_total->GetNbinsX(), "x");
    h_cov_total->SetNdivisions(h_cov_total->GetNbinsY(), "y");
    h_cov_total->GetYaxis()->SetTitleOffset(0.6);
    h_cov_total->GetZaxis()->SetTitleOffset(0.8);

    float tot_median = 0.5 * (h_cov_total->GetMaximum() - h_cov_total->GetMinimum()) + h_cov_total->GetMinimum();
    h_cov_total->DrawClone("COLZ");
    h_cov_total->GetZaxis()->SetRangeUser(h_cov_total->GetMinimum(), tot_median);
    h_cov_total->SetMarkerColor(kWhite);
    h_cov_total->DrawClone("TEXT SAME");
    h_cov_total->GetZaxis()->SetRangeUser(tot_median, h_cov_total->GetMaximum());
    h_cov_total->SetMarkerColor(kBlack);
    h_cov_total->DrawClone("TEXT SAME");

    c_tot->AutoExec();
    c_tot->SaveAs(".pdf");


    // Write to file
    c_tot->Write();
    c_stat->Write();
    c_syst->Write();

    outFile->Close();
    cout << "Wrote canvases to " << outName << endl;
}
