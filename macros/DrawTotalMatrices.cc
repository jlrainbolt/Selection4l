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

    TH2 *muID[H], *elID[H], *reco[H], *other[H], *syst[H], *unf[H], *tot[H];

    for (unsigned h = 0; h < H; h++)
    {
        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_MuonID", muID[h]);
        muID[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_ElecID", elID[h]);
        elID[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_ElecReco", reco[h]);
        reco[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_Others", other[h]);
        other[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_syst", syst[h]);
        syst[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_cov_unf", unf[h]);
        unf[h]->SetDirectory(0);

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



        // Muon ID covariance matrix
        TString title = "Muon ID covariance matrix";
        TCanvas *c_muID = new TCanvas("comb_" + hnames[h] + "_MuonID_cov", "", 
                lCanvasSize, lCanvasSize);
        c_muID->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_muID->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_muID->cd();
        c_muID->AddExec("t_muID", "gStyle->SetPaintTextFormat(\".1e\")");
        c_muID->AddExec("p_muID", "gStyle->SetPalette(kBird)");
        muID[h]->SetTitle(title);
        muID[h]->SetStats(0);
        muID[h]->SetMarkerSize(2);
        muID[h]->SetZTitle("\\sigma^{2}_{ij} (kEv/bin)^{2}");
        muID[h]->SetTitleSize(0.05, "xyz");
        muID[h]->SetTickLength(0, "xy");
        muID[h]->SetNdivisions(muID[h]->GetNbinsX(), "x");
        muID[h]->SetNdivisions(muID[h]->GetNbinsY(), "y");
        muID[h]->GetYaxis()->SetTitleOffset(0.6);
        muID[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            muID[h]->SetMarkerSize(1);
        else if (hnames[h].EqualTo("b_l1p"))
            muID[h]->SetMarkerSize(1.75);
        else if (hnames[h].EqualTo("b_ttm") || hnames[h].EqualTo("b_z2m"))
            muID[h]->SetMarkerSize(1.5);

        float muID_median = 0.5 * (muID[h]->GetMaximum() - muID[h]->GetMinimum()) + muID[h]->GetMinimum();
        muID[h]->DrawClone("COLZ");
        muID[h]->GetZaxis()->SetRangeUser(muID[h]->GetMinimum(), muID_median);
        muID[h]->SetMarkerColor(kWhite);
        muID[h]->DrawClone("TEXT SAME");
        muID[h]->GetZaxis()->SetRangeUser(muID_median, muID[h]->GetMaximum());
        muID[h]->SetMarkerColor(kBlack);
        muID[h]->DrawClone("TEXT SAME");

        c_muID->AutoExec();
        c_muID->SaveAs(".pdf");


        // Electron ID covariance matrix
        title = "Electron ID covariance matrix";
        TCanvas *c_elID = new TCanvas("comb_" + hnames[h] + "_ElecID_cov", "", 
                lCanvasSize, lCanvasSize);
        c_elID->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_elID->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_elID->cd();
        c_elID->AddExec("t_elID", "gStyle->SetPaintTextFormat(\".1e\")");
        c_elID->AddExec("p_elID", "gStyle->SetPalette(kBird)");
        elID[h]->SetTitle(title);
        elID[h]->SetStats(0);
        elID[h]->SetMarkerSize(2);
        elID[h]->SetZTitle("\\sigma^{2}_{ij} (kEv/bin)^{2}");
        elID[h]->SetTitleSize(0.05, "xyz");
        elID[h]->SetTickLength(0, "xy");
        elID[h]->SetNdivisions(elID[h]->GetNbinsX(), "x");
        elID[h]->SetNdivisions(elID[h]->GetNbinsY(), "y");
        elID[h]->GetYaxis()->SetTitleOffset(0.6);
        elID[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            elID[h]->SetMarkerSize(1);
        else if (hnames[h].EqualTo("b_l1p"))
            elID[h]->SetMarkerSize(1.75);
        else if (hnames[h].EqualTo("b_ttm") || hnames[h].EqualTo("b_z2m"))
            elID[h]->SetMarkerSize(1.5);

        float elID_median = 0.5 * (elID[h]->GetMaximum() - elID[h]->GetMinimum()) + elID[h]->GetMinimum();
        elID[h]->DrawClone("COLZ");
        elID[h]->GetZaxis()->SetRangeUser(elID[h]->GetMinimum(), elID_median);
        elID[h]->SetMarkerColor(kWhite);
        elID[h]->DrawClone("TEXT SAME");
        elID[h]->GetZaxis()->SetRangeUser(elID_median, elID[h]->GetMaximum());
        elID[h]->SetMarkerColor(kBlack);
        elID[h]->DrawClone("TEXT SAME");

        c_elID->AutoExec();
        c_elID->SaveAs(".pdf");


        // Electron reco covariance matrix
        title = "Electron reco covariance matrix";
        TCanvas *c_reco = new TCanvas("comb_" + hnames[h] + "_ElecReco_cov", "", 
                lCanvasSize, lCanvasSize);
        c_reco->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_reco->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_reco->cd();
        c_reco->AddExec("t_reco", "gStyle->SetPaintTextFormat(\".1e\")");
        c_reco->AddExec("p_reco", "gStyle->SetPalette(kBird)");
        reco[h]->SetTitle(title);
        reco[h]->SetStats(0);
        reco[h]->SetMarkerSize(2);
        reco[h]->SetZTitle("\\sigma^{2}_{ij} (kEv/bin)^{2}");
        reco[h]->SetTitleSize(0.05, "xyz");
        reco[h]->SetTickLength(0, "xy");
        reco[h]->SetNdivisions(reco[h]->GetNbinsX(), "x");
        reco[h]->SetNdivisions(reco[h]->GetNbinsY(), "y");
        reco[h]->GetYaxis()->SetTitleOffset(0.6);
        reco[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            reco[h]->SetMarkerSize(1);
        else if (hnames[h].EqualTo("b_l1p"))
            reco[h]->SetMarkerSize(1.75);
        else if (hnames[h].EqualTo("b_ttm") || hnames[h].EqualTo("b_z2m"))
            reco[h]->SetMarkerSize(1.5);

        float reco_median = 0.5 * (reco[h]->GetMaximum() - reco[h]->GetMinimum()) + reco[h]->GetMinimum();
        reco[h]->DrawClone("COLZ");
        reco[h]->GetZaxis()->SetRangeUser(reco[h]->GetMinimum(), reco_median);
        reco[h]->SetMarkerColor(kWhite);
        reco[h]->DrawClone("TEXT SAME");
        reco[h]->GetZaxis()->SetRangeUser(reco_median, reco[h]->GetMaximum());
        reco[h]->SetMarkerColor(kBlack);
        reco[h]->DrawClone("TEXT SAME");

        c_reco->AutoExec();
        c_reco->SaveAs(".pdf");


        // Other systematics covariance matrix
        title = "Other systematics covariance matrix";
        TCanvas *c_other = new TCanvas("comb_" + hnames[h] + "_other_cov", "", 
                lCanvasSize, lCanvasSize);
        c_other->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_other->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_other->cd();
        c_other->AddExec("t_other", "gStyle->SetPaintTextFormat(\".1e\")");
        c_other->AddExec("p_other", "gStyle->SetPalette(kBird)");
        other[h]->SetTitle(title);
        other[h]->SetStats(0);
        other[h]->SetMarkerSize(2);
        other[h]->SetZTitle("\\sigma^{2}_{ij} (kEv/bin)^{2}");
        other[h]->SetTitleSize(0.05, "xyz");
        other[h]->SetTickLength(0, "xy");
        other[h]->SetNdivisions(other[h]->GetNbinsX(), "x");
        other[h]->SetNdivisions(other[h]->GetNbinsY(), "y");
        other[h]->GetYaxis()->SetTitleOffset(0.6);
        other[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            other[h]->SetMarkerSize(1);
        else if (hnames[h].EqualTo("b_l1p"))
            other[h]->SetMarkerSize(1.75);
        else if (hnames[h].EqualTo("b_ttm") || hnames[h].EqualTo("b_z2m"))
            other[h]->SetMarkerSize(1.5);

        float other_median = 0.5 * (other[h]->GetMaximum() - other[h]->GetMinimum()) + other[h]->GetMinimum();
        other[h]->DrawClone("COLZ");
        other[h]->GetZaxis()->SetRangeUser(other[h]->GetMinimum(), other_median);
        other[h]->SetMarkerColor(kWhite);
        other[h]->DrawClone("TEXT SAME");
        other[h]->GetZaxis()->SetRangeUser(other_median, other[h]->GetMaximum());
        other[h]->SetMarkerColor(kBlack);
        other[h]->DrawClone("TEXT SAME");

        c_other->AutoExec();
        c_other->SaveAs(".pdf");


        // Systematic covariance matrix
        title = "Systematic covariance matrix";
        TCanvas *c_syst = new TCanvas("comb_" + hnames[h] + "_syst_cov", "", 
                lCanvasSize, lCanvasSize);
        c_syst->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_syst->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_syst->cd();
        c_syst->AddExec("t_syst", "gStyle->SetPaintTextFormat(\".1e\")");
        c_syst->AddExec("p_syst", "gStyle->SetPalette(kBird)");
        syst[h]->SetTitle(title);
        syst[h]->SetStats(0);
        syst[h]->SetMarkerSize(2);
        syst[h]->SetZTitle("\\sigma^{2}_{ij} (kEv/bin)^{2}");
        syst[h]->SetTitleSize(0.05, "xyz");
        syst[h]->SetTickLength(0, "xy");
        syst[h]->SetNdivisions(syst[h]->GetNbinsX(), "x");
        syst[h]->SetNdivisions(syst[h]->GetNbinsY(), "y");
        syst[h]->GetYaxis()->SetTitleOffset(0.6);
        syst[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            syst[h]->SetMarkerSize(1);
        else if (hnames[h].EqualTo("b_l1p"))
            syst[h]->SetMarkerSize(1.75);
        else if (hnames[h].EqualTo("b_ttm") || hnames[h].EqualTo("b_z2m"))
            syst[h]->SetMarkerSize(1.5);

        float syst_median = 0.5 * (syst[h]->GetMaximum() - syst[h]->GetMinimum()) + syst[h]->GetMinimum();
        syst[h]->DrawClone("COLZ");
        syst[h]->GetZaxis()->SetRangeUser(syst[h]->GetMinimum(), syst_median);
        syst[h]->SetMarkerColor(kWhite);
        syst[h]->DrawClone("TEXT SAME");
        syst[h]->GetZaxis()->SetRangeUser(syst_median, syst[h]->GetMaximum());
        syst[h]->SetMarkerColor(kBlack);
        syst[h]->DrawClone("TEXT SAME");

        c_syst->AutoExec();
        c_syst->SaveAs(".pdf");


        // Unfolding covariance matrix
        title = "Unfolding covariance matrix";
        TCanvas *c_unf = new TCanvas("comb_" + hnames[h] + "_unf_cov", "", 
                lCanvasSize, lCanvasSize);
        c_unf->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_unf->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_unf->cd();
        c_unf->AddExec("t_unf", "gStyle->SetPaintTextFormat(\".1e\")");
        c_unf->AddExec("p_unf", "gStyle->SetPalette(kBird)");
        unf[h]->SetTitle(title);
        unf[h]->SetStats(0);
        unf[h]->SetMarkerSize(2);
        unf[h]->SetZTitle("\\sigma^{2}_{ij} (kEv/bin)^{2}");
        unf[h]->SetTitleSize(0.05, "xyz");
        unf[h]->SetTickLength(0, "xy");
        unf[h]->SetNdivisions(syst[h]->GetNbinsX(), "x");
        unf[h]->SetNdivisions(syst[h]->GetNbinsY(), "y");
        unf[h]->GetYaxis()->SetTitleOffset(0.6);
        unf[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            unf[h]->SetMarkerSize(1);
        else if (hnames[h].EqualTo("b_l1p"))
            unf[h]->SetMarkerSize(1.75);
        else if (hnames[h].EqualTo("b_ttm") || hnames[h].EqualTo("b_z2m"))
            unf[h]->SetMarkerSize(1.5);

        float unf_median = 0.5 * (unf[h]->GetMaximum() - unf[h]->GetMinimum()) + unf[h]->GetMinimum();
        unf[h]->DrawClone("COLZ");
        unf[h]->GetZaxis()->SetRangeUser(unf[h]->GetMinimum(), unf_median);
        unf[h]->SetMarkerColor(kWhite);
        unf[h]->DrawClone("TEXT SAME");
        unf[h]->GetZaxis()->SetRangeUser(unf_median, unf[h]->GetMaximum());
        unf[h]->SetMarkerColor(kBlack);
        unf[h]->DrawClone("TEXT SAME");

        c_unf->AutoExec();
        c_unf->SaveAs(".pdf");


        // Total covariance matrix
        title = "Total covariance matrix";
        TCanvas *c_tot = new TCanvas("comb_" + hnames[h] + "_tot_cov", "", 
                lCanvasSize, lCanvasSize);
        c_tot->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
        c_tot->SetMargin(0.6*lCanvasMargin, 1.1*lCanvasMargin, lCanvasMargin, lCanvasMargin);

        c_tot->cd();
        c_tot->AddExec("t_tot", "gStyle->SetPaintTextFormat(\".1e\")");
        c_tot->AddExec("p_tot", "gStyle->SetPalette(kBird)");
        tot[h]->SetTitle(title);
        tot[h]->SetStats(0);
        tot[h]->SetMarkerSize(2);
        tot[h]->SetZTitle("\\sigma^{2}_{ij} (kEv/bin)^{2}");
        tot[h]->SetTitleSize(0.05, "xyz");
        tot[h]->SetTickLength(0, "xy");
        tot[h]->SetNdivisions(tot[h]->GetNbinsX(), "x");
        tot[h]->SetNdivisions(tot[h]->GetNbinsY(), "y");
        tot[h]->GetYaxis()->SetTitleOffset(0.6);
        tot[h]->GetZaxis()->SetTitleOffset(0.8);
        if (hnames[h].EqualTo("sin_phi"))
            tot[h]->SetMarkerSize(1);
        else if (hnames[h].EqualTo("b_l1p"))
            tot[h]->SetMarkerSize(1.75);
        else if (hnames[h].EqualTo("b_ttm") || hnames[h].EqualTo("b_z2m"))
            tot[h]->SetMarkerSize(1.5);

        float tot_median = 0.5 * (tot[h]->GetMaximum() - tot[h]->GetMinimum()) + tot[h]->GetMinimum();
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
        c_unf->Write();
        c_syst->Write();

        c_muID->Write();
        c_elID->Write();
        c_reco->Write();
        c_other->Write();
    }

    outFile->Close();
}
