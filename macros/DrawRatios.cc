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
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TMathText.h"

// Custom
//#include "Cuts2016.hh"
#include "Cuts2017.hh"

using namespace std;


/*
**  DrawRatios
**
**  Draws ratio plots from unfolding studies
*/ 

void DrawRatios()
{

    //
    //  INPUT FILES
    //

    TString inName  = "unfolding_" + YEAR_STR + ".root";
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

    TString zzName = "4l_" + YEAR_STR + "_zz_4l.root";
    TFile   *zzFile = TFile::Open(zzName);
    cout << endl << endl << "Opened " << zzName << endl;

    TString psName = "4l_" + YEAR_STR + "_phase_space.root";
    TFile   *psFile = TFile::Open(psName);
    cout << endl << "Opened " << psName << endl << endl;



    //
    //  GET HISTOGRAMS
    //

    // Get histograms
    TH1 *data[H], *gen[H], *reco[H], *sel[H], *ps[H], *res[H];

    for (unsigned h = 0; h < H; h++)
    {
        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_data", data[h]);
        data[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_bkg", res[h]);
        res[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_gen", gen[h]);
        gen[h]->SetDirectory(0);

        inFile->GetObject(hnames[h] + "/" + hnames[h] + "_reco", reco[h]);
        reco[h]->SetDirectory(0);

        zzFile->GetObject("4l/" + hnames[h] + "_zz_4l", sel[h]);
        sel[h]->SetDirectory(0);

        psFile->GetObject("4l/" + hnames[h] + "_phase_space", ps[h]);
        ps[h]->SetDirectory(0);
    }
    cout << "Got histograms" << endl;

    inFile->Close();
    zzFile->Close();
    psFile->Close();



    //
    //  OUTPUT FILE
    //

    TString outName = "ratios_" + YEAR_STR + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");


    // Draw
    for (unsigned h = 0; h < H; h++)
    {
        TString title = data[h]->GetYaxis()->GetTitle();
        TString width; width.Form("%g", data[h]->GetBinWidth(1));
        title.Prepend("\\mbox{Events/}" + width + "\\ ");

        
        // Reco vs. gen

        TCanvas *c_reco_gen = new TCanvas(YEAR_STR + "_" + hnames[h] + "_reco_vs_gen",
                "", lCanvasSize, lCanvasSize);
        Facelift(c_reco_gen);

        gen[h]->SetMinimum(0);
        gen[h]->SetStats(0);
        gen[h]->SetLineColor(lRed);
        gen[h]->SetLineWidth(4);

        reco[h]->SetMinimum(0);
        reco[h]->SetMaximum(1.25 * gen[h]->GetMaximum());
        reco[h]->SetStats(0);
        reco[h]->SetLineColor(lBlue);
        reco[h]->SetLineWidth(4);

        TRatioPlot *r_reco_gen = new TRatioPlot(reco[h], gen[h], "divsym");
        r_reco_gen->SetH1DrawOpt("E");
        r_reco_gen->SetH2DrawOpt("E");
        r_reco_gen->SetSeparationMargin(0.01);

        float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
        float RightPosition = 1,        RightMargin = -lLegendMargin;
        float TopPosition = 1,          TopMargin = -lLegendMargin;
        float BottomPosition = TopPosition - 0.085 * 3;
        float BottomMargin = 2. * lCanvasMargin - lLegendMargin;
        TLegend *l_reco_gen = new TLegend(LeftPosition + LeftMargin, BottomPosition - TopMargin,
                TopPosition + TopMargin, TopPosition + TopMargin);
        l_reco_gen->AddEntry(gen[h], "Gen MC", "L");
        l_reco_gen->AddEntry(reco[h], "Reco MC", "L");
        Facelift(l_reco_gen);

        c_reco_gen->cd();
        r_reco_gen->Draw();

        Facelift(r_reco_gen->GetUpperRefYaxis());
        r_reco_gen->GetUpperRefYaxis()->SetTitle(title);
        r_reco_gen->GetUpperRefYaxis()->SetTitleOffset(lTitleOffsetY);
        r_reco_gen->GetUpperPad()->SetLeftMargin(1.2 * lCanvasMargin);
        r_reco_gen->GetUpperPad()->Modified();

        Facelift(r_reco_gen->GetLowerRefXaxis());
        Facelift(r_reco_gen->GetLowerRefYaxis());
        r_reco_gen->GetLowerRefGraph()->SetMinimum(0.8);
        r_reco_gen->GetLowerRefGraph()->SetMaximum(1.2);
        r_reco_gen->GetLowerRefYaxis()->SetTitle("Reco/Gen");
        r_reco_gen->GetLowerPad()->SetBottomMargin(3 * lCanvasMargin);
        r_reco_gen->GetLowerPad()->SetLeftMargin(1.2 * lCanvasMargin);
        r_reco_gen->GetLowerPad()->Modified();

        r_reco_gen->GetUpperPad()->cd();
        l_reco_gen->Draw();

        
        // Selected vs. phase space

        TCanvas *c_axe = new TCanvas(YEAR_STR + "_" + hnames[h] + "_axe",
                "", lCanvasSize, lCanvasSize);
        Facelift(c_axe);
        c_axe->SetLogy();

        ps[h]->SetMinimum(100);
        ps[h]->SetStats(0);
        ps[h]->SetLineColor(lPurple);
        ps[h]->SetLineWidth(4);

        sel[h]->SetMinimum(100);
        sel[h]->SetMaximum(10 * ps[h]->GetMaximum());
        sel[h]->SetStats(0);
        sel[h]->SetLineColor(lBlue);
        sel[h]->SetLineWidth(4);

        TRatioPlot *r_axe = new TRatioPlot(sel[h], ps[h], "divsym");
        r_axe->SetH1DrawOpt("E");
        r_axe->SetH2DrawOpt("E");
        r_axe->SetSeparationMargin(0.01);

        TLegend *l_axe = new TLegend(LeftPosition + 2*lLegendMargin, BottomPosition - TopMargin,
                TopPosition + TopMargin, TopPosition + TopMargin);
        l_axe->AddEntry(ps[h], "Phase Space", "L");
        l_axe->AddEntry(sel[h], "Selected", "L");
        Facelift(l_axe);

        c_axe->cd();
        r_axe->Draw();

        Facelift(r_axe->GetUpperRefYaxis());
        r_axe->GetUpperRefYaxis()->SetTitle(title);
        r_axe->GetUpperRefYaxis()->SetTitleOffset(lTitleOffsetY);
        r_axe->GetUpperPad()->SetLeftMargin(1.2 * lCanvasMargin);
        r_axe->GetUpperPad()->Modified();

        Facelift(r_axe->GetLowerRefXaxis());
        Facelift(r_axe->GetLowerRefYaxis());
        r_axe->GetLowerRefYaxis()->SetTitle("(A \\times \\epsilon)");
        r_axe->GetLowerPad()->SetBottomMargin(3 * lCanvasMargin);
        r_axe->GetLowerPad()->SetLeftMargin(1.2 * lCanvasMargin);
        r_axe->GetLowerPad()->Modified();

        r_axe->GetUpperPad()->cd();
        l_axe->Draw();

        
        // Folded vs. unfolded

//      TCanvas *c_fol_unf = new TCanvas(YEAR_STR + "_" + hnames[h] + "_fol_unf",
        TCanvas *c_fol_unf = new TCanvas(YEAR_STR + "_" + hnames[h] + "_bkg_data",
                "", lCanvasSize, lCanvasSize);
        Facelift(c_fol_unf);
        c_fol_unf->SetLogy();

//      res[h]->SetMinimum(0);
        res[h]->SetMinimum(1);
//      res[h]->SetMaximum(1.5 * data[h]->GetMaximum());
        res[h]->SetMaximum(10 * data[h]->GetMaximum());
        res[h]->SetStats(0);
//      res[h]->SetLineColor(lOrange);
        res[h]->SetLineColor(lGreen);
        res[h]->SetLineWidth(2);
//      res[h]->SetMarkerColor(lOrange);
        res[h]->SetMarkerColor(lGreen);
        res[h]->SetMarkerStyle(22);
        res[h]->SetMarkerSize(2);

//      data[h]->SetMinimum(0);
        data[h]->SetMinimum(1);
//      data[h]->SetMaximum(1.5 * data[h]->GetMaximum());
        data[h]->SetMaximum(10 * data[h]->GetMaximum());
        data[h]->SetStats(0);
        data[h]->SetLineColor(kBlack);
        data[h]->SetLineWidth(2);
        data[h]->SetMarkerColor(kBlack);
        data[h]->SetMarkerStyle(20);
        data[h]->SetMarkerSize(2);

        TRatioPlot *r_fol_unf = new TRatioPlot(res[h], data[h], "divsym");
        r_fol_unf->SetH1DrawOpt("E1");
        r_fol_unf->SetH2DrawOpt("E1");
        r_fol_unf->SetSeparationMargin(0.01);

//      TLegend *l_fol_unf = new TLegend(LeftPosition + LeftMargin, BottomPosition - TopMargin,
        TLegend *l_fol_unf = new TLegend(LeftPosition + 3*lLegendMargin, BottomPosition - TopMargin,
                TopPosition + TopMargin, TopPosition + TopMargin);
        l_fol_unf->AddEntry(data[h], "Data", "LP");
        l_fol_unf->AddEntry(res[h], "Background", "LP");
        Facelift(l_fol_unf);

        c_fol_unf->cd();
        r_fol_unf->Draw();

        Facelift(r_fol_unf->GetUpperRefYaxis());
        r_fol_unf->GetUpperRefYaxis()->SetTitle(title);
        r_fol_unf->GetUpperRefYaxis()->SetTitleOffset(lTitleOffsetY);
        r_fol_unf->GetUpperPad()->SetLeftMargin(1.2 * lCanvasMargin);
        r_fol_unf->GetUpperPad()->Modified();

        Facelift(r_fol_unf->GetLowerRefXaxis());
        Facelift(r_fol_unf->GetLowerRefYaxis());
        r_fol_unf->GetLowerRefYaxis()->SetTitleOffset(lTitleOffsetY);
//      r_fol_unf->GetLowerRefGraph()->SetMinimum(0.5);
//      r_fol_unf->GetLowerRefGraph()->SetMaximum(1.5);
        r_fol_unf->GetLowerRefGraph()->SetMinimum(0);
        r_fol_unf->GetLowerRefGraph()->SetMaximum(0.3);
//      r_fol_unf->GetLowerRefGraph()->SetLineColor(lOrange);
        r_fol_unf->GetLowerRefYaxis()->SetTitle("Bkg/Data");
        r_fol_unf->GetLowerPad()->SetBottomMargin(3 * lCanvasMargin);
        r_fol_unf->GetLowerPad()->SetLeftMargin(1.2 * lCanvasMargin);
        r_fol_unf->GetLowerPad()->Modified();

        r_fol_unf->GetUpperPad()->cd();
        l_fol_unf->Draw();


        // Write to file
        outFile->mkdir(hnames[h]);
        outFile->cd(hnames[h]);
        c_reco_gen->Write();
        c_fol_unf->Write();
        c_axe->Write();
    }

    outFile->Close();
}
