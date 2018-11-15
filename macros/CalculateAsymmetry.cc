// STL
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "PlotUtils.hh"

// Cuts
#include "Cuts2017.hh"
#include "PlotUtils.hh"

using namespace std;


/*
**  CalculateAsymmetry
**
**  Calculates triple product asymmetry observable using sin(phi) from unweighted4l histograms
*/

void CalculateAsymmetry(bool useDY = kFALSE)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;   // Channel indices
    unsigned                L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;
    TString selection[N] = {"4l",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N] = { 5,      6,      7,      8,      9};
    TString lepChan[N]  = {_4l,     _4mu,   _2mu2e, "",     _4e};

    // Index 0 corresponds to full binning, index 1 corresponds to 2 bins
    TH1 *hObserved[2][N], *hExpected[2][N], *hBackground[2][N], *hObsMinusBg[2][N];

    // Index 0 is raw, index 1 is scaled to data
    TH1 *hTotal[2][N][2];
    TString hname[2] = {"sin_phi", "sin_phi_2"},    htype[2] = {"raw", "scaled"};






    ////
    ////
    ////    FILL HISTOGRAMS
    ////
    ////


//  TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/";
    TString inPath = "";
    TString prefix = "unscaled4l";



    //
    //  DATA
    //

    // Muon file
    TString muName = inPath + prefix + "_" + MU_SUFF + ".root";
    TFile *muFile = TFile::Open(muName);
    cout << "Opened " << muName << endl;

    // Electron file
    TString elName = inPath + prefix + "_" + EL_SUFF + ".root";
    TFile *elFile = TFile::Open(elName);
    cout << "Opened " << elName << endl;

    for (unsigned i = 1; i < N; i++)
    {
        TH1 *hist[2];

        for (unsigned h = 0; h < 2; h++)
        {
            if      ((i == M4) || (i == ME))
                muFile->GetObject(selection[i] + "/" + hname[h] + "_" + MU_SUFF, hist[h]);
            else if ((i == E4) || (i == EM))
                elFile->GetObject(selection[i] + "/" + hname[h] + "_" + EL_SUFF, hist[h]);

            hist[h]->SetName(hname[h] + "_data");
            hist[h]->SetDirectory(0);
            hist[h]->Sumw2();
            hObserved[h][i] = hist[h];
        }
    }
    muFile->Close();
    elFile->Close();



    //
    //  MONTE CARLO
    //

    TH1 *mcHist[2][N][N_MC];
    for (unsigned j = 0; j < N_MC; j++) // sample loop
    {
        TString inName = inPath + prefix + "_" + MC_SUFF[j] + ".root";
        TFile *inFile = TFile::Open(inName);

        cout << "Opened " << inName << endl;

        for (unsigned i = 1; i < N; i++)
        {
            // Scale by appropriate lumi, etc.
            float LUMI;
            if      ((i == M4) || (i == ME))
                LUMI = MUON_TRIG_LUMI;
            else if ((i == E4) || (i == EM))
                LUMI = ELEC_TRIG_LUMI;
            float sf = LUMI * 1000 * XSEC[j] / NGEN[j];

            TH1 *hist[2];
            for (unsigned h = 0; h < 2; h++)
            {
                inFile->GetObject(selection[i] + "/" + hname[h] + "_" + MC_SUFF[j], hist[h]);

                hist[h]->SetDirectory(0);
                hist[h]->Sumw2();
                hist[h]->Scale(sf);
                mcHist[h][i][j] = hist[h];
            }
        }
        inFile->Close();
    }






    ////
    ////
    ////    CALCULATIONS
    ////
    ////

 
    //
    //  SUMS
    //

    // 4l
    for (unsigned h = 0; h < 2; h++)
    {
        hObserved[h][0] = (TH1*) hObserved[h][1]->Clone();
        for (unsigned i = 2; i < N; i++)
            hObserved[h][0]->Add(hObserved[h][i]);
        for (unsigned j = 0; j < N_MC; j++)
        {
            mcHist[h][0][j] = (TH1*) mcHist[h][1][j]->Clone();
            for (unsigned i = 2; i < N; i++)
                mcHist[h][0][j]->Add(mcHist[h][i][j]);
        }

        // 2m2e
        hObserved[h][ME]->Add(hObserved[h][EM]);
        for (unsigned j = 0; j < N_MC; j++)
            mcHist[h][ME][j]->Add(mcHist[h][EM][j]);

        // Expected, total
        for (unsigned i = 0; i < N; i++)
        {
            hExpected[h][i] = (TH1*) mcHist[h][i][0]->Clone("sin_phi_exp");
            for (unsigned j = 1; j < N_MC; j++)
                hExpected[h][i]->Add(mcHist[h][i][j]);

            for (unsigned k = 0; k < 2; k++)
            {
                if (useDY)
                    hTotal[h][i][k] = (TH1*) hExpected[h][i]->Clone();
                else
                {
                    hTotal[h][i][k] = (TH1*) mcHist[h][i][0]->Clone("sin_phi_tot");
                    for (unsigned j = 2; j < N_MC; j++)
                        hTotal[h][i][k]->Add(mcHist[h][i][j]);
                }
            }
        }

        // Background
        for (unsigned i = 0; i < N; i++)
        {
            hBackground[h][i] = (TH1*) mcHist[h][i][2]->Clone("sin_phi_bkg");

            if (useDY)
                hBackground[h][i]->Add(mcHist[h][i][DY]);

            for (unsigned j = 3; j < N_MC; j++)
                hBackground[h][i]->Add(mcHist[h][i][j]);

            // Observed minus background
            hObsMinusBg[h][i] = (TH1*) hObserved[h][i]->Clone("sin_phi_omb");
            hObsMinusBg[h][i]->Add(hBackground[h][i], -1);
        }
    }



    //
    //  PRINT
    //

    cout << endl << endl;
    for (unsigned b = 1; b <= 2; b++)
    {
        if (b == 1)
            cout << "NEGATIVE" << endl;
        else
            cout << "POSITIVE" << endl;
        cout << "\t\t" << "Observed" << "\t" << "Expected" << "\t" << "Signal  " << "\t";
        cout << "Background" << "\t" << "Obs - Bkg" << endl;
        for (unsigned i = 0; i < N; i++)
        {
            if (i == EM)
                continue;

            cout << selection[i] << "\t\t";
            cout << setw(8) << hObserved[1][i]->GetBinContent(b) << "\t";
            cout << setw(8) << hExpected[1][i]->GetBinContent(b) << "\t";
            cout << setw(8) << mcHist[1][i][0]->GetBinContent(b) << "\t";
            cout << setw(8) << hBackground[1][i]->GetBinContent(b) << "\t";
            cout << setw(8) << hObsMinusBg[1][i]->GetBinContent(b) << endl;
        }
    }






    ////
    ////
    ////    ASYMMETRIES
    ////
    ////


    //
    //  TOTAL
    //

    float asymmetry[N], unc[N], fracUnc[N], sig[N];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)
            continue;

        float neg = hObserved[1][i]->GetBinContent(1), pos = hObserved[1][i]->GetBinContent(2);

        asymmetry[i] = (pos - neg) / (pos + neg);

        unc[i] = 1 / sqrt(neg + pos);
        fracUnc[i] = fabs(unc[i] / asymmetry[i]);
        sig[i] = fabs(asymmetry[i] / unc[i]);
    }

    cout << endl << endl;
    cout << "ASYMMETRIES" << endl;
    cout << "\t\t" << "Value" << "\t\t\t\t\t" << "Frac. Unc.\t\tSignificance" << endl;
    for (unsigned i = L4; i < N; i++)
    {
        if (i == EM)
            continue;

        cout << selection[i] << "\t\t" << setw(6) << asymmetry[i];
        cout << " +- " << "\t" << unc[i] << "\t\t";
        cout << fracUnc[i] << "\t\t" << sig[i] << endl;
    }

    // Graph of results
    float yval[4] = {5, 3, 2, 1},   xerr[4] = {unc[L4], unc[M4], unc[ME], unc[E4]};
    float xval[4] = {asymmetry[L4], asymmetry[M4], asymmetry[ME], asymmetry[E4]};
    TGraphErrors *graph = new TGraphErrors(4, xval, yval, xerr, 0);


    float yedge[6] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
    graph->GetYaxis()->Set(5, yedge);
    graph->GetYaxis()->SetBinLabel(5, _4l);
    graph->GetYaxis()->SetBinLabel(3, _4mu);
    graph->GetYaxis()->SetBinLabel(2, _2mu2e);
    graph->GetYaxis()->SetBinLabel(1, _4e);
//  graph->GetYaxis()->SetTicks("-");

    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("A_{" + lepChan[0] + "}");
    graph->GetXaxis()->SetTitleOffset(0.7 * lTitleOffsetY);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerSize(2);
    graph->SetLineWidth(2);

    TCanvas *asyCanvas = new TCanvas("asymmetry", "", lCanvasSize, lCanvasSize);
    asyCanvas->cd();
    Facelift(asyCanvas);
    asyCanvas->SetCanvasSize(lCanvasSize, 0.5*lCanvasSize);
    asyCanvas->SetMargin(1.3*lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);
    graph->Draw("AP");

    Facelift(graph->GetXaxis());
    Facelift(graph->GetYaxis());
    graph->GetYaxis()->SetLabelSize(lLarge);
    graph->GetYaxis()->SetLabelOffset(0.01 * lTitleOffsetY);
    graph->GetYaxis()->SetTickLength(0);
    graph->GetYaxis()->SetRangeUser(yedge[0] - 0.5, yedge[5] + 0.5);
    graph->GetXaxis()->SetLimits(graph->GetXaxis()->GetXmin(), fabs(graph->GetXaxis()->GetXmin()));
    gPad->Modified();

    TLine *line = new TLine(0, yedge[0] - 0.5, 0, yedge[5] + 0.5);
    line->SetLineColor(kRed);
    line->Draw();
    graph->Draw("P");
    gPad->RedrawAxis();




    ////
    ////
    ////    WITH CUT
    ////
    ////


    const unsigned M = 10, X = 7;
    float phi_min[M] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};
    float phi_max[X] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

    float min_cut[N][M], min_unc[N][M], max_cut[N][M], max_unc[N][M];
    TGraphErrors *min_graph[N], *max_graph[N];
    TCanvas *minCanvas[N], *maxCanvas[N];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)
            continue;

        for (unsigned j = 0; j < M; j++)
        {
            // Minimum cut
            float neg = hObserved[0][i]->Integral(1, 1 + j);
            float pos = hObserved[0][i]->Integral(20 - j, 20);

            min_cut[i][j] = (pos - neg) / (pos + neg);
            min_unc[i][j] = 1 / sqrt(neg + pos);


            // Maximum cut
            neg = hObserved[0][i]->Integral(7 - j, 10);
            pos = hObserved[0][i]->Integral(11, 14 + j);

            max_cut[i][j] = (pos - neg) / (pos + neg);
            max_unc[i][j] = 1 / sqrt(neg + pos);
        }
        min_graph[i] = new TGraphErrors(M, phi_min, min_cut[i], 0, min_unc[i]);
        min_graph[i]->SetTitle("");
        min_graph[i]->GetXaxis()->SetTitle("\\min\\,|\\sin\\phi|");
        min_graph[i]->GetYaxis()->SetTitle("A_{" + lepChan[i] + "}");
        min_graph[i]->GetYaxis()->SetTitleOffset(0.7 * lTitleOffsetY);
        min_graph[i]->SetMarkerStyle(kFullCircle);
        min_graph[i]->SetMarkerSize(2);
        min_graph[i]->SetLineWidth(2);

        minCanvas[i] = new TCanvas("abs_sinphi_min_" + selection[i], "", lCanvasSize, lCanvasSize);
        minCanvas[i]->cd();
        Facelift(minCanvas[i]);
//      minCanvas[i]->SetLeftMargin(1.2 * lCanvasMargin);
        minCanvas[i]->SetCanvasSize(lCanvasSize, 0.5 * lCanvasSize);
        minCanvas[i]->SetMargin(1.3*lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);
        min_graph[i]->Draw("AP");

        Facelift(min_graph[i]->GetXaxis());
        Facelift(min_graph[i]->GetYaxis());
        if (min_graph[i]->GetYaxis()->GetXmax() < 0)
            min_graph[i]->GetYaxis()->SetRangeUser(min_graph[i]->GetYaxis()->GetXmin(), 0);
        min_graph[i]->GetXaxis()->SetLimits(-0.05, 0.95);
        gPad->Modified();



        max_graph[i] = new TGraphErrors(X, phi_max, max_cut[i], 0, max_unc[i]);
        max_graph[i]->SetTitle("");
        max_graph[i]->GetXaxis()->SetTitle("\\max\\,|\\sin\\phi|");
        max_graph[i]->GetYaxis()->SetTitle("A_{" + lepChan[i] + "}");
        max_graph[i]->GetYaxis()->SetTitleOffset(0.7 * lTitleOffsetY);
        max_graph[i]->SetMarkerStyle(kFullCircle);
        max_graph[i]->SetMarkerSize(2);
        max_graph[i]->SetLineWidth(2);

        maxCanvas[i] = new TCanvas("abs_sinphi_max_" + selection[i], "", lCanvasSize, lCanvasSize);
        maxCanvas[i]->cd();
        Facelift(maxCanvas[i]);
//      maxCanvas[i]->SetLeftMargin(1.2 * lCanvasMargin);
        maxCanvas[i]->SetCanvasSize(lCanvasSize, 0.5 * lCanvasSize);
        maxCanvas[i]->SetMargin(1.3*lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);
        max_graph[i]->Draw("AP");

        Facelift(max_graph[i]->GetXaxis());
        Facelift(max_graph[i]->GetYaxis());
        if (max_graph[i]->GetYaxis()->GetXmax() < 0)
            max_graph[i]->GetYaxis()->SetRangeUser(max_graph[i]->GetYaxis()->GetXmax(), 0);
//      max_graph[i]->GetXaxis()->SetLimits(-0.1, 1);
        gPad->Modified();
    }



    //
    //  DIFFERENCE
    //

    float min_diff[N][M], min_dunc[N][M], max_diff[N][M], max_dunc[N][M];
    TGraphErrors *min_diff_graph[N], *max_diff_graph[N];
    TCanvas *minDiffCanvas[N], *maxDiffCanvas[N];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)
            continue;

        float min_neg0 = hObserved[0][i]->Integral(1, 1);
        float min_pos0 = hObserved[0][i]->Integral(20, 20);
        float min_asy0 = (min_pos0 - min_neg0) / (min_pos0 + min_neg0);
        float min_unc0 = 1 / sqrt(min_neg0 + min_neg0);

        float max_neg0 = hObserved[0][i]->Integral(7, 10);
        float max_pos0 = hObserved[0][i]->Integral(11, 14);
        float max_asy0 = (max_pos0 - max_neg0) / (max_pos0 + max_neg0);
        float max_unc0 = 1 / sqrt(max_neg0 + max_neg0);

        for (unsigned j = 0; j < M; j++)
        {
            // Minimum cut difference
            float neg = hObserved[0][i]->Integral(1, 1 + j);
            float pos = hObserved[0][i]->Integral(20 - j, 20);

            float asy = (pos - neg) / (pos + neg);
            min_diff[i][j] = asy - min_asy0;

            float term1 = (4 * min_neg0 * min_pos0) / pow(min_neg0 + min_pos0, 3);
            float term2 = (4 * neg * neg) / pow(neg + pos, 3);
            float term3 = (4 * neg) / pow(neg + pos, 2);
            float term4 = (8 * min_neg0 * min_pos0) / (pow(min_neg0 + min_pos0, 2) * (neg + pos));

            min_dunc[i][j] = sqrt(term1 - term2 + term3 - term4);
//          cout << term1 << "\t" << term2 << "\t" << term3 << "\t" << term4 << endl;


            // Maximum cut difference
            neg = hObserved[0][i]->Integral(7 - j, 10);
            pos = hObserved[0][i]->Integral(11, 14 + j);

            asy = (pos - neg) / (pos + neg);
            max_diff[i][j] = asy - max_asy0;

            term1 = (4 * max_neg0 * max_pos0) / pow(max_neg0 + max_pos0, 3);
            term2 = (4 * neg * neg) / pow(neg + pos, 3);
            term3 = (4 * neg) / pow(neg + pos, 2);
            term4 = (8 * max_neg0 * max_pos0) / (pow(max_neg0 + max_pos0, 2) * (neg + pos));

            max_dunc[i][j] = sqrt(term1 - term2 + term3 - term4);
//          cout << term1 << "\t" << term2 << "\t" << term3 << "\t" << term4 << endl;
        }

        min_diff_graph[i] = new TGraphErrors(M, phi_min, min_diff[i], 0, min_dunc[i]);
        min_diff_graph[i]->SetTitle("");
        min_diff_graph[i]->GetXaxis()->SetTitle("\\min\\,|\\sin\\phi|");
        min_diff_graph[i]->GetYaxis()->SetTitle("\\Delta A_{" + lepChan[i] + "}");
        min_diff_graph[i]->GetYaxis()->SetTitleOffset(0.7 * lTitleOffsetY);
        min_diff_graph[i]->SetMarkerStyle(kFullCircle);
        min_diff_graph[i]->SetMarkerSize(2);
        min_diff_graph[i]->SetLineWidth(2);

        minDiffCanvas[i] = new TCanvas("diff_min_sinphi" + selection[i], "", lCanvasSize, lCanvasSize);
        minDiffCanvas[i]->cd();
        Facelift(minDiffCanvas[i]);
//      minDiffCanvas[i]->SetLeftMargin(1.2 * lCanvasMargin);
        minDiffCanvas[i]->SetCanvasSize(lCanvasSize, 0.5 * lCanvasSize);
        minDiffCanvas[i]->SetMargin(1.3*lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);
        min_diff_graph[i]->Draw("AP");

        min_diff_graph[i]->GetXaxis()->SetLimits(-0.05, 0.95);
        min_diff_graph[i]->GetYaxis()->SetRangeUser(-0.17, 0.17);
        Facelift(min_diff_graph[i]->GetXaxis());
        Facelift(min_diff_graph[i]->GetYaxis());


        float xfake[2] = {-0.5, 1.5}, yfake[2] = {0, 0}, errfake[2] = {min_unc0, min_unc0};
        TGraphErrors *min_unc_graph = new TGraphErrors(2, xfake, yfake, 0, errfake);
        min_unc_graph->SetFillColor(lPurple);
        min_unc_graph->SetFillStyle(3003);
        min_unc_graph->Draw("3 SAME");
        min_diff_graph[i]->Draw("P");
        gPad->Modified();






        max_diff_graph[i] = new TGraphErrors(X, phi_max, max_diff[i], 0, max_dunc[i]);
        max_diff_graph[i]->SetTitle("");
        max_diff_graph[i]->GetXaxis()->SetTitle("\\max\\,|\\sin\\phi|");
        max_diff_graph[i]->GetYaxis()->SetTitle("\\Delta A_{" + lepChan[i] + "}");
        max_diff_graph[i]->GetYaxis()->SetTitleOffset(0.7 * lTitleOffsetY);
        max_diff_graph[i]->SetMarkerStyle(kFullCircle);
        max_diff_graph[i]->SetMarkerSize(2);
        max_diff_graph[i]->SetLineWidth(2);

        maxDiffCanvas[i] = new TCanvas("diff_max_sinphi" + selection[i], "", lCanvasSize, lCanvasSize);
        maxDiffCanvas[i]->cd();
        Facelift(maxDiffCanvas[i]);
//      maxDiffCanvas[i]->SetLeftMargin(1.2 * lCanvasMargin);
        maxDiffCanvas[i]->SetCanvasSize(lCanvasSize, 0.5 * lCanvasSize);
        maxDiffCanvas[i]->SetMargin(1.3*lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);
        max_diff_graph[i]->Draw("AP");

        max_diff_graph[i]->GetYaxis()->SetRangeUser(-0.17, 0.17);
        Facelift(max_diff_graph[i]->GetXaxis());
        Facelift(max_diff_graph[i]->GetYaxis());


        errfake[0] = max_unc0; errfake[1] = max_unc0;
        TGraphErrors *max_unc_graph = new TGraphErrors(2, xfake, yfake, 0, errfake);
        max_unc_graph->SetFillColor(lRed);
        max_unc_graph->SetFillStyle(3003);
        max_unc_graph->Draw("3 SAME");
        max_diff_graph[i]->Draw("P");
        gPad->Modified();
    }






    ////
    ////
    ////    MAKE STACKS
    ////
    ////


    THStack *stack[2][N][2];
    float nObserved[2], nTotal[2];

    for (unsigned i = 0; i < N; i++)    // channel loop
    {
        if (i == EM)    // skip "2e2m" since we added those (gross)
            continue;

        nObserved[i] = hObserved[1][i]->Integral();
        nTotal[i] = hTotal[1][i][1]->Integral();

        for (unsigned h = 0; h < 2; h++)    // distribution loop
        {
            hObserved[h][i]->SetTitle("");
            hObserved[h][i]->SetStats(0);
            hObserved[h][i]->SetMarkerColor(kBlack);
            hObserved[h][i]->SetMarkerStyle(kFullCircle);
            hObserved[h][i]->SetMarkerSize(2);
            hObserved[h][i]->SetLineWidth(2);
            hObserved[h][i]->SetLineColor(kBlack);


            // Create and fill stack object
            hTotal[h][i][1]->Scale(nObserved[i] / nTotal[i]);

            for (unsigned k = 0; k < 2; k++)
            {
                stack[h][i][k] = new THStack(hname[h] + "_" + htype[k] + "_" + selection[i], "");

                for (unsigned j = 0; j < N_MC; j++) // sample loop
                {
                    if (!useDY && j == DY)
                        continue;

                    if (k == 1)
                    {
                        TH1 *mcClone = (TH1*) mcHist[h][i][j]->Clone();
                        mcClone->Scale(nObserved[i] / nTotal[i]);

                        stack[h][i][k]->Add(mcClone);
                    }
                    else
                    {
                        mcHist[h][i][j]->SetTitle("");
                        mcHist[h][i][j]->SetStats(0);
                        mcHist[h][i][j]->SetFillColor(COLOR[j]);
                        mcHist[h][i][j]->SetLineColor(COLOR[j]);

                        stack[h][i][k]->Add(mcHist[h][i][j]);
                    }
                }
            }
        }
    }



    //
    //  LEGEND
    //

    float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
    float RightPosition = 1,        RightMargin = -lLegendMargin;
    float TopPosition = 1,          TopMargin = -lLegendMargin;
    float BottomPosition = TopPosition - 0.085 * 5;//0.065 * 6.;
    float BottomMargin = 2. * lCanvasMargin - lLegendMargin;
    TLegend *legend = new TLegend(LeftPosition + LeftMargin, BottomPosition - TopMargin,
            TopPosition + TopMargin, TopPosition + TopMargin);

    TString dentry = _sp+"\\mbox{Data}";
    TString lentry[5] = {_sp+_ZZ+_to+_4l, _sp+_Z+_to+_ll, _sp+_H, _sp+_ttbar, _sp+_V+_V};
    Int_t lfill[5] = {lLightBlue, lYellow, lPurple, lGreen, lRed};

    TH1D* dummy = new TH1D(dentry, "", 1, 0, 1);
    dummy->SetMarkerColor(kBlack);
    dummy->SetMarkerStyle(kFullCircle);
    dummy->SetMarkerSize(2);
    dummy->SetLineWidth(2);
    dummy->SetLineColor(kBlack);
    legend->AddEntry(dummy, dentry, "LP");

    for (unsigned h = 0; h < 5; h++)
    {
        if (!useDY && (h == 1))
            continue;

        TH1D* hist = new TH1D(lentry[h], "", 1, 0, 1);
        hist->SetFillColor(lfill[h]);
        hist->SetLineColor(lfill[h]);
        legend->AddEntry(hist, lentry[h], "F");
    }
    Facelift(legend);



    //
    //  OUTPUT FILE
    //

    TString typeName = useDY ? "incDY" : "noDY";
    TString outName = "asymmetry_" + typeName + "_2017.root";
    TFile *outFile  = new TFile(outName, "RECREATE");
    outFile->cd();

    asyCanvas->Write();



    //
    //  DRAW CANVASES
    //

    TCanvas *canvas[2][N][2];
    TRatioPlot *ratio[2][N][2];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)    // skip "2e2m" since we added those (gross)
            continue;

        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        for (unsigned h = 0; h < 2; h++)
        {
            for (unsigned k = 0; k < 2; k++)
            {
                canvas[h][i][k] = new TCanvas(hname[h] + "_" + htype[k] + "_" + selection[i],
                        "", lCanvasSize, lCanvasSize);

                Facelift(canvas[h][i][k]);
                canvas[h][i][k]->cd();

                hObserved[h][i]->GetXaxis()->SetTitle("\\sin\\phi");
                hTotal[h][i][k]->SetLineColor(0);

                hObserved[h][i]->SetMinimum(0);
                hTotal[h][i][k]->SetMinimum(0);

                ratio[h][i][k] = new TRatioPlot(hObserved[h][i], hTotal[h][i][k], "divsym");
                //          ratio[h][i][k]->SetGraphDrawOpt("B");
                ratio[h][i][k]->SetH1DrawOpt("E");
                ratio[h][i][k]->SetH2DrawOpt("E");
                ratio[h][i][k]->SetSeparationMargin(0.0);
                ratio[h][i][k]->Draw();

                TPad *upper = ratio[h][i][k]->GetUpperPad(), *lower = ratio[h][i][k]->GetLowerPad();
                upper->cd();

                stack[h][i][k]->Draw("HIST SAME");
                Facelift(stack[h][i][k]);
                hObserved[h][i]->Draw("E SAME");
                //          upper->RedrawAxis();
                upper->Modified();

                Facelift(ratio[h][i][k]->GetLowerRefXaxis());
                Facelift(ratio[h][i][k]->GetLowerRefYaxis());
                lower->SetBottomMargin(3 * lCanvasMargin);
                lower->Modified();

                legend->Draw();

                canvas[h][i][k]->Write();
            }
        }
        minCanvas[i]->Write();
        minDiffCanvas[i]->Write();
        maxCanvas[i]->Write();
        maxDiffCanvas[i]->Write();
    }
    cout << "done!" << endl << endl;
    outFile->Close();

    cout << "Wrote canvases to " << outName << endl << endl << endl;
}
