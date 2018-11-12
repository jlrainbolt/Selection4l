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
    cout << "BEFORE CORRECTIONS" << endl;
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
    float xbin[4] = {1, 3, 4, 5},   yerr[4] = {unc[L4], unc[M4], unc[ME], unc[E4]};
    float yval[4] = {asymmetry[L4], asymmetry[M4], asymmetry[ME], asymmetry[E4]};
    TGraphErrors *graph = new TGraphErrors(4, xbin, yval, 0, yerr);
    graph->SetNameTitle("asymmetry", "");
    float xedge[6] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
    graph->GetXaxis()->Set(5, xedge);
    graph->GetXaxis()->SetBinLabel(1, _4l);
    graph->GetXaxis()->SetBinLabel(3, _4mu);
    graph->GetXaxis()->SetBinLabel(4, _2mu2e);
    graph->GetXaxis()->SetBinLabel(5, _4e);
    graph->GetXaxis()->SetTicks("-");
    graph->SetMarkerStyle(8);
//  graph->SetLineWidth(0);
//  graph->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, _4l);



    //
    //  WITH CUT
    //

    const unsigned M = 10;
    float asymmetry_cut[N][M], unc_cut[N][M];
    float phi_min[M] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};
    TGraphErrors *asymm_func[M];
    TCanvas *fCanvas[N];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)
            continue;

        for (unsigned j = 0; j < M; j++)
        {
            float neg = hObserved[0][i]->Integral(1, 1 + j);
            float pos = hObsMinusBg[0][i]->Integral(20 - j, 20);
            float half = (pos + neg) / 2;

            asymmetry_cut[i][j] = (pos - neg) / (pos + neg);
            unc_cut[i][j] = sqrt(4 * pow(half, 2) / pow(neg + pos, 3));
        }
        asymm_func[i] = new TGraphErrors(M, phi_min, asymmetry_cut[i], 0, unc_cut[i]);
        asymm_func[i]->SetTitle("");
        asymm_func[i]->GetXaxis()->SetTitle("\\min\\,|\\sin\\phi|");
        asymm_func[i]->GetYaxis()->SetTitle("A_{" + lepChan[i] + "}");
        asymm_func[i]->GetYaxis()->SetTitleOffset(0.7 * lTitleOffsetY);
        asymm_func[i]->SetMarkerStyle(kFullCircle);
        asymm_func[i]->SetMarkerSize(2);
        asymm_func[i]->SetLineWidth(2);

        fCanvas[i] = new TCanvas("abs_sinphi_min_" + selection[i], "", lCanvasSize, lCanvasSize);
        fCanvas[i]->cd();
        Facelift(fCanvas[i]);
//      fCanvas[i]->SetLeftMargin(1.2 * lCanvasMargin);
        fCanvas[i]->SetCanvasSize(lCanvasSize, 0.5 * lCanvasSize);
        fCanvas[i]->SetMargin(1.3*lCanvasMargin, lCanvasMargin/2, 1.8*lCanvasMargin, lCanvasMargin);
        asymm_func[i]->Draw("AP");

        Facelift(asymm_func[i]->GetXaxis());
        Facelift(asymm_func[i]->GetYaxis());
        if (asymm_func[i]->GetYaxis()->GetXmax() < 0)
            asymm_func[i]->GetYaxis()->SetRangeUser(asymm_func[i]->GetYaxis()->GetXmin(), 0);
        asymm_func[i]->GetXaxis()->SetLimits(-0.1, 1);
        gPad->Modified();
    }



    //
    //  DIFFERENCE
    //

    float var_A0[N], negObs0[N], posObs0[N];
    float unc_deltaA[N][M];
    TGraph *deltaA_unc[N];

    for (unsigned i = 0; i < N; i++)
    {
        if (i == EM)
            continue;

        negObs0[i] = hObserved[0][i]->Integral(7, 10);
        posObs0[i] = hObserved[0][i]->Integral(11, 14);

        var_A0[i] = 4 * negObs0[i] * posObs0[i] / pow(negObs0[i] + posObs0[i], 3);

        for (unsigned j = 0; j < M; j++)
        {
            float negObs = hObserved[0][i]->Integral(1, 2 + j);
            float posObs = hObserved[0][i]->Integral(19 - j, 20);

            float term1 = var_A0[i];
            float term2 = 4 * pow(posObs, 2) / pow(posObs + negObs, 3);
            float term3 = 4 * posObs / pow(posObs + negObs, 2);
            float term4 = 8 * negObs0[i] * posObs0[i] / pow(posObs0[i] + negObs0[i], 2);
            term4 /= (posObs + negObs);

            unc_deltaA[i][j] = sqrt(term1 - term2 + term3 - term4);
        }

        deltaA_unc[i] = new TGraph(M, phi_min, unc_deltaA[i]);
        deltaA_unc[i]->SetNameTitle("unc_deltaA_" + selection[i], "unc_deltaA");
        deltaA_unc[i]->GetXaxis()->SetTitle("max |sin phi|");
        deltaA_unc[i]->GetYaxis()->SetTitle("unc. in A - A(x = 0.4)");
        deltaA_unc[i]->SetMarkerStyle(8);
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

    TCanvas *gCanvas = new TCanvas("asymmetry_canvas", "", lCanvasSize, lCanvasSize);
    gCanvas->cd();
    Facelift(gCanvas);
    graph->Draw();
    Facelift(graph->GetXaxis());
    Facelift(graph->GetYaxis());
    gPad->Modified();

    graph->Write();
    gCanvas->Write();



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
        fCanvas[i]->Write();
        deltaA_unc[i]->Write();
    }
    cout << "done!" << endl << endl;
    outFile->Close();

    cout << "Wrote canvases, graphs to " << outName << endl << endl << endl;
}
