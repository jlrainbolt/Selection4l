#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TUnfoldDensity.h"
#include "TSpline.h"

using namespace std;

void UnfoldingExercise_Landau()
{

    unsigned N = 10;


    //
    //  CREATE HISTOGRAMS
    // 

    int bins_gen = 10, bins_reco = 20;
    float xmin = 0, xmax = 60;

    // "MC"
    TH1D *h_gen = new TH1D("h_gen", "", bins_gen, xmin, xmax);
    TH1D *h_reco = new TH1D("h_reco", "", bins_reco, xmin, xmax);
    TH1D *h_reco_coarse = new TH1D("h_reco_coarse", "", bins_gen, xmin, xmax);

    // Migration matrix: "true" on x axis, "expected" on y axis
    TH2D *h_A = new TH2D("h_A", "A", bins_reco, xmin, xmax, bins_gen, xmin, xmax);
    h_A->GetXaxis()->SetTitle("reco");
    h_A->GetYaxis()->SetTitle("gen");

    // "Data"
    TH1D *h_data[N+1], *h_data_coarse[N+1];
    for (unsigned i = 0; i < N; i++)
    {
        TString index;
        index.Form("%i", i);

        h_data[i] = new TH1D("h_data_" + index, "", bins_reco, xmin, xmax);
        h_data_coarse[i] = new TH1D("h_data_coarse_" + index, "", bins_gen, xmin, xmax);
    }
    TH1 *h_unfolded[N+1];



    //
    //  FILL HISTOGRAMS
    //

    float mu_gen = 20, sigma_gen = 5;
    float sigma_smear = 2;

    TRandom3 rng(27);
    int n_mc = 57500, n_data = 430;
//  float scale = float(n_mc) / float(n_data);
    float scale = float(n_data) / float(n_mc);


    // MC: create gen histogram and smear to get reco histogram
    for (unsigned j = 0; j < n_mc; j++)
    {
        float gen_j = rng.Landau(mu_gen, sigma_gen);
        float reco_j = gen_j + rng.Gaus(0, sigma_smear);

        h_gen->Fill(gen_j);
        h_reco->Fill(reco_j);
        h_reco_coarse->Fill(reco_j);

        h_A->Fill(reco_j, gen_j);
    }

    // Create 10 iterations of the data histograms
    for (unsigned i = 0; i < N; i++)
    {
        for (unsigned j = 0; j < n_data; j++)
        {
            float gen_j = rng.Landau(mu_gen, sigma_gen);
            float reco_j = gen_j + rng.Gaus(0, sigma_smear);

            h_data[i]->Fill(reco_j);
            h_data_coarse[i]->Fill(reco_j);
        }
//      h_data[i]->Scale(scale);
//      h_data_coarse[i]->Scale(scale);
    }



    //
    //  UNFOLD
    //

    TUnfoldDensity unfolder(h_A, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone);

    for (unsigned i = 0; i < N; i++)
    {
        cout << "Unfolding distribution " << i << "..." << endl;

        TString index;
        index.Form("%i", i);

//      unfolder.SetInput(h_data[i], scale);
        unfolder.SetInput(h_reco);
//      unfolder.DoUnfold(pow(10, -0.875));     // kRegModeCurvature
//      unfolder.DoUnfold(pow(10, -0.875));     // kRegModeDerivative
        unfolder.DoUnfold(0);                   // kRegModeNone
        h_unfolded[i] = unfolder.GetOutput("Unfolded");
        h_unfolded[i]->SetName("h_unfolded_" + index);

        cout << endl;
    }

    TSpline *spline;
    unfolder.ScanTau(20, 1e-7, 1, &spline);
    spline->Draw();


    // Find mean of data histograms
    h_data[N] = (TH1D*) h_data[0]->Clone();
    h_data_coarse[N] = (TH1D*) h_data_coarse[0]->Clone();
    h_unfolded[N] = (TH1*) h_unfolded[0]->Clone();
    for (unsigned i = 1; i < N; i++)
    {
        h_data[N]->Add(h_data[i]);
        h_data_coarse[N]->Add(h_data_coarse[i]);
        h_unfolded[N]->Add(h_unfolded[i]);
    }
    h_data[N]->Scale(pow(float(N), -1));
    h_data[N]->SetName("h_data_mean");
    h_data_coarse[N]->Scale(pow(float(N), -1));
    h_data_coarse[N]->SetName("h_data_coarse_mean");
    h_unfolded[N]->Scale(pow(float(N), -1));
    h_unfolded[N]->SetName("h_unfolded_mean");


//  TUnfoldDensity unfolder(h_A, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone);
    // kRegModeCurvature: k = 2 (second derivative)
//  TH1D* h_bias = (TH1D*) h_gen->Clone();
//  h_bias->Scale(scale);
//  unfolder.SetBias(h_bias);
//  unfolder.DoUnfold(pow(10, -2.769796));


/*
    //
    //  SCAN TAU
    //
  
*/


    //
    //  DRAW
    //

    TCanvas *c_mc = new TCanvas("mc", "mc", 800, 600);
    h_gen->SetLineColor(kBlack);
    h_gen->SetLineWidth(2);
    h_gen->SetMarkerStyle(20);
    h_gen->SetMarkerSize(1);
    h_reco_coarse->SetLineColor(kViolet);
    h_reco_coarse->SetLineWidth(2);

    TLegend *l_mc = new TLegend(0.65, 0.65, 0.95, 0.95);
    l_mc->AddEntry(h_gen, "Gen", "LP");
    l_mc->AddEntry(h_reco_coarse, "Reco (rebinned)", "L");

    c_mc->cd();
    h_gen->DrawCopy();
    h_reco_coarse->Draw("SAME");
    l_mc->Draw();

//  h_gen->Scale(scale);


    TCanvas *c_data[N+1];
    for (unsigned i = 0; i <= N; i++)
    {
        TString index;
        index.Form("%i", i);

        c_data[i] = new TCanvas("data_" + index, "data", 800, 600);
        h_unfolded[i]->SetLineColor(kRed);
        h_unfolded[i]->SetLineWidth(2);
        h_data_coarse[i]->SetLineColor(kBlue);
        h_data_coarse[i]->SetLineWidth(2);

        TLegend *l_data = new TLegend(0.65, 0.65, 0.95, 0.95);
        l_data->AddEntry(h_gen, "Gen MC", "LP");
        l_data->AddEntry(h_data_coarse[i], "Data (rebinned)", "L");
        l_data->AddEntry(h_unfolded[i], "Unfolded data", "L");

        c_data[i]->cd();
        h_gen->Draw();
//      h_data_coarse[i]->Draw();
//      h_gen->Draw("SAME");
//      h_data_coarse[i]->Draw("SAME");
        h_unfolded[i]->Draw("SAME");
        l_data->Draw();
    }
    c_data[N]->SetName("data_mean");


/*
    //
    //  TEST
    //

    double chi2_before = h_gen->Chi2Test(h_data_coarse, "UU NORM");
    double chi2_after = h_gen->Chi2Test(h_unfolded, "UW");

    cout << endl;
    cout << "chi2 before unfolding: " << chi2_before << endl;
    cout << "chi2 after unfolding: " << chi2_after << endl;
*/


    //
    //  WRITE
    //

    TFile *outFile = new TFile("unfolding_landau.root", "RECREATE");

    outFile->mkdir("Histograms");
    outFile->cd("Histograms");

    h_gen->Write();
    h_reco->Write();
    h_A->Write();
    for (unsigned i = 0; i <= N; i++)
    {
        h_data[i]->Write();
        h_unfolded[i]->Write();
    }

    outFile->mkdir("Canvases");
    outFile->cd("Canvases");

    c_mc->Write();
    for (unsigned i = 0; i <= N; i++)
        c_data[i]->Write();

    outFile->Close();
}
