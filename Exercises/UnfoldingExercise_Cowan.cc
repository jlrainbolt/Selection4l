#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TUnfoldDensity.h"
#include "TSpline.h"

using namespace std;

void UnfoldingExercise_Cowan()
{


    //
    //  CREATE DISTRIBUTIONS
    // 

    int bins_gen = 20, bins_reco = 40;

    TH1D *h_gen = new TH1D("h_gen", "gen", bins_gen, 0, 1);     // "true" histogram
    TH1D *h_reco = new TH1D("h_reco", "reco", bins_reco, 0, 1); // "expectation value" histogram
    TH1D *h_data = new TH1D("h_data", "data", bins_reco, 0, 1); // "observed data" histogram

    // Migration matrix: "true" on x axis, "expected" on y axis
    TH2D *h_A = new TH2D("h_A", "A", bins_reco, 0, 1, bins_gen, 0, 1);
    h_A->GetXaxis()->SetTitle("reco");
    h_A->GetYaxis()->SetTitle("gen");

    float mean_gen[2] = {0.3, 0.7};
    float sigma_gen = 0.12, sigma_smear = 0.0375;

    TRandom3 rng(27);
    int n_entries = 12500;


    // Create gen histogram and smear to get reco histogram
    for (unsigned i = 0; i < 2; i++)
    {
        for (unsigned j = 0; j < n_entries; j++)
        {
            float gen_j = rng.Gaus(mean_gen[i], sigma_gen);
            float reco_j = gen_j + rng.Gaus(0, sigma_smear);

            h_gen->Fill(gen_j);
            h_reco->Fill(reco_j);

            h_A->Fill(reco_j, gen_j);
        }
    }

    for (unsigned k = 0; k <= bins_reco + 1; k++)
    {
        float mu_k = h_reco->GetBinContent(k);
        float n_k = rng.Poisson(mu_k);

        h_data->SetBinContent(k, n_k);
        h_data->SetBinError(k, sqrt(n_k));
    } 



    //
    //  UNFOLD
    //

    // kRegModeCurvature: k = 2 (second derivative)
//  TUnfoldDensity unfolder(h_A, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone);
    TUnfoldDensity unfolder(h_A, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature);
    unfolder.SetInput(h_data);
    unfolder.DoUnfold(pow(10, -3.881395));

    TH1 *h_unfolded = unfolder.GetOutput("Unfolded");



    //
    //  SCAN TAU
    //
  
    TSpline *spline;

    unfolder.ScanTau(20, 1e-7, 1, &spline);
    spline->Draw();
  


    //
    //  DRAW
    //

    TCanvas *canvas = new TCanvas("overlay", "overlay", 800, 600);
    h_gen->SetLineColor(kBlack);
    h_unfolded->SetLineColor(kRed);

    canvas->cd();
    h_unfolded->Draw();
    h_gen->Draw("SAME");
    h_unfolded->Draw("SAME");

 

    //
    //  WRITE
    //

    TFile *outFile = new TFile("unfolding_cowan.root", "RECREATE");

    outFile->mkdir("Peaks");
    outFile->cd("Peaks");

    h_gen->Write();
    h_reco->Write();
    h_data->Write();

    h_A->Write();
    h_unfolded->Write();

    canvas->Write();

    outFile->Close();
}
