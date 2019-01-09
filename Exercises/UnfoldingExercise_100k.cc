#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TUnfoldDensity.h"
#include "TSpline.h"

using namespace std;

void UnfoldingExercise_100k()
{

    unsigned N = 10;


    //
    //  CREATE HISTOGRAMS
    // 

    // Binning, etc.
    int bins_gen[N], bins_reco[N];
    for (unsigned i = 0; i < N; i++)
    {
        bins_gen[i] = 5 * (i + 1);
        bins_reco[i] = 2 * bins_gen[i];
    }
    float xmin = 0, xmax = 60;

    // "MC"
    TH1D *h_gen[N], *h_reco[N], *h_reco_coarse[N];
    for (unsigned i = 0; i < N; i++)
    {
        TString index;
        index.Form("%i", bins_gen[i]);

        h_gen[i] = new TH1D("h_gen_" + index, "", bins_gen[i], xmin, xmax);
        h_reco[i] = new TH1D("h_reco_" + index, "", bins_reco[i], xmin, xmax);
        h_reco_coarse[i] = new TH1D("h_reco_coarse_" + index, "", bins_gen[i], xmin, xmax);
    }

    // Migration matrix: "true" on x axis, "expected" on y axis
    TH2D *m_A[N];
    for (unsigned i = 0; i < N; i++)
    {
        TString index;
        index.Form("%i", bins_gen[i]);

        m_A[i] = new TH2D("m_A_" + index, "", bins_gen[i], xmin, xmax, bins_reco[i], xmin, xmax);
        m_A[i]->GetXaxis()->SetTitle("gen");
        m_A[i]->GetYaxis()->SetTitle("reco");
        m_A[i]->SetStats(0);
    }

    // "Data"
    TH1D *h_data[N], *h_data_coarse[N];
    for (unsigned i = 0; i < N; i++)
    {
        TString index;
        index.Form("%i", bins_gen[i]);

        h_data[i] = new TH1D("h_data_" + index, "", bins_reco[i], xmin, xmax);
        h_data_coarse[i] = new TH1D("h_data_coarse_" + index, "", bins_gen[i], xmin, xmax);
    }



    //
    //  FILL HISTOGRAMS
    //

    float mu_gen = 20, sigma_gen = 5;
    float sigma_smear = 2;

    TRandom3 rng(27);
    int n_mc = 100000, n_data = 100000;
//  float scale = float(n_data) / float(n_mc);

    // MC: create gen histogram and smear to get reco histogram
    cout << "Generating " << n_mc << " MC events..." << endl << endl;

    for (unsigned j = 0; j < n_mc; j++)
    {
        float gen_j = rng.Landau(mu_gen, sigma_gen);
        float reco_j = gen_j + rng.Gaus(0, sigma_smear);

        for (unsigned i = 0; i < N; i++)
        {
            h_gen[i]->Fill(gen_j);
            h_reco[i]->Fill(reco_j);
            h_reco_coarse[i]->Fill(reco_j);

            m_A[i]->Fill(gen_j,reco_j);
        }
    }

    // Create data histograms
    cout << "Generating " << n_data << " data events..." << endl << endl;

    for (unsigned j = 0; j < n_data; j++)
    {
        float gen_j = rng.Landau(mu_gen, sigma_gen);
        float reco_j = gen_j + rng.Gaus(0, sigma_smear);

        for (unsigned i = 0; i < N; i++)
        {
            h_data[i]->Fill(reco_j);
            h_data_coarse[i]->Fill(reco_j);
        }
    }



    //
    //  UNFOLD
    //

    TH1 *h_result[N];
    TH2 *m_response[N], *m_unfolding[N];

    for (unsigned i = 0; i < N; i++)
    {
        cout << "Unfolding distribution " << i << "..." << endl;

        TUnfoldDensity unfolder(m_A[i], TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone);

        TString index;
        index.Form("%i", bins_gen[i]);

        unfolder.SetInput(h_reco[i]);
        unfolder.DoUnfold(0);
        h_result[i] = unfolder.GetOutput("Unfolded");
        h_result[i]->SetName("h_result_" + index);

        m_response[i] = unfolder.GetProbabilityMatrix("m_response_" + index);
        m_response[i]->SetStats(0);
        m_unfolding[i] = unfolder.GetRhoIJtotal("m_unfolding_" + index);
        m_unfolding[i]->SetStats(0);

        cout << endl;
    }



    //
    //  DRAW
    //

    cout << "Drawing distributions..." << endl << endl;

    TCanvas *c[N];
    for (unsigned i = 0; i < N; i++)
    {
        TString index;
        index.Form("%i", bins_gen[i]);

        c[i] = new TCanvas("c_" + index, index + " bins", 800, 600);

        h_gen[i]->SetLineColor(kBlack);
        h_gen[i]->SetLineWidth(2);
        h_gen[i]->SetMarkerStyle(20);
        h_gen[i]->SetMarkerSize(2);
        h_gen[i]->SetTitle(index + " bins");

        h_reco_coarse[i]->SetLineColor(kViolet);
        h_reco_coarse[i]->SetLineWidth(2);

        h_result[i]->SetLineColor(kRed);
        h_result[i]->SetLineWidth(2);

        h_data_coarse[i]->SetLineColor(kBlue);
        h_data_coarse[i]->SetLineWidth(2);

        TLegend *l = new TLegend(0.65, 0.65, 0.95, 0.95);
        l->AddEntry(h_gen[i], "Gen", "LP");
        l->AddEntry(h_reco_coarse[i], "Reco (rebinned)", "L");
        l->AddEntry(h_data_coarse[i], "Data (rebinned)", "L");
        l->AddEntry(h_result[i], "Unfolded data", "L");

        c[i]->cd();
        h_gen[i]->Draw();
        h_reco_coarse[i]->Draw("SAME");
        h_data_coarse[i]->Draw("SAME");
        h_result[i]->Draw("SAME");
        l->Draw();
    }



    //
    //  WRITE
    //

    TFile *outFile = new TFile("unfolding_100k.root", "RECREATE");

    for (unsigned i = 0; i < N; i++)
    {
        TString index;
        index.Form("%i", bins_gen[i]);

        outFile->mkdir(index + " bins");
        outFile->cd(index + " bins");

        c[i]->Write();
        m_A[i]->Write();
        m_response[i]->Write();
        m_unfolding[i]->Write();

        h_gen[i]->Write();
        h_reco[i]->Write();
        h_reco_coarse[i]->Write();
        h_data[i]->Write();
        h_data_coarse[i]->Write();
        h_result[i]->Write();
    }

    outFile->Close();
}
