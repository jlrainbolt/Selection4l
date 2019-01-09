// STL
#include <iostream>
#include <vector>
#include <tuple>

//  ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TUnfoldDensity.h"

// Cuts
#include "Cuts2017.hh"

using namespace std;


/*
 **  Unfold4l
 **
 **  Unfolds all differential distributions.  Takes signal input from "migration_" file.
 */

void Unfold4l()
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                    L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = { "4l",   "4m",   "2m2e", "2e2m", "4e"};

    // Signal scale factor
    float sf[N];
    for (unsigned i = 1; i < N; i++)
    {
        float LUMI;
        if      ((i == M4) || (i == ME))
            LUMI = MUON_TRIG_LUMI;
        else if ((i == E4) || (i == EM))
            LUMI = ELEC_TRIG_LUMI * ELEC_TRIG_SF;
        sf[i] = LUMI * 1000 * XSEC[ZZ] / NGEN[ZZ];
    }



    //
    //  GET HISTOGRAM KEYS
    //

    vector<TString> hname;
    TString inName = "migration_zz_4l.root";
    TFile *inFile = TFile::Open(inName);
    cout << "Opened " << inName << endl;

    TDirectory *keyDir = inFile->GetDirectory("/" + selection[M4] + "/2d", kTRUE, "GetDirectory");
    TKey *histKey;
    TIter next(keyDir->GetListOfKeys());
    while ((histKey = (TKey*) next()))
    {
        TString hname_ = histKey->GetName();
        hname_.Resize(hname_.Length() - (3));    // truncate before suffix
        hname.push_back(hname_);
    }
    cout << "Got histogram keys" << endl;

    const unsigned H = hname.size();



    //
    //  SIGNAL MONTE CARLO
    // 

    TH1 *h_gen[N][H], *h_reco[N][H], *h_gen_sc[N][H], *h_reco_sc[N][H];
    TH2 *m_A[N][H];

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        TH1 *h_gen_, *h_reco_;
        TH2 *m_A_;

        for (unsigned i = M4; i < N; i++)    // channel loop
        {
            // Get histograms

            inFile->GetObject(selection[i] + "/gen/" + hname[h] + "_gen", h_gen_);
            h_gen_->SetDirectory(0);
            h_gen_->SetStats(0);
            h_gen[i][h] = (TH1*) h_gen_->Clone();
            h_gen_->Scale(sf[i]);
            h_gen_sc[i][h] = (TH1*) h_gen_->Clone();

            inFile->GetObject(selection[i] + "/reco/" + hname[h] + "_reco", h_reco_);
            h_reco_->SetDirectory(0);
            h_reco_->SetStats(0);
            h_reco[i][h] = (TH1*) h_reco_->Clone();
            h_reco_->Scale(sf[i]);
            h_reco_sc[i][h] = (TH1*) h_reco_->Clone();

            inFile->GetObject(selection[i] + "/2d/" + hname[h] + "_2d", m_A_);
            m_A_->SetDirectory(0);
            m_A_->SetStats(0);
            m_A[i][h] = m_A_;
        }
    }
    cout << "Got signal MC histograms and matrices" << endl << endl;




    //
    //  DATA
    //


    TString prefix  = "unscaled4l";

    // Muon file
    TString muName = prefix + "_" + MU_SUFF + ".root";
    TFile *muFile = TFile::Open(muName);
    cout << "Opened " << muName << endl;

    // Electron file
    TString elName = prefix + "_" + EL_SUFF + ".root";
    TFile *elFile = TFile::Open(elName);
    cout << "Opened " << elName << endl;


    // Now get the data histograms
    TH1 *h_data[N][H];

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {   
        TH1 *hist;

        for (unsigned i = M4; i < N; i++)    // channel loop
        {   
            if      ((i == M4) || (i == ME))
                muFile->GetObject(selection[i] + "/" + hname[h] + "_" + MU_SUFF, hist);
            else if ((i == E4) || (i == EM))
                elFile->GetObject(selection[i] + "/" + hname[h] + "_" + EL_SUFF, hist);

            hist->SetDirectory(0);
            h_data[i][h] = hist;
        }
    }
    muFile->Close();
    elFile->Close();

    cout << "Got data histograms" << endl << endl;
    


    //
    //  ADD CHANNELS
    //

    cout << endl << "Combining..." << flush;
    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        h_data[L4][h] = (TH1*) h_data[1][h]->Clone();
        h_gen[L4][h] = (TH1*) h_gen[1][h]->Clone();
        h_gen_sc[L4][h] = (TH1*) h_gen_sc[1][h]->Clone();
        h_reco[L4][h] = (TH1*) h_reco[1][h]->Clone();
        h_reco_sc[L4][h] = (TH1*) h_reco_sc[1][h]->Clone();
        m_A[L4][h] = (TH2*) m_A[1][h]->Clone();

        for (unsigned i = 2; i < N; i++)
        {
            h_data[L4][h]->Add(h_data[i][h]);
            h_gen[L4][h]->Add(h_gen[i][h]);
            h_gen_sc[L4][h]->Add(h_gen_sc[i][h]);
            h_reco[L4][h]->Add(h_reco[i][h]);
            h_reco_sc[L4][h]->Add(h_reco_sc[i][h]);
            m_A[L4][h]->Add(m_A[i][h]);
        }
    }

    // Calculate overall scale factor
    sf[L4] = h_reco_sc[L4][0]->Integral() / h_reco[L4][0]->Integral();



    //
    //  UNFOLD
    //

    TH1 *h_result[N][H];
    TH2 *m_response[N][H], *m_unfolding[N][H];

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        cout << "Unfolding " << hname[h] << " distribution..." << endl;

        for (unsigned i = 0; i < N; i++)
        {
            TUnfoldDensity unfolder(m_A[i][h], TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone);

//          unfolder.SetInput(h_data[i][h], sf[i]);
            unfolder.SetBias(h_gen_sc[i][h]);
            unfolder.SetInput(h_data[i][h]);
            unfolder.DoUnfold(0);

            h_result[i][h] = unfolder.GetOutput("h_result_" + hname[h] + "_" + selection[i], "");
            h_result[i][h]->SetStats(0);
            h_result[i][h]->GetXaxis()->SetTitle(h_gen[i][h]->GetXaxis()->GetTitle());

            m_response[i][h] = unfolder.GetProbabilityMatrix("m_resp_"+hname[h]+"_"+selection[i]);
            m_response[i][h]->SetStats(0);

            m_unfolding[i][h] = unfolder.GetRhoIJtotal("m_unf_" + hname[h] + "_" + selection[i]);
            m_unfolding[i][h]->SetStats(0);

            cout << endl;
        }
    }



    //
    //  DRAW
    //

    cout << "Drawing distributions..." << endl << endl;

    TCanvas *c[N][H];
    for (unsigned h = 0; h < H; h++)
    {
        for (unsigned i = 0; i < N; i++)
        {
            c[i][h] = new TCanvas("c_" + hname[h] + "_" + selection[i], "", 800, 600);

            h_data[i][h]->SetLineColor(kBlack);
            h_data[i][h]->SetLineWidth(2);
            h_data[i][h]->SetMarkerStyle(20);
            h_data[i][h]->SetMarkerSize(2);

            h_reco_sc[i][h]->SetLineColor(kViolet);
            h_reco_sc[i][h]->SetLineWidth(2);

            h_result[i][h]->SetLineColor(kRed);
            h_result[i][h]->SetLineWidth(2);
            h_result[i][h]->SetMarkerColor(kRed);
            h_result[i][h]->SetMarkerStyle(20);
            h_result[i][h]->SetMarkerSize(2);

            h_gen_sc[i][h]->SetLineColor(kBlue);
            h_gen_sc[i][h]->SetLineWidth(2);

            TLegend *l = new TLegend(0.8, 0.8, 1, 1);
            l->AddEntry(h_gen_sc[i][h], "Gen (scaled)", "L");
            l->AddEntry(h_reco_sc[i][h], "Reco (scaled)", "L");
            l->AddEntry(h_data[i][h], "Data", "LP");
            l->AddEntry(h_result[i][h], "Unfolded data", "LP");

            c[i][h]->cd();
            h_result[i][h]->Draw();
            h_gen_sc[i][h]->Draw("SAME");
            h_result[i][h]->Draw("SAME");
            h_reco_sc[i][h]->Draw("SAME");
            h_data[i][h]->Draw("SAME");
            l->Draw();
        }
    }



    //
    //  WRITE
    //

    TString outName = "unfolding_" + YEAR_STR + ".root";
    TFile *outFile = new TFile(outName, "RECREATE");

    for (unsigned i = 0; i < N; i++)
    {
        outFile->mkdir(selection[i]);

        for (unsigned h = 0; h < H; h++)
        {
            outFile->mkdir(selection[i] + "/" + hname[h]);
            outFile->cd(selection[i] + "/" + hname[h]);

            c[i][h]->Write();

            m_A[i][h]->SetName("m_A_" + hname[h] + "_" + selection[i]); 
            m_A[i][h]->Write();
            m_response[i][h]->Write();
            m_unfolding[i][h]->Write();

            h_gen[i][h]->SetName("h_gen_" + hname[h] + "_" + selection[i]); 
            h_gen[i][h]->Write();
            h_reco[i][h]->SetName("h_reco_" + hname[h] + "_" + selection[i]); 
            h_reco[i][h]->Write();
            h_data[i][h]->SetName("h_data_" + hname[h] + "_" + selection[i]); 
            h_data[i][h]->Write();
            h_result[i][h]->Write();
        }
    }
    outFile->Close();

    cout << "Wrote results to " << outName << endl << endl;
}
