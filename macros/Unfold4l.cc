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
#include "TMatrixTBase.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldResponse.h"

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
    }



    //
    //  GET HISTOGRAM KEYS
    //

    vector<TString> hname = {"b_l1p"};
//  vector<TString> hname;
    TString inName = "migration_" + YEAR_STR + "_zz_4l.root";
    TFile *inFile = TFile::Open(inName);
    cout << "Opened " << inName << endl;

//  TDirectory *keyDir = inFile->GetDirectory("/" + selection[M4] + "/2d", kTRUE, "GetDirectory");
//  TKey *histKey;
//  TIter next(keyDir->GetListOfKeys());
//  while ((histKey = (TKey*) next()))
//  {
//      TString hname_ = histKey->GetName();
//      hname_.Resize(hname_.Length() - (3));    // truncate before suffix
//      hname.push_back(hname_);
//  }
//  cout << "Got histogram keys" << endl;

    const unsigned H = hname.size();



    //
    //  BUILD RESPONSE MATRIX
    // 

    TH1 *h_gen[N][H], *h_reco[N][H];
    TH2 *m_A[N][H];

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        TH1 *h_gen_, *h_reco_;
        TH2 *m_A_;

        for (unsigned i = M4; i < N; i++)    // channel loop
        {
            float LUMI;
            if      ((i == M4) || (i == ME))
                LUMI = MUON_TRIG_LUMI;
            else if ((i == E4) || (i == EM))
                LUMI = ELEC_TRIG_LUMI * ELEC_TRIG_SF;
            float sf = LUMI * 1000 * XSEC[ZZ] / NGEN[ZZ];

            // Get histograms
            inFile->GetObject(selection[i] + "/" + hname[h] + "_gen", h_gen_);
            h_gen_->SetDirectory(0);
            h_gen_->SetStats(0);
            h_gen_->Scale(sf);
            h_gen[i][h] = h_gen_;

            inFile->GetObject(selection[i] + "/" + hname[h] + "_reco", h_reco_);
            h_reco_->SetDirectory(0);
            h_reco_->SetStats(0);
            h_reco_->Scale(sf);
            h_reco[i][h] = h_reco_;

            inFile->GetObject(selection[i] + "/" + hname[h] + "_2d", m_A_);
            m_A_->SetDirectory(0);
            m_A_->SetStats(0);
            m_A_->Scale(sf);
            m_A[i][h] = m_A_;
        }
    }
    cout << "Got migration histograms and matrices" << endl << endl;




    //
    //  DATA
    //


    TString prefix  = "4l_" + YEAR_STR;

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
        h_data[L4][h]   = (TH1*) h_data[1][h]->Clone();
        h_gen[L4][h]    = (TH1*) h_gen[1][h]->Clone();
        h_reco[L4][h]   = (TH1*) h_reco[1][h]->Clone();
        m_A[L4][h]      = (TH2*) m_A[1][h]->Clone();

        for (unsigned i = 2; i < N; i++)
        {
            h_data[L4][h]->Add(h_data[i][h]);
            h_gen[L4][h]->Add(h_gen[i][h]);
            h_reco[L4][h]->Add(h_reco[i][h]);
            m_A[L4][h]->Add(m_A[i][h]);
        }
    }





    ////
    ////
    ////    UNFOLD
    ////
    ////


    //
    //  RESPONSE OBJECTS
    //

    TH1 *h_result[N][H];
    TH2D *m_cov[N][H], *m_unf[N][H];

    for (unsigned h = 0; h < H; h++)    // distribution loop
    {
        cout << "Unfolding " << hname[h] << " distribution..." << endl;

        for (unsigned i = L4; i < 1; i++)
        {
            // Don't use gen since it doesn't include inefficiences
//          RooUnfoldResponse response(h_reco[i][h], 0, m_A[i][h]);
            RooUnfoldResponse response;
            response.UseOverflow();
            response.Setup(0, 0, m_A[i][h]);

            RooUnfoldBayes unfolder(&response, h_data[i][h], 4);

            h_result[i][h] = unfolder.Hreco();
            h_result[i][h]->SetStats(0);
            h_result[i][h]->GetXaxis()->SetTitle(h_gen[i][h]->GetXaxis()->GetTitle());

            TMatrixD covariance = unfolder.Ereco();
            m_cov[i][h] = new TH2D(covariance);

            TMatrixD unfolding = unfolder.UnfoldingMatrix();
            m_unf[i][h] = new TH2D(unfolding);

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
        for (unsigned i = 0; i < 1; i++)
        {
            c[i][h] = new TCanvas("c_" + hname[h] + "_" + selection[i], "", 800, 600);

            h_data[i][h]->SetLineColor(kBlack);
            h_data[i][h]->SetLineWidth(2);
            h_data[i][h]->SetMarkerStyle(20);
            h_data[i][h]->SetMarkerSize(2);

            h_reco[i][h]->SetLineColor(8);
            h_reco[i][h]->SetLineWidth(2);

            h_result[i][h]->SetLineColor(kRed);
            h_result[i][h]->SetLineWidth(2);
            h_result[i][h]->SetMarkerColor(kRed);
            h_result[i][h]->SetMarkerStyle(22);
            h_result[i][h]->SetMarkerSize(2);

            h_gen[i][h]->SetLineColor(kBlue);
            h_gen[i][h]->SetLineWidth(2);

            TLegend *l = new TLegend(0.8, 0.8, 1, 1);
            l->AddEntry(h_gen[i][h], "Gen", "L");
            l->AddEntry(h_reco[i][h], "Reco", "L");
            l->AddEntry(h_data[i][h], "Data", "LP");
            l->AddEntry(h_result[i][h], "Unfolded data", "LP");

            c[i][h]->cd();
            h_result[i][h]->Draw();
            h_gen[i][h]->Draw("SAME");
            h_reco[i][h]->Draw("SAME");
            h_data[i][h]->Draw("SAME");
            h_result[i][h]->Draw("SAME");
            l->Draw();
        }
    }



    //
    //  WRITE
    //

    TString outName = "unfolding_roo_" + YEAR_STR + ".root";
    TFile *outFile = new TFile(outName, "RECREATE");

    for (unsigned i = 0; i < 1; i++)
    {
        outFile->mkdir(selection[i]);

        for (unsigned h = 0; h < H; h++)
        {
            outFile->mkdir(selection[i] + "/" + hname[h]);
            outFile->cd(selection[i] + "/" + hname[h]);

            c[i][h]->Write();

            m_A[i][h]->SetName("m_A_" + hname[h] + "_" + selection[i]); 
            m_A[i][h]->Write();
            m_cov[i][h]->SetName("m_cov_" + hname[h] + "_" + selection[i]); 
            m_cov[i][h]->Write();
            m_unf[i][h]->SetName("m_unf_" + hname[h] + "_" + selection[i]); 
            m_unf[i][h]->Write();

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
