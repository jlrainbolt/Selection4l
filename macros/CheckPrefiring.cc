// STL
#include <vector>
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"

// Cuts
//#include "Cuts2017.hh"
#include "Cuts2016.hh"

using namespace std;



void CheckPrefiring(const TString suffix, const TString year)
{
    if (!year.EqualTo(YEAR_STR))
    {
        cout << "Wrong year in header file!" << endl;
        return;
    }


    //
    //  OUTPUT FILE
    //

    TString prefix  = "prefiring";
    TString outName = prefix + "_" + year + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    const unsigned N = 7;
    float minEta[N] = {1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

    TH1D *h_uncorr[N], *h_corr[N];
    for (unsigned i = 0; i < N; i++)
    {
        TString etaVal = TString::Format("%.1f", minEta[i]);

        h_uncorr[i] = new TH1D("eta_" + etaVal + "_unweighted", "#eta > " + etaVal,
                5, -0.5, 4.5);
        h_uncorr[i]->SetXTitle("Number of electrons");

        h_corr[i] = new TH1D("eta_" + etaVal + "_weighted", "#eta > " + etaVal,
                5, -0.5, 4.5);
        h_corr[i]->SetXTitle("Number of electrons");
    }



    //
    //  EVENT LOOP
    //

    TString inName  = "selected_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + year + "_new/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl << endl;

    TTreeReader reader("4e_zz_4l", inFile);

    TTreeReaderValue    <Float_t>               weight_         (reader,    "weight");
    TTreeReaderValue    <Float_t>               ecalWeight_     (reader,    "ecalWeight");
    TTreeReaderValue    <Bool_t>                hasTauDecay_    (reader,    "hasTauDecay");
    TTreeReaderValue    <TLorentzVector>        l1p4_           (reader,    "l1p4");
    TTreeReaderValue    <TLorentzVector>        l2p4_           (reader,    "l2p4");
    TTreeReaderValue    <TLorentzVector>        l3p4_           (reader,    "l3p4");
    TTreeReaderValue    <TLorentzVector>        l4p4_           (reader,    "l4p4");

    cout << "Loaded branches" << endl;


    // Event loop

    int nEvents = reader.GetEntries(kTRUE);

    cout << "Running over " << nEvents << " total 4e events" << endl;
    cout << endl;

    while (reader.Next())
    {
        if (*hasTauDecay_)
            continue;


        // Get electron momenta
        vector<TLorentzVector> p4;
        p4.push_back(*l1p4_);
        p4.push_back(*l2p4_);
        p4.push_back(*l3p4_);
        p4.push_back(*l4p4_);


        // Get weights
        float sf = INT_LUMI * 1000. * XSEC_ZZ_4L / NGEN_ZZ_4L;
        float weight = sf * (*weight_);
        float unweight = sf * (*weight_) / (*ecalWeight_);


        for (unsigned i = 0; i < N; i++)    // loop over minimum eta values
        {
            int nElecs = 0;     // number of electrons in prefiring region

            for (unsigned l = 0; l < p4.size(); l++)
            {
                if (p4[l].Eta() > minEta[i])
                    nElecs++;
            }

            h_uncorr[i]->Fill(nElecs, unweight);
            h_corr[i]->Fill(nElecs, weight);
        }
    }
    inFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  GRAPHS
    //

    float y_corr[N], y_uncorr[N];
    for (unsigned i = 0; i < N; i++)
    {
        y_corr[i] = h_corr[i]->Integral(3, 5) / h_corr[i]->Integral(0, N+1);
        y_uncorr[i] = h_uncorr[i]->Integral(3, 5) / h_uncorr[i]->Integral(0, N+1);
    }

    TGraph *g_corr = new TGraph(N, minEta, y_corr);
    TGraph *g_uncorr = new TGraph(N, minEta, y_uncorr);

    TCanvas *c_corr = new TCanvas("graph_weighted", "", 800, 600);
    c_corr->cd();
    g_corr->SetMarkerColor(kBlue);
    g_corr->Draw("AP*");
    g_corr->SetTitle("");
    g_corr->GetXaxis()->SetTitle("#eta_{min}");
    g_corr->GetYaxis()->SetTitle("Fraction of events with #geq 2 electrons with #eta > #eta_{min}");

    TCanvas *c_uncorr = new TCanvas("graph_unweighted", "", 800, 600);
    c_uncorr->cd();
    g_uncorr->SetMarkerColor(kBlue);
    g_uncorr->Draw("AP*");
    g_uncorr->SetTitle("");
    g_uncorr->GetXaxis()->SetTitle("#eta_{min}");
    g_uncorr->GetYaxis()->SetTitle("Fraction of events with #geq 2 electrons with #eta > #eta_{min}");



    //
    //  WRITE
    //

    outFile->cd();
    for (unsigned i = 0; i < N; i++)
    {
        h_corr[i]->Write();
        h_uncorr[i]->Write();
    }
    g_corr->Write();
    g_uncorr->Write();
    c_corr->Write();
    c_uncorr->Write();

    cout << "Wrote histograms to " << outName << endl;
}
