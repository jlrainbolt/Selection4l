// STL
#include <vector>
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TMath.h"

// Cuts
//#include "Cuts2017.hh"
#include "Cuts2016.hh"

using namespace std;



void CalculateLeptonReco(const TString flavor, const TString type, const TString suff = "")
{


    //
    //  INPUT
    //

    unsigned PDG;

    if      (flavor.EqualTo("muon"))
        PDG = 13;
    else if (flavor.EqualTo("electron"))
        PDG = 11;



    //
    //  CONTAINERS
    //

    const unsigned  N = 300;    // Number of Gaussian iterations

    const unsigned  M4L = 4,            M2L = 2;    // Number of final state trees
    TString sel_4l[M4L]   = { "4l",   "4m",   "2m2e", "4e"};
    TString sel_2l[M2L]     = { "mumu", "ee"};

    float   sel_4l_gen[M4L],        sel_2l_gen[M2L];
    float   sel_4l_nom[M4L],        sel_2l_nom[M2L];
    float   sel_4l_var[M4L][N],     sel_2l_var[M2L][N];

    for (unsigned i = 0; i < M4L; i++)
    {
        sel_4l_nom[i] = 0;
        sel_4l_gen[i] = 0;

        for (unsigned j = 0; j < N; j++)
            sel_4l_var[i][j] = 0;
    }

    for (unsigned i = 0; i < M2L; i++)
    {
        sel_2l_nom[i] = 0;
        sel_2l_gen[i] = 0;

        for (unsigned j = 0; j < N; j++)
            sel_2l_var[i][j] = 0;
    }



    //
    //  HISTOGRAMS
    //

    TString idName, idName_suff;
    idName = flavor + "_" + type + "_smear_" + YEAR_STR + ".root";
    idName_suff = flavor + "_" + type + "_smear_" + YEAR_STR + "_" + suff + ".root";

    TH2     *h_nom, *h_nom_suff, *h_smr[N], *h_smr_suff[N];

    TString inPath  = BLT_PATH + "/BLTAnalysis/data/" + idName;
    TFile   *idFile = TFile::Open(inPath);

    TString inPath_suff  = BLT_PATH + "/BLTAnalysis/data/" + idName_suff;
    TFile   *idFile_suff = TFile::Open(inPath_suff);

    cout << "Opened " << inPath << endl;

    idFile->GetObject("FINAL", h_nom);
    h_nom->SetDirectory(0);
    idFile_suff->GetObject("FINAL", h_nom_suff);
    h_nom_suff->SetDirectory(0);

    float PT_MIN = h_nom_suff->GetYaxis()->GetXmin();
    float PT_THRESH = h_nom->GetYaxis()->GetXmin();
    float PT_MAX = h_nom->GetYaxis()->GetXmax();

    cout << "Limits: " << PT_MIN << ", " << PT_THRESH << ", " << PT_MAX << endl << endl;

    for (unsigned j = 0; j < N; j++)
    {
        TString id = TString::Format("%i", j);
        idFile->GetObject("SMEAR" + id, h_smr[j]);
        h_smr[j]->SetDirectory(0);
        idFile_suff->GetObject("SMEAR" + id, h_smr_suff[j]);
        h_smr_suff[j]->SetDirectory(0);
    }

    idFile->Close();
    cout << "Closed " << inPath << endl;



    //
    //  FOUR-LEPTON
    //

    inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "_new/selected_zz_4l.root";
    TFile   *zzFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    unsigned i = 0;

    for (unsigned i = 0; i < M4L; i++)
    {
        TTreeReader reader(sel_4l[i] + "_zz_4l", zzFile);

        TTreeReaderValue    <Float_t>               genWeight_      (reader,    "genWeight");
        TTreeReaderValue    <TLorentzVector>        l1p4_           (reader,    "l1p4");
        TTreeReaderValue    <Short_t>               l1pdg_          (reader,    "l1pdg");
        TTreeReaderValue    <TLorentzVector>        l2p4_           (reader,    "l2p4");
        TTreeReaderValue    <Short_t>               l2pdg_          (reader,    "l2pdg");
        TTreeReaderValue    <TLorentzVector>        l3p4_           (reader,    "l3p4");
        TTreeReaderValue    <Short_t>               l3pdg_          (reader,    "l3pdg");
        TTreeReaderValue    <TLorentzVector>        l4p4_           (reader,    "l4p4");
        TTreeReaderValue    <Short_t>               l4pdg_          (reader,    "l4pdg");

        cout << "Loaded branches" << endl;


        // Event loop

        int nEvents = reader.GetEntries(kTRUE);

        cout << "Running over " << nEvents << " total " << sel_4l[i] << " events" << endl;
        cout << endl;

        while (reader.Next())
        {
            // Get lepton momenta for the flavor we care about

            vector<TLorentzVector> p4;

            if (abs(*l1pdg_) == PDG)
                p4.push_back(*l1p4_);
            if (abs(*l2pdg_) == PDG)
                p4.push_back(*l2p4_);
            if (abs(*l3pdg_) == PDG)
                p4.push_back(*l3p4_);
            if (abs(*l4pdg_) == PDG)
                p4.push_back(*l4p4_);


            // Get bin numbers for leptons we care about

            vector<int> bins;
            vector<bool> isSuff;
            for (unsigned l = 0; l < p4.size(); l++)
            {
                float pt = p4[l].Pt();

                // all but low-et reco histograms are inclusive of higher pt
                if (pt > PT_MAX)
                    pt = 0.99 * PT_MAX;
                if (pt < PT_MIN)
                    pt = 1.01 * PT_MIN;

                // in any case, don't include leptons below pt range
                if (pt > PT_THRESH)
                {
                    isSuff.push_back(kFALSE);
                    bins.push_back(h_nom->FindBin(p4[l].Eta(), pt));
                }
                else
                {
                    isSuff.push_back(kTRUE);
                    bins.push_back(h_nom_suff->FindBin(p4[l].Eta(), pt));
                }
            }


            // Nominal weight

            float genWeight = (*genWeight_);
            sel_4l_gen[i] += genWeight;

            float nomWeight = genWeight;

            for (unsigned l = 0; l < bins.size(); l++)
            {
                if (isSuff[l])
                    nomWeight *= h_nom_suff->GetBinContent(bins[l]);
                else
                    nomWeight *= h_nom->GetBinContent(bins[l]);
            }

            sel_4l_nom[i] += nomWeight;


            // Variations

            for (unsigned j = 0; j < N; j++)
            {
                float varWeight = genWeight;

                for (unsigned l = 0; l < bins.size(); l++)
                {
                    if (isSuff[l])
                        varWeight *= h_smr_suff[j]->GetBinContent(bins[l]);
                    else
                        varWeight *= h_smr[j]->GetBinContent(bins[l]);
                }

                sel_4l_var[i][j] += varWeight;
            }
        }
    }
    zzFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  DILEPTON
    //

    // Open one file

    inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "_new/selected_zjets_m-50.root";
    TFile   *dyFile = TFile::Open(inPath);

    float mod = 100;

    cout << "Opened " << inPath << endl;

    for (unsigned i = 0; i < M2L; i++)
    {
        TTreeReader reader(sel_2l[i] + "_zjets_m-50", dyFile);

        TTreeReaderValue    <Float_t>               genWeight_      (reader,    "genWeight");
        TTreeReaderValue    <Float_t>               qtWeight_       (reader,    "qtWeight");
        TTreeReaderValue    <TLorentzVector>        l1p4_           (reader,    "l1p4");
        TTreeReaderValue    <Short_t>               l1pdg_          (reader,    "l1pdg");
        TTreeReaderValue    <TLorentzVector>        l2p4_           (reader,    "l2p4");

        cout << "Loaded branches" << endl;


        // Event loop

        int nEvents = reader.GetEntries(kTRUE) / mod;

        cout << endl;
        cout << "Running over " << nEvents << " total " << sel_2l[i] << " events" << endl;
        cout << endl;

        while (reader.Next() && (reader.GetCurrentEntry() < nEvents))
        {
            if ((reader.GetCurrentEntry() % 100000 == 0) && (reader.GetCurrentEntry() > 0))
            {
                cout << "Processed " << reader.GetCurrentEntry() << " of " << nEvents << " events";
                cout << endl;
            }


            // Gen weight

            float genWeight = (*genWeight_) * (*qtWeight_);
            sel_2l_gen[i] += genWeight;

            float nomWeight = genWeight;


            if (abs(*l1pdg_) == PDG)
            {
                // Get lepton momenta for the flavor we care about
                vector<TLorentzVector> p4 = {*l1p4_, *l2p4_};

                vector<int> bins;
                vector<bool> isSuff;
                for (unsigned l = 0; l < p4.size(); l++)
                {
                    float pt = p4[l].Pt();

                    // all but low-et reco histograms are inclusive of higher pt
                    if (pt > PT_MAX)
                        pt = 0.99 * PT_MAX;
                    if (pt < PT_MIN)
                        pt = 1.01 * PT_MIN;

                    // in any case, don't include leptons below pt range
                    if (pt > PT_THRESH)
                    {
                        isSuff.push_back(kFALSE);
                        bins.push_back(h_nom->FindBin(p4[l].Eta(), pt));
                    }
                    else
                    {
                        isSuff.push_back(kTRUE);
                        bins.push_back(h_nom_suff->FindBin(p4[l].Eta(), pt));
                    }
                }


                // Nominal
                for (unsigned l = 0; l < bins.size(); l++)
                {
                    if (isSuff[l])
                        nomWeight *= h_nom_suff->GetBinContent(bins[l]);
                    else
                        nomWeight *= h_nom->GetBinContent(bins[l]);
                }


                // Variations

                for (unsigned j = 0; j < N; j++)
                {
                    float varWeight = genWeight;

                    for (unsigned l = 0; l < bins.size(); l++)
                    {
                        if (isSuff[l])
                            varWeight *= h_smr_suff[j]->GetBinContent(bins[l]);
                        else
                            varWeight *= h_smr[j]->GetBinContent(bins[l]);
                    }

                    sel_2l_var[i][j] += varWeight;
                }
            }
            else
            {
                for (unsigned j = 0; j < N; j++)
                    sel_2l_var[i][j] += nomWeight;
            }

            sel_2l_nom[i] += nomWeight;
        }
    }
    dyFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  PRINTOUT
    //

    cout << setprecision(16) << endl;

    cout << "%\t" << "Efficiency" << "\t\t\t" << "4l" << "\t\t\t" << "4m";
    cout << "\t\t\t" << "2m2e" << "\t\t\t" << "4e" << endl << endl;

    cout << "eff_4l_nom = [" << endl;
    for (unsigned i = 0; i < M4L; i++)
        cout << "\t" << sel_4l_gen[i] / sel_4l_nom[i];
    cout << endl << "\t" << "];" << endl << endl;

    cout << "eff_4l_var = [" << endl;
    for (unsigned j = 0; j < N; j++)
    {
        for (unsigned i = 0; i < M4L; i++)
            cout << "\t" << sel_4l_gen[i] / sel_4l_var[i][j];
        cout << endl;
    }
    cout << "\t" << "];" << endl << endl << endl;


    cout << "%\t" << "Efficiency" << "\t\t\t" << "mumu" << "\t\t\t" << "ee" << endl << endl;

    cout << "eff_2l_nom = [" << endl;
    cout << "\t" << (sel_2l_gen[0] + sel_2l_gen[1]) / (sel_2l_nom[0] + sel_2l_nom[1]);
    for (unsigned i = 0; i < M2L; i++)
        cout << "\t" << sel_2l_gen[i] / sel_2l_nom[i];
    cout << endl << "\t" << "];" << endl << endl;

    cout << "eff_2l_var = [" << endl;
    for (unsigned j = 0; j < N; j++)
    {
        cout << "\t" << (sel_2l_gen[0] + sel_2l_gen[1]) / (sel_2l_var[0][j] + sel_2l_var[1][j]);

        for (unsigned i = 0; i < M2L; i++)
            cout << "\t" << sel_2l_gen[i] / sel_2l_var[i][j];
        cout << endl;
    }
    cout << "\t" << "];" << endl;
}
