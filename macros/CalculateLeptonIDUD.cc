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

// Cuts
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2012.hh"

using namespace std;



void CalculateLeptonIDUD(const TString flavor, const TString type, const TString suff = "")
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

    const unsigned  N4L = 4,            N2L = 3;    // Number of final state trees
    TString sel_4l[N4L+1]   = { "4l",   "4m",   "2m2e", "2e2m", "4e"};
    TString sel_2l[N2L]     = { "ll",   "mumu", "ee"};

    float   sel_4l_gen[N4L],    sel_4l_nom[N4L],    sel_4l_up[N4L],     sel_4l_dn[N4L];
    float   sel_2l_gen[N2L],    sel_2l_nom[N2L],    sel_2l_up[N2L],     sel_2l_dn[N2L];

    for (unsigned i = 0; i < N4L; i++)
    {
        sel_4l_nom[i] = 0;  sel_4l_gen[i] = 0;  sel_4l_up[i] = 0;   sel_4l_dn[i] = 0;
    }

    for (unsigned i = 0; i < N2L; i++)
    {
        sel_2l_nom[i] = 0;  sel_2l_gen[i] = 0;  sel_2l_up[i] = 0;   sel_2l_dn[i] = 0;
    }



    //
    //  HISTOGRAMS
    //

    TString idName;
    if (suff.IsNull())
        idName = flavor + "_" + type + "_smear_" + YEAR_STR + ".root";
    else
        idName = flavor + "_" + type + "_smear_" + YEAR_STR + "_" + suff + ".root";

    TH2     *h_nom, *h_err;

    TString inPath  = BLT_PATH + "/BLTAnalysis/data/" + idName;
    TFile   *idFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    idFile->GetObject("FINAL", h_nom);
    idFile->GetObject("ERROR", h_err);
    h_nom->SetDirectory(0);
    h_err->SetDirectory(0);

    const float PT_MIN = h_nom->GetYaxis()->GetXmin();
    const float PT_MAX = h_nom->GetYaxis()->GetXmax();

    cout << "Limits: " << PT_MIN << ", " << PT_MAX << endl << endl;

    idFile->Close();
    cout << "Closed " << inPath << endl;



    //
    //  FOUR-LEPTON
    //

    inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/selected_zz_4l.root";
    TFile   *zzFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    unsigned i = 0;

    for (unsigned i_ = 0; i_ < N4L+1; i_++)
    {
        TTreeReader reader(sel_4l[i_] + "_zz_4l", zzFile);

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

        cout << "Running over " << nEvents << " total " << sel_4l[i_] << " events" << endl;
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
            for (unsigned l = 0; l < p4.size(); l++)
            {
                float pt = p4[l].Pt();

                // all but low-et reco histograms are inclusive of higher pt
                if (pt > PT_MAX && !suff.EqualTo("lowEt"))
                    pt = 0.99 * PT_MAX;

                // in any case, don't include leptons below pt range
                else if (pt < PT_MIN)
                    pt = 1.01 * PT_MIN;

                bins.push_back(h_nom->FindBin(p4[l].Eta(), pt));
            }


            // Calculate event weight variations

            float genWeight = *genWeight_;
            float nomWeight = genWeight,    upWeight = genWeight,     dnWeight = genWeight;

            for (unsigned l = 0; l < bins.size(); l++)
            {
                float sf = h_nom->GetBinContent(bins[l]);
                float err = h_err->GetBinContent(bins[l]);

                nomWeight *= sf;
                upWeight *= sf + err;
                dnWeight *= sf - err;
            }

            sel_4l_gen[i] += genWeight;     sel_4l_nom[i] += nomWeight;
            sel_4l_up[i] += upWeight;       sel_4l_dn[i] += dnWeight;
        }

        // This adds "2m2e" and "2e2m"
        if (!sel_4l[i_].EqualTo("2m2e"))
            i++;
    }
    zzFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  DILEPTON
    //

    // Open one file

    inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/selected_zjets_m-50.root";
    TFile   *dyFile = TFile::Open(inPath);

    const float MOD = 100;

    cout << "Opened " << inPath << endl;

    for (unsigned i = 1; i < N2L; i++)
    {
        TTreeReader reader(sel_2l[i] + "_zjets_m-50", dyFile);

        TTreeReaderValue    <Float_t>               genWeight_      (reader,    "genWeight");
        TTreeReaderValue    <Float_t>               qtWeight_       (reader,    "qtWeight");
        TTreeReaderValue    <TLorentzVector>        l1p4_           (reader,    "l1p4");
        TTreeReaderValue    <Short_t>               l1pdg_          (reader,    "l1pdg");
        TTreeReaderValue    <TLorentzVector>        l2p4_           (reader,    "l2p4");

        cout << "Loaded branches" << endl;


        // Event loop

        int nEvents = reader.GetEntries(kTRUE) / MOD;

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


            // Nominal weight

            float genWeight = (*genWeight_) * (*qtWeight_);
            float nomWeight = genWeight,    upWeight = genWeight,     dnWeight = genWeight;

            if (abs(*l1pdg_) == PDG)
            {
                // Get lepton momenta for the flavor we care about
                vector<TLorentzVector> p4 = {*l1p4_, *l2p4_};

                vector<int> bins;
                for (unsigned l = 0; l < p4.size(); l++)
                {
                    float pt = p4[l].Pt();

                    // all but low-et reco histograms are inclusive of higher pt
                    if (pt > PT_MAX && !suff.EqualTo("lowEt"))
                        pt = 0.99 * PT_MAX;

                    // in any case, don't include leptons below pt range
                    else if (pt < PT_MIN)
                        continue;

                    bins.push_back(h_nom->FindBin(p4[l].Eta(), pt));
                }


                // Calculate event weight variations

                for (unsigned l = 0; l < bins.size(); l++)
                {
                    float sf = h_nom->GetBinContent(bins[l]);
                    float err = h_err->GetBinContent(bins[l]);

                    nomWeight *= sf;
                    upWeight *= sf + err;
                    dnWeight *= sf - err;
                }
            }

            sel_2l_gen[i] += genWeight;     sel_2l_nom[i] += nomWeight;
            sel_2l_up[i] += upWeight;       sel_2l_dn[i] += dnWeight;

            sel_2l_gen[0] += genWeight;     sel_2l_nom[0] += nomWeight;
            sel_2l_up[0] += upWeight;       sel_2l_dn[0] += dnWeight;
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
    for (unsigned i = 0; i < N4L; i++)
        cout << "\t" << sel_4l_gen[i] / sel_4l_nom[i];
    cout << endl << "\t" << "];" << endl << endl;

    cout << "eff_4l_up = [" << endl;
    for (unsigned i = 0; i < N4L; i++)
        cout << "\t" << sel_4l_gen[i] / sel_4l_up[i];
    cout << endl << "\t" << "];" << endl << endl;

    cout << "eff_4l_dn = [" << endl;
    for (unsigned i = 0; i < N4L; i++)
        cout << "\t" << sel_4l_gen[i] / sel_4l_dn[i];
    cout << endl << "\t" << "];" << endl << endl;


    cout << "%\t" << "Efficiency" << "\t\t\t" << "ll" << "\t\t\t";
    cout << "mumu" << "\t\t\t" << "ee" << endl << endl;

    cout << "eff_2l_nom = [" << endl;
    for (unsigned i = 0; i < N2L; i++)
        cout << "\t" << sel_2l_gen[i] / sel_2l_nom[i];
    cout << endl << "\t" << "];" << endl << endl;

    cout << "eff_2l_up = [" << endl;
    for (unsigned i = 0; i < N2L; i++)
        cout << "\t" << sel_2l_gen[i] / sel_2l_up[i];
    cout << endl << "\t" << "];" << endl << endl;

    cout << "eff_2l_dn = [" << endl;
    for (unsigned i = 0; i < N2L; i++)
        cout << "\t" << sel_2l_gen[i] / sel_2l_dn[i];
    cout << endl << "\t" << "];" << endl << endl;
}
