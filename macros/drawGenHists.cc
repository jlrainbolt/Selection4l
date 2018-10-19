#include <tuple>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"


using namespace std;



void drawGenHists(const TString suffix)
{
    // Selection
    const unsigned N = 4;
    unsigned                M4 = 0, ME = 1, EM = 2, E4 = 3;     // Indices
    TString selection[N] = {"4m",   "2m2e", "2e2m", "4e"};
    bool muPair1[N]      = {kTRUE,  kTRUE,  kFALSE, kFALSE};


    // Level
    // (phase space, fiducial region, gen-level matched selected events)
    const unsigned M = 4;
    unsigned                PS = 0,         FR = 1,         GS = 2,         RS = 3;
    TString inPrefix[M]  = {"gen",          "gen",          "match",        "match"};
    TString outPrefix[M] = {"PhaseSpace",   "Fiducial",     "GenSelected",  "RecoSelected"};
    int color[M]         = {46,             30,             9,              lLightBlue};




    //--- SCALING INFO ---//

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections

    // Only valid for ZZTo4L!!!
    Float_t xsec = 1.212;



    //--- HISTOGRAMS ---//

    vector<tuple<TString, TString, Double_t, Double_t>> val2l, val4l;
    Int_t bins = 100;
    //                           name           member              xmin        xmax
    // Four lepton
    val4l.push_back(make_tuple( "zzm",          "zzp4.M()",         80,         100));
    val4l.push_back(make_tuple( "zzpt",         "zzp4.Pt()",        0,          100));

    val4l.push_back(make_tuple( "z1m",          "z1p4.M()",         0,          100));
    val4l.push_back(make_tuple( "z1pt",         "z1p4.Pt()",        0,          100));
    val4l.push_back(make_tuple( "z1pdg",        "z1pdg",            8.5,        15.5));

    val4l.push_back(make_tuple( "z2m",          "z2p4.M()",         0,          50));
    val4l.push_back(make_tuple( "z2pt",         "z2p4.Pt()",        0,          60));
    val4l.push_back(make_tuple( "z2pdg",        "z2pdg",            8.5,        15.5));

    val4l.push_back(make_tuple( "l1pt",         "l1p4.Pt()",        0,          100));
    val4l.push_back(make_tuple( "l1eta",        "l1p4.Eta()",       -6,         6));
    val4l.push_back(make_tuple( "l1pdg",        "l1pdg",            -15.5,      15.5));

    val4l.push_back(make_tuple( "l2pt",         "l2p4.Pt()",        0,          50));
    val4l.push_back(make_tuple( "l2eta",        "l2p4.Eta()",       -6,         6));
    val4l.push_back(make_tuple( "l2pdg",        "l2pdg",            -15.5,      15.5));

    val4l.push_back(make_tuple( "l3pt",         "l3p4.Pt()",        0,          30));
    val4l.push_back(make_tuple( "l3eta",        "l3p4.Eta()",       -6,         6));
    val4l.push_back(make_tuple( "l3pdg",        "l3pdg",            -15.5,      15.5));

    val4l.push_back(make_tuple( "l4pt",         "l4p4.Pt()",        0,          20));
    val4l.push_back(make_tuple( "l4eta",        "l4p4.Eta()",       -6,         6));
    val4l.push_back(make_tuple( "l4pdg",        "l4pdg",            -15.5,      15.5));


    // Differential distributions
    val4l.push_back(make_tuple( "b_z1p",        "b_z1p4.P()",       0,          50));
    val4l.push_back(make_tuple( "b_z2p",        "b_z2p4.P()",       0,          50));
    val4l.push_back(make_tuple( "b_ttm",        "b_ttp4.M()",       0,          70));

    val4l.push_back(make_tuple( "b_l1p",        "b_l1p4.P()",       20,         50));
    val4l.push_back(make_tuple( "b_l2p",        "b_l2p4.P()",       10,         50));
    val4l.push_back(make_tuple( "b_l3p",        "b_l3p4.P()",       0,          30));
    val4l.push_back(make_tuple( "b_l4p",        "b_l4p4.P()",       0,          20));

    val4l.push_back(make_tuple( "b_theta",      "b_theta",          0,          3.142));
    val4l.push_back(make_tuple( "b_phi",        "b_phi",            0,          3.142));
    val4l.push_back(make_tuple( "b_z1alpha",    "b_z1alpha",        0,          3.142));
    val4l.push_back(make_tuple( "b_z2alpha",    "b_z2alpha",        0,          3.142));
    val4l.push_back(make_tuple( "bb_z1theta",   "bb_z1theta",       0,          3.142));
    val4l.push_back(make_tuple( "bb_z2theta",   "bb_z2theta",       0,          3.142));



    // Draw histograms and write to file
    for (unsigned h = 0; h < M; h++)
    {
        TFile *inFile = TFile::Open("trees_" + inPrefix[h] + "_" + suffix + ".root");
        TFile *outFile = new TFile("hists_" + outPrefix[h] + "_" + suffix + ".root", "RECREATE");
  
        TH1D *hTotalEvents, *hPhaseSpaceEvents, *hFiducialEvents;   //, *hSelectedEvents;
        inFile->GetObject("TotalEvents_" + suffix,      hTotalEvents);
        inFile->GetObject("PhaseSpaceEvents_" + suffix, hPhaseSpaceEvents);
        inFile->GetObject("FiducialEvents_" + suffix,   hFiducialEvents);
//      inFile->GetObject("SelectedEvents_" + suffix,   hSelectedEvents);

        Float_t nGen    = hTotalEvents->GetBinContent(1) - 2 * hTotalEvents->GetBinContent(10);


        for (unsigned i = 0; i < N; i++)
        {
            cout << "Drawing " << selection[i] << " histograms..." << endl;



            //--- DIRECTORY, TREE ---//

            outFile->mkdir(selection[i]);
            outFile->cd(selection[i]);

            TTree *tree;
            inFile->GetObject(selection[i] + "_" + suffix, tree);



            //--- SCALING ---//

            Float_t lumi    = muPair1[i] ? 36.735 : 41.529;
            Float_t scale   = lumi * 1000. * xsec / nGen;



            //--- DRAW HISTS ---//

            vector<tuple<TString, TString, Double_t, Double_t>> val = val4l;

            for (unsigned j = 0; j < val.size(); j++)
            {
                TString hname = get<0>(val[j]), quantity = get<1>(val[j]);
                Float_t xmin = get<2>(val[j]),  xmax = get<3>(val[j]);
                TString weight = "weight";

                if (h == GS)
                    quantity = "gen_" + quantity;
                if (h == FR)
                    weight = weight + "*isFid";

                TString varexp;
                varexp.Form(quantity + ">>" + hname + "(%i,%g,%g)", bins, xmin, xmax);
                tree->Draw(varexp, weight);


                TH1F* hist;
                gDirectory->GetObject(hname, hist);

                if (hist)
                    hist->SetName(hname);
                else
                    hist = new TH1F(hname, quantity+" {"+weight+"}", bins, xmin, xmax);


                hist->Scale(scale);
                hist->SetLineColor(color[h]);
                hist->SetFillColor(color[h]);
                hist->Sumw2(kTRUE);

                hist->Write();
            }
        }

        outFile->cd();
/*
        hTotalEvents->Write();
        hPhaseSpaceEvents->Write();
        hFiducialEvents->Write();
        hSelectedEvents->Write();
*/
        outFile->Close();
        inFile->Close();

        cout << "Created file hists_" + outPrefix[h] + "_" + suffix + ".root" << endl;
        cout << endl;
    }
}
