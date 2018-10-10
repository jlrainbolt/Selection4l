#include <tuple>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"


using namespace std;



void drawHists(const TString suffix)
{
    // Selection
    const unsigned N = 6;
    unsigned                MM = 0, EE = 1, M4 = 2, ME = 3, EM = 4, E4 = 5; // Indices
    TString selection[N] = {"mumu", "ee",   "4m",   "2m2e", "2e2m", "4e"};
    bool muPair1[N]      = {kTRUE,  kFALSE, kTRUE,  kTRUE,  kFALSE, kFALSE};



    //--- SCALING INFO ---//

    Bool_t isData;
    Int_t color;


    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections
    Float_t xsec = 0;

    if (suffix.Contains("2017"))
    {
        isData = kTRUE;
        color = kBlack;
    }
    else
    {
        isData = kFALSE;
        if (suffix == "zjets_m-50")
        {
            xsec = 5765.4;
            color = lYellow;
        }
        else if (suffix == "ggH_zz_4l")
        {
            xsec = 0.01212;
            color = lPurple;
        }
        else if (suffix == "vbfH_zz_4l")
        {
            xsec = 0.001034;
            color = lPurple;
        }
        else if (suffix == "ttbar")
        {
            xsec = 831.76;
            color = lGreen;
        }
        else if (suffix == "ww_2l2nu")
        {
            xsec = 12.178;
            color = lOrange;
        }
        else if (suffix == "wz_2l2q")
        {
            xsec = 5.595;
            color = lOrange;
        }
        else if (suffix == "wz_3lnu")
        {
            xsec = 4.42965;
            color = lOrange;
        }
        else if (suffix == "zz_2l2q")
        {
            xsec = 3.22;
            color = lOrange;
        }
        else if (suffix == "zz_4l")
        {
            xsec = 1.212;
            color = lLightBlue;
        }
    }



    //--- HISTOGRAMS ---//

    vector<tuple<TString, TString, Int_t, Double_t, Double_t>> val2l, val4l;
    //                           name       member          bins    xmin    xmax
    // Dilepton
    val2l.push_back(make_tuple("nPV",     "nPV",          51,     -0.5,   50.5));
    val2l.push_back(make_tuple("met",     "met",          100,    0,      100));

    val2l.push_back(make_tuple("z1m",     "z1p4.M()",     40,    80,     100));
    val2l.push_back(make_tuple("z1pt",    "z1p4.Pt()",    40,    0,      80));

    val2l.push_back(make_tuple("l1pt",    "l1p4.Pt()",    40,    20,     100));
    val2l.push_back(make_tuple("l1eta",   "l1p4.Eta()",   40,    -2.5,   2.5));
    val2l.push_back(make_tuple("l1iso",   "l1iso",        100,    0,      0.35));
    val2l.push_back(make_tuple("l1pdg",   "l1pdg",        31,     -15.5,  15.5));

    val2l.push_back(make_tuple("l2pt",    "l2p4.Pt()",    40,     0,      80));
    val2l.push_back(make_tuple("l2eta",   "l2p4.Eta()",   40,    -2.5,   2.5));
    val2l.push_back(make_tuple("l2iso",   "l2iso",        100,    0,      0.35));
    val2l.push_back(make_tuple("l2pdg",   "l2pdg",        31,     -15.5,  15.5));


    // Four lepton

    val4l.push_back(make_tuple("nPV",     "nPV",          25,     0,      50));
    val4l.push_back(make_tuple("met",     "met",          20,     0,      100));

    val4l.push_back(make_tuple("zzm",     "zzp4.M()",     20,     80,     100));
    val4l.push_back(make_tuple("zzpt",    "zzp4.Pt()",    20,     0,      80));

    val4l.push_back(make_tuple("z1m",     "z1p4.M()",     20,     10,      90));
    val4l.push_back(make_tuple("z1pt",    "z1p4.Pt()",    20,     0,      80));

    val4l.push_back(make_tuple("z2m",     "z2p4.M()",     20,     0,      40));
    val4l.push_back(make_tuple("z2pt",    "z2p4.Pt()",    20,     0,      80));

    val4l.push_back(make_tuple("l1pt",    "l1p4.Pt()",    20,     20,      100));
    val4l.push_back(make_tuple("l1eta",   "l1p4.Eta()",   25,     -2.5,   2.5));
    val4l.push_back(make_tuple("l1iso",   "l1iso",        35,     0.001,      0.35));
    val4l.push_back(make_tuple("l1pdg",   "l1pdg",        31,     -15.5,  15.5));

    val4l.push_back(make_tuple("l2pt",    "l2p4.Pt()",    20,     10,      90));
    val4l.push_back(make_tuple("l2eta",   "l2p4.Eta()",   25,     -2.5,   2.5));
    val4l.push_back(make_tuple("l2iso",   "l2iso",        35,     0.001,      0.35));
    val4l.push_back(make_tuple("l2pdg",   "l2pdg",        31,     -15.5,  15.5));

    val4l.push_back(make_tuple("l3pt",    "l3p4.Pt()",    20,     5,      45));
    val4l.push_back(make_tuple("l3eta",   "l3p4.Eta()",   25,     -2.5,   2.5));
    val4l.push_back(make_tuple("l3iso",   "l3iso",        35,     0.001,      0.35));
    val4l.push_back(make_tuple("l3pdg",   "l3pdg",        31,     -15.5,  15.5));

    val4l.push_back(make_tuple("l4pt",    "l4p4.Pt()",    20,     5,      25));
    val4l.push_back(make_tuple("l4eta",   "l4p4.Eta()",   25,     -2.5,   2.5));
    val4l.push_back(make_tuple("l4iso",   "l4iso",        35,     0.001,      0.35));
    val4l.push_back(make_tuple("l4pdg",   "l4pdg",        31,     -15.5,  15.5));


    // Differential distributions
    val4l.push_back(make_tuple("b_z1p", "b_z1p4.P()",     15,     0,      45));
    val4l.push_back(make_tuple("b_z2p", "b_z2p4.P()",     15,     0,      45));
    val4l.push_back(make_tuple("b_ttm", "b_ttp4.M()",     15,     5,      65));

    val4l.push_back(make_tuple("b_l1p", "b_l1p4.P()",    15,     20,     65));
    val4l.push_back(make_tuple("b_l2p", "b_l2p4.P()",    15,     15,      45));
    val4l.push_back(make_tuple("b_l3p", "b_l3p4.P()",    15,     0,      30));
    val4l.push_back(make_tuple("b_l4p", "b_l4p4.P()",    15,     0,      30));

    val4l.push_back(make_tuple("b_theta", "b_theta",      16,      0,      3.2));
    val4l.push_back(make_tuple("b_phi", "b_phi",      16,      0,      3.2));
    val4l.push_back(make_tuple("b_sin_phi", "sin(b_phi)",      20,      0,      1));
    val4l.push_back(make_tuple("b_z1alpha", "b_z1alpha",    16,      0,      3.2));
    val4l.push_back(make_tuple("b_z2alpha", "b_z2alpha",    16,      0,      3.2));
    val4l.push_back(make_tuple("bb_z1theta", "bb_z1theta",    16,      0,      3.2));
    val4l.push_back(make_tuple("bb_z2theta", "bb_z2theta",    16,      0,      3.2));



    // Draw histograms and write to file
    TFile *inFile = TFile::Open("trees_" + suffix + ".root");
    TFile *outFile = new TFile("hists_" + suffix + ".root", "RECREATE");

    TH1D *hTotalEvents, *hSelectedEvents;
    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    inFile->GetObject("SelectedEvents_" + suffix, hSelectedEvents);

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

        vector<tuple<TString, TString, Int_t, Double_t, Double_t>> val = i < M4 ? val2l : val4l;

        for (unsigned j = 0; j < val.size(); j++)
        {
            TString varexp, hname = get<0>(val[j]);
            varexp.Form(get<1>(val[j]) + ">>" + hname + "(%i,%g,%g)",
                        get<2>(val[j]), get<3>(val[j]), get<4>(val[j]));
            tree->Draw(varexp, "weight");


            TH1F* hist;
            gDirectory->GetObject(hname, hist);

            if (hist)
                hist->SetName(hname);
            else
                hist = new TH1F(hname, get<1>(val[j])+" {weight}",
                                get<2>(val[j]), get<3>(val[j]), get<4>(val[j]));


            hist->SetLineColor(color);
            if (isData)
            {
                hist->SetMarkerColor(color);
                hist->SetMarkerStyle(kFullCircle);
                hist->SetMarkerSize(2);
                hist->SetLineWidth(2);
            }
            else
            {
                hist->Scale(scale);
                hist->SetFillColor(color);
            }

            hist->Write();
        }
    }

    outFile->cd();
    hTotalEvents->Write();
    hSelectedEvents->Write();

    outFile->Close();
    inFile->Close();
    
    cout << "Created file hists_" + suffix + ".root" << endl;
}
