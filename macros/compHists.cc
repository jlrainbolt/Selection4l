#include <tuple>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"


using namespace std;



void compHists(const TString suffix)
{
    // Selection
    // (only considering 4l cases...)
    const unsigned N = 4;
    unsigned                M4 = 0, ME = 1, EM = 2, E4 = 3;     // Indices
    TString selection[N] = {"4m",   "2m2e", "2e2m", "4e"};
    bool muPair1[N]      = {kTRUE,  kTRUE,  kFALSE, kFALSE};



    //--- SCALING INFO ---//

    Int_t color;

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections
    Float_t xsec = 0;
    if (suffix.Contains("zjets_m-50"))
    {
        xsec = 5765.4;
        color = lYellow;
    }
    else if (suffix.Contains("ggH_zz_4l"))
    {
        xsec = 0.01212;
        color = lPurple;
    }
    else if (suffix.Contains("vbfH_zz_4l"))
    {
        xsec = 0.001034;
        color = lPurple;
    }
    else if (suffix.Contains("ttbar"))
    {
        xsec = 831.76;
        color = lGreen;
    }
    else if (suffix.Contains("ww_2l2nu"))
    {
        xsec = 12.178;
        color = lOrange;
    }
    else if (suffix.Contains("wz_2l2q"))
    {
        xsec = 5.595;
        color = lOrange;
    }
    else if (suffix.Contains("wz_3lnu"))
    {
        xsec = 4.42965;
        color = lOrange;
    }
    else if (suffix.Contains("zz_2l2q"))
    {
        xsec = 3.22;
        color = lOrange;
    }
    else if (suffix.Contains("zz_4l"))
    {
        xsec = 1.212;
        color = lLightBlue;
    }

    color = 9;



    //--- HISTOGRAMS ---//

    // (Right now I am not using the bin/min/max info)
    vector<tuple<TString, TString, Int_t, Double_t, Double_t>> val2l, val4l;

    //                          name            member                              bins    xmin    xmax
    // Kinematic quantities
    val4l.push_back(make_tuple( "zzm",          "gen_zzp4.M()-zzp4.M()",        20,     80,     100));
    val4l.push_back(make_tuple( "zzpt",         "gen_zzp4.Pt()-zzp4.Pt()",      20,     0,      80));

    val4l.push_back(make_tuple( "z1m",          "gen_z1p4.M()-z1p4.M()",        20,     10,     90));
    val4l.push_back(make_tuple( "z1pt",         "gen_z1p4.Pt()-z1p4.Pt()",      20,     0,      80));

    val4l.push_back(make_tuple( "z2m",          "gen_z2p4.M()-z2p4.M()",        20,     0,      40));
    val4l.push_back(make_tuple( "z2pt",         "gen_z2p4.Pt()-z2p4.Pt()",      20,     0,      80));

    val4l.push_back(make_tuple( "l1pt",         "gen_l1p4.Pt()-l1p4.Pt()",      20,     20,     100));
    val4l.push_back(make_tuple( "l1eta",        "gen_l1p4.Eta()-l1p4.Eta()",    25,     -2.5,   2.5));
    val4l.push_back(make_tuple( "l1pdg",        "gen_l1pdg-l1pdg",              31,     -15.5,  15.5));

    val4l.push_back(make_tuple( "l2pt",         "gen_l2p4.Pt()-l2p4.Pt()",      20,     10,     90));
    val4l.push_back(make_tuple( "l2eta",        "gen_l2p4.Eta()-l2p4.Eta()",    25,     -2.5,   2.5));
    val4l.push_back(make_tuple( "l2pdg",        "gen_l2pdg-l2pdg",              31,     -15.5,  15.5));

    val4l.push_back(make_tuple( "l3pt",         "gen_l3p4.Pt()-l3p4.Pt()",      20,     5,      45));
    val4l.push_back(make_tuple( "l3eta",        "gen_l3p4.Eta()-l3p4.Eta()",    25,     -2.5,   2.5));
    val4l.push_back(make_tuple( "l3pdg",        "gen_l3pdg-l3pdg",              31,     -15.5,  15.5));

    val4l.push_back(make_tuple( "l4pt",         "gen_l4p4.Pt()-l4p4.Pt()",      20,     5,      25));
    val4l.push_back(make_tuple( "l4eta",        "gen_l4p4.Eta()-l4p4.Eta()",    25,     -2.5,   2.5));
    val4l.push_back(make_tuple( "l4pdg",        "gen_l4pdg-l4pdg",              31,     -15.5,  15.5));


    // Differential distributions
//  val4l.push_back(make_tuple( "b_z1p",        "gen_b_z1p4.P()-b_z1p4.P()",    15,     0,      45));
//  val4l.push_back(make_tuple( "b_z2p",        "gen_b_z2p4.P()-b_z2p4.P()",    15,     0,      45));
    val4l.push_back(make_tuple( "b_ttm",        "gen_b_ttp4.M()-b_ttp4.M()",    15,     5,      65));

    val4l.push_back(make_tuple( "b_l1p",        "gen_b_l1p4.P()-b_l1p4.P()",    15,     20,     65));
    val4l.push_back(make_tuple( "b_l1pdg",      "gen_b_l1pdg-b_l1pdg",          31,     -15.5,  15.5));

    val4l.push_back(make_tuple( "b_l2p",        "gen_b_l2p4.P()-b_l2p4.P()",    15,     15,     45));
    val4l.push_back(make_tuple( "b_l2pdg",      "gen_b_l2pdg-b_l2pdg",          31,     -15.5,  15.5));

    val4l.push_back(make_tuple( "b_l3p",        "gen_b_l3p4.P()-b_l3p4.P()",    15,     0,      30));
    val4l.push_back(make_tuple( "b_l3pdg",      "gen_b_l3pdg-b_l3pdg",          31,     -15.5,  15.5));

    val4l.push_back(make_tuple( "b_l4p",        "gen_b_l4p4.P()-b_l4p4.P()",    15,     0,      30));
    val4l.push_back(make_tuple( "b_l4pdg",      "gen_b_l4pdg-b_l4pdg",          31,     -15.5,  15.5));

    val4l.push_back(make_tuple( "b_theta",      "gen_b_theta-b_theta",          16,     0,      3.2));
    val4l.push_back(make_tuple( "b_phi",        "gen_b_phi-b_phi",              16,     0,      3.2));
    val4l.push_back(make_tuple( "b_sin_phi",    "sin(gen_b_phi)-sin(b_phi)",    20,     0,      1));

    val4l.push_back(make_tuple( "b_z1alpha",    "gen_b_z1alpha-b_z1alpha",      16,     0,      3.2));
    val4l.push_back(make_tuple( "b_z2alpha",    "gen_b_z2alpha-b_z2alpha",      16,     0,      3.2));

    val4l.push_back(make_tuple( "bb_z1theta",   "gen_bb_z1theta-bb_z1theta",    16,     0,      3.2));
    val4l.push_back(make_tuple( "bb_z2theta",   "gen_bb_z2theta-bb_z2theta",    16,     0,      3.2));



    //--- DRAW, WRITE ---//

    TFile *inFile = TFile::Open("trees_match_" + suffix + ".root");
    TFile *outFile = new TFile("hists_match_" + suffix + ".root", "RECREATE");


    // Copy over other histograms

    TH1D *hTotalEvents, *hPhaseSpaceEvents, *hFiducialEvents, *hSelectedEvents;

    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);

    inFile->GetObject("PhaseSpaceEvents_" + suffix, hPhaseSpaceEvents);
    hPhaseSpaceEvents->SetDirectory(outFile);

    inFile->GetObject("FiducialEvents_" + suffix, hFiducialEvents);
    hFiducialEvents->SetDirectory(outFile);

    inFile->GetObject("SelectedEvents_" + suffix, hSelectedEvents);
    hSelectedEvents->SetDirectory(outFile);


    // Get number of generator events for scaling
    Float_t nGen = hTotalEvents->GetBinContent(1) - 2 * hTotalEvents->GetBinContent(10);


    // Loop over selection directories
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
            varexp = get<1>(val[j]) + ">>" + hname;
//          varexp.Form(get<1>(val[j]) + ">>" + hname + "(%i,%g,%g)", get<2>(val[j]), get<3>(val[j]), get<4>(val[j]));
            tree->Draw(varexp, "weight");


            TH1F* hist;
            gDirectory->GetObject(hname, hist);

            if (hist)
                hist->SetName(hname);
            else    // this uses the binning but I guess it doesn't really matter
                hist = new TH1F(hname, get<1>(val[j])+" {weight}",
                                get<2>(val[j]), get<3>(val[j]), get<4>(val[j]));


            hist->Sumw2(kTRUE);
            hist->SetLineColor(color);
            hist->Scale(scale);
            hist->SetFillColor(color);

            hist->Write();
        }
    }

    
    outFile->cd();
    hTotalEvents->Write();
    hPhaseSpaceEvents->Write();
    hFiducialEvents->Write();
    hSelectedEvents->Write();



    //--- CLOSE ---//

    outFile->Close();
    inFile->Close();
    
    cout << "Created file hists_match_" + suffix + ".root" << endl;
}
