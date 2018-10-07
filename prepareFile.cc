#include <tuple>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"


using namespace std;



void prepareFile(const TString path, const TString selection)
{
    //--- PARSING ---//


    // Selection
    Bool_t sel2l = kFALSE, sel4l = kFALSE;
    Bool_t muPair1, muPair2;    // TRUE for muon, FALSE for electron
    if (selection == "mumu" || selection == "2mu" || selection == "2m")
    {
        sel2l = kTRUE;      muPair1 = kTRUE;    muPair2 = muPair1;
    }
    else if (selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;      muPair1 = kFALSE;   muPair2 = muPair1;
    }
    else if (selection == "4mu" || selection == "4m")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;    muPair2 = kTRUE;
    }
    else if (selection == "2mu2e" || selection == "2m2e")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;    muPair2 = kFALSE;
    }
    else if (selection == "2e2mu" || selection == "2e2m")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;   muPair2 = kTRUE;
    }
    else if (selection == "4e")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;   muPair2 = kFALSE;
    }


    // Output name (in future pass outSuff as arg?)
    int year = 2017;    TString yearStr;    yearStr.Form("%i", year);
    TString output = selection + "_" + yearStr + ".root";



    //--- SCALING INFO ---//
    // Lepton names
    TString Lep12, Lep34, Lep, lep12, lep34, lep;
    Lep12   = muPair1 ? "Muon" : "Electron";        lep12   = muPair1 ? "muon" : "electron";
    Lep34   = muPair2 ? "Muon" : "Electron";        lep34   = muPair2 ? "muon" : "electron";
    Lep     = (muPair1 == muPair2) ? Lep12 : "Lepton";
    lep     = (muPair1 == muPair2) ? lep12 : "lepton";

    Float_t lumi = muPair1 ? 36.735 : 41.529;


    // Sample indices       // Suffixes                         // Cross sections 
    const unsigned N = 10;  TString suffix[N];                  Float_t xs[N];
    unsigned RD = 0;        suffix[RD] = lep12 + "_" + yearStr;
    unsigned DY = 1;        suffix[DY] = "zjets_m-50";          xs[DY] = 5765.4;
    unsigned GH = 2;        suffix[GH] = "ggH_zz_4l";           xs[GH] = 0.01212;
    unsigned VH = 3;        suffix[VH] = "vbfH_zz_4l";          xs[VH] = 0.001034;
    unsigned TT = 4;        suffix[TT] = "ttbar";               xs[TT] = 831.76;
    unsigned WW = 5;        suffix[WW] = "ww_2l2nu";            xs[WW] = 12.178;
    unsigned W2 = 6;        suffix[W2] = "wz_2l2q";             xs[W2] = 5.595;
    unsigned W3 = 7;        suffix[W3] = "wz_3lnu";             xs[W3] = 4.42965;
    unsigned Z2 = 8;        suffix[Z2] = "zz_2l2q";             xs[Z2] = 3.22;
    unsigned Z4 = 9;        suffix[Z4] = "zz_4l";               xs[Z4] = 1.212;
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2016#Samples_Cross_sections


    // Colors                   // Data tag             // Signal tag
    Int_t col[N];               Bool_t data[N],         signal[N];
    col[RD] = kBlack;           data[RD] = kTRUE;       signal[RD] = kFALSE;
    col[DY] = lYellow;          data[DY] = kFALSE;      signal[DY] = sel2l ? kTRUE : kFALSE;
    col[GH] = lPurple;          data[GH] = kFALSE;      signal[GH] = kFALSE;
    col[VH] = lPurple;          data[VH] = kFALSE;      signal[VH] = kFALSE;
    col[TT] = lGreen;           data[TT] = kFALSE;      signal[TT] = kFALSE;
    col[WW] = lOrange;          data[WW] = kFALSE;      signal[WW] = kFALSE;
    col[W2] = lOrange;          data[W2] = kFALSE;      signal[W2] = kFALSE;
    col[W3] = lOrange;          data[W3] = kFALSE;      signal[W3] = kFALSE;
    col[Z2] = lOrange;          data[Z2] = kFALSE;      signal[Z2] = kFALSE;
    col[Z4] = lLightBlue;       data[Z4] = kFALSE;      signal[Z4] = sel4l ? kTRUE : kFALSE;


    // Convert to TParameters
    TParameter<Float_t> *luminosity = new TParameter<Float_t>("lumi", lumi);
    TParameter<Float_t> *xsec[N];
    TParameter<Int_t> *color[N];
    TParameter<Bool_t> *isData[N], *isSignal[N];
    for (unsigned i = 0; i < N; i++)
    {
        xsec[i] = new TParameter<Float_t>("xsec", xs[i]);
        isData[i] = new TParameter<Bool_t>("isData", data[i]);
        isSignal[i] = new TParameter<Bool_t>("isSignal", signal[i]);
        color[i] = new TParameter<Int_t>("color", col[i]);
    }

    // Draw validation histograms
    vector<tuple<TString, TString, Int_t, Double_t, Double_t>> val;
    //                           name       member          bins    xmin    xmax
    if (sel2l)
    {
        val.push_back(make_tuple("nPV",     "nPV",          51,     -0.5,   50.5));
        val.push_back(make_tuple("met",     "met",          100,    0,      100));

        val.push_back(make_tuple("z1m",     "z1p4.M()",     40,    80,     100));
        val.push_back(make_tuple("z1pt",    "z1p4.Pt()",    40,    0,      80));

        val.push_back(make_tuple("l1pt",    "l1p4.Pt()",    40,    20,     100));
        val.push_back(make_tuple("l1eta",   "l1p4.Eta()",   40,    -2.5,   2.5));
        val.push_back(make_tuple("l1iso",   "l1iso",        100,    0,      0.35));
        val.push_back(make_tuple("l1pdg",   "l1pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l2pt",    "l2p4.Pt()",    40,     0,      80));
        val.push_back(make_tuple("l2eta",   "l2p4.Eta()",   40,    -2.5,   2.5));
        val.push_back(make_tuple("l2iso",   "l2iso",        100,    0,      0.35));
        val.push_back(make_tuple("l2pdg",   "l2pdg",        31,     -15.5,  15.5));
    }
    else if (sel4l)
    {
        val.push_back(make_tuple("nPV",     "nPV",          25,     0,      50));
        val.push_back(make_tuple("met",     "met",          20,     0,      100));

        val.push_back(make_tuple("zzm",     "zzp4.M()",     20,     80,     100));
        val.push_back(make_tuple("zzpt",    "zzp4.Pt()",    20,     0,      80));

        val.push_back(make_tuple("z1m",     "z1p4.M()",     20,     10,      90));
        val.push_back(make_tuple("z1pt",    "z1p4.Pt()",    20,     0,      80));

        val.push_back(make_tuple("z2m",     "z2p4.M()",     20,     0,      40));
        val.push_back(make_tuple("z2pt",    "z2p4.Pt()",    20,     0,      80));

        val.push_back(make_tuple("l1pt",    "l1p4.Pt()",    20,     20,      100));
        val.push_back(make_tuple("l1eta",   "l1p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l1iso",   "l1iso",        35,     0.001,      0.35));
        val.push_back(make_tuple("l1pdg",   "l1pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l2pt",    "l2p4.Pt()",    20,     10,      90));
        val.push_back(make_tuple("l2eta",   "l2p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l2iso",   "l2iso",        35,     0.001,      0.35));
        val.push_back(make_tuple("l2pdg",   "l2pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l3pt",    "l3p4.Pt()",    20,     5,      45));
        val.push_back(make_tuple("l3eta",   "l3p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l3iso",   "l3iso",        35,     0.001,      0.35));
        val.push_back(make_tuple("l3pdg",   "l3pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l4pt",    "l4p4.Pt()",    20,     5,      25));
        val.push_back(make_tuple("l4eta",   "l4p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l4iso",   "l4iso",        35,     0.001,      0.35));
        val.push_back(make_tuple("l4pdg",   "l4pdg",        31,     -15.5,  15.5));


        // Differential distributions
        // FIXME to include overflow?
//      val.push_back(make_tuple("l1pt_d",  "l1p4.Pt()",    10,      20,     70));
//      val.push_back(make_tuple("z2m_d",   "z2p4.M()",     8,      0,      40));
//      val.push_back(make_tuple("m3l_d",   "tlp4.M()",     8,      0,      80));
//      val.push_back(make_tuple("angle_d", "angle",        8,      0,      3.2));
//      val.push_back(make_tuple("cos_d",   "cos(angle)",   8,      -1,     1));
//      val.push_back(make_tuple("cos2_d",  "cos(angle)*cos(angle)", 10, 0, 1));
    }


    // Draw histograms and write to file
    TFile *outFile = new TFile(output, "RECREATE");
    TDirectory *histDir = outFile->mkdir("Histograms");
    TDirectory *histSubdir[N];

    for (unsigned i = 0; i < N; i++)
    {
        cout << "Creating " << suffix[i] << " histograms..." << endl;

        TFile *inFile = TFile::Open(path + selection + "_" + suffix[i] + ".root");
        TTree *tree;
        inFile->GetObject("tree_" + suffix[i], tree);

        histSubdir[i] = histDir->mkdir(suffix[i]);
        histSubdir[i]->cd();
        for (unsigned j = 0; j < val.size(); j++)
        {
            TString varexp, hname = get<0>(val[j]);
            varexp.Form(get<1>(val[j]) + ">>" + hname + "(%i,%g,%g)",
                        get<2>(val[j]), get<3>(val[j]), get<4>(val[j]));
            tree->Draw(varexp, "weight");
            TH1F* hist;
            gDirectory->GetObject(hname, hist);

            if (hist)
                hist->SetName(hname + "_" + suffix[i]);
            else
                hist = new TH1F(hname + "_" + suffix[i], get<1>(val[j])+" {weight}",
                                get<2>(val[j]), get<3>(val[j]), get<4>(val[j]));

            hist->Write();
        }
        TH1D *hTotalEvents, *hAcceptedEvents;
        inFile->GetObject("TotalEvents_" + suffix[i], hTotalEvents);
        hTotalEvents->Write();

        inFile->GetObject("AcceptedEvents_" + suffix[i], hAcceptedEvents);
        hAcceptedEvents->Write();
        inFile->Close();

        isData[i]->Write();
        isSignal[i]->Write();
        color[i]->Write();
        if (i == 0)
            luminosity->Write();
        else
            xsec[i]->Write();
    }
    outFile->Close();
    
    cout << "Created file " + output + "." << endl;


    // Delete things
    delete outFile;
    delete luminosity;
    for (unsigned i = 0; i < N; i++)
    {
        delete xsec[i];
        delete isData[i];
        delete color[i];
        delete isSignal[i];
    }
}
