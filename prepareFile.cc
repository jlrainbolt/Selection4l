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
    int year = 2016;    TString yearStr;    yearStr.Form("%i", year);
    TString output = selection + "_" + yearStr + ".root";



    //--- SCALING INFO ---//
    // Lepton names
    TString Lep12, Lep34, Lep, lep12, lep34, lep;
    Lep12   = muPair1 ? "Muon" : "Electron";        lep12   = muPair1 ? "muon" : "electron";
    Lep34   = muPair2 ? "Muon" : "Electron";        lep34   = muPair2 ? "muon" : "electron";
    Lep     = (muPair1 == muPair2) ? Lep12 : "Lepton";
    lep     = (muPair1 == muPair2) ? lep12 : "lepton";


    // Sample indices       // Suffixes                         // Cross sections 
    const unsigned N = 18;  TString suffix[N];                  Float_t xs[N];
    unsigned RD = 0;        suffix[RD] = lep12 + "_2016";       Float_t lumi = 35.9;
    unsigned DI = 1;        suffix[DI] = "dy_m-10to50";         xs[DI] = 18610.;
    unsigned YI = 2;        suffix[YI] = "dy_m-50";             xs[YI] = 5765.4;
    unsigned D1 = 3;        suffix[D1] = "z1jets_m-10to50";     xs[D1] = 730.3;     // 855.5;
    unsigned Y1 = 4;        suffix[Y1] = "z1jets_m-50";         xs[Y1] = 1012.0;    // 1198.88;
    unsigned D2 = 5;        suffix[D2] = "z2jets_m-10to50";     xs[D2] = 387.4;     // 466.1;
    unsigned Y2 = 6;        suffix[Y2] = "z2jets_m-50";         xs[Y2] = 334.7;     // 390.58;
    unsigned D3 = 7;        suffix[D3] = "z3jets_m-10to50";     xs[D3] = 95.02;     // 114.46;
    unsigned Y3 = 8;        suffix[Y3] = "z3jets_m-50";         xs[Y3] = 102.3;     // 113.28;
    unsigned D4 = 9;        suffix[D4] = "z4jets_m-10to50";     xs[D4] = 36.71;     // 36.4;
    unsigned Y4 = 10;       suffix[Y4] = "z4jets_m-50";         xs[Y4] = 54.52;     // 60.18;
    unsigned GH = 11;       suffix[GH] = "ggH_zz_4l";           xs[GH] = 0.01218;
    unsigned QH = 12;       suffix[QH] = "H_zz_4l";             xs[QH] = 0.001044;
    unsigned TT = 13;       suffix[TT] = "ttbar";               xs[TT] = 831.76;
    unsigned TZ = 14;       suffix[TZ] = "ttz_2l2nu";           xs[TZ] = 0.2529;
    unsigned WW = 15;       suffix[WW] = "ww_2l2nu";            xs[WW] = 12.178;
    unsigned WZ = 16;       suffix[WZ] = "wz_3lnu";             xs[WZ] = 4.42965;
    unsigned ZZ = 17;       suffix[ZZ] = "zz_4l";               xs[ZZ] = 1.212;

/*    
    // Create colormap (Matlab "lines")
    unsigned lBlack = 0, lBlue = 1, lOrange = 2, lYellow = 3, lPurple = 4, lGreen = 5, lLtBlue = 6,
             lRed = 7;
    const unsigned L = 8;   Int_t colIdx[L];
    for (unsigned i = 0; i < L; i++)
        colIdx[i] = 1179 + (Int_t) i;
*/

    // Colors                   // Data tag             // Signal tag
    Int_t col[N];               Bool_t data[N],         signal[N];
    col[RD] = kBlack;           data[RD] = kTRUE;       signal[RD] = kFALSE;
    col[DI] = lOrange;          data[DI] = kFALSE;      signal[DI] = sel2l ? kTRUE : kFALSE;
    col[YI] = lOrange;          data[YI] = kFALSE;      signal[YI] = sel2l ? kTRUE : kFALSE;
    col[D1] = lRed;             data[D1] = kFALSE;      signal[D1] = sel2l ? kTRUE : kFALSE;
    col[Y1] = lRed;             data[Y1] = kFALSE;      signal[Y1] = sel2l ? kTRUE : kFALSE;
    col[D2] = lRed;             data[D2] = kFALSE;      signal[D2] = sel2l ? kTRUE : kFALSE;
    col[Y2] = lRed;             data[Y2] = kFALSE;      signal[Y2] = sel2l ? kTRUE : kFALSE;
    col[D3] = lRed;             data[D3] = kFALSE;      signal[D3] = sel2l ? kTRUE : kFALSE;
    col[Y3] = lRed;             data[Y3] = kFALSE;      signal[Y3] = sel2l ? kTRUE : kFALSE;
    col[D4] = lRed;             data[D4] = kFALSE;      signal[D4] = sel2l ? kTRUE : kFALSE;
    col[Y4] = lRed;             data[Y4] = kFALSE;      signal[Y4] = sel2l ? kTRUE : kFALSE;
    col[GH] = lPurple;          data[GH] = kFALSE;      signal[GH] = kFALSE;
    col[QH] = lPurple;          data[QH] = kFALSE;      signal[QH] = kFALSE;
    col[TT] = lBlue;            data[TT] = kFALSE;      signal[TT] = kFALSE;
    col[TZ] = lBlue;            data[TZ] = kFALSE;      signal[TZ] = kFALSE;
    col[WW] = lLightBlue;       data[WW] = kFALSE;      signal[WW] = kFALSE;
    col[WZ] = lLightBlue;       data[WZ] = kFALSE;      signal[WZ] = kFALSE;
    col[ZZ] = lGreen;           data[ZZ] = kFALSE;      signal[ZZ] = sel4l ? kTRUE : kFALSE;


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

        val.push_back(make_tuple("z1m",     "z1p4.M()",     100,    80,     100));
        val.push_back(make_tuple("z1pt",    "z1p4.Pt()",    100,    0,      200));

        val.push_back(make_tuple("l1pt",    "l1p4.Pt()",    100,    0,      200));
        val.push_back(make_tuple("l1eta",   "l1p4.Eta()",   100,    -2.5,   2.5));
        val.push_back(make_tuple("l1iso",   "l1iso",        100,    0,      0.35));
        val.push_back(make_tuple("l1pdg",   "l1pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l2pt",    "l2p4.Pt()",    100,    0,      200));
        val.push_back(make_tuple("l2eta",   "l2p4.Eta()",   100,    -2.5,   2.5));
        val.push_back(make_tuple("l2iso",   "l2iso",        100,    0,      0.35));
        val.push_back(make_tuple("l2pdg",   "l2pdg",        31,     -15.5,  15.5));
    }
    else if (sel4l)
    {
//      int bins = (muPair1 && muPair2) ?  
        val.push_back(make_tuple("nPV",     "nPV",          25,     0,      50));
        val.push_back(make_tuple("met",     "met",          20,     0,      100));

        val.push_back(make_tuple("zzm",     "zzp4.M()",     20,     80,     100));
        val.push_back(make_tuple("zzpt",    "zzp4.Pt()",    20,     0,      200));

        val.push_back(make_tuple("z1m",     "z1p4.M()",     20,     0,      100));
        val.push_back(make_tuple("z1pt",    "z1p4.Pt()",    20,     0,      200));

        val.push_back(make_tuple("z2m",     "z2p4.M()",     20,     0,      50));
        val.push_back(make_tuple("z2pt",    "z2p4.Pt()",    20,     0,      100));

        val.push_back(make_tuple("l1pt",    "l1p4.Pt()",    20,     0,      200));
        val.push_back(make_tuple("l1eta",   "l1p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l1iso",   "l1iso",        35,     0,      0.35));
        val.push_back(make_tuple("l1pdg",   "l1pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l2pt",    "l2p4.Pt()",    20,     0,      100));
        val.push_back(make_tuple("l2eta",   "l2p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l2iso",   "l2iso",        35,     0,      0.35));
        val.push_back(make_tuple("l2pdg",   "l2pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l3pt",    "l3p4.Pt()",    20,     0,      50));
        val.push_back(make_tuple("l3eta",   "l3p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l3iso",   "l3iso",        35,     0,      0.35));
        val.push_back(make_tuple("l3pdg",   "l3pdg",        31,     -15.5,  15.5));

        val.push_back(make_tuple("l4pt",    "l4p4.Pt()",    20,     0,      50));
        val.push_back(make_tuple("l4eta",   "l4p4.Eta()",   25,     -2.5,   2.5));
        val.push_back(make_tuple("l4iso",   "l4iso",        35,     0,      0.35));
        val.push_back(make_tuple("l4pdg",   "l4pdg",        31,     -15.5,  15.5));


        // Differential distributions
        // FIXME to include overflow?
        val.push_back(make_tuple("l1pt_d",  "l1p4.Pt()",    8,      20,     100));
        val.push_back(make_tuple("z2m_d",   "z2p4.M()",     8,      0,      40));
        val.push_back(make_tuple("m3l_d",   "tlp4.M()",     8,      0,      80));
        val.push_back(make_tuple("angle_d", "angle",        8,      0,      3.2));
        val.push_back(make_tuple("cos_d",   "cos(angle)",   8,      -1,     1));
        val.push_back(make_tuple("cos2_d",  "cos(angle)*cos(angle)", 10, 0, 1));
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
