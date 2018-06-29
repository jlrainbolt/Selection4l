#include "TString.h"
#include "TFile.h"
#include "TParameter.h"


using namespace std;

void prepareFile(const TString selection)
{
    /* FILE INFO */
    // Parse selection
    Bool_t sel2l = kFALSE, sel4l = kFALSE;
    Bool_t muPair1, muPair2;    // TRUE for muon, FALSE for electron
    if (selection == "mumu" || selection == "2m")
    {   
        sel2l = kTRUE;      muPair1 = kTRUE;
    }
    else if (selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;      muPair1 = kFALSE;
    }
    else if (selection == "4mu" || selection == "4m")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;    muPair2 = kTRUE;
    }
    else if (selection == "4e")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;   muPair2 = kFALSE;
    }


    // Name of output file
    TString year   = "_2016";
    TString output = selection + year + ".root";



    /* SCALING INFO */
    // Lepton names
    TString Lep12, Lep34, Lep, lep12, lep34, lep;
    Lep12   = muPair1 ? "Muon" : "Electron";        lep12   = muPair1 ? "muon" : "electron";
    Lep34   = muPair2 ? "Muon" : "Electron";        lep34   = muPair2 ? "muon" : "electron";
    Lep     = (muPair1 == muPair2) ? Lep12 : "Lepton";
    lep     = (muPair1 == muPair2) ? lep12 : "lepton";


    // Sample indices       // Suffixes                         // Cross sections 
    const unsigned N = 10;  TString suffix[N];                  Float_t xs[N];
    unsigned RD = 0;        suffix[RD] = lep12 + "_2016";       Float_t lumi = 35.9;
    unsigned D1 = 1;        suffix[D1] = "dy_m-10to50";         xs[D1] = 18610.;
    unsigned D5 = 2;        suffix[D5] = "dy_m-50";             xs[D5] = 5765.4;
    unsigned GH = 3;        suffix[GH] = "ggH_zz_4l";           xs[GH] = 0.01218;
    unsigned QH = 4;        suffix[QH] = "H_zz_4l";             xs[QH] = 0.001044;
    unsigned TT = 5;        suffix[TT] = "ttbar";               xs[TT] = 831.76;
    unsigned TZ = 6;        suffix[TZ] = "ttz_2l2nu";           xs[TZ] = 0.2529;
    unsigned WW = 7;        suffix[WW] = "ww_2l2nu";            xs[WW] = 12.178;
    unsigned WZ = 8;        suffix[WZ] = "wz_3lnu";             xs[WZ] = 4.42965;
    unsigned ZZ = 9;        suffix[ZZ] = "zz_4l";               xs[ZZ] = 1.212;

    // Colors                   // Data tag
    Int_t col[N];               Bool_t data[N];
    col[RD] = kBlack;           data[RD] = kTRUE;
    col[D1] = kRed + 2;         data[D1] = kFALSE;
    col[D5] = kRed;             data[D5] = kFALSE;
    col[GH] = kMagenta + 2;     data[GH] = kFALSE;
    col[QH] = kMagenta;         data[QH] = kFALSE;
    col[TT] = kBlue;            data[TT] = kFALSE;
    col[TZ] = kBlue + 2;        data[TZ] = kFALSE;
    col[WW] = kCyan;            data[WW] = kFALSE;
    col[WZ] = kCyan + 2;;       data[WZ] = kFALSE;
    col[ZZ] = kGreen;           data[ZZ] = kFALSE;


    // Convert to TParameters
    TParameter<Float_t> *luminosity = new TParameter<Float_t>("lumi", lumi);
    TParameter<Float_t> *xsec[N];   TParameter<Int_t> *color[N];    TParameter<Bool_t> *isData[N];
    for (unsigned i = 0; i < N; i++)
    {
        xsec[i] = new TParameter<Float_t>("xsec", xs[i]);
        color[i] = new TParameter<Int_t>("color", col[i]);
        isData[i] = new TParameter<Bool_t>("isData", data[i]);
    }


    // Write to file
    TFile *outFile = new TFile(output, "NEW");
    TDirectory *dir = outFile->mkdir("Histograms");
    dir->cd();
    TDirectory *subdir[N];
    for (unsigned i = 0; i < N; i++)
    {
        subdir[i] = dir->mkdir(suffix[i]);
        subdir[i]->cd();
        isData[i]->Write();
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
        delete color[i];
        delete isData[i];
    }
}
