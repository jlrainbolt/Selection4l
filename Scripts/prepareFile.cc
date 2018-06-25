#include "TString.h"
#include "TFile.h"
#include "TParameter.h"


using namespace std;

void prepareFile(const TString selection)
{
    // Choose electron or muon
    Bool_t sel2l = kFALSE, sel4l = kFALSE;
    Bool_t sel2m = kFALSE, sel2e = kFALSE;
    Bool_t sel4e = kFALSE, sel4m = kFALSE;
//  Bool_t sel2e2m = kFALSE, sel2m2e = kFALSE;
    if (selection == "mumu" || selection == "2m")
    {
        sel2l = kTRUE;
        sel2m = kTRUE;
    }
    else if (selection == "ee" || selection == "2e")
    {
        sel2l = kTRUE;
        sel2e = kTRUE;
    }
    else if (selection == "4e")
    {
        sel4l = kTRUE;
        sel4e = kTRUE;
    }
    else if (selection == "4mu" || selection == "4m")
    {
        sel4l = kTRUE;
        sel4m = kTRUE;
    }


    // Name of output file
    TString year   = "_2016";
    TString output = selection + year + ".root";
    TString lepton, Lepton;
    if (sel2m)
    {
        lepton = "muon";        Lepton = "Muon";
    }
    else if (sel2e)
    {
        lepton = "electron";    Lepton = "Electron";
    }
    else if (sel4e)
    {
        lepton = "electron";    Lepton = "Electron";
    }
    else if (sel4m)
    {
        lepton = "muon";        Lepton = "muon";
    }


    // Sample indices           // Suffixes                         // Cross sections 
    const unsigned N  = 8;      TString suffix[N];                  Float_t xs[N];
    const unsigned RD = 0;      suffix[RD] = lepton + "_2016";      Float_t lumi = 35.9;
    const unsigned D1 = 1;      suffix[D1] = "dy_m-10to50";         xs[D1] = 18610.;
    const unsigned D5 = 2;      suffix[D5] = "dy_m-50";             xs[D5] = 5765.4;
//  const unsigned D1 = 1;      suffix[D1] = "zjets_m-10to50";      xs[D1] = 18610.; //21959.8;
//  const unsigned D5 = 2;      suffix[D5] = "zjets_m-50";          xs[D5] = 5765.4; //6803.172;
    const unsigned TT = 3;      suffix[TT] = "ttbar";               xs[TT] = 831.76;
    const unsigned TZ = 4;      suffix[TZ] = "ttz_2l2nu";           xs[TZ] = 0.2529;
    const unsigned WW = 5;      suffix[WW] = "ww_2l2nu";            xs[WW] = 12.178;
    const unsigned WZ = 6;      suffix[WZ] = "wz_3lnu";             xs[WZ] = 4.42965;
    const unsigned ZZ = 7;      suffix[ZZ] = "zz_4l";               xs[ZZ] = 1.212;


    // Colors                   // Data tag
    Int_t col[N];               Bool_t data[N];
    col[RD] = kBlack;           data[RD] = kTRUE;
    col[D1] = kYellow;          data[D1] = kFALSE;
    col[D5] = kRed;             data[D5] = kFALSE;
    col[TT] = kBlue;            data[TT] = kFALSE;
    col[TZ] = kMagenta;         data[TZ] = kFALSE;
    col[WW] = kCyan;            data[WW] = kFALSE;
    col[WZ] = kOrange;          data[WZ] = kFALSE;
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
