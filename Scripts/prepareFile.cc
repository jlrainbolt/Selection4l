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
    else if (selection == "2mu2e" || selection == "2m2e")
    {
        sel4l = kTRUE;      muPair1 = kTRUE;    muPair2 = kFALSE;
    }
    else if (selection == "2e2mu" || selection == "2e2m")
    {
        sel4l = kTRUE;      muPair1 = kFALSE;    muPair2 = kTRUE;
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
    const unsigned N = 18;  TString suffix[N];                  Float_t xs[N];
    unsigned RD = 0;        suffix[RD] = lep12 + "_2016";       Float_t lumi = 35.9;
    unsigned D1 = 1;        suffix[D1] = "dy_m-10to50";         xs[D1] = 18610.;
    unsigned D5 = 2;        suffix[D5] = "dy_m-50";             xs[D5] = 5765.4;
    unsigned Z1 = 3;        suffix[Z1] = "z1jets_m-10to50";     xs[Z1] = 730.3;     // 855.5;
    unsigned Y1 = 4;        suffix[Y1] = "z1jets_m-50";         xs[Y1] = 1012.0;    // 1198.88;
    unsigned Z2 = 5;        suffix[Z2] = "z2jets_m-10to50";     xs[Z2] = 387.4;     // 466.1;
    unsigned Y2 = 6;        suffix[Y2] = "z2jets_m-50";         xs[Y2] = 334.7;     // 390.58;
    unsigned Z3 = 7;        suffix[Z3] = "z3jets_m-10to50";     xs[Z3] = 95.02;     // 114.46;
    unsigned Y3 = 8;        suffix[Y3] = "z3jets_m-50";         xs[Y3] = 102.3;     // 113.28;
    unsigned Z4 = 9;        suffix[Z4] = "z4jets_m-10to50";     xs[Z4] = 36.71;     // 36.4;
    unsigned Y4 = 10;       suffix[Y4] = "z4jets_m-50";         xs[Y4] = 54.52;     // 60.18;
    unsigned GH = 11;       suffix[GH] = "ggH_zz_4l";           xs[GH] = 0.01218;
    unsigned QH = 12;       suffix[QH] = "H_zz_4l";             xs[QH] = 0.001044;
    unsigned TT = 13;       suffix[TT] = "ttbar";               xs[TT] = 831.76;
    unsigned TZ = 14;       suffix[TZ] = "ttz_2l2nu";           xs[TZ] = 0.2529;
    unsigned WW = 15;       suffix[WW] = "ww_2l2nu";            xs[WW] = 12.178;
    unsigned WZ = 16;       suffix[WZ] = "wz_3lnu";             xs[WZ] = 4.42965;
    unsigned ZZ = 17;       suffix[ZZ] = "zz_4l";               xs[ZZ] = 1.212;

    // Colors                   // Data tag             // Signal tag
    Int_t col[N];               Bool_t data[N],         signal[N];
    col[RD] = kBlack;           data[RD] = kTRUE;       signal[RD] = kFALSE;
    col[D1] = kRed + 2;         data[D1] = kFALSE;      signal[D1] = sel2l ? kTRUE : kFALSE;
    col[D5] = kRed;             data[D5] = kFALSE;      signal[D5] = sel2l ? kTRUE : kFALSE;
    col[Z1] = kRed - 3;         data[Z1] = kFALSE;      signal[Z1] = sel2l ? kTRUE : kFALSE;
    col[Y1] = kRed - 4;         data[Y1] = kFALSE;      signal[Y1] = sel2l ? kTRUE : kFALSE;
    col[Z2] = kRed - 5;         data[Z2] = kFALSE;      signal[Z2] = sel2l ? kTRUE : kFALSE;
    col[Y2] = kRed - 6;         data[Y2] = kFALSE;      signal[Y2] = sel2l ? kTRUE : kFALSE;
    col[Z3] = kRed - 7;         data[Z3] = kFALSE;      signal[Z3] = sel2l ? kTRUE : kFALSE;
    col[Y3] = kRed - 8;         data[Y3] = kFALSE;      signal[Y3] = sel2l ? kTRUE : kFALSE;
    col[Z4] = kRed - 9;         data[Z4] = kFALSE;      signal[Z4] = sel2l ? kTRUE : kFALSE;
    col[Y4] = kRed - 10;        data[Y4] = kFALSE;      signal[Y4] = sel2l ? kTRUE : kFALSE;
    col[GH] = kMagenta + 2;     data[GH] = kFALSE;      signal[GH] = kFALSE;
    col[QH] = kMagenta;         data[QH] = kFALSE;      signal[QH] = kFALSE;
    col[TT] = kBlue;            data[TT] = kFALSE;      signal[TT] = kFALSE;
    col[TZ] = kBlue + 2;        data[TZ] = kFALSE;      signal[TZ] = kFALSE;
    col[WW] = kCyan;            data[WW] = kFALSE;      signal[WW] = kFALSE;
    col[WZ] = kCyan + 2;;       data[WZ] = kFALSE;      signal[WZ] = kFALSE;
    col[ZZ] = kGreen;           data[ZZ] = kFALSE;      signal[ZZ] = sel4l ? kTRUE : kFALSE;


    // Convert to TParameters
    TParameter<Float_t> *luminosity = new TParameter<Float_t>("lumi", lumi);
    TParameter<Float_t> *xsec[N];   TParameter<Int_t> *color[N];
    TParameter<Bool_t> *isData[N], *isSignal[N];
    for (unsigned i = 0; i < N; i++)
    {
        xsec[i] = new TParameter<Float_t>("xsec", xs[i]);
        color[i] = new TParameter<Int_t>("color", col[i]);
        isData[i] = new TParameter<Bool_t>("isData", data[i]);
        isSignal[i] = new TParameter<Bool_t>("isSignal", signal[i]);
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
        delete color[i];
        delete isData[i];
        delete isSignal[i];
    }
}
