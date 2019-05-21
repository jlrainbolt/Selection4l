//
// A sample fit for measuring lepton efficiencies.
//
// Michael Schmitt     Northwestern University      May 20, 2019
//

#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TString.h"
#include "TExec.h"

#include "Cuts2017.hh"

double quad_sum(double a, double b) { return sqrt(a * a + b * b); }


//
//  PARAMETERS
//

TString selection = "mumu";

TFile *inFile, *outFile;
TH1 *passed, *failed;

const unsigned  N_PT_MM = 18,   N_ETA_MM = 4;
const unsigned  N_PT_EE = 17,   N_ETA_EE = 4;

float   mumuPt[N_PT_MM+1] = {5, 8,  10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 40, 45, 50, 60, 75, 100, 200 };
float   mumuEta[N_ETA_MM+1]={-2.4,  -1.2,       0,  1.2,    2.4 };

float   eePt[N_PT_EE+1] =   {7, 10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 40, 45, 50, 60, 75, 100, 200 };
float   eeEta[N_ETA_EE+1] = {-2.5,  -1.479,     0,  1.479,  2.5 };

const unsigned  N_PT = selection.Contains("mumu") ? N_PT_MM : N_PT_EE;
float   *binsPt = selection.Contains("mumu") ? mumuPt : eePt;
float   *binsEta = selection.Contains("mumu") ? mumuEta : eeEta;

const unsigned  N_M = 60;
const double    M_MIN_TNP = 60, M_MAX_TNP = 120;

const unsigned  N_PAR = 10;


//
//  GET HISTOGRAMS
//

void GetHistograms(TString dirName, TString histName)
{
    inFile->GetObject(dirName + "/passed_" + histName, passed);
    inFile->GetObject(dirName + "/failed_" + histName, failed);

    passed->SetDirectory(0);
    failed->SetDirectory(0);
}



//
//  FIT FUNCTIONS
//

// Signal is the sum of two modified Gaussians
double sig_func(double x, double mean, double amp1, double rms1, double del1,
               double amp2, double rms2, double del2)
{
    double gaus1 = amp1 * TMath::Gaus(x, mean, rms1 * (1 + del1 * x), 1);
    double gaus2 = amp2 * TMath::Gaus(x, mean, rms2 * (1 + del2 * x), 1);
    
    return gaus1 + gaus2;
}

// Background is a quadratic function
double bkg_func(double x, double b0, double b1, double b2)
{
    return b2 * pow(x, 2) + b1 * x + b0;
}

// Fit the total
double fit_func(double x, double mean, double amp1, double rms1, double del1,
       double amp2, double rms2, double del2, double b0, double b1, double b2)
{
    double sig = sig_func(x, mean, amp1, rms1, del1, amp2, rms2, del2);
    double bkg = bkg_func(x, b0, b1, b2);

    return sig + bkg;
}

double fit_func(double *V, double *par)     // overloaded version for hist fitter
{
    return fit_func(V[0],
            par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]);
}



//
//  HISTOGRAM FIT
//

void HistogramFit(double ptMin, double ptMax)
{

    // Set up fit function

    TF1 *func = new TF1("fit", fit_func, M_MIN_TNP, M_MAX_TNP, N_PAR);
    func->SetNpx(1000);     // number of points for drawing the function

    // Parameter names and initial values
    func->SetParNames(  "mean", "amp1", "rms1", "del1", "amp2", "rms2", "del2",
                                "b_0",  "b_1",  "b_2");
    int mean = func->GetParNumber("mean");
    int amp1 = func->GetParNumber("amp1"),      amp2 = func->GetParNumber("amp2");
    int rms1 = func->GetParNumber("rms1"),      rms2 = func->GetParNumber("rms2");
    int del1 = func->GetParNumber("del1"),      del2 = func->GetParNumber("del2");

    double pos_inf = 1e7;
    func->SetParLimits(mean, 88.5, 93.5);   // Make sure the Gaussians are used for the signal
    func->SetParLimits(rms1, 0, 1);         // And put them on a diet...
    if (ptMin > 45)
        func->SetParLimits(rms2, 0, 1);
    else
        func->SetParLimits(rms2, 0, 100);
    func->SetParLimits(del1, -1, 1);
    func->SetParLimits(del2, -1, 1);


    // Perform fit
    double  mean_0 = 91;
    double  rms1_0 = 0.5,   del1_0 = 0,     amp1_0 = passed->GetEntries();
    double  rms2_0 = 1,     del2_0 = 0,     amp2_0 = 0;
    double  b2_0 = 0,       b1_0 = 0,       b0_0 = 0.01 * passed->GetEntries();


    // Passed
    func->SetParameters(mean_0, amp1_0, rms1_0, del1_0, amp2_0, rms2_0, del2_0, b0_0, b1_0, b2_0);
    func->SetParLimits(amp1, 0, amp1_0);
    func->SetParLimits(amp2, 0, amp1_0);

    passed->Fit(func, "QBM", "", M_MIN_TNP, M_MAX_TNP);   // "Q" = quiet mode
    double P  = func->GetParameter(amp1) + func->GetParameter(amp2);
    double uP = quad_sum(func->GetParError(amp1), func->GetParError(amp2));

    // Failed
    amp1_0 = failed->GetEntries();
    b0_0 = 0.01 * failed->GetEntries();
    func->SetParameters(mean_0, amp1_0, rms1_0, del1_0, amp2_0, rms2_0, del2_0, b0_0, b1_0, b2_0);
    func->SetParLimits(amp1, 0, amp1_0);

    if (ptMin > 25)
        func->SetParLimits(amp2, 0, amp1_0);
    else
        func->SetParLimits(amp2, 0, 1e-7);

    failed->Fit(func, "QBM", "", M_MIN_TNP, M_MAX_TNP);
    double F  = func->GetParameter(amp1) + func->GetParameter(amp2);
    double uF = quad_sum(func->GetParError(amp1), func->GetParError(amp2));


    // Print results

    double eff = P/(F+P);
    double effUnc = eff * sqrt( pow( (1./P-1./(F+P))*uP, 2 ) + pow( -1./(F+P)*uF, 2 ) );

    cout << fixed << setprecision(6) << eff << " +/- " << effUnc << endl;
}

void PlotHistograms(TString histName)
{
    outFile->cd();
//  gROOT->SetStyle("Plain");
//  gStyle->SetOptStat(10);
//  gStyle->SetOptFit(11);  
//  gStyle->SetStatFontSize(0.04);

    TCanvas *CP = new TCanvas("passed_" + histName, "", 800, 500);
    CP->SetGridx(1);
    CP->SetGridy(1);
    CP->AddExec("stat", "gStyle->SetOptStat(10)");
    CP->AddExec("fit", "gStyle->SetOptFit(11)");
    passed->SetMinimum(0);
    passed->Draw("E");
    CP->AutoExec();
//  CP->Print(fitType + "_" + ptRange + "_passed.pdf");
    CP->Write();

    TCanvas *CF = new TCanvas("failed_" + histName, "", 800, 500);
    CF->SetGridx(1);
    CF->SetGridy(1);
    CP->AddExec("stat", "gStyle->SetOptStat(10)");
    CP->AddExec("fit", "gStyle->SetOptFit(11)");
    failed->SetMinimum(0);
    failed->Draw("E");
    CF->AutoExec();
//  CF->Print(fitType + "_" + ptRange + "_failed.pdf");
    CF->Write();
}





void sampleFit2()
{
    // Open file
    TString prefix = "mumu", suffix = "muon_2017", tag = "iso35";
    TString inName = prefix + "_" + tag + "_" + suffix + ".root";
    TString inPath = EOS_PATH + "/TagAndProbe/" + YEAR_STR + "/";

    inFile = TFile::Open(inPath + inName);
    cout << endl << endl << "Opened " << inPath + inName << endl;

    TString dirName = "hists_" + suffix;


    // Create output file
    TString outName = inName;
    outFile = new TFile(outName, "RECREATE");

    cout << "Min pt" << "\t" << "Max pt" << "\t" << "Efficiency" << endl;

    // Loop
    for (unsigned i = 0; i <= N_PT; i++)
    {
        double ptMin = binsPt[i], ptMax = binsPt[i+1];

        TString ptString;
        if (i < N_PT)
            ptString.Form("pt(%g,%g)", ptMin, ptMax);
        else
            ptString.Form("pt(%g+)", ptMin);
        TString histName = ptString + "_eta(-2.4,2.4)";

        cout << fixed << setprecision(0) << ptMin << "\t" << ptMax << "\t" << flush;

        GetHistograms(dirName, histName);
        HistogramFit(ptMin, ptMax);
        PlotHistograms(histName);
    }

    outFile->Write();
    outFile->Close();
    cout << endl << "Wrote canvases to " << outName << endl;

    inFile->Close();
    cout << "Closed " << inPath + inName << endl;
}
