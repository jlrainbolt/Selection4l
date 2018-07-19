#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"


using namespace std;


void GetQuantities(const TString rootFile,
                   Double_t *nEvents, Double_t *nUnc, Double_t *nFid, Double_t *nAcc);


void calculateBF(const TString filePath, const TString fileSuffix)
{
    // Constants
    Double_t Zto2l = 0.03366, f_nr = 0.96;


    // Storage for quantities                                       // Not actually used
    const unsigned N = 6;           TString prefix[N];              Double_t frac[N];
    unsigned MM = 0;                prefix[MM] = "mumu";            frac[MM] = 0.5;
    unsigned EE = 1;                prefix[EE] = "ee";              frac[EE] = 0.5;
    unsigned M4 = 2;                prefix[M4] = "4m";              frac[M4] = 0.2655;
    unsigned ME = 3;                prefix[ME] = "2m2e";            frac[ME] = 0.4690 * 0.5;
    unsigned EM = 4;                prefix[EM] = "2e2m";            frac[EM] = 0.4690 * 0.5;
    unsigned E4 = 5;                prefix[E4] = "4e";              frac[E4] = 0.2655;

    Double_t nEvents[N], nUnc[N], nFid[N], nAcc[N];




    // Table of yields
    cout << endl << endl << " YIELDS" << endl;
    cout << "==================================================";
    cout << "==================================================" << endl;
    cout << "           DATA\t\t\t\tALL MC\t\t\t\tOBS - BG" << endl;
    cout << "--------------------------------------------------";
    cout << "--------------------------------------------------" << endl;


    // Get quantities from files
    for (unsigned i = 0; i < N; i++)
    {
        cout << setw(5) << prefix[i] << "\t   ";
        GetQuantities(filePath + prefix[i] + "_" + fileSuffix + ".root", 
                      &nEvents[i], &nUnc[i], &nFid[i], &nAcc[i]);
        cout << endl;

        if (i == EE)
            cout << endl;
    }
    cout << "==================================================";
    cout << "==================================================" << endl;





    // Table of gen-level acceptance for sanity check
    cout << endl << endl << " ACCEPTANCE" << endl;
    cout << "==================================================";
    cout << "==================================================" << endl;
    cout << "           FIDUCIAL (GEN)\t\tACCEPTED (RECO)\t\t\tSCALE FACTOR" << endl;
    cout << "--------------------------------------------------";
    cout << "--------------------------------------------------" << endl;


    // Find scale factors
    Double_t scale[N], s_err[N];
    for (unsigned i = 0; i < N; i++)
    {
        scale[i] = nFid[i] / nAcc[i];
        s_err[i] = sqrt(1./nFid[i] + 1./nAcc[i]) * scale[i];

        cout << setw(5) << prefix[i] << "\t   ";
        cout << setw(11) << nFid[i] << "  +-  " << setw(7) << sqrt(nFid[i]) << "\t";
        cout << setw(11) << nAcc[i] << "  +-  " << setw(7) << sqrt(nAcc[i]) << "\t";
        cout << setw(11) << scale[i] << "  +-  " << setw(7) << s_err[i] << endl;

        if (i == EE)
            cout << endl;
    }
    cout << "==================================================";
    cout << "==================================================" << endl;





    // Table of scaled yields to compare to expectation
    cout << endl << endl << " SCALED YIELDS" << endl;
    cout << "==================================================";
    cout << "==================================================" << endl;
    cout << "           OBS - BG\t\t\tFRACTION\t\t\tEXPECTED (?)" << endl;
    cout << "--------------------------------------------------";
    cout << "--------------------------------------------------" << endl;


    // Correct yields for acceptance/efficiency
    Double_t yield[N], error[N];
    for (unsigned i = 0; i < N; i++)
    {
        yield[i] = nEvents[i] * scale[i];
//      error[i] = nUnc[i] * scale[i];
        error[i] = (nUnc[i] * nUnc[i]) / (nEvents[i] * nEvents[i]);
        error[i] += 1./nFid[i] + 1./nAcc[i];
        error[i] = sqrt(error[i]) * yield[i];

        cout << setw(5) << prefix[i] << "\t   ";
        cout << setw(11) << yield[i] << "  +-  " << setw(7) << error[i] << "\t";
        cout << "\t\t\t\t";
        cout << setw(11) << frac[i] << endl;

        if (i == EE)
            cout << endl;
    }
    cout << "==================================================";
    cout << "==================================================" << endl;




    // Calculate total BF (no lepton discrimination)
    unsigned T = 0;
    Double_t result[N], uncf[N], unc[N];

    Double_t yield_2l = 0, yield_4l = 0, var_2l = 0, var_4l = 0;
    for (unsigned i = 0; i < M4; i++)
    {
        yield_2l += yield[i];
        var_2l += error[i] * error[i];
    }
    for (unsigned i = M4; i < N; i++)
    {
        yield_4l += yield[i];
        var_4l += error[i] * error[i];
    }
    result[T] = 2. * Zto2l * yield_4l / yield_2l;
    uncf[T] = var_2l / (yield_2l * yield_2l);
    uncf[T] += var_4l / (yield_4l * yield_4l);
    uncf[T] = sqrt(uncf[T]);
    unc[T] = uncf[T] * result[T];


    // Individual BFs for each 4l channel
    for (unsigned i = M4; i < N; i++)
    {
        // Index of DY selection
        unsigned j = i/2 - 1;

        result[i] = Zto2l * yield[i] / yield[j]; 
        uncf[i] += (error[i] * error[i]) / (yield[i] * yield[i]);
        uncf[i] += (error[j] * error[j]) / (yield[j] * yield[j]);
        uncf[i] = sqrt(uncf[i]);
        unc[i] = uncf[i] * result[i];
    }


    // Sum of individual channels
    unsigned S = 1;
    result[S] = 0;
    for (unsigned i = M4; i < N; i++)
    {
        result[S] += result[i];
        unc[S] += unc[i] * unc[i];
    }
    unc[S] = sqrt(unc[S]);
    uncf[S] = unc[S] / result[S];


    // For printing
    TString prefix2[N],             suffix[N];
    prefix2[T] = "4l";              suffix[T] = "(TOTAL)\t";
    prefix2[S] = "4l";              suffix[S] = "(SUM)\t";
    prefix2[M4] = "4m";             suffix[M4] = "(MUMU)\t";
    prefix2[ME] = "2m2e";           suffix[ME] = "(MUMU)\t";
    prefix2[EM] = "2e2m";           suffix[EM] = "(EE)\t";
    prefix2[E4] = "4e";             suffix[E4] = "(EE)\t";




    // Print table of BFs
    cout << endl << endl << " BRANCHING FRACTIONS" << endl;
    cout << "==================================================";
    cout << "==================================================" << endl;
    cout << " CHANNEL\t\t\tRESULT\t\t\t\t\tFRAC UNC" << endl;
    cout << "--------------------------------------------------";
    cout << "--------------------------------------------------" << endl;
    for (unsigned i = 0; i < N; i++)
    {
        cout << " Z -> " << prefix2[i] << "\t" << suffix[i] << "\t";
        cout << setw(11) << result[i] << "  +-  " << setw(8) << unc[i] << "\t\t";
        cout << uncf[i] << endl;

        if (i == S)
            cout << endl;
    }
    cout << "--------------------------------------------------";
    cout << "--------------------------------------------------" << endl << endl << endl;
}






void GetQuantities(const TString rootFile,
                   Double_t *nEvents, Double_t *nUnc, Double_t *nFid, Double_t *nAcc)
{
    TFile *file = TFile::Open(rootFile, "READ");
    TDirectory *dir = file->GetDirectory("/Calculation", kTRUE, "GetDirectory");
    TH1 *hAccepted;
    dir->GetObject("AcceptedEvents", hAccepted);
    *nFid = hAccepted->GetBinContent(2);
    *nAcc = hAccepted->GetBinContent(3);

    TH1 *hData;
    TDirectory *data_subdir = dir->GetDirectory("All Data");
    data_subdir->GetObject("TotalEvents", hData);
    cout << setw(11) << hData->GetBinContent(7) << "  +-  " << setw(7) << hData->GetBinError(7);
    cout << "\t"; 

    TH1 *hAllMC;
    TDirectory *allMC_subdir = dir->GetDirectory("All MC");
    allMC_subdir->GetObject("TotalEvents", hAllMC);
    cout << setw(11) << hAllMC->GetBinContent(7) << "  +-  " << setw(7) << hAllMC->GetBinError(7);
    cout << "\t";

    TH1 *hBgMC;
    TDirectory *bgMC_subdir = dir->GetDirectory("Background MC");
    bgMC_subdir->GetObject("TotalEvents", hBgMC);
    
    TH1 *hDiff = (TH1*) hData->Clone();     hDiff->Add(hBgMC, -1);
    *nEvents = hDiff->GetBinContent(7);
    *nUnc = hDiff->GetBinError(7);
    cout << setw(11) << *nEvents << "  +-  " << setw(7) << *nUnc;


    file->Close();
}
