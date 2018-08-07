#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"


using namespace std;


void GetQuantities(const TString rootFile,
                   Double_t *nEvents, Double_t *nUnc, Double_t *nSpc, Double_t *nSel);


void calculateMuonBF(const TString filePath, const TString fileSuffix)
{
    // Constants
    Double_t Zto2l = 0.03366, f_nr = 0.96;


    // Storage for quantities                                       // Not actually used
    const unsigned N = 2;           TString prefix[N];              Double_t frac[N];
    unsigned MM = 0;                prefix[MM] = "mumu";            frac[MM] = 0.5;
    unsigned M4 = 1;                prefix[M4] = "4m";              frac[M4] = 0.2655;

    Double_t nEvents[N], nUnc[N], nSpc[N], nSel[N];




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
                      &nEvents[i], &nUnc[i], &nSpc[i], &nSel[i]);
        cout << endl;
    }
    cout << "==================================================";
    cout << "==================================================" << endl;




    // Table of scaled yields for sanity check
    cout << "           SCALED\t\t\tPHASE SPACE\t\t\tSELECTED" << endl;
    cout << "--------------------------------------------------";
    cout << "--------------------------------------------------" << endl;


    // Correct yields for acceptance/efficiency
    Double_t scale[N], yield[N], error[N];
    for (unsigned i = 0; i < N; i++)
    {
        scale[i] = nSpc[i] / nSel[i];
        yield[i] = nEvents[i] * scale[i];
        error[i] = nUnc[i] * scale[i];

        cout << setw(5) << prefix[i] << "\t   ";
        cout << setw(11) << yield[i] << "  +-  " << setw(7) << error[i] << "\t";
        cout << nSpc[i] << "\t\t\t" << nSel[i] << endl;
    }
    cout << "==================================================";
    cout << "==================================================" << endl << endl << endl;


    // Individual BFs for each 4l channel
    Double_t result[N], uncf[N], unc[N];
    for (unsigned i = M4; i < N; i++)
    {
        // Index of DY selection
        unsigned j = 0;

        result[i] = Zto2l * yield[i] / yield[j]; 
        uncf[i] += (error[i] * error[i]) / (yield[i] * yield[i]);
        uncf[i] += (error[j] * error[j]) / (yield[j] * yield[j]);
        uncf[i] = sqrt(uncf[i]);
        unc[i] = uncf[i] * result[i];
    }




    // Print table of BFs
    cout << endl << endl << " BRANCHING FRACTIONS" << endl;
    cout << "==================================================";
    cout << "==================================================" << endl;
    cout << " CHANNEL\t\tRESULT\t\t\t\t\tFRAC UNC" << endl;
    cout << "--------------------------------------------------";
    cout << "--------------------------------------------------" << endl;

    for (unsigned i = M4; i < N; i++)
    {
        cout << " Z -> " << prefix[i] << "\t\t";
        cout << setw(11) << result[i] << "  +-  " << setw(8) << unc[i] << "\t\t";
        cout << uncf[i] << endl;
    }
    cout << "==================================================";
    cout << "==================================================" << endl << endl << endl;
}






void GetQuantities(const TString rootFile,
                   Double_t *nEvents, Double_t *nUnc, Double_t *nSpc, Double_t *nSel)
{
    TFile *file = TFile::Open(rootFile, "READ");
    TDirectory *dir = file->GetDirectory("/Calculation", kTRUE, "GetDirectory");
    TH1 *hAccepted;
    dir->GetObject("AcceptedEvents", hAccepted);
    *nSpc = hAccepted->GetBinContent(2);
    *nSel = hAccepted->GetBinContent(3);

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
