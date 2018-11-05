// STL
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"

// Cuts
#include "Cuts2017.hh"

using namespace std;


/*
**  CalculateAccEff2l
**
**  Calculates acceptance * efficiency for input file name
*/

void CalculateAccEff2l(const TString inName)
{

    //
    //  SAMPLE INFO
    //

    TString suffix = "zjets_m-50";
    const unsigned N = 3;   // Channel indices
    unsigned                LL = 0, MM = 1, EE = 2;
    unsigned chanIdx[N] = { 2,      3,      4};
    float nSelected[N], nPhaseSpace[N],xAccEff[N];
    TString selection[N] = {"ll", "mumu", "ee"};



    //
    //  INPUT FILE
    //

    TFile *inFile = TFile::Open(inName);

    cout << "Opened " << inName << endl;

    TH1 *hSelected;
    inFile->GetObject("SelectedEvents_" + suffix, hSelected);
    hSelected->SetDirectory(0);
    hSelected->Sumw2();

    inFile->Close();



    //
    //  GEN LEVEL
    //

    TString psPath = EOS_PATH + "/BLT/" + YEAR_STR + "/";
    TString prefix = "genHardProc";

    TString psName = psPath + "gen_" + suffix + "/" + prefix + "_" + suffix + ".root";
    TFile *psFile = TFile::Open(psName);

    cout << "Opened " << psName << endl;

    TH1 *hPhaseSpace;
    psFile->GetObject("PhaseSpaceEvents_" + suffix, hPhaseSpace);
    hPhaseSpace->SetDirectory(0);
    hPhaseSpace->Sumw2();

    psFile->Close();



    //
    //  CALCULATE
    //

    nSelected[1] = hSelected->GetBinContent(chanIdx[MM]);
    nSelected[2] = hSelected->GetBinContent(chanIdx[EE]);
    nSelected[0] = nSelected[1] + nSelected[2];

    nPhaseSpace[1] = hPhaseSpace->GetBinContent(chanIdx[MM]);
    nPhaseSpace[2] = hPhaseSpace->GetBinContent(chanIdx[EE]);
    nPhaseSpace[0] = nPhaseSpace[1] + nPhaseSpace[2];

    for (unsigned i = 0; i < N; i++)
        xAccEff[i] = nSelected[i] / nPhaseSpace[i];



    //
    //  PRINT
    //

    cout << endl << endl;
    cout << "\t\t" << "Selected" << "\t" << "Phase Space" << "\t" << "Acc * Eff" << endl;
    for (unsigned i = 0; i < N; i++)
    {
        cout << selection[i] << "\t\t" << setw(8) << nSelected[i] << "\t";
        cout << setw(8) << nPhaseSpace[i] << "\t" << setw(8) << xAccEff[i] << endl;
    }
}
