// STL
#include <iostream>
#include <fstream>

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

void CalculateAccEff2l(const TString suffix, const TString systematics, const TString hID)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 3;   // Channel indices
    unsigned                LL = 0, MM = 1, EE = 2;
    unsigned chanIdx[N] = { 2,      3,      4};
    float nSelected[N], nPhaseSpace[N],xAccEff[N];
    TString selection[N] = {"ll", "mumu", "ee"};



    //
    //  INPUT FILE
    //

    TString inPath = EOS_PATH + "/Selected/" + YEAR_STR + "/";
    TString inName = inPath + "selected_" + suffix + ".root";
//  TString inPath = "output/";
//  TString inName = systematics + hID + "_" + suffix + ".root";
//  TString inPath = "";
//  TString inName = inPath + "rescaled_" + suffix + ".root";
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
    TString psPrefix = "genHardProc";

    TString psName = psPath + "gen_" + suffix + "/" + psPrefix + "_" + suffix + ".root";
//  TString psName = "rescaled_phase_space.root";
    TFile *psFile = TFile::Open(psName);

    cout << "Opened " << psName << endl;

    TH1 *hPhaseSpace;
    psFile->GetObject("PhaseSpaceEvents_" + suffix, hPhaseSpace);
//  psFile->GetObject("PhaseSpaceEvents_phase_space", hPhaseSpace);
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



    //
    //  WRITE OUT
    //

    TString ofsName = systematics + "_" + suffix + "_" + hID + ".txt";
    ofstream ofs;
    ofs.open(inPath + "/" + ofsName);

    ofs << hID;
    for (unsigned i = 0; i < N; i++)
        ofs << "\t" << xAccEff[i];

    ofs << endl;
    ofs.close();
}
