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
**  CalculateAccEff4l
**
**  Calculates acceptance * efficiency for input file name
*/

void CalculateAccEff4l(const TString suffix, const TString systematics, const TString hID)
{

    //
    //  SAMPLE INFO
    //

    unsigned                L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;
    unsigned chanIdx[5] = { 5,      6,      7,      8,      9};

    // Reduced to calculate acc * eff
    const unsigned N = 4;   // Channel indices
    float nSelected[N], nPhaseSpace[N],xAccEff[N];
    TString selection[N] = {"4l", "4m", "2m2e", "4e"};



    //
    //  INPUT FILE
    //

    TString inPath = "output";
    TString inName = systematics + hID + "_" + suffix + ".root";
    TFile *inFile = TFile::Open(inPath + "/" + inName);

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

    nSelected[1] = hSelected->GetBinContent(chanIdx[M4]);
    nSelected[2] = hSelected->GetBinContent(chanIdx[ME]) + hSelected->GetBinContent(chanIdx[EM]);
    nSelected[3] = hSelected->GetBinContent(chanIdx[E4]);
    nSelected[0] = nSelected[1] + nSelected[2] + nSelected[3];

    nPhaseSpace[1] = hPhaseSpace->GetBinContent(chanIdx[M4]);
    nPhaseSpace[2] = hPhaseSpace->GetBinContent(chanIdx[ME]);
    nPhaseSpace[2] += hPhaseSpace->GetBinContent(chanIdx[EM]);
    nPhaseSpace[3] = hPhaseSpace->GetBinContent(chanIdx[E4]);
    nPhaseSpace[0] = nPhaseSpace[1] + nPhaseSpace[2] + nPhaseSpace[3];

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
