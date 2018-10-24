#include <vector>
#include <sstream>
#include <tuple>

#include "TString.h"
#include "TFile.h"
#include "TParameter.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMathText.h"


using namespace std;


/*
**  DivideDists
**
**  Creates acceptance/efficiency ratio histograms for distributions in "unscaled_"
**  and maybe also overlays them on canvases if I feel up to it someday
*/

void DivideDists()
{

    //
    //  OPTIONS
    //

    const int nRebin = 5;



    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                    L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = { "4l",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = { 5,      6,      7,      8,      9};


    const unsigned M = 3;
    unsigned                    PS = 0,         FR = 1,         RS = 2;
    TString suffix[M]       = { "phase_space",  "fiducial",     "zz_4l"};
    TString histName[M]     = { "PhaseSpace",   "Fiducial",     "Selected"};


    const unsigned P = 3;
    unsigned                    ACC = 0,        EFF = 1,        AXE = 2;
    TString prefix[P]       = { "acceptance",   "efficiency",   "acc_x_eff"};
    unsigned numerator[P]   = { FR,             RS,             RS};
    unsigned denominator[P] = { PS,             FR,             PS};



    //
    //  INPUT FILES
    //

    TFile *inFile[M];

    cout << endl << endl;
    for (unsigned j = 0; j < M; j++)
    {
        TString inName  = "unscaled_" + suffix[j] + ".root";
        TString inPath  = inName; //HOME_PATH + "/Boosted/" + YEAR_STR + "/" + inName;
        inFile[j] = new TFile(inPath);

        cout << "Opened " << inPath << endl;
    }
    cout << endl;



    //
    //  GET HISTOGRAMS
    //
    
    // Get histogram names from one of the files
    vector<TString> hname;

    TDirectory *keyDir = inFile[PS]->GetDirectory("/" + selection[L4], kTRUE, "GetDirectory");
    TKey *histKey;
    TIter next(keyDir->GetListOfKeys());
    while ((histKey = (TKey*) next()))
    {
        TString hname_ = histKey->GetName();
        hname_.Resize(hname_.Length() - (1 + suffix[PS].Length()));     // truncate before suffix
        hname.push_back(hname_);
    }


    // Now get all of the histograms
    const unsigned H = hname.size();
    TH1* hist[M][N][H];

    cout << "Loading " << H * N << " histograms from each file..." << flush;
    for (unsigned j = 0; j < M; j++)    // file loop
    {
        for (unsigned i = 0; i < N; i++)    // directory (selection) loop
        {
            for (unsigned h = 0; h < H; h++)    // histogram loop
            {
                TH1 *hist_;

                inFile[j]->GetObject(selection[i] + "/" + hname[h] + "_" + suffix[j], hist_);
                hist_->SetDirectory(0);

                if (hist_->GetNbinsX() > 99)
                    hist_->Rebin(nRebin);

                hist[j][i][h] = hist_;
            }
        }
        inFile[j]->Close();
    }
    cout << "done!" << endl << endl;



    //
    //  CALCULATE RATIOS
    //
    //  &
    //
    //  WRITE FILES
    //


    cout << "Dividing histograms..." << endl;

    for (unsigned p = 0; p < P; p++)    // file loop
    {
        TString outName = prefix[p] + ".root";
        TFile *outFile = new TFile(outName, "RECREATE");

        for (unsigned i = 0; i < N; i++)    // directory (selection) loop
        {
            outFile->mkdir(selection[i]);
            outFile->cd(selection[i]);

            for (unsigned h = 0; h < H; h++)    // histogram loop
            {
                TH1* ratio = (TH1*) hist[numerator[p]][i][h]->Clone();
                ratio->Divide(hist[denominator[p]][i][h]);

                ratio->SetName(hname[h]);
                ratio->Write();
            }
        }
        outFile->Close();

        cout << "Wrote acceptance histograms to " << outName << endl;
    }
}
