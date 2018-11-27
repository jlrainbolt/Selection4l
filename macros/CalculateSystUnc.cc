// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"

// Cuts
#include "Cuts2017.hh"

using namespace std;


/*
**  CalculateSystUnc
**
**  Calculates acceptance * efficiency variance for input systematic uncertainty source
*/

void CalculateSystUnc(const TString systematics, const TString suffix4l, const TString suffix2l)
{

    //
    //  SAMPLE INFO
    //

    const unsigned N = 4;   // Channel indices
    unsigned                L4 = 0, M4 = 1, ME = 2, E4 = 4;
    TString selection[N] = {"4l", "4m", "2m2e", "4e"};

    float axe4l[N][N_HISTS], axe2l[N_HISTS], axeTotal[N][N_HISTS], fracDiff[N][N_HISTS];



    //
    //  READ IN
    //

    TString inPath = "output";

    // Four-lepton
    TString inName4l = systematics + "_" + suffix4l + ".txt";
    ifstream inFile4l;
    inFile4l.open(inPath + "/" + inName4l);
    cout << "Opened " << inName4l << endl;

    // Dilepton
    TString inName2l = systematics + "_" + suffix2l + ".txt";
    ifstream inFile2l;
    inFile2l.open(inPath + "/" + inName2l);
    cout << "Opened " << inName2l << endl;

    // Loop over smear histogram iterations
    string line4l, line2l;
    for (unsigned j = 0; j < N_HISTS; j++)
    {
        unsigned idx4l, idx2l;

        getline(inFile4l, line4l);
        stringstream stream4l(line4l);
        stream4l >> idx4l;
        for (unsigned i = 0; i < N; i++)
            stream4l >> axe4l[i][idx4l];

        getline(inFile2l, line2l);
        stringstream stream2l(line2l);
        stream2l >> idx2l;
        stream2l >> axe2l[idx2l];
    }

    inFile4l.close();
    inFile2l.close();



    //
    //  CALCULATE TOTAL
    //

    // Without smearing
    float axeTrue[N];
    for (unsigned i = 0; i < N; i++)
        axeTrue[i] = AXE_2L[0] / AXE_4L[i];
    cout << "\t" << AXE_2L[0] << "\t" << AXE_4L[0] << "\t" << axeTrue[0] << endl; 

    // With smearing
    float variance[N] = {0, 0, 0, 0};
    for (unsigned j = 0; j < N_HISTS; j++)
    {
        for (unsigned i = 0; i < N; i++)
        {
            axeTotal[i][j] = axe2l[j] / axe4l[i][j];
            fracDiff[i][j] = fabs(axeTotal[i][j] - axeTrue[i]) / axeTrue[i];
//          variance[i] += pow(fracDiff[i][j], 2);
            variance[i] += pow(axeTotal[i][j] - axeTrue[i], 2);
        }

        cout << j << "\t" << axe2l[j] << "\t" << axe4l[0][j] << "\t" << axeTotal[0][j] << "\t";
        cout << fracDiff[0][j] << endl;
    }
    
    float fracUnc[N];
    for (unsigned i = 0; i < N; i++)
    {
        variance[i] /= (float) N_HISTS;
        fracUnc[i] = sqrt(variance[i]) / axeTrue[i];
    }
    cout << endl << "Variance:\t" << variance[0] << endl;




    //
    //  PRINT
    //

    cout << endl << endl;
    cout << "\t\t" << "Fractional uncertainty" << endl;
    for (unsigned i = 0; i < N; i++)
        cout << selection[i] << "\t" << fracUnc[i] << endl;



/*
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
*/
}
