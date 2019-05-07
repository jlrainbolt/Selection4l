// STL
#include <vector>
#include <iostream>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// Cuts
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2012.hh"

using namespace std;



void CalculateScales()
{


    //
    //  CONTAINERS
    //

    const unsigned  Q = 6,      P = 100;    // QCD, PDF
    unsigned        id_qcd[Q],  id_pdf[P];

    float   fid_4l_nom = 0,     gen_4l_nom = 0,     fid_2l_nom = 0,     gen_2l_nom = 0;
    float   fid_4l_qcd[Q],      gen_4l_qcd[Q],      fid_2l_qcd[Q],      gen_2l_qcd[Q];
    float   fid_4l_pdf[P],      gen_4l_pdf[P],      fid_2l_pdf[P],      gen_2l_pdf[P];

    for (unsigned i = 0; i < Q; i++)
    {
        fid_4l_qcd[i] = 0;
        gen_4l_qcd[i] = 0;
        fid_2l_qcd[i] = 0;
        gen_2l_qcd[i] = 0;
    }

    for (unsigned j = 0; j < P; j++)
    {
        fid_4l_pdf[j] = 0;
        gen_4l_pdf[j] = 0;
        fid_2l_pdf[j] = 0;
        gen_2l_pdf[j] = 0;
    }



    //
    //  FOUR-LEPTON
    //

    TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/gen_zz_4l_0.root";
    TFile   *zzFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    TTreeReader reader("tree_zz_4l", zzFile);

    TTreeReaderValue    <Bool_t>                isFiducial_     (reader,    "isFiducial");
    TTreeReaderValue    <vector<UShort_t>>      qcdID_          (reader,    "qcdID");
    TTreeReaderValue    <vector<Float_t>>       qcdWeight_      (reader,    "qcdWeight");
    TTreeReaderValue    <vector<UShort_t>>      pdfID_          (reader,    "pdfID");
    TTreeReaderValue    <vector<Float_t>>       pdfWeight_      (reader,    "pdfWeight");

    cout << "Loaded branches" << endl;


    // Event loop

    int nEvents = reader.GetEntries(kTRUE);

    cout << endl;
    cout << "Running over " << nEvents << " total ZZTo4L events" << endl;
    cout << endl;

    while (reader.Next())
    {
        if ((reader.GetCurrentEntry() % 100000 == 0) && (reader.GetCurrentEntry() > 0))
        {
            cout << "Processed " << reader.GetCurrentEntry() << " of " << nEvents << " events";
            cout << endl;
        }


        // Nominal weight

        gen_4l_nom += (*qcdWeight_)[0];

        if (*isFiducial_)
            fid_4l_nom += (*qcdWeight_)[0];


        // QCD scales

        unsigned i = 0;
        for (unsigned i_ = 1; i_ < (*qcdID_).size(); i_++)
        {
            // Skip 0.5, 2 combos   (https://twiki.cern.ch/twiki/bin/view/CMS/LHEReaderCMSSW)
            if ((i_ == 5) || (i_ == 7))
                continue;

            gen_4l_qcd[i] += (*qcdWeight_)[i_];

            if (*isFiducial_)
                fid_4l_qcd[i] += (*qcdWeight_)[i_];

            if (reader.GetCurrentEntry() == 0)
                id_qcd[i] = (*qcdID_)[i_];
            
            i++;
        }


        // PDF scales

        for (unsigned j = 0; j < P; j++)
        {
            unsigned j_ = j + 1;

            gen_4l_pdf[j] += (*pdfWeight_)[j_];

            if (*isFiducial_)
                fid_4l_pdf[j] += (*pdfWeight_)[j_];

            if (reader.GetCurrentEntry() == 0)
                id_pdf[j] = (*pdfID_)[j_];
        }
    }
    zzFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  DILEPTON
    //
/*
    // This doesn't even work...

    TChain  *dyTree = new TChain("tree_zjets_m-50");
    for (unsigned i = 0; i < N_DY; i++)
    {
        TString index = TString::Format("%i", i);
        inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/gen_zjets_m-50_" + index + ".root";
        dyTree->Add(inPath);
        cout << "Opened " << inPath << endl;
    }
    reader.SetTree(dyTree);
*/

    // Open one file

    inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/gen_zjets_m-50_0.root";
    TFile   *dyFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;

    reader.SetTree("tree_zjets_m-50", dyFile);

    reader.Restart();
    cout << "Loaded branches" << endl;


    // Event loop

    nEvents = reader.GetEntries(kTRUE);

    cout << endl;
    cout << "Running over " << nEvents << " total DYJetsToLL events" << endl;
    cout << endl;

    while (reader.Next())
    {
        if ((reader.GetCurrentEntry() % 100000 == 0) && (reader.GetCurrentEntry() > 0))
        {
            cout << "Processed " << reader.GetCurrentEntry() << " of " << nEvents << " events";
            cout << endl;
        }



        // Nominal weight
        gen_2l_nom += (*qcdWeight_)[0];

        if (*isFiducial_)
            fid_2l_nom += (*qcdWeight_)[0];


        // QCD scales

        unsigned i = 0;
        for (unsigned i_ = 1; i_ < (*qcdID_).size(); i_++)
        {
            // Skip 0.5, 2 combos   (https://twiki.cern.ch/twiki/bin/view/CMS/LHEReaderCMSSW)
            if ((i_ == 5) || (i_ == 7))
                continue;

            gen_2l_qcd[i] += (*qcdWeight_)[i_];

            if (*isFiducial_)
                fid_2l_qcd[i] += (*qcdWeight_)[i_];
            
            i++;
        }


        // PDF scales

        for (unsigned j = 0; j < P; j++)
        {
            unsigned j_ = j + 1;

            gen_2l_pdf[j] += (*pdfWeight_)[j_];

            if (*isFiducial_)
                fid_2l_pdf[j] += (*pdfWeight_)[j_];
        }
    }
    dyFile->Close();
    cout << "Closed " << inPath << endl << endl;



    //
    //  PRINTOUT
    //

    cout << setprecision(16) << endl;

    float   acc_4l_nom = fid_4l_nom / gen_4l_nom,   acc_2l_nom = fid_2l_nom / gen_2l_nom;
    float   acc_4l_qcd[Q],  acc_2l_qcd[Q],          acc_4l_pdf[P],  acc_2l_pdf[P];
    float   rat_nom = acc_2l_nom / acc_4l_nom,      rat_qcd[Q],     rat_pdf[P];

    cout << "%\t" << "Ratio" << "\t\t\t" << "2l" << "\t\t\t" << "4l" << endl << endl;
    cout << "x_nom = [" << endl;
    cout << "\t" << rat_nom << "\t" << acc_2l_nom << "\t" << acc_4l_nom << endl;
    cout << "\t" << "];" << endl;

    cout << "x_qcd = [" << endl;
    for (unsigned i = 0; i < Q; i++)
    {
        acc_4l_qcd[i] = fid_4l_qcd[i] / gen_4l_qcd[i];
        acc_2l_qcd[i] = fid_2l_qcd[i] / gen_2l_qcd[i];

        rat_qcd[i] = acc_2l_qcd[i] / acc_4l_qcd[i];

        cout << "\t" << rat_qcd[i] << "\t" << acc_2l_qcd[i] << "\t" << acc_4l_qcd[i] << endl;
    }
    cout << "\t" << "];" << endl << endl;

    cout << "x_pdf = [" << endl;
    for (unsigned j = 0; j < P; j++)
    {
        acc_4l_pdf[j] = fid_4l_pdf[j] / gen_4l_pdf[j];
        acc_2l_pdf[j] = fid_2l_pdf[j] / gen_2l_pdf[j];

        rat_pdf[j] = acc_2l_pdf[j] / acc_4l_pdf[j];

        cout << "\t" << rat_pdf[j] << "\t" << acc_2l_pdf[j] << "\t" << acc_4l_pdf[j] << endl;
    }
    cout << "\t" << "];" << endl << endl;

    cout << "id_qcd = [" << endl;
    for (unsigned i = 0; i < Q; i++)
        cout << "\t" << id_qcd[i] << endl;
    cout << "\t" << "];" << endl << endl;

    cout << "id_pdf = [" << endl;
    for (unsigned j = 0; j < P; j++)
        cout << "\t" << id_pdf[j] << endl;
    cout << "\t" << "];" << endl << endl;
}
