// STL
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

// Cuts
#include "Cuts2017.hh"

using namespace std;



/*
**  GenTests
**
**  Some checks on the "signal" MC sample, and possibly others.  Reads from post-BLT analyzer 
**  (PhaseSpaceAnalyzer) ntuples.
*/

void GenTests()
{

    //
    //    OPTIONS
    //

    // fidOnly = kTRUE  => only events which fall into the fiducial region are written out
    //         = kFALSE => all input (phase space region) events are written out


    bool debug = kFALSE;
    int selectEvents = 1000;
    int printEvery = 30000;



    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                   L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = {"4l",   "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = {5,      6,      7,      8,      9};



    //
    //  INPUT FILE
    //

    TString inDir   = "gen_zz_4l";
    TString inName  = "genHardProc_zz_4l.root";
    TString inPath  = EOS_PATH + "/BLT/" + YEAR_STR + "/" + inDir + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << "Opened " << inPath << endl;



    //
    //  INPUT BRANCHES
    //

    TTreeReader reader("tree_zz_4l", inFile);

    TTreeReaderValue    <Float_t>               genWeight_  (reader,    "genWeight");
    TTreeReaderValue    <UShort_t>              channel_    (reader,    "decayChannel");
    TTreeReaderValue    <UShort_t>              nMuons_     (reader,    "nHardProcMuons");
    TTreeReaderValue    <UShort_t>              nElecs_     (reader,    "nHardProcElectrons");
    TTreeReaderValue    <UShort_t>              nLeps_      (reader,    "nHardProcLeptons");
    TTreeReaderArray    <TLorentzVector>        muonP4_     (reader,    "hardProcMuonP4");
    TTreeReaderValue    <vector<Short_t>>       muonQ_      (reader,    "hardProcMuonQ");
    TTreeReaderValue    <vector<UShort_t>>      muonZ_      (reader,    "hardProcMuonZIndex");
    TTreeReaderArray    <TLorentzVector>        elecP4_     (reader,    "hardProcElectronP4");
    TTreeReaderValue    <vector<Short_t>>       elecQ_      (reader,    "hardProcElectronQ");
    TTreeReaderValue    <vector<UShort_t>>      elecZ_      (reader,    "hardProcElectronZIndex");
//  TTreeReaderValue    <TLorentzVector>        lepsP4_     (reader,    "hardProcLeptonsP4");

    cout << "Loaded branches" << endl;



    //
    //  OUTPUT FILE
    //

    TString prefix  = "tests";
    TString suffix  = "phase_space"; 
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    TH1D *hMassSFOS, *hMassSFSS, *hMassDFOS, *hMassDFSS;

    hMassSFOS = new TH1D("MassSFOS", "MassSFOS", 100, 0, 100);
    hMassSFOS->Sumw2();
    hMassSFSS = new TH1D("MassSFSS", "MassSFSS", 100, 0, 100);
    hMassSFSS->Sumw2();
    hMassDFOS = new TH1D("MassDFOS", "MassDFOS", 100, 0, 100);
    hMassDFOS->Sumw2();
    hMassDFSS = new TH1D("MassDFSS", "MassDFSS", 100, 0, 100);
    hMassDFSS->Sumw2();

    TH1D *hDeltaR = new TH1D("DeltaR", "DeltaR", 100, 0, 1);
    hDeltaR->Sumw2();






    ////
    ////
    ////    EVENT LOOP
    ////
    ////


    int nEvents = reader.GetEntries(kTRUE);

    cout << endl;
    cout << "Running over " << nEvents << " total phase space events" << endl;
    if (debug)  cout << "Will select " << selectEvents << " events" << endl;
    cout << endl;

    long count = 0;     // track total number of selected events

    while (reader.Next())
    {

        //
        //  DEBUG PRINTOUT
        //

        const long currentEntry = reader.GetCurrentEntry();
        bool print = kFALSE;

        if (debug)
        {
            if (count >= selectEvents)
                break;

            if (currentEntry % printEvery == 0)
            {
                print = kTRUE;
                cout << endl << count << endl;
            }
        }
        else if (currentEntry % 10000 == 0)
        {
            cout << "Processed " << reader.GetCurrentEntry() << " of " << nEvents << " events";
            cout << endl;
        }



        //
        //  EVENT INFO
        //                

        // Quantities copied directly to output tree
        float       weight  = *genWeight_;
        unsigned    channel = *channel_;

        // Quantities used in analysis, but not written out
        unsigned    nMuons  = *nMuons_,         nElecs  = *nElecs_,         nLeps   = *nLeps_;



        //
        //  LEPTONS
        //

        vector<Lepton> muons, elecs;
//      TVector3 zz_boost = zzp4.BoostVector();

        // Muons
        for (unsigned i = 0; i < nMuons; i++)
        {
            Lepton      muon;

            muon.p4     = muonP4_.At(i);
            muon.q      = (*muonQ_)[i];
            muon.pdg    = -13 * muon.q;
            muon.mother = (*muonZ_)[i];

//          muon.SetBoostedP4(zz_boost);

            muons.push_back(muon);
        }
        sort(muons.begin(), muons.end(), DecreasingPt);


        // Electrons
        for (unsigned i = 0; i < nElecs; i++)
        {
            Lepton      elec;

            elec.p4     = elecP4_.At(i);
            elec.q      = (*elecQ_)[i];
            elec.pdg    = -11 * elec.q;
            elec.mother = (*elecZ_)[i];

//          elec.SetBoostedP4(zz_boost);

            elecs.push_back(elec);
        }
        sort(elecs.begin(), elecs.end(), DecreasingPt);


        // All leptons
        vector<Lepton> leps = muons;
        leps.insert(leps.end(), elecs.begin(), elecs.end());

        if (leps.size() != 4)
        {
            cout << "Event does not have four leptons" << endl;
            continue;
        }



        //
        //  FILL HISTOGRAMS
        //

        for (unsigned j = 1; j < nLeps;j++)
        {
            for (unsigned i = 0; i < j; i++)
            {
                TLorentzVector dilep = leps[i].p4 + leps[j].p4;

                if (abs(leps[i].pdg) == abs(leps[j].pdg))   // same flavor
                {
                    if (leps[i].q == leps[j].q)     // same sign
                        hMassSFSS->Fill(dilep.M(), weight);
                    else                            // opposite sign
                        hMassSFOS->Fill(dilep.M(), weight);
                }

                else                                        // different flavor
                {
                    if (leps[i].q == leps[j].q)     // same sign
                        hMassDFSS->Fill(dilep.M(), weight);
                    else                            // opposite sign
                        hMassDFOS->Fill(dilep.M(), weight);
                }

                hDeltaR->Fill(leps[i].p4.DeltaR(leps[j].p4), weight);
            }
        }


        count++;

    } // END event loop



    //
    //  WRITE FILE
    //

    outFile->cd();

    hMassSFOS->Write();
    hMassSFSS->Write();
    hMassDFOS->Write();
    hMassDFSS->Write();
    hDeltaR->Write();

    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
