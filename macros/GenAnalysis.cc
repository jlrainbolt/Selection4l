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
**  GenAnalysis
**
**
**  All "selected" events are within the phase space region.  Events which also fall into the
**  fiducial region are flagged (isFiducial) based on cuts from header file (CutsXXXX.hh).  The
**  input argument determines whether all or only fiducial events are written out.
**
**  Reads from post-BLT analyzer (PhaseSpaceAnalyzer) ntuples and write out gen-level files that
**  are identical to the BoostedAnalysis output format ("boosted_").
**
**  Currently only implemented for signal (zz_4l) events...and only uses hard leptons.
*/

void GenAnalysis(const bool fidOnly)
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
    //  OUTPUT FILE
    //

    TString prefix  = "boosted";
    TString suffix  = fidOnly ? "fiducial" : "phase_space"; 
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);



    //
    //  OUTPUT BRANCHES
    //

    // Event info
    Int_t               runNum,     evtNum,     lumiSec;
    Float_t             weight;
    UInt_t              channel;
    Bool_t              isFiducial;

    // Lab frame objects
    TLorentzVector      z1p4,       z2p4,       zzp4;
    Short_t             z1pdg,      z2pdg;

    TLorentzVector      l1p4,       l2p4,       l3p4,       l4p4;
    Short_t             l1pdg,      l2pdg,      l3pdg,      l4pdg;
    UShort_t            l1z,        l2z,        l3z,        l4z;


    // Z rest frame objects
    TLorentzVector      b_z1p4,     b_z2p4,     b_ttp4;

    TVector3            b_l1v3,     b_l2v3,     b_l3v3,     b_l4v3;
    Short_t             b_l1pdg,    b_l2pdg,    b_l3pdg,    b_l4pdg;
    UShort_t            b_l1z,      b_l2z,      b_l3z,      b_l4z;

    // Observables
    Float_t             psi,                    sin_phi;
    Float_t             theta_z1,               theta_z2;
    Float_t             angle_z1leps,           angle_z2leps;
    Float_t             angle_z1l2_z2;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("runNum",   &runNum);               tree[i]->Branch("evtNum",   &evtNum);
        tree[i]->Branch("lumiSec",  &lumiSec);
        tree[i]->Branch("weight",   &weight);               tree[i]->Branch("channel",  &channel);
        tree[i]->Branch("isFiducial", &isFiducial);

        tree[i]->Branch("psi", &psi);                       tree[i]->Branch("sin_phi",  &sin_phi);
        tree[i]->Branch("theta_z1", &theta_z1);             tree[i]->Branch("theta_z2", &theta_z2);
        tree[i]->Branch("angle_z1leps", &angle_z1leps);
        tree[i]->Branch("angle_z2leps", &angle_z2leps);
        tree[i]->Branch("angle_z1l2_z2",  &angle_z1l2_z2);

        tree[i]->Branch("b_z1p4",   &b_z1p4);               tree[i]->Branch("b_z2p4",   &b_z2p4);
        tree[i]->Branch("b_ttp4",   &b_ttp4);

        tree[i]->Branch("b_l1v3",   &b_l1v3);               tree[i]->Branch("b_l1pdg",  &b_l1pdg);
        tree[i]->Branch("b_l1z",    &b_l1z);
        tree[i]->Branch("b_l2v3",   &b_l2v3);               tree[i]->Branch("b_l2pdg",  &b_l2pdg);
        tree[i]->Branch("b_l2z",    &b_l2z);
        tree[i]->Branch("b_l3v3",   &b_l3v3);               tree[i]->Branch("b_l3pdg",  &b_l3pdg);
        tree[i]->Branch("b_l3z",    &b_l3z);
        tree[i]->Branch("b_l4v3",   &b_l4v3);               tree[i]->Branch("b_l4pdg",  &b_l4pdg);
        tree[i]->Branch("b_l4z",    &b_l4z);

        tree[i]->Branch("zzp4",     &zzp4);
        tree[i]->Branch("z1p4",     &z1p4);                 tree[i]->Branch("z1pdg",    &z1pdg);
        tree[i]->Branch("z2p4",     &z2p4);                 tree[i]->Branch("z2pdg",    &z2pdg);
        tree[i]->Branch("l1p4",     &l1p4);                 tree[i]->Branch("l1pdg",    &l1pdg);
        tree[i]->Branch("l1z",      &l1z);
        tree[i]->Branch("l2p4",     &l2p4);                 tree[i]->Branch("l2pdg",    &l2pdg);
        tree[i]->Branch("l2z",      &l2z);
        tree[i]->Branch("l3p4",     &l3p4);                 tree[i]->Branch("l3pdg",    &l3pdg);
        tree[i]->Branch("l3z",      &l3z);
        tree[i]->Branch("l4p4",     &l4p4);                 tree[i]->Branch("l4pdg",    &l4pdg);
        tree[i]->Branch("l4z",      &l4z);
    }



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

    TTreeReaderValue    <Int_t>                 runNum_     (reader,    "runNumber");
    TTreeReaderValue    <Int_t>                 evtNum_     (reader,    "evtNumber.eventNumber");
    TTreeReaderValue    <Int_t>                 lumiSec_    (reader,    "lumiSection");
    TTreeReaderValue    <Float_t>               genWeight_  (reader,    "genWeight");
    TTreeReaderValue    <UInt_t>                channel_    (reader,    "decayChannel");
    TTreeReaderValue    <UShort_t>              nMuons_     (reader,    "nHardProcMuons");
    TTreeReaderValue    <UShort_t>              nElecs_     (reader,    "nHardProcElectrons");
    TTreeReaderArray    <TLorentzVector>        muonP4_     (reader,    "hardProcMuonP4");
    TTreeReaderValue    <vector<Short_t>>       muonQ_      (reader,    "hardProcMuonQ");
    TTreeReaderValue    <vector<UShort_t>>      muonZ_      (reader,    "hardProcMuonZIndex");
    TTreeReaderArray    <TLorentzVector>        elecP4_     (reader,    "hardProcElectronP4");
    TTreeReaderValue    <vector<Short_t>>       elecQ_      (reader,    "hardProcElectronQ");
    TTreeReaderValue    <vector<UShort_t>>      elecZ_      (reader,    "hardProcElectronZIndex");
    TTreeReaderValue    <TLorentzVector>        lepsP4_     (reader,    "hardProcLeptonsP4");

    cout << "Loaded branches" << endl;



    //
    //  HISTOGRAMS
    //

    TH1D *hPhaseSpaceEvents, *hFiducialEvents;

    inFile->GetObject("PhaseSpaceEvents_zz_4l", hPhaseSpaceEvents);
    hPhaseSpaceEvents->SetDirectory(outFile);
    hPhaseSpaceEvents->SetName("PhaseSpaceEvents_" + suffix);
    hPhaseSpaceEvents->Sumw2();
    // forgot to write bin 1 in BLT analyzer :(
    hPhaseSpaceEvents->SetBinContent(chanIdx[L4], hPhaseSpaceEvents->GetBinContent(chanIdx[M4])
                                                + hPhaseSpaceEvents->GetBinContent(chanIdx[ME])
                                                + hPhaseSpaceEvents->GetBinContent(chanIdx[EM])
                                                + hPhaseSpaceEvents->GetBinContent(chanIdx[E4]));

    hFiducialEvents = new TH1D("FiducialEvents_" + suffix, "FiducialEvents", 10, 0.5, 10.5);
    hFiducialEvents->SetDirectory(outFile);
    hFiducialEvents->Sumw2();






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
        runNum  = *runNum_;         evtNum  = *evtNum_;         lumiSec     = *lumiSec_;
        weight  = *genWeight_;      channel = *channel_;        zzp4    = *lepsP4_;

        // Quantities used in analysis, but not written out
        unsigned    nMuons  = *nMuons_,         nElecs  = *nElecs_;



        //
        //  LEPTONS
        //

        vector<Lepton> muons, elecs;
        TVector3 zz_boost = zzp4.BoostVector();

        // Muons
        for (unsigned i = 0; i < nMuons; i++)
        {
            Lepton      muon;

            muon.p4     = muonP4_.At(i);
            muon.q      = (*muonQ_)[i];
            muon.pdg    = -13 * muon.q;
            muon.mother = (*muonZ_)[i];

            muon.SetBoostedP4(zz_boost);

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

            elec.SetBoostedP4(zz_boost);

            elecs.push_back(elec);
        }
        sort(elecs.begin(), elecs.end(), DecreasingPt);


        // All leptons
        vector<Lepton> tmp_leps = muons;
        tmp_leps.insert(tmp_leps.end(), elecs.begin(), elecs.end());

        if (tmp_leps.size() != 4)
        {
            cout << "Event does not have four leptons" << endl;
            continue;
        }



        //
        //  CATEGORIZE EVENT
        //

        unsigned C;     // index for filling trees by channel
        LeptonPair z1, z2;
        

        if (print)
            cout << channel << endl;


        if      (channel == 6)                      // 4m
        {
            C = M4;
            MakePairsFromMother(muons, &z1, &z2);
        }
        else if (channel == 7)                      // 2m2e
        {
            C = ME;
            z1.SetMembers(muons[0], muons[1]);
            z2.SetMembers(elecs[0], elecs[1]);
        }
        else if (channel == 8)                      // 2e2m
        {
            C = EM;
            z1.SetMembers(elecs[0], elecs[1]);
            z2.SetMembers(muons[0], muons[1]);
        }
        else if (channel == 9)                      // 4e
        {
            C = E4;
            MakePairsFromMother(elecs, &z1, &z2);
        }
        else
        {
            cout << "Invalid decay channel" << endl;
            continue;
        }

        z1.SetMothers(1);
        z2.SetMothers(2);

        vector<Lepton> leps = z1.GetMembers(), z2_leps = z2.GetMembers();
        leps.insert(leps.end(), z2_leps.begin(), z2_leps.end());
        sort(leps.begin(), leps.end(), DecreasingPt);

        if (print)
        {
            cout << "Lab Pt:\t";
            for (unsigned i = 0; i < leps.size(); i++)
                cout << leps[i].p4.Pt() << "\t";
            cout << endl;
        }



        //
        //  FIDUCIAL CHECK
        //

        isFiducial = kTRUE;     // innocent until proven guilty?

        // Pt requirement
        if (leps[0].p4.Pt() < FID_PT1_MIN)
            isFiducial = kFALSE;

        if (leps[1].p4.Pt() < FID_PT2_MIN)
            isFiducial = kFALSE;

        if (leps[2].p4.Pt() < FID_PT_MIN)
            isFiducial = kFALSE;

        if (leps[3].p4.Pt() < FID_PT_MIN)
            isFiducial = kFALSE;


        // Eta requirement
        for (unsigned i = 0; i < leps.size(); i++)
        {
            if (fabs(leps[i].p4.Eta()) > FID_ETA_MAX)
                isFiducial = kFALSE;
        }


        hFiducialEvents->Fill(1, weight);

        if (isFiducial)
            hFiducialEvents->Fill(chanIdx[C], weight);
            hFiducialEvents->Fill(chanIdx[L4], weight);

        if (fidOnly && !isFiducial)
            continue;






        ////
        ////
        ////    OBSERVABLES
        ////
        ////


        TVector3    z1_plus = z1.Plus().b_v3,           z1_minus = z1.Minus().b_v3;
        TVector3    z2_plus = z2.Plus().b_v3,           z2_minus = z2.Minus().b_v3;

        // Normals to z1, z2 decay planes
        TVector3    N_z1 = z1_plus.Cross(z1_minus),     N_z2 = z2_plus.Cross(z2_minus);
        TVector3    n_z1 = N_z1.Unit(),                 n_z2 = N_z2.Unit();

        // Sclar triple product
        psi = z2_plus.Dot(N_z1);

        // Angle between decay planes
        TVector3    n_cross_n = n_z1.Cross(n_z2);
        sin_phi = n_cross_n.Dot(z1.b_v3.Unit());

        // Angles between paired leptons
        angle_z1leps = z1_plus.Angle(z1_minus);         angle_z2leps = z2_plus.Angle(z2_minus);

        // "theta": angle between trailing pair 1 lepton and low-mass pair
        TVector3    z1_low = z1.BSecond().b_v3;
        angle_z1l2_z2 = z2.b_p4.Angle(z1_low);


        // Boosted lepton-pair angles

        TVector3    z1_boost = z1.p4.BoostVector(),     z2_boost = z2.p4.BoostVector();
        LeptonPair  b1_z1 = z1,     b1_z2 = z2,         b2_z1 = z1,     b2_z2 = z2;

        b1_z1.SetBoostedP4(z1_boost);                   b1_z2.SetBoostedP4(z1_boost);
        b2_z1.SetBoostedP4(z2_boost);                   b2_z2.SetBoostedP4(z2_boost);

        // "theta_zX": angle between positive pair X lepton and Y pair in X pair CM frame
        theta_z1  = b1_z2.b_v3.Angle(b1_z1.Plus().b_v3);
        theta_z2  = b2_z1.b_v3.Angle(b2_z2.Plus().b_v3);



        //
        //  FILL TREE
        //

        // Sort all leptons by P in Z CM frame
        vector<Lepton> b_leps = leps;
        sort(b_leps.begin(), b_leps.end(), DecreasingBoostedP);

        z1p4    = z1.p4;                z1pdg   = z1.pdg;
        z2p4    = z2.p4;                z2pdg   = z2.pdg;

        l1p4    = leps[0].p4;           l1pdg   = leps[0].pdg;          l1z = leps[0].mother;
        l2p4    = leps[1].p4;           l2pdg   = leps[1].pdg;          l2z = leps[1].mother;
        l3p4    = leps[2].p4;           l3pdg   = leps[2].pdg;          l3z = leps[2].mother;
        l4p4    = leps[3].p4;           l4pdg   = leps[3].pdg;          l4z = leps[3].mother;


        b_ttp4  = b_leps[1].b_p4 + b_leps[2].b_p4 + b_leps[3].b_p4;
        b_z1p4  = z1.b_p4;              b_z2p4  = z2.b_p4;

        b_l1v3  = b_leps[0].b_v3;       b_l1pdg = b_leps[0].pdg;        b_l1z = b_leps[0].mother;
        b_l2v3  = b_leps[1].b_v3;       b_l2pdg = b_leps[1].pdg;        b_l2z = b_leps[1].mother;
        b_l3v3  = b_leps[2].b_v3;       b_l3pdg = b_leps[2].pdg;        b_l3z = b_leps[2].mother;
        b_l4v3  = b_leps[3].b_v3;       b_l4pdg = b_leps[3].pdg;        b_l4z = b_leps[3].mother;


        tree[C]->Fill();
        tree[L4]->Fill();

        count++;

    } // END event loop



    //
    //  PRINT RESULTS
    //

    int nFiducial = tree[L4]->GetEntries("isFiducial");

    cout << endl << endl;
    cout << "Done!" << endl;
    cout << "Found " << nFiducial << "/" << reader.GetCurrentEntry() << " fiducial events" << endl;

    if (!fidOnly)
    {
        cout << endl << endl;
        cout << "Yields for phase space region:" << endl;
        for (unsigned i = 0; i < N; i++)
        {
            if (i == L4)
                cout << endl;
            else
                cout << "\t\t";

            cout << selection[i] << ":\t";
            cout << tree[i]->GetEntries() << endl;
        }
    }

    cout << endl << endl;
    cout << "Yields for fiducial region:" << endl;
    for (unsigned i = 0; i < N; i++)
    {
        if (i == L4)
            cout << endl;
        else
            cout << "\t\t";

        cout << selection[i] << ":\t";
        cout << tree[i]->GetEntries("isFiducial") << endl;
    }
    cout << endl << endl;



    //
    //  WRITE FILE
    //

    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();

    hPhaseSpaceEvents->Write();
    hFiducialEvents->Write();

    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
