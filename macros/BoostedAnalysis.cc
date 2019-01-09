// STL
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

// Cuts
#include "Cuts2017.hh"

using namespace std;



/*
**  BoostedAnalysis
**
**  Reads 4-lepton events from a "Selected" level tree.  Boosts objects into Z ("zz") rest frame
**  and calculates observables for differential distributions.
**
**  Also calculates event weights per source using scale factors from external histograms.
*/

void BoostedAnalysis(const TString suffix, bool doWeights = kFALSE, bool doRescale = kFALSE)
{

    //
    //  SAMPLE INFO
    //

    const bool isData = suffix.Contains(YEAR_STR);

    const unsigned N = 5;
    unsigned                    L4 = 0, M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = { "4l",   "4m",   "2m2e", "2e2m", "4e"    };
    TString selection2l[N]  = { "",     "mumu", "mumu", "ee",   "ee"    };
    unsigned chanIdx[N]     = { 5,      6,      7,      8,      9       };



    //
    //  OUTPUT FILE
    //

    TString prefix  = "boosted";
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
    Float_t             met,        weight,     rescale;
    UInt_t              channel,    nPV;


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
    Float_t             cos_theta_z1,           cos_theta_z2;
    Float_t             cos_zeta_z1,            cos_zeta_z2;
    Float_t             angle_z1leps,           angle_z2leps;
    Float_t             angle_z1l2_z2;


    // Weights
    Float_t             gen_weight,             pileup_weight;  //,          trigger_weight;
    Float_t             muonID_weight,          elecID_weight,          elecReco_weight;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("runNum",   &runNum);               tree[i]->Branch("evtNum",   &evtNum);
        tree[i]->Branch("lumiSec",  &lumiSec);              tree[i]->Branch("nPV",      &nPV);
        tree[i]->Branch("met",      &met);                  tree[i]->Branch("channel",  &channel);
        tree[i]->Branch("weight",   &weight);

        if (doRescale)  tree[i]->Branch("rescale",  &rescale);

        tree[i]->Branch("psi",              &psi);
        tree[i]->Branch("sin_phi",          &sin_phi);
        tree[i]->Branch("cos_theta_z1",     &cos_theta_z1);
        tree[i]->Branch("cos_theta_z2",     &cos_theta_z2);
        tree[i]->Branch("cos_zeta_z1",      &cos_zeta_z1);
        tree[i]->Branch("cos_zeta_z2",      &cos_zeta_z2);
        tree[i]->Branch("angle_z1leps",     &angle_z1leps);
        tree[i]->Branch("angle_z2leps",     &angle_z2leps);
        tree[i]->Branch("angle_z1l2_z2",    &angle_z1l2_z2);

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

        if (doWeights)
        {
            tree[i]->Branch("gen_weight",       &gen_weight);
            tree[i]->Branch("pileup_weight",    &pileup_weight);
//          tree[i]->Branch("trigger_weight",   &trigger_weight);
            tree[i]->Branch("muonID_weight",    &muonID_weight);
            tree[i]->Branch("elecID_weight",    &elecID_weight);
            tree[i]->Branch("elecReco_weight",  &elecReco_weight);
        }
    }



    //
    //  INPUT FILE
    //

    TString inName  = "selected_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl;



    //
    //  HISTOGRAMS
    //

    TH1D *hTotalEvents, *hSelectedEvents;

    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);
    hTotalEvents->Sumw2();

    if (doRescale)
        hSelectedEvents = new TH1D("SelectedEvents_" + suffix, "SelectedEvents", 10, 0.5, 10.5);
    else
        inFile->GetObject("SelectedEvents_" + suffix, hSelectedEvents);

    hSelectedEvents->SetDirectory(outFile);
    hSelectedEvents->Sumw2();



    //
    //  WEIGHT UTILS
    //

    TString sfPath = "../data/",    sfSuffix = "sf";

    TGraph  *gPileup;
    TH1     *hRescale[N];

    if (doWeights)
    {
        TString puName = sfPath + "pileup_" + sfSuffix + "_" + YEAR_STR + ".root";  
        TFile*  puFile = new TFile(puName, "OPEN");
        puFile->GetObject("pileup_sf", gPileup);
    }

    if (doRescale)
    {
        TString sfName = "../data/z1pt_rescale_" + YEAR_STR + ".root";
        TFile*  sfFile = new TFile(sfName, "OPEN");

        for (unsigned i = 1; i < N; i++)
        {
            sfFile->GetObject(selection2l[i] + "/hist_z1pt_" + selection2l[i], hRescale[i]);
            hRescale[i]->SetDirectory(0);
        }
    }



    //
    //  INPUT BRANCHES
    //

    for (unsigned i = 1; i < N; i++)
    {
        TTreeReader reader(selection[i] + "_" + suffix, inFile);

        TTreeReaderValue    <Int_t>                 runNum_     (reader,    "runNum");
        TTreeReaderValue    <Int_t>                 evtNum_     (reader,    "evtNum");
        TTreeReaderValue    <Int_t>                 lumiSec_    (reader,    "lumiSec");
        TTreeReaderValue    <UShort_t>              nPV_        (reader,    "nPV");
        TTreeReaderValue    <Float_t>               met_        (reader,    "met");
        TTreeReaderValue    <Float_t>               weight_     (reader,    "weight");
        TTreeReaderValue    <UInt_t>                channel_    (reader,    "channel");
        TTreeReaderValue    <TLorentzVector>        zzp4_       (reader,    "zzp4");
        TTreeReaderValue    <TLorentzVector>        z1p4_       (reader,    "z1p4");
        TTreeReaderValue    <Short_t>               z1pdg_      (reader,    "z1pdg");
        TTreeReaderValue    <TLorentzVector>        z2p4_       (reader,    "z2p4");
        TTreeReaderValue    <Short_t>               z2pdg_      (reader,    "z2pdg");
        TTreeReaderValue    <TLorentzVector>        l1p4_       (reader,    "l1p4");
        TTreeReaderValue    <Short_t>               l1pdg_      (reader,    "l1pdg");
        TTreeReaderValue    <UShort_t>              l1z_        (reader,    "l1z");
        TTreeReaderValue    <TLorentzVector>        l2p4_       (reader,    "l2p4");
        TTreeReaderValue    <Short_t>               l2pdg_      (reader,    "l2pdg");
        TTreeReaderValue    <UShort_t>              l2z_        (reader,    "l2z");
        TTreeReaderValue    <TLorentzVector>        l3p4_       (reader,    "l3p4");
        TTreeReaderValue    <Short_t>               l3pdg_      (reader,    "l3pdg");
        TTreeReaderValue    <UShort_t>              l3z_        (reader,    "l3z");
        TTreeReaderValue    <TLorentzVector>        l4p4_       (reader,    "l4p4");
        TTreeReaderValue    <Short_t>               l4pdg_      (reader,    "l4pdg");
        TTreeReaderValue    <UShort_t>              l4z_        (reader,    "l4z");






        ////
        ////
        ////    EVENT LOOP
        ////
        ////


        int nEvents = reader.GetEntries(kTRUE);

        cout << endl;
        cout << "Running over " << nEvents << " selected " + selection[i] + " events..." << flush;

        while (reader.Next())
        {

            //
            //  EVENT INFO
            //                

            // Quantities copied directly to output tree
            runNum  = *runNum_;         evtNum  = *evtNum_;         lumiSec = *lumiSec_;
            weight  = *weight_;         nPV     = *nPV_;            met     = *met_;
            channel = *channel_;
            zzp4    = *zzp4_;
            z1p4    = *z1p4_;           z1pdg   = *z1pdg_;
            z2p4    = *z2p4_;           z2pdg   = *z2pdg_;
            l1p4    = *l1p4_;           l1pdg   = *l1pdg_;          l1z     = *l1z_;
            l2p4    = *l2p4_;           l2pdg   = *l2pdg_;          l2z     = *l2z_;
            l3p4    = *l3p4_;           l3pdg   = *l3pdg_;          l3z     = *l3z_;
            l4p4    = *l4p4_;           l4pdg   = *l4pdg_;          l4z     = *l4z_;

            // Initialize weights
            gen_weight      = 0;        pileup_weight   = 0;        //trigger_weight  = 1;
            muonID_weight   = 1;        elecID_weight   = 1;        elecReco_weight = 1;



            //
            //  LEPTONS
            //

            vector<Lepton> leps(4);

            leps[0].p4  = l1p4;         leps[0].pdg = l1pdg;        leps[0].mother  = l1z;
            leps[1].p4  = l2p4;         leps[1].pdg = l2pdg;        leps[1].mother  = l2z;
            leps[2].p4  = l3p4;         leps[2].pdg = l3pdg;        leps[2].mother  = l3z;
            leps[3].p4  = l4p4;         leps[3].pdg = l4pdg;        leps[3].mother  = l4z;

            for (unsigned i = 0; i < leps.size(); i++)
                leps[i].q = -1 * copysign(1, leps[i].pdg);



            //
            //  BOOST
            //

            // "ZZ" rest frame
            TVector3 zz_boost = zzp4.BoostVector();
            for (unsigned i = 0; i < leps.size(); i++)
                leps[i].SetBoostedP4(zz_boost);

            sort(leps.begin(), leps.end(), DecreasingBoostedP);



            //
            //  PAIRS
            //

            LeptonPair z1, z2;
            MakePairsFromMother(leps, &z1, &z2);






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

            // "beta": angle between trailing pair 1 lepton and low-mass pair
            TVector3    z1_low = z1.BSecond().b_v3;
            angle_z1l2_z2 = z2.b_p4.Angle(z1_low);


            // Boosted lepton-pair angles
            TVector3    z1_boost = z1.p4.BoostVector(),     z2_boost = z2.p4.BoostVector();
            LeptonPair  b1_z1 = z1,     b1_z2 = z2,         b2_z1 = z1,     b2_z2 = z2;

            b1_z1.SetBoostedP4(z1_boost);                   b1_z2.SetBoostedP4(z1_boost);
            b2_z1.SetBoostedP4(z2_boost);                   b2_z2.SetBoostedP4(z2_boost);


            // "theta_zX": angle between positive pair X lepton and Y pair in X pair CM frame
            TVector3    u_b1_z2 = b1_z2.b_v3.Unit(),        u_b2_z1 = b2_z1.b_v3.Unit();
            TVector3    u_b1_z1_plus = b1_z1.Plus().b_v3.Unit();
            TVector3    u_b2_z2_plus = b2_z2.Plus().b_v3.Unit();

            cos_theta_z1 = u_b1_z2.Dot(u_b1_z1_plus);
            cos_theta_z2 = u_b2_z1.Dot(u_b2_z2_plus);


            // "zeta_zX": polarization angle for positive pair X lepton and X pair in Z CM frame
            TVector3    u_z1 = z1.b_v3.Unit(),              u_z2 = z2.b_v3.Unit();

            cos_zeta_z1 = u_z1.Dot(u_b1_z1_plus);
            cos_zeta_z2 = u_z2.Dot(u_b2_z2_plus);






            ////
            ////
            ////    WEIGHTS
            ////
            ////


            // "Gen" weight: whether event weight is positively or negatively weighted
            gen_weight = copysign(1, weight);

            // Pileup weight
            
            // Rescale

            if (doRescale)
            {
                if (!isData)
                    rescale = hRescale[i]->GetBinContent(hRescale[i]->FindBin(zzp4.Pt()));
                else
                    rescale = 1;

                hSelectedEvents->Fill(chanIdx[i], weight * rescale);
                hSelectedEvents->Fill(chanIdx[L4], weight * rescale);
            }



            //
            //  FILL TREE
            //

            b_ttp4  = leps[1].b_p4 + leps[2].b_p4 + leps[3].b_p4;
            b_z1p4  = z1.b_p4;          b_z2p4  = z2.b_p4;

            b_l1v3  = leps[0].b_v3;     b_l1pdg = leps[0].pdg;      b_l1z   = leps[0].mother;
            b_l2v3  = leps[1].b_v3;     b_l2pdg = leps[1].pdg;      b_l2z   = leps[1].mother;
            b_l3v3  = leps[2].b_v3;     b_l3pdg = leps[2].pdg;      b_l3z   = leps[2].mother;
            b_l4v3  = leps[3].b_v3;     b_l4pdg = leps[3].pdg;      b_l4z   = leps[3].mother;


            tree[i]->Fill();
            tree[L4]->Fill();

        } // END event loop

        cout << "done!" << endl;

    } // END tree loop


    cout << endl << endl;
    cout << "Processed all trees in " << inName << endl;



    //
    //  WRITE FILE
    //

    outFile->cd();

    for (unsigned i = 0; i < N; i++)
        tree[i]->Write();
    hTotalEvents->Write();
    hSelectedEvents->Write();

    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
