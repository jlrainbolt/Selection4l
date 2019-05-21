// STL
#include <vector>
#include <cmath>
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
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

// Cuts
//#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
#include "Cuts2012.hh"

using namespace std;



/*
**  MatchedAnalysis
**
**  Reads signal events from "Selected" level tree.  Performs lepton matching and writes out all
**  "Boosted" quantities.
*/

void MatchedAnalysis()
{

    //
    //    OPTIONS
    //

    bool debug = kFALSE;
    int analyzeEvents = 1000;
    int printEvery = 100;



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

    TString prefix  = "matched",    suffix = "zz_4l";
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    TTree *tree[N];
    for (unsigned i = 0; i < N; i++)
        tree[i] = new TTree(selection[i] + "_" + suffix, suffix);



    //
    //  OUTPUT BRANCHES
    //

    // Event info
    Bool_t              isMatched;
    Int_t               runNum,     evtNum,     lumiSec;
    UShort_t            nPV;
    Float_t             weight,     genWeight,  qtWeight,   puWeight,   ecalWeight;
    Float_t             trigWeight, idWeight,   recoWeight;
    UInt_t              channel;


    // Lab frame objects
    TLorentzVector      z1p4,           z2p4,           zzp4;
    TLorentzVector      gen_z1p4,       gen_z2p4,       gen_zzp4;
    Short_t             z1pdg,          z2pdg;

    TLorentzVector      l1p4,           l2p4,           l3p4,           l4p4;
    TLorentzVector      gen_l1p4,       gen_l2p4,       gen_l3p4,       gen_l4p4;
    Float_t             l1dr,           l2dr,           l3dr,           l4dr;
    Short_t             l1pdg,          l2pdg,          l3pdg,          l4pdg;
    UShort_t            l1z,            l2z,            l3z,            l4z;


    // Z rest frame objects
    TLorentzVector      b_z1p4,         b_z2p4,         b_ttp4;
    TLorentzVector      gen_b_z1p4,     gen_b_z2p4,     gen_b_ttp4;

    TVector3            b_l1v3,         b_l2v3,         b_l3v3,         b_l4v3;
    TVector3            gen_b_l1v3,     gen_b_l2v3,     gen_b_l3v3,     gen_b_l4v3;
    Short_t             b_l1pdg,        b_l2pdg,        b_l3pdg,        b_l4pdg;
    UShort_t            b_l1z,          b_l2z,          b_l3z,          b_l4z;

    // Observables
    Float_t     psi,            gen_psi,                sin_phi,        gen_sin_phi;
    Float_t     cos_theta_z1,   gen_cos_theta_z1,       cos_theta_z2,   gen_cos_theta_z2;
    Float_t     cos_zeta_z1,    gen_cos_zeta_z1,        cos_zeta_z2,    gen_cos_zeta_z2;
    Float_t     angle_z1leps,   gen_angle_z1leps,       angle_z2leps,   gen_angle_z2leps;
    Float_t     angle_z1l2_z2,  gen_angle_z1l2_z2;

    // Matched quantities


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("isMatched",    &isMatched);

        tree[i]->Branch("runNum",       &runNum);       tree[i]->Branch("evtNum",       &evtNum);
        tree[i]->Branch("lumiSec",      &lumiSec);      tree[i]->Branch("nPV",          &nPV);
        tree[i]->Branch("weight",       &weight);       tree[i]->Branch("genWeight",    &genWeight);
        tree[i]->Branch("qtWeight",     &qtWeight);     tree[i]->Branch("puWeight",     &puWeight);
        tree[i]->Branch("ecalWeight",   &ecalWeight);   tree[i]->Branch("trigWeight",   &trigWeight);
        tree[i]->Branch("idWeight",     &idWeight);     tree[i]->Branch("recoWeight",   &recoWeight);
        tree[i]->Branch("channel",      &channel);

        tree[i]->Branch("psi",      &psi);          tree[i]->Branch("gen_psi",      &gen_psi);                       
        tree[i]->Branch("sin_phi",  &sin_phi);      tree[i]->Branch("gen_sin_phi",  &gen_sin_phi);
        tree[i]->Branch("cos_theta_z1",     &cos_theta_z1);
        tree[i]->Branch("gen_cos_theta_z1", &gen_cos_theta_z1);             
        tree[i]->Branch("cos_theta_z2",     &cos_theta_z2);
        tree[i]->Branch("gen_cos_theta_z2", &gen_cos_theta_z2);
        tree[i]->Branch("cos_zeta_z1",      &cos_zeta_z1);
        tree[i]->Branch("gen_cos_zeta_z1",  &gen_cos_zeta_z1);             
        tree[i]->Branch("cos_zeta_z2",      &cos_zeta_z2);
        tree[i]->Branch("gen_cos_zeta_z2",  &gen_cos_zeta_z2);
        tree[i]->Branch("angle_z1leps",     &angle_z1leps);
        tree[i]->Branch("gen_angle_z1leps", &gen_angle_z1leps);
        tree[i]->Branch("angle_z2leps",     &angle_z2leps);
        tree[i]->Branch("gen_angle_z2leps", &gen_angle_z2leps);
        tree[i]->Branch("angle_z1l2_z2",    &angle_z1l2_z2);
        tree[i]->Branch("gen_angle_z1l2_z2",&gen_angle_z1l2_z2);

        tree[i]->Branch("b_z1p4",   &b_z1p4);       tree[i]->Branch("gen_b_z1p4",   &gen_b_z1p4);
        tree[i]->Branch("b_z2p4",   &b_z2p4);       tree[i]->Branch("gen_b_z2p4",   &gen_b_z2p4);
        tree[i]->Branch("b_ttp4",   &b_ttp4);       tree[i]->Branch("gen_b_ttp4",   &gen_b_ttp4);
        tree[i]->Branch("b_l1v3",   &b_l1v3);       tree[i]->Branch("gen_b_l1v3",   &gen_b_l1v3);
        tree[i]->Branch("b_l1pdg",  &b_l1pdg);      tree[i]->Branch("b_l1z",        &b_l1z);
        tree[i]->Branch("b_l2v3",   &b_l2v3);       tree[i]->Branch("gen_b_l2v3",   &gen_b_l2v3);
        tree[i]->Branch("b_l2pdg",  &b_l2pdg);      tree[i]->Branch("b_l2z",        &b_l2z);
        tree[i]->Branch("b_l3v3",   &b_l3v3);       tree[i]->Branch("gen_b_l3v3",   &gen_b_l3v3);
        tree[i]->Branch("b_l3pdg",  &b_l3pdg);      tree[i]->Branch("b_l3z",        &b_l3z);
        tree[i]->Branch("b_l4v3",   &b_l4v3);       tree[i]->Branch("gen_b_l4v3",   &gen_b_l4v3);
        tree[i]->Branch("b_l4pdg",  &b_l4pdg);      tree[i]->Branch("b_l4z",        &b_l4z);

        tree[i]->Branch("zzp4",     &zzp4);         tree[i]->Branch("gen_zzp4",     &gen_zzp4);
        tree[i]->Branch("z1p4",     &z1p4);         tree[i]->Branch("gen_z1p4",     &gen_z1p4);
        tree[i]->Branch("z1pdg",    &z1pdg);
        tree[i]->Branch("z2p4",     &z2p4);         tree[i]->Branch("gen_z2p4",     &gen_z2p4); 
        tree[i]->Branch("z2pdg",    &z2pdg);

        tree[i]->Branch("l1p4",         &l1p4);
        tree[i]->Branch("gen_l1p4",     &gen_l1p4);         tree[i]->Branch("l1dr",     &l1dr);
        tree[i]->Branch("l1pdg",        &l1pdg);            tree[i]->Branch("l1z",      &l1z);
        tree[i]->Branch("l2p4",         &l2p4);                                         
        tree[i]->Branch("gen_l2p4",     &gen_l2p4);         tree[i]->Branch("l2dr",     &l2dr);
        tree[i]->Branch("l2pdg",        &l2pdg);            tree[i]->Branch("l2z",      &l2z);
        tree[i]->Branch("l3p4",         &l3p4);                                         
        tree[i]->Branch("gen_l3p4",     &gen_l3p4);         tree[i]->Branch("l3dr",     &l3dr);
        tree[i]->Branch("l3pdg",        &l3pdg);            tree[i]->Branch("l3z",      &l3z);
        tree[i]->Branch("l4p4",         &l4p4);                                         
        tree[i]->Branch("gen_l4p4",     &gen_l4p4);         tree[i]->Branch("l4dr",     &l4dr);
        tree[i]->Branch("l4pdg",        &l4pdg);            tree[i]->Branch("l4z",      &l4z);
    }



    //
    //    INPUT FILE
    //

    TString inName  = "selected_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl;



    //
    //  HISTOGRAMS
    //

    TH1D *hTotalEvents, *hSelectedEvents, *hMatchedEvents;

    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);

    inFile->GetObject("SelectedEvents_" + suffix, hSelectedEvents);
    hSelectedEvents->SetDirectory(outFile);

    hMatchedEvents = new TH1D("MatchedEvents_" + suffix, "MatchedEvents", 10, 0.5, 10.5);
    hMatchedEvents->SetDirectory(outFile);



    //
    //  INPUT BRANCHES
    //

    for (unsigned i = 1; i < N; i++)
    {
        TTreeReader reader(selection[i] + "_" + suffix, inFile);

        TTreeReaderValue    <Int_t>                 runNum_         (reader,    "runNum");
        TTreeReaderValue    <Int_t>                 evtNum_         (reader,    "evtNum");
        TTreeReaderValue    <Int_t>                 lumiSec_        (reader,    "lumiSec");
        TTreeReaderValue    <UShort_t>              nPV_            (reader,    "nPV");
        TTreeReaderValue    <Float_t>               weight_         (reader,    "weight");
        TTreeReaderValue    <Float_t>               genWeight_      (reader,    "genWeight");
        TTreeReaderValue    <Float_t>               qtWeight_       (reader,    "qtWeight");
        TTreeReaderValue    <Float_t>               puWeight_       (reader,    "puWeight");
        TTreeReaderValue    <Float_t>               ecalWeight_     (reader,    "ecalWeight");
        TTreeReaderValue    <Float_t>               trigWeight_     (reader,    "trigWeight");
        TTreeReaderValue    <Float_t>               idWeight_       (reader,    "idWeight");
        TTreeReaderValue    <Float_t>               recoWeight_     (reader,    "recoWeight");
        TTreeReaderValue    <UInt_t>                channel_        (reader,    "channel");
        TTreeReaderValue    <TLorentzVector>        zzp4_           (reader,    "zzp4");
        TTreeReaderValue    <TLorentzVector>        z1p4_           (reader,    "z1p4");
        TTreeReaderValue    <Short_t>               z1pdg_          (reader,    "z1pdg");
        TTreeReaderValue    <TLorentzVector>        z2p4_           (reader,    "z2p4");
        TTreeReaderValue    <Short_t>               z2pdg_          (reader,    "z2pdg");
        TTreeReaderValue    <TLorentzVector>        l1p4_           (reader,    "l1p4");
        TTreeReaderValue    <Short_t>               l1pdg_          (reader,    "l1pdg");
        TTreeReaderValue    <UShort_t>              l1z_            (reader,    "l1z");
        TTreeReaderValue    <TLorentzVector>        l2p4_           (reader,    "l2p4");
        TTreeReaderValue    <Short_t>               l2pdg_          (reader,    "l2pdg");
        TTreeReaderValue    <UShort_t>              l2z_            (reader,    "l2z");
        TTreeReaderValue    <TLorentzVector>        l3p4_           (reader,    "l3p4");
        TTreeReaderValue    <Short_t>               l3pdg_          (reader,    "l3pdg");
        TTreeReaderValue    <UShort_t>              l3z_            (reader,    "l3z");
        TTreeReaderValue    <TLorentzVector>        l4p4_           (reader,    "l4p4");
        TTreeReaderValue    <Short_t>               l4pdg_          (reader,    "l4pdg");
        TTreeReaderValue    <UShort_t>              l4z_            (reader,    "l4z");
                                                                
        TTreeReaderValue    <UShort_t>          nGenMuons_      (reader,    "nFinalStateMuons");
        TTreeReaderValue    <UShort_t>          nGenElecs_      (reader,    "nFinalStateElectrons");
        TTreeReaderArray    <TLorentzVector>    genMuonP4_      (reader,    "finalStateMuonP4");
        TTreeReaderValue    <vector<Short_t>>   genMuonQ_       (reader,    "finalStateMuonQ");
        TTreeReaderArray    <TLorentzVector>    genElecP4_      (reader,    "finalStateElectronP4");
        TTreeReaderValue    <vector<Short_t>>   genElecQ_       (reader,    "finalStateElectronQ");
        TTreeReaderValue    <TLorentzVector>    genLepsP4_      (reader,    "finalStateLeptonsP4");
        TTreeReaderValue    <TLorentzVector>    u_l1p4_         (reader,    "uncorr_l1p4");
        TTreeReaderValue    <TLorentzVector>    u_l2p4_         (reader,    "uncorr_l2p4");
        TTreeReaderValue    <TLorentzVector>    u_l3p4_         (reader,    "uncorr_l3p4");
        TTreeReaderValue    <TLorentzVector>    u_l4p4_         (reader,    "uncorr_l4p4");






        ////
        ////
        ////    EVENT LOOP
        ////
        ////


        int nEvents = reader.GetEntries(kTRUE);

        cout << endl;
        cout << "Running over " << nEvents << " selected " + selection[i] + " events..." << flush;
        if (debug)  cout << "Will analyze " << analyzeEvents << " events" << endl;
        cout << endl;

        long count = 0;     // track total number of analyzed events

        while (reader.Next())
        {

            //
            //  DEBUG PRINTOUT
            //

            const long currentEntry = reader.GetCurrentEntry();
            bool print = kFALSE;

            if (debug)
            {
                if (count >= analyzeEvents)
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

            // Innocent until proven guilty
            isMatched   = kTRUE;

            // Quantities copied directly to output tree
            runNum      = *runNum_;     evtNum      = *evtNum_;     lumiSec     = *lumiSec_;
            nPV         = *nPV_;        weight      = *weight_;     genWeight   = *genWeight_;
            qtWeight    = *qtWeight_;   puWeight    = *puWeight_;   ecalWeight  = *ecalWeight_;
            trigWeight  = *trigWeight_; idWeight    = *idWeight_;   recoWeight  = *recoWeight_;
            channel     = *channel_;

            zzp4    = *zzp4_;           gen_zzp4    = *genLepsP4_;
            z1p4    = *z1p4_;           z1pdg   = *z1pdg_;
            z2p4    = *z2p4_;           z2pdg   = *z2pdg_;
            l1p4    = *l1p4_;           l1pdg   = *l1pdg_;          l1z     = *l1z_;
            l2p4    = *l2p4_;           l2pdg   = *l2pdg_;          l2z     = *l2z_;
            l3p4    = *l3p4_;           l3pdg   = *l3pdg_;          l3z     = *l3z_;
            l4p4    = *l4p4_;           l4pdg   = *l4pdg_;          l4z     = *l4z_;

            // Quantities used in analysis, but not written out
            unsigned    nGenMuons   = *nGenMuons_,          nGenElecs   = *nGenElecs_;

            // Reset gen quantities in case event is fake
            gen_b_z1p4.Delete();        gen_b_z2p4.Delete();        gen_b_ttp4.Delete();
            gen_b_l1v3.Delete();        gen_b_l2v3.Delete();
            gen_b_l3v3.Delete();        gen_b_l4v3.Delete();
            gen_psi = NAN;              gen_sin_phi = NAN;
            gen_cos_theta_z1 = NAN;     gen_cos_theta_z2 = NAN;
            gen_cos_zeta_z1 = NAN;      gen_cos_zeta_z2 = NAN;
            gen_angle_z1leps = NAN;     gen_angle_z2leps = NAN;     gen_angle_z1l2_z2 = NAN;



            //
            //  RECO LEPTONS
            //

            vector<Lepton> leps(4);

            leps[0].p4  = l1p4;         leps[0].pdg = l1pdg;        leps[0].mother  = l1z;
            leps[1].p4  = l2p4;         leps[1].pdg = l2pdg;        leps[1].mother  = l2z;
            leps[2].p4  = l3p4;         leps[2].pdg = l3pdg;        leps[2].mother  = l3z;
            leps[3].p4  = l4p4;         leps[3].pdg = l4pdg;        leps[3].mother  = l4z;

            // Fill uncorrected p4
            leps[0].u_p4     = *u_l1p4_;                leps[1].u_p4     = *u_l2p4_;
            leps[2].u_p4     = *u_l3p4_;                leps[3].u_p4     = *u_l4p4_;

            for (unsigned i = 0; i < leps.size(); i++)
                leps[i].q = -1 * copysign(1, leps[i].pdg);



            //
            //  PRESELECTION
            //

            // Make sure there are enough lepons available for matching
            // (Remember: there *can* be extras!

            hMatchedEvents->Fill(1, weight);

            if (nGenMuons + nGenElecs != 4)
            {
                if (print)
                    cout << "Wrong number of gen leptons" << endl;

                isMatched = kFALSE;

                tree[i]->Fill();
                tree[L4]->Fill();

                continue;
            }



            //
            //  GEN LEPTONS
            //

            vector<Lepton> gen_muons, gen_elecs;
            TVector3 gen_zz_boost = gen_zzp4.BoostVector();

            for (unsigned i = 0; i < nGenMuons; i++)
            {
                Lepton muon;

                muon.p4     = genMuonP4_.At(i);
                muon.q      = (*genMuonQ_)[i];
                muon.pdg    = -13 * muon.q;

                muon.SetBoostedP4(gen_zz_boost);

                gen_muons.push_back(muon);
            }

            for (unsigned i = 0; i < nGenElecs; i++)
            {
                Lepton elec;

                elec.p4     = genElecP4_.At(i);
                elec.q      = (*genElecQ_)[i];
                elec.pdg    = -11 * elec.q;

                elec.SetBoostedP4(gen_zz_boost);

                gen_elecs.push_back(elec);
            }

            vector<Lepton> gen_leps = gen_muons;
            gen_leps.insert(gen_leps.end(), gen_elecs.begin(), gen_elecs.end());






            ////
            ////
            ////    MATCHING
            ////
            ////


            bool failedDRCut = kFALSE;

            for (unsigned i = 0; i < leps.size(); i++)  // loop over reco muons
            {
                // Find DeltaR between reco lep and each gen lep
                vector<Float_t> deltaR(gen_leps.size());

                for (unsigned j = 0; j < gen_leps.size(); j++)
                {
                    deltaR[j] = leps[i].u_p4.DeltaR(gen_leps[j].p4);

                    if (print)
                        cout << deltaR[j] << "\t";
                }
                if (print)  cout << endl;


                // Find index of minimum DeltaR
                unsigned m = min_element(deltaR.begin(), deltaR.end()) - deltaR.begin();
                if (deltaR[m] > MATCH_DR_MAX)
                {
                    failedDRCut = kTRUE;
                    break;
                }

                leps[i].SetMatch(gen_leps[m]);

                if (print)
                {
                    cout << "Lepton " << i << ":\t" << leps[i].dr << "\t";
                    cout << leps[i].pdg << ", " << gen_leps[m].pdg << endl;
                }
                
                gen_leps.erase(gen_leps.begin() + m);

            }
            if (failedDRCut)
            {
                if (print)
                    cout << "Could not find a match!" << endl;

                isMatched = kFALSE;

                tree[i]->Fill();
                tree[L4]->Fill();

                continue;
            }


            // Now we for sure have a match!
            hMatchedEvents->Fill(chanIdx[i], weight);
            hMatchedEvents->Fill(chanIdx[L4], weight);



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



            //
            //  FILL TREE
            //

            gen_z1p4    = z1.m_p4;                  gen_z2p4    = z2.m_p4;
            gen_l1p4    = leps[0].m_p4;             l1dr        = leps[0].dr;
            gen_l2p4    = leps[1].m_p4;             l2dr        = leps[1].dr;
            gen_l3p4    = leps[2].m_p4;             l3dr        = leps[2].dr;
            gen_l4p4    = leps[3].m_p4;             l4dr        = leps[3].dr;






            ////
            ////
            ////    OBSERVABLES
            ////
            ////


            TVector3    z1_plus = z1.Plus().b_v3,           m_z1_plus = z1.Plus().m_b_v3;
            TVector3    z1_minus = z1.Minus().b_v3,         m_z1_minus = z1.Minus().m_b_v3;
            TVector3    z2_plus = z2.Plus().b_v3,           m_z2_plus = z2.Plus().m_b_v3;
            TVector3    z2_minus = z2.Minus().b_v3,         m_z2_minus = z2.Minus().m_b_v3;

            // Normals to z1, z2 decay planes
            TVector3    N_z1 = z1_plus.Cross(z1_minus),     m_N_z1 = m_z1_plus.Cross(m_z1_minus);
            TVector3    N_z2 = z2_plus.Cross(z2_minus),     m_N_z2 = m_z2_plus.Cross(m_z2_minus);
            TVector3    n_z1 = N_z1.Unit(),                 m_n_z1 = m_N_z1.Unit();
            TVector3    n_z2 = N_z2.Unit(),                 m_n_z2 = m_N_z2.Unit();

            // Sclar triple product
            psi = z2_plus.Dot(N_z1);                        gen_psi = m_z2_plus.Dot(m_N_z1);

            // Angle between decay planes
            TVector3    n_cross_n = n_z1.Cross(n_z2),       m_n_cross_n = m_n_z1.Cross(m_n_z2);
            sin_phi     = n_cross_n.Dot(z1.b_v3.Unit());
            gen_sin_phi = m_n_cross_n.Dot(z1.m_b_v3.Unit());

            // Angles between paired leptons
            angle_z1leps = z1_plus.Angle(z1_minus); gen_angle_z1leps = m_z1_plus.Angle(m_z1_minus);
            angle_z2leps = z2_plus.Angle(z2_minus); gen_angle_z2leps = m_z2_plus.Angle(m_z2_minus);

            // "theta": angle between trailing pair 1 lepton and low-mass pair
            TVector3    z1_low = z1.BSecond().b_v3, m_z1_low = z1.BSecond().m_b_v3;
            angle_z1l2_z2 = z2.b_p4.Angle(z1_low);  gen_angle_z1l2_z2 = z2.m_b_p4.Angle(m_z1_low);


            // Boosted lepton-pair angles
            TVector3    z1_boost = z1.p4.BoostVector(),     m_z1_boost = z1.m_p4.BoostVector();
            TVector3    z2_boost = z2.p4.BoostVector(),     m_z2_boost = z2.m_p4.BoostVector();
            LeptonPair  b1_z1 = z1,     b1_z2 = z2,         b2_z1 = z1,     b2_z2 = z2;

            b1_z1.SetBoostedP4(z1_boost, m_z1_boost);   b1_z2.SetBoostedP4(z1_boost, m_z1_boost);
            b2_z1.SetBoostedP4(z2_boost, m_z2_boost);   b2_z2.SetBoostedP4(z2_boost, m_z2_boost);


            // "theta_zX": angle between positive pair X lepton and Y pair in X pair CM frame
            TVector3    u_b1_z2 = b1_z2.b_v3.Unit(),        u_b2_z1 = b2_z1.b_v3.Unit();
            TVector3    u_b1_z1_plus = b1_z1.Plus().b_v3.Unit();
            TVector3    u_b2_z2_plus = b2_z2.Plus().b_v3.Unit();
            cos_theta_z1 = u_b1_z2.Dot(u_b1_z1_plus);
            cos_theta_z2 = u_b2_z1.Dot(u_b2_z2_plus);

            TVector3    m_u_b1_z2 = b1_z2.m_b_v3.Unit(),    m_u_b2_z1 = b2_z1.m_b_v3.Unit();
            TVector3    m_u_b1_z1_plus = b1_z1.Plus().m_b_v3.Unit();
            TVector3    m_u_b2_z2_plus = b2_z2.Plus().m_b_v3.Unit();
            gen_cos_theta_z1 = m_u_b1_z2.Dot(m_u_b1_z1_plus);
            gen_cos_theta_z2 = m_u_b2_z1.Dot(m_u_b2_z2_plus);


            // "zeta_zX": polarization angle for positive pair X lepton and X pair in Z CM frame
            TVector3    u_z1 = z1.b_v3.Unit(),              u_z2 = z2.b_v3.Unit();
            cos_zeta_z1 = u_z1.Dot(u_b1_z1_plus);
            cos_zeta_z2 = u_z2.Dot(u_b2_z2_plus);

            TVector3    m_u_z1 = z1.m_b_v3.Unit(),          m_u_z2 = z2.m_b_v3.Unit();
            gen_cos_zeta_z1 = m_u_z1.Dot(m_u_b1_z1_plus);
            gen_cos_zeta_z2 = m_u_z2.Dot(m_u_b2_z2_plus);



            //
            //  FILL TREE
            //

            b_ttp4      = leps[1].b_p4 + leps[2].b_p4 + leps[3].b_p4;
            gen_b_ttp4  = leps[1].m_b_p4 + leps[2].m_b_p4 + leps[3].m_b_p4;
            b_z1p4  = z1.b_p4;          gen_b_z1p4  = z1.m_b_p4;
            b_z2p4  = z2.b_p4;          gen_b_z2p4  = z2.m_b_p4;

            b_l1v3  = leps[0].b_v3;     b_l1pdg = leps[0].pdg;      b_l1z   = leps[0].mother;
            b_l2v3  = leps[1].b_v3;     b_l2pdg = leps[1].pdg;      b_l2z   = leps[1].mother;
            b_l3v3  = leps[2].b_v3;     b_l3pdg = leps[2].pdg;      b_l3z   = leps[2].mother;
            b_l4v3  = leps[3].b_v3;     b_l4pdg = leps[3].pdg;      b_l4z   = leps[3].mother;

            gen_b_l1v3  = leps[0].m_b_v3;       gen_b_l2v3  = leps[1].m_b_v3;
            gen_b_l3v3  = leps[2].m_b_v3;       gen_b_l4v3  = leps[3].m_b_v3;


            tree[i]->Fill();
            tree[L4]->Fill();

            count++;

        } // END event loop

            cout << "done!" << endl;

            cout << "Matched " << tree[i]->GetEntries("isMatched") << "/";
            cout << reader.GetCurrentEntry() << " events" << endl;
            cout << endl << endl;

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
    hMatchedEvents->Write();

    outFile->Purge();
    outFile->Close();
    inFile->Close();

    cout << "Wrote trees to " << outName << endl << endl << endl;
}
