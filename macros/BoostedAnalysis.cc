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
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

// Cuts
#include "Cuts2018.hh"
//#include "Cuts2017.hh"
//#include "Cuts2016.hh"
//#include "Cuts2012.hh"

using namespace std;



/*
**  BoostedAnalysis
**
**  Reads 4-lepton events from a "Selected" level tree.  Boosts objects into Z ("zz") rest frame
**  and calculates observables for differential distributions.
**
**  Also calculates event weights per source using scale factors from external histograms.
*/

void BoostedAnalysis(const TString suffix, const bool isBkg = kFALSE)
{

    //
    //  SAMPLE INFO
    //

    const bool isData       = suffix.Contains(YEAR_STR);
    const bool isSignal     = suffix.EqualTo("zz_4l");
    const bool isDrellYan   = suffix.EqualTo("zjets_m-50");

    const unsigned N = 4;
    unsigned                    L4 = 0, M4 = 1, ME = 2, E4 = 3;     // Indices
    TString selection[N]    = { "4l",   "4m",   "2m2e", "4e"    };
    TString selection2l[N]  = { "",     "mumu", "mumu", "ee"    };
    unsigned chanIdx[N]     = { 5,      6,      7,      9       };

    TRandom3 *rng = new TRandom3(RNG_SEED);


/*
    //
    //  SF HISTOGRAMS
    //


    // Muon ID

    TH2     *mu_err;
    TString idName = "muon_id_smear_" + YEAR_STR + ".root";
    TString idPath = BLT_PATH + "/BLTAnalysis/data/" + idName;
    TFile   *idFile = TFile::Open(idPath);

    cout << "Opened " << idPath << endl;

    idFile->GetObject("ERROR", mu_err);
    mu_err->SetDirectory(0);

    idFile->Close();
    cout << "Closed " << idPath << endl;

    float MU_PT_MIN = mu_err->GetYaxis()->GetXmin();
    float MU_PT_MAX = mu_err->GetYaxis()->GetXmax();

    if (YEAR_STR.EqualTo("2012"))
    {
        MU_PT_MIN = mu_err->GetXaxis()->GetXmin();
        MU_PT_MAX = mu_err->GetXaxis()->GetXmax();
    }

    cout << "Limits: " << MU_PT_MIN << ", " << MU_PT_MAX << endl << endl;


    // Electron ID

    TH2     *el_err;
    idName = "electron_id_smear_" + YEAR_STR + ".root";
    idPath = BLT_PATH + "/BLTAnalysis/data/" + idName;
    idFile = TFile::Open(idPath);

    cout << "Opened " << idPath << endl;

    idFile->GetObject("ERROR", el_err);
    el_err->SetDirectory(0);

    idFile->Close();
    cout << "Closed " << idPath << endl;

    float EL_PT_MIN = el_err->GetYaxis()->GetXmin();
    float EL_PT_MAX = el_err->GetYaxis()->GetXmax();

    if (YEAR_STR.EqualTo("2012"))
    {
        EL_PT_MIN = el_err->GetXaxis()->GetXmin();
        EL_PT_MAX = el_err->GetXaxis()->GetXmax();
    }

    cout << "Limits: " << EL_PT_MIN << ", " << EL_PT_MAX << endl << endl;


    // Electron reco

    float   RECO_PT_MIN, RECO_PT_THRESH, RECO_PT_MAX;
    TH2     *reco_err, *reco_err_lowEt;

    if (!YEAR_STR.EqualTo("2012"))
    {
        idName = "electron_reco_smear_" + YEAR_STR + ".root";
        idPath = BLT_PATH + "/BLTAnalysis/data/" + idName;
        idFile = TFile::Open(idPath);

        cout << "Opened " << idPath << endl;

        idFile->GetObject("ERROR", reco_err);
        reco_err->SetDirectory(0);

        idFile->Close();
        cout << "Closed " << idPath << endl;

        if (YEAR_STR.EqualTo("2018"))
            RECO_PT_MIN = reco_err->GetYaxis()->GetXmin();
        else
        {
            idName = "electron_reco_smear_" + YEAR_STR + "_lowEt.root";
            idPath = BLT_PATH + "/BLTAnalysis/data/" + idName;
            idFile = TFile::Open(idPath);

            cout << "Opened " << idPath << endl;

            idFile->GetObject("ERROR", reco_err_lowEt);
            reco_err_lowEt->SetDirectory(0);

            idFile->Close();
            cout << "Closed " << idPath << endl;

            RECO_PT_MIN = reco_err_lowEt->GetYaxis()->GetXmin();
        }

        RECO_PT_THRESH = reco_err->GetYaxis()->GetXmin();
        RECO_PT_MAX = reco_err->GetYaxis()->GetXmax();

        cout << "Limits: " << RECO_PT_MIN << ", " << RECO_PT_THRESH << ", " << RECO_PT_MAX << endl << endl;
    }
*/


    //
    //  OUTPUT FILE
    //

    TString prefix  = "boosted";
    if (isBkg)
        prefix = prefix + "_bkg";
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
    UShort_t            nPV;
    Float_t             weight,     genWeight,  qtWeight,   puWeight,   ecalWeight;
    Float_t             trigWeight, idWeight;
    Float_t             wtMuonIDUp, wtMuonIDDn, wtElecIDUp, wtElecIDDn, wtElecRecoUp,wtElecRecoDn;
    UInt_t              channel;
    Bool_t              hasTauDecay;


    // Lab frame objects
    TLorentzVector      z1p4,       z2p4,       zzp4;
    UShort_t            z1pdg,      z2pdg;

    TLorentzVector      l1p4,       l2p4,       l3p4,       l4p4;
    Short_t             l1pdg,      l2pdg,      l3pdg,      l4pdg;
    Float_t             l1iso,      l2iso,      l3iso,      l4iso;
    UShort_t            l1z,        l2z,        l3z,        l4z;


    // Z rest frame objects
    TLorentzVector      b_z1p4,     b_z2p4,     b_ttp4;

    TVector3            b_l1v3,     b_l2v3,     b_l3v3,     b_l4v3;
    Short_t             b_l1pdg,    b_l2pdg,    b_l3pdg,    b_l4pdg;
    Float_t             b_l1iso,    b_l2iso,    b_l3iso,    b_l4iso;
    UShort_t            b_l1z,      b_l2z,      b_l3z,      b_l4z;


    // Observables
    Float_t             psi;
    Float_t             phi,                    sin_phi,                cos_phi;
    Float_t             cos_theta_z1,           cos_theta_z2;
    Float_t             angle_z1leps,           angle_z2leps,           angle_z1l2_z2;


    for (unsigned i = 0; i < N; i++)
    {
        tree[i]->Branch("runNum",       &runNum);       tree[i]->Branch("evtNum",       &evtNum);
        tree[i]->Branch("lumiSec",      &lumiSec);      tree[i]->Branch("nPV",          &nPV);
        tree[i]->Branch("weight",       &weight);       tree[i]->Branch("genWeight",    &genWeight);
        tree[i]->Branch("qtWeight",     &qtWeight);     tree[i]->Branch("puWeight",     &puWeight);
        tree[i]->Branch("ecalWeight",   &ecalWeight);   tree[i]->Branch("trigWeight",   &trigWeight);
        tree[i]->Branch("idWeight",     &idWeight);
        tree[i]->Branch("channel",      &channel);      tree[i]->Branch("hasTauDecay",  &hasTauDecay);
        tree[i]->Branch("wtMuonIDUp",   &wtMuonIDUp);   tree[i]->Branch("wtMuonIDDn",   &wtMuonIDDn);
        tree[i]->Branch("wtElecIDUp",   &wtElecIDUp);   tree[i]->Branch("wtElecIDDn",   &wtElecIDDn);
        tree[i]->Branch("wtElecRecoUp", &wtElecRecoUp); tree[i]->Branch("wtElecRecoDn", &wtElecRecoDn);

        tree[i]->Branch("psi",              &psi);      tree[i]->Branch("phi",          &phi);
        tree[i]->Branch("sin_phi",          &sin_phi);  tree[i]->Branch("cos_phi",      &cos_phi);
        tree[i]->Branch("cos_theta_z1",     &cos_theta_z1);
        tree[i]->Branch("cos_theta_z2",     &cos_theta_z2);
        tree[i]->Branch("angle_z1leps",     &angle_z1leps);
        tree[i]->Branch("angle_z2leps",     &angle_z2leps);
        tree[i]->Branch("angle_z1l2_z2",    &angle_z1l2_z2);

        tree[i]->Branch("b_z1p4",   &b_z1p4);           tree[i]->Branch("b_z2p4",   &b_z2p4);
        tree[i]->Branch("b_ttp4",   &b_ttp4);

        tree[i]->Branch("b_l1v3",   &b_l1v3);           tree[i]->Branch("b_l1pdg",  &b_l1pdg);
        tree[i]->Branch("b_l1z",    &b_l1z);            tree[i]->Branch("b_l1iso",  &b_l1iso);
        tree[i]->Branch("b_l2v3",   &b_l2v3);           tree[i]->Branch("b_l2pdg",  &b_l2pdg);
        tree[i]->Branch("b_l2z",    &b_l2z);            tree[i]->Branch("b_l2iso",  &b_l2iso);
        tree[i]->Branch("b_l3v3",   &b_l3v3);           tree[i]->Branch("b_l3pdg",  &b_l3pdg);
        tree[i]->Branch("b_l3z",    &b_l3z);            tree[i]->Branch("b_l3iso",  &b_l3iso);
        tree[i]->Branch("b_l4v3",   &b_l4v3);           tree[i]->Branch("b_l4pdg",  &b_l4pdg);
        tree[i]->Branch("b_l4z",    &b_l4z);            tree[i]->Branch("b_l4iso",  &b_l4iso);

        tree[i]->Branch("zzp4",     &zzp4);
        tree[i]->Branch("z1p4",     &z1p4);             tree[i]->Branch("z1pdg",    &z1pdg);
        tree[i]->Branch("z2p4",     &z2p4);             tree[i]->Branch("z2pdg",    &z2pdg);
        tree[i]->Branch("l1p4",     &l1p4);             tree[i]->Branch("l1pdg",    &l1pdg);
        tree[i]->Branch("l1z",      &l1z);              tree[i]->Branch("l1iso",    &l1iso);
        tree[i]->Branch("l2p4",     &l2p4);             tree[i]->Branch("l2pdg",    &l2pdg);
        tree[i]->Branch("l2z",      &l2z);              tree[i]->Branch("l2iso",    &l2iso);
        tree[i]->Branch("l3p4",     &l3p4);             tree[i]->Branch("l3pdg",    &l3pdg);
        tree[i]->Branch("l3z",      &l3z);              tree[i]->Branch("l3iso",    &l3iso);
        tree[i]->Branch("l4p4",     &l4p4);             tree[i]->Branch("l4pdg",    &l4pdg);
        tree[i]->Branch("l4z",      &l4z);              tree[i]->Branch("l4iso",    &l4iso);
    }



    //
    //  INPUT FILE
    //

    TString inName  = "selected_" + suffix + ".root";
    if (isBkg)
        inName = "background_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "_v1/" + inName;
    TFile   *inFile = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl;



    //
    //  HISTOGRAMS
    //

    TH1D *hTotalEvents, *hSelectedEvents;

    inFile->GetObject("TotalEvents_" + suffix, hTotalEvents);
    hTotalEvents->SetDirectory(outFile);
    hTotalEvents->Sumw2();

    inFile->GetObject("SelectedEvents_" + suffix, hSelectedEvents);
    hSelectedEvents->SetDirectory(outFile);
    hSelectedEvents->Sumw2();



    //
    //  INPUT BRANCHES
    //

    for (unsigned i = 1; i < N; i++)
    {
        TTreeReader reader(selection[i] + "_" + suffix, inFile);

        TTreeReaderValue    <Int_t>             runNum_         (reader,    "runNum");
        TTreeReaderValue    <Int_t>             evtNum_         (reader,    "evtNum");
        TTreeReaderValue    <Int_t>             lumiSec_        (reader,    "lumiSec");
        TTreeReaderValue    <UShort_t>          nPV_            (reader,    "nPV");
        TTreeReaderValue    <Float_t>           weight_         (reader,    "weight");
        TTreeReaderValue    <Float_t>           genWeight_      (reader,    "genWeight");
        TTreeReaderValue    <Float_t>           qtWeight_       (reader,    "qtWeight");
        TTreeReaderValue    <Float_t>           puWeight_       (reader,    "puWeight");
        TTreeReaderValue    <Float_t>           ecalWeight_     (reader,    "ecalWeight");
        TTreeReaderValue    <Float_t>           trigWeight_     (reader,    "trigWeight");
        TTreeReaderValue    <Float_t>           idWeight_       (reader,    "idWeight");
        TTreeReaderValue    <UShort_t>          channel_        (reader,    "channel");
        TTreeReaderValue    <Bool_t>            hasTauDecay_    (reader,    "hasTauDecay");

        TTreeReaderValue    <TLorentzVector>    zzp4_           (reader,    "zzp4");
        TTreeReaderValue    <TLorentzVector>    z1p4_           (reader,    "z1p4");
        TTreeReaderValue    <UShort_t>          z1pdg_          (reader,    "z1pdg");
        TTreeReaderValue    <TLorentzVector>    z2p4_           (reader,    "z2p4");
        TTreeReaderValue    <UShort_t>          z2pdg_          (reader,    "z2pdg");
        TTreeReaderValue    <TLorentzVector>    l1p4_           (reader,    "l1p4");
        TTreeReaderValue    <Short_t>           l1pdg_          (reader,    "l1pdg");
        TTreeReaderValue    <Float_t>           l1iso_          (reader,    "l1iso");
        TTreeReaderValue    <UShort_t>          l1z_            (reader,    "l1z");
        TTreeReaderValue    <TLorentzVector>    l2p4_           (reader,    "l2p4");
        TTreeReaderValue    <Short_t>           l2pdg_          (reader,    "l2pdg");
        TTreeReaderValue    <Float_t>           l2iso_          (reader,    "l2iso");
        TTreeReaderValue    <UShort_t>          l2z_            (reader,    "l2z");
        TTreeReaderValue    <TLorentzVector>    l3p4_           (reader,    "l3p4");
        TTreeReaderValue    <Short_t>           l3pdg_          (reader,    "l3pdg");
        TTreeReaderValue    <Float_t>           l3iso_          (reader,    "l3iso");
        TTreeReaderValue    <UShort_t>          l3z_            (reader,    "l3z");
        TTreeReaderValue    <TLorentzVector>    l4p4_           (reader,    "l4p4");
        TTreeReaderValue    <Short_t>           l4pdg_          (reader,    "l4pdg");
        TTreeReaderValue    <Float_t>           l4iso_          (reader,    "l4iso");
        TTreeReaderValue    <UShort_t>          l4z_            (reader,    "l4z");





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
            runNum      = *runNum_;     evtNum      = *evtNum_;     lumiSec     = *lumiSec_;
            nPV         = *nPV_;        weight      = *weight_;     genWeight   = *genWeight_;
            qtWeight    = *qtWeight_;   puWeight    = *puWeight_;   ecalWeight  = *ecalWeight_;
            trigWeight  = *trigWeight_; idWeight    = *idWeight_;
            channel     = *channel_;    hasTauDecay = *hasTauDecay_;

            zzp4    = *zzp4_;
            z1p4    = *z1p4_;       z1pdg   = *z1pdg_;
            z2p4    = *z2p4_;       z2pdg   = *z2pdg_;
            l1p4    = *l1p4_;       l1pdg   = *l1pdg_;      l1z = *l1z_;        l1pdg = *l1pdg_;
            l2p4    = *l2p4_;       l2pdg   = *l2pdg_;      l2z = *l2z_;        l2pdg = *l2pdg_;
            l3p4    = *l3p4_;       l3pdg   = *l3pdg_;      l3z = *l3z_;        l3pdg = *l3pdg_;
            l4p4    = *l4p4_;       l4pdg   = *l4pdg_;      l4z = *l4z_;        l4pdg = *l4pdg_;

            wtMuonIDUp  = weight;       wtElecIDUp  = weight;       wtElecRecoUp    = weight;
            wtMuonIDDn  = weight;       wtElecIDDn  = weight;       wtElecRecoDn    = weight;



            //
            //  LEPTONS
            //

            vector<Lepton> leps(4);

            leps[0].p4 = l1p4; leps[0].pdg = l1pdg; leps[0].mother  = l1z; leps[0].iso = l1iso;
            leps[1].p4 = l2p4; leps[1].pdg = l2pdg; leps[1].mother  = l2z; leps[1].iso = l2iso;
            leps[2].p4 = l3p4; leps[2].pdg = l3pdg; leps[2].mother  = l3z; leps[2].iso = l3iso;
            leps[3].p4 = l4p4; leps[3].pdg = l4pdg; leps[3].mother  = l4z; leps[3].iso = l4iso;

            for (unsigned i = 0; i < leps.size(); i++)
                leps[i].q = -1 * copysign(1, leps[i].pdg);



            //
            //  ID/RECO WEIGHT
            //

            for (unsigned i = 0; i < leps.size(); i++)
            {
                // Muons
                if (abs(leps[i].pdg) == 13)
                {
                    wtMuonIDUp *= (1 + 0.01);
                    wtMuonIDDn *= (1 - 0.01);
                }

                // Electrons
                else if (abs(leps[i].pdg) == 11)
                {
                    wtElecIDUp *= (1 + 0.01);
                    wtElecIDDn *= (1 - 0.01);

                    wtElecRecoUp *= (1 + 0.01);
                    wtElecRecoDn *= (1 - 0.01);
                }
            }



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



            // Blinding
            if      (isBkg)
                z2.BlindCharges(rng->Rndm());       // Randomly assign z2 charges
            else if (isData)
            {
                z1.BlindCharges(rng->Rndm());       // 50% chance to flip z1 charges
                z2.BlindCharges(rng->Rndm());       // 50% chance to flip z2 charges
            }






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

            // Scalar triple product
            psi = z2_plus.Dot(N_z1); 

            // Angle between decay planes
            TVector3    n_cross_n = n_z1.Cross(n_z2);
            sin_phi = n_cross_n.Dot(z1.b_v3.Unit());
            cos_phi = n_z1.Dot(n_z2);
            phi     = atan2(sin_phi, cos_phi);


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



            //
            //  FILL TREE
            //

            b_ttp4  = leps[1].b_p4 + leps[2].b_p4 + leps[3].b_p4;
            b_z1p4  = z1.b_p4;          b_z2p4  = z2.b_p4;

            b_l1v3 = leps[0].b_v3; b_l1pdg = leps[0].pdg; b_l1z = leps[0].mother; b_l1iso = leps[0].iso;
            b_l2v3 = leps[1].b_v3; b_l2pdg = leps[1].pdg; b_l2z = leps[1].mother; b_l2iso = leps[1].iso;
            b_l3v3 = leps[2].b_v3; b_l3pdg = leps[2].pdg; b_l3z = leps[2].mother; b_l3iso = leps[2].iso;
            b_l4v3 = leps[3].b_v3; b_l4pdg = leps[3].pdg; b_l4z = leps[3].mother; b_l4iso = leps[3].iso;


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
