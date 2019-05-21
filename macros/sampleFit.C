//
// A sample fit for measuring lepton efficiencies.
//
// Michael Schmitt     Northwestern University      May 20, 2019
//
#include <TF1.h>
#include <TRandom2.h>


TH1D *tnpPass, *tnpFail;
const Double_t mLo = 60.; const Double_t mHi = 120.;  const Int_t nBins = 60;
const Double_t effTrue = 0.9;
const Int_t nE = 150000;
const Double_t bgFrac = 0.25;
const Double_t peakMean = 91.; const Double_t peakRMS = 3.;

const Int_t nP = 5;
const Int_t nVMax = 100;
Int_t nV;
Double_t XV[nVMax], YV[nVMax], EV[nVMax];


void getHistograms() {

  tnpPass = new TH1D("tnpPass","PASS;mass;entries",nBins,mLo,mHi);
  tnpFail = new TH1D("tnpFail","FAIL;mass;entries",nBins,mLo,mHi);  

  TRandom2 *RAN = new TRandom2();

  for (Int_t i=0; i<nE; ++i) {
    Double_t bgTest = RAN->Uniform();
    if (bgTest < bgFrac) {
      Double_t x = mLo + (mHi-mLo)*RAN->Uniform();
      Double_t whichHist = RAN->Uniform();
      if (whichHist > 0.5) {
	tnpPass->Fill(x);
      } else {
	tnpFail->Fill(x);
      }
    } else {
      Double_t x = RAN->Gaus(peakMean,peakRMS);
      Double_t passfail = RAN->Uniform();
      if (passfail < effTrue) {
	tnpPass->Fill(x);
      } else {
	tnpFail->Fill(x);
      }
    }
  }
}

//
// Signal is a Gaussian
// Background is a linear function
//
Double_t chisqFun(Double_t AG, Double_t mG, Double_t sG, Double_t aB, Double_t bB) {
  Double_t sum = 0;
  for (Int_t i=0; i<nV; ++i) {
    Double_t val = AG*TMath::Gaus(XV[i],mG,sG,1) + aB + bB*(XV[i]-mLo);
    Double_t term = (YV[i] - val) / EV[i];
    sum += pow( term, 2 );
  }
  return sum;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &thefun, Double_t *par, Int_t iflag) {
  Double_t AG = par[0];
  Double_t mG = par[1];
  Double_t sG = par[2];
  Double_t aB = par[3];
  Double_t bB = par[4];
  Double_t val = chisqFun( AG, mG, sG, aB, bB );
  thefun = val;
}



void fitHistograms() {
  //
  // Set up fit
  //
  TMinuit *gMinuit = new TMinuit(nP);
  gMinuit->SetFCN(fcn);
  gMinuit->SetPrintLevel(-1); // Set to +1 to see the fit in action
  gMinuit->SetErrorDef(1.);  // This is a chi-squared fit
  //
  Double_t vstart[nP] = { 10000., 91.,  3.   , 300., 0.};
  Double_t step[nP]   = { 100.  , 0.1,  0.1  ,  10., 0.01};
  Int_t ierflg = 0;
  gMinuit->mnparm(0, "A",    vstart[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "mean", vstart[1], step[1], 0, 0, ierflg);
  gMinuit->mnparm(2, "rms",  vstart[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "a",    vstart[3], step[3], 0, 0, ierflg);
  gMinuit->mnparm(4, "b",    vstart[4], step[4], 0, 0, ierflg);    
  //
  // extract data from PASS histogram
  //
  nV = 0;
  for (Int_t k=0; k<nBins; ++k) {
    if ( tnpPass->GetBinContent(k+1) > 0. ) {
      XV[k] = tnpPass->GetBinCenter(k+1);
      YV[k] = tnpPass->GetBinContent(k+1);
      EV[k] = tnpPass->GetBinError(k+1);
      nV++;
    }
  }
  //
  // fit the PASS histogram
  //
  Double_t arglist[10];
  arglist[0] = 1500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  Double_t vpar1[nP], vunc1[nP];
  for (Int_t k=0; k<nP; ++k) {gMinuit->GetParameter( k, vpar1[k], vunc1[k] );}
  //
  //
  // extract data from FAIL histogram
  //
  nV = 0;
  for (Int_t k=0; k<nBins; ++k) {
    if ( tnpFail->GetBinContent(k+1) > 0. ) {
      XV[k] = tnpFail->GetBinCenter(k+1);
      YV[k] = tnpFail->GetBinContent(k+1);
      EV[k] = tnpFail->GetBinError(k+1);
      nV++;
    }
  }
  //
  // fit the FAIL histogram
  //
  vstart[0] = 1000.;
  gMinuit->mnparm(0, "A",    vstart[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "mean", vstart[1], step[1], 0, 0, ierflg);
  gMinuit->mnparm(2, "rms",  vstart[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "a",    vstart[3], step[3], 0, 0, ierflg);
  gMinuit->mnparm(4, "b",    vstart[4], step[4], 0, 0, ierflg);    
  //
  arglist[0] = 1500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  Double_t vpar2[nP], vunc2[nP];
  for (Int_t k=0; k<nP; ++k) {gMinuit->GetParameter( k, vpar2[k], vunc2[k] );}
  //
  // Process results
  //
  cout << "\nEfficiency results:\n\n";
  //
  Double_t P = vpar1[0]; Double_t uP = vunc1[0];
  Double_t F = vpar2[0]; Double_t uF = vunc2[0];  
  Double_t eff = P/(F+P);
  Double_t effUnc = eff * sqrt( pow( (1./P-1./(F+P))*uP, 2 ) + pow( -1./(F+P)*uF, 2 ) );
  printf("\teff= %8.4f +/- %7.4f \n", eff, effUnc );
  
}




void plotHistograms() {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(10);
  //  gStyle->SetOptFit(0011);  

  TCanvas *CF = new TCanvas("CF","Fail",10,10,800,500);
  CF->SetGridx(1);
  CF->SetGridy(1);
  tnpFail->SetMinimum(0);
  tnpFail->Draw("e");
  CF->Print("plot_tnpFail.pdf");

  TCanvas *CP = new TCanvas("CP","Pass",810,10,800,500);
  CP->SetGridx(1);
  CP->SetGridy(1);
  tnpPass->SetMinimum(0);
  tnpPass->Draw("e");
  CP->Print("plot_tnpPass.pdf");
}





void sampleFit() {
  cout << "...........A sample fit...........\n";

  getHistograms();

  fitHistograms();

  plotHistograms();

}
