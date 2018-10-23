// STL
#include <vector>
#include <iostream>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "AnaUtils.hh"


using namespace std;



//
//  BOOSTING
//

// Get TLorentzVector boosted in frame of "system" where beta = system.BoostVector()
TLorentzVector  BoostP4(    const TLorentzVector &p4,  const TVector3 &beta)
{
    TLorentzVector _p4(p4);
    _p4.Boost(-beta);
    return _p4;
}

// Get TVector3 in frame of "system" where beta = system.BoostVector()
TVector3        BoostP3(    const TLorentzVector &p4,  const TVector3 &beta)
{
    TLorentzVector _p4 = BoostP4(p4, beta);
    return _p4.Vect();
}



//
//  SYSTEMATICS
//
/*
// GetBinContentPtEta
double GetBinContentPtEta(const TH2 *hist, const TLorentzVector &p4)
{
    int xbin = GetXbin(hist, p4.Eta());
    int ybin = GetYbin(hist, p4.Pt());

    return hist->GetBinContent(xbin, ybin);
}
*/