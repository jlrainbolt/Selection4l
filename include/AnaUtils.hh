#ifndef ANAUTILS_HH
#define ANAUTILS_HH


// STL
#include <vector>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"

// Custom
#include "Lepton.hh"

using namespace std;



//
//  BOOSTING
//

// Get TLorentzVector boosted into TVector3 frame
TLorentzVector  BoostP4(        const TLorentzVector&,  const TVector3&);

// Get TVector3 boosted into TVector3 frame
TVector3        BoostP3(        const TLorentzVector&,  const TVector3&);



//
//  SYSTEMATICS
//
/*
// Get TH2 bin from particle Pt, Eta
double GetBinContentPtEta(const TH2 *hist, const TLorentzVector &p4);
*/

#endif  // ANAUTILS_HH
