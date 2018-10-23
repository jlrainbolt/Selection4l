#ifndef SELECTIONTOOLS_HH
#define SELECTIONTOOLS_HH


// STL
#include <vector>

// ROOT

// Custom
#include "Lepton.hh"
#include "LeptonPair.hh"

using namespace std;



//
//  PAIRING
//

// Make pairs using mother indices
bool MakePairsFromMother(const vector<Lepton>&, LeptonPair*, LeptonPair*);

// Make pairs by maximizing mass difference of opposite-sign pairs
bool MakePairsMaxDiff(const vector<Lepton>&, LeptonPair*, LeptonPair*);



#endif  // SELECTIONTOOLS_HH
