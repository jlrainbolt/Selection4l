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

// Make pairs by maximizing mass of one same-flavor, opposite-sign pair
bool MakePairsMaxZ1(const vector<Lepton>&, LeptonPair*, LeptonPair*);

// Make pairs by naively looking for opposite sign
bool MakePairs6l(const vector<Lepton>&, LeptonPair*, LeptonPair*, LeptonPair*);



//
//  SCALE FACTORS
//

//float GetTriggerWeight(const vector<Lepton>&);



#endif  // SELECTIONTOOLS_HH
