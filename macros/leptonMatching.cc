/*
vector<tuple<TLorentzVector, Short_t, Float_t>> SortGenLeps(
        vector<tuple<TLorentzVector, Short_t, Float_t>> _genLeps,
        vector<tuple<TLorentzVector, Short_t, Float_t>> recoLeps,
        const TString &comp)
{
    //  vector<tuple<TLorentzVector, Short_t, Float_t>> genLeps = _genLeps;

    vector<tuple<TLorentzVector, Short_t, Float_t>> genLeps;

    // Make sure there are enough gen leptons
    if (recoLeps.size() > _genLeps.size())
        return genLeps;


    // Make sure both input vectors are sorted
    if      (comp == "Pt")
    {
        sort(recoLeps.begin(), recoLeps.end(), SortDecPt);
        sort(_genLeps.begin(), _genLeps.end(), SortDecPt);
    }
    else if (comp == "P")
    {
        sort(recoLeps.begin(), recoLeps.end(), SortDecP);
        sort(_genLeps.begin(), _genLeps.end(), SortDecP);
    }


    // Loop over each reco lepton
    for (unsigned i = 0; i < recoLeps.size(); i++)
    {
        TLorentzVector recoP4 = get<0>(recoLeps[i]);

        // Find DeltaR between reco lep and each gen lep
        vector<Float_t> deltaR;
        for (unsigned j = 0; j < _genLeps.size(); j++)
        {
            TLorentzVector genP4 = get<0>(_genLeps[j]);
            deltaR.push_back(recoP4.DeltaR(genP4));
        }

        // Find index of minimum DeltaR
        unsigned m = min_element(deltaR.begin(), deltaR.end()) - deltaR.begin();
        genLeps.push_back(make_tuple(get<0>(_genLeps[m]), get<1>(_genLeps[m]), deltaR[m]));

        // Remove gen lepton from input list
        _genLeps.erase(_genLeps.begin() + m);
    }

    return genLeps;
}
*/
