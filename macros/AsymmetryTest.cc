void AsymmetryTest()
{
    TH1D *asy = new TH1D("asymmetry", "asymmetry", 250, -1, 1);

    TRandom3 *rng = new TRandom3();

    for (unsigned i = 0; i < 1000; i++)
    {
        float pos = rng->PoissonD(186);
        float neg = rng->PoissonD(237);
        float asy_ = (pos - neg) / (pos + neg);

        asy->Fill(asy_);
    }

    asy->Draw();
}
