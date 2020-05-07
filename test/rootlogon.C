//// rootlogon()
{
    TH1::SetDefaultSumw2();


    // Load source files
    gInterpreter->AddIncludePath("../include");
    gROOT->SetMacroPath("../src");

    TSystemDirectory src("src", "../src");
    TIter next(src.GetListOfFiles());
    TSystemFile *_srcFile;
    while ((_srcFile = (TSystemFile*) next()))
    {
        TString srcFile = _srcFile->GetName();
        if (srcFile.Contains(".cc") && !srcFile.Contains(".sw"))
            gROOT->ProcessLine(".L " + srcFile);
    }
    cout << "Loaded files in ../src" << endl;


    // Link macros
    gROOT->SetMacroPath("../macros");
    cout << "Linked macros in ../macros" << endl;


    // Load custom colors
    gROOT->ProcessLine(".x LoadColors.cc");


    // Set plotting options
//  gStyle->SetErrorX(0);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1100);
/*
//  gStyle->SetPaintTextFormat(".4f");


    // Load RooUnfold
    gInterpreter->AddIncludePath("~/RooUnfold-trunk/src");
    gSystem->Load("~/RooUnfold-trunk/libRooUnfold");
    cout << "Loaded RooUnfold libraries" << endl;
*/

    cout << endl;
}
