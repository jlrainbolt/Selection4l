//// rootlogon()
{
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


    // Load plot utilities
    gROOT->ProcessLine(".x loadPlotSetup.cc");


    cout << endl;
}
