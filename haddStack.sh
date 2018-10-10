#! /bin/sh

hadd hists_data.root hists_muon_2017.root hists_electron_2017.root

# rm hists_muon_2017.root hists_electron_2017.root

root.exe -q -b "stackHists.cc()"
