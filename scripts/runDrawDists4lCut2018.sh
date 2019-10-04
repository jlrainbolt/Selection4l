#! /bin/bash

##
##  runDrawDists4l2018.sh
##
##  Runs DrawDists4l.cc on all 2018 (data and MC) samples
##

for suffix in   "electron_2018" "muon_2018"     "zjets_m-50"    "ttbar"     "tt_2l2nu"  "ttz_2l2nu" \
                "ww_2l2nu"      "wz_2l2q"       "wz_3lnu"       "zz_2l2nu"  "zz_2l2q"   "zz_4l"     \
                "ggH_zz_4l"     "vbfH_zz_4l"    "wwz_4l2nu"     "wzz_4l2nu" "zzz_4l2nu" "zzg_4l2nu"
do
    root.exe -q -b "../macros/DrawDists4lCut.cc(\"$suffix\", \"2018\", \"!singleMuTrig\")"
    root.exe -q -b "../macros/DrawDists4lCut.cc(\"$suffix\", \"2018\", \"!doubleMuTrig\")"
done

for suffix in   "electron_2018" "muon_2018"
do
    root.exe -q -b "../macros/DrawDists4lCut.cc(\"$suffix\", \"2018\", \"!singleMuTrig\", \"kTRUE\")"
    root.exe -q -b "../macros/DrawDists4lCut.cc(\"$suffix\", \"2018\", \"!doubleMuTrig\", \"kTRUE\")"
done
