#! /bin/bash

##
##  runDrawDists2l2017.sh
##
##  Runs DrawDists2l.cc on all 2017 (data and MC) samples
##

for suffix in   "electron_2017" "muon_2017"     "zjets_m-50"    "ttbar"     "tt_2l2nu"  "ttz_2l2nu" \
                "ww_2l2nu"      "wz_2l2q"       "wz_3lnu"       "zz_2l2nu"  "zz_2l2q"   "zz_4l"     \
                "ggH_zz_4l"     "vbfH_zz_4l"    "wwz_4l2nu"     "wzz_4l2nu" "zzz_4l2nu" "zzg_4l2nu"
do
    root.exe -q -b "../macros/DrawBackground.cc(\"$suffix\", \"2017\")"
done
