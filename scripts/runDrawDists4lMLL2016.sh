#! /bin/bash

##
##  runDrawDists4l2016.sh
##
##  Runs DrawDists4l.cc on all 2016 (data and MC) samples
##

for suffix in   "electron_2016" "muon_2016"     "zjets_m-50"    "ttbar"     "tt_2l2nu"  "ttz_2l2nu" \
                "ww_2l2nu"      "wz_2l2q"       "wz_3lnu"       "zz_2l2nu"  "zz_2l2q"   "zz_4l_m-1" \
                "ggH_zz_4l"     "vbfH_zz_4l"    "wwz_4l2nu"     "wzz_4l2nu" "zzz_4l2nu" "zzg_4l2nu"

do
    root.exe -q -b "../macros/DrawDists4lMLL.cc(\"$suffix\", \"2016\")"
done

for suffix in   "electron_2016" "muon_2016"
do
    root.exe -q -b "../macros/DrawDists4lMLL.cc(\"$suffix\", \"2016\", \"kTRUE\")"
done
