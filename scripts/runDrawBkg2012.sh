#! /bin/bash

##
##  runDrawBkg2012.sh
##
##  Runs DrawBackground.cc on all 2012 (data and MC) samples
##

for suffix in   "electron_2012" "muon_2012"     "zjets_m-50"    "ttbar"     "ttz_2l2nu" \
                "ww_2l2nu"      "wz_2l2q"       "wz_3lnu"       "zz_2l2nu"  "zz_2l2q"   "zz_4l"
do
    root.exe -q -b "../macros/DrawBackground.cc(\"$suffix\", \"2012\")"
    root.exe -q -b "../macros/DrawBackground.cc(\"$suffix\", \"2012\", kTRUE)"
done
