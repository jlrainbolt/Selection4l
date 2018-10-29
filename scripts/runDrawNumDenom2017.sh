#! /bin/bash

##
##  runDrawNumDenom.sh
##
##  Runs DrawNumDenom.cc on all 2017 (data and MC) samples
##

for suffix in "electron_2017" "muon_2017" "zjets_m-50" "ttbar" "ww_2l2nu" "wz_2l2q" "wz_3lnu" "zz_2l2q" "zz_4l" "ggH_zz_4l" "vbfH_zz_4l"
do
    root.exe -q -b "../macros/DrawNumDenom.cc(\"$suffix\")"
done

for suffix in "zjets_m-50" "zz_4l"
do
    root.exe -q -b "../macros/DrawNumDenom.cc(\"$suffix\")"
done
