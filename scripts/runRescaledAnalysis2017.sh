#! /bin/bash

##
##  runRescaledAnalysis2016.sh
##
##  Runs RescaledAnalysis.cc on all 2016 (data and MC) samples
##

for suffix in   "zjets_m-50" "ttbar" "ww_2l2nu" "wz_2l2q" \
                "wz_3lnu" "zz_2l2q" "zz_4l" "ggH_zz_4l" "vbfH_zz_4l"
do
    root.exe -q -b "../macros/RescaledAnalysis.cc(\"$suffix\")"
done
