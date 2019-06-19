#! /bin/bash

##
##  runDrawDists4l2017.sh
##
##  Runs DrawDists4l.cc on all 2017 (data and MC) samples
##

for suffix in   "electron_2016" "muon_2016" "zjets_m-50" "ttbar" "tt_2l2nu" "ttz_2l2nu" \
                "ww_2l2nu" "wz_2l2q" "wz_3lnu" "zz_2l2q" "zz_2l2nu" "zz_4l" "phase_space" \
                "ggH_zz_4l" "vbfH_zz_4l"
do
    root.exe -q -b "../macros/DrawDists4l.cc(\"$suffix\", \"2016\")"
done

for suffix in   "electron_2016" "muon_2016"
do
    root.exe -q -b "../macros/DrawDists4l.cc(\"$suffix\", \"2016\", \"kTRUE\")"
done
