#! /bin/bash

##
##  runDrawDists2l2016.sh
##
##  Runs DrawDists2l.cc on all 2016 (data and MC) samples
##

for suffix in   "electron_2016" "muon_2016" "zjets_m-50" "ttbar" "tt_2l2nu" "ttz_2l2nu" \
                "ww_2l2nu" "wz_2l2q" "wz_3lnu" "zz_2l2q" "zz_2l2nu" "zz_4l" \
                "ggH_zz_4l" "vbfH_zz_4l"
do
    root.exe -q -b "../macros/DrawDists2lExt.cc(\"$suffix\", \"2016\")"
done
