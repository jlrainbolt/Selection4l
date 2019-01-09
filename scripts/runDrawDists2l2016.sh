#! /bin/bash

##
##  runDrawDists2l2016.sh
##
##  Runs DrawDists2l.cc on all 2016 (data and MC) samples
##

for suffix in   "electron_2016" "muon_2016" "dy_m-50" "dy_m-10to50" "ttbar" "ttz_2l2nu" "ww_2l2nu" \
                "wz_3lnu" "zz_4l" "ggH_zz_4l" "H_zz_4l"
do
    root.exe -q -b "../macros/DrawDists2l.cc(\"$suffix\")"
done
