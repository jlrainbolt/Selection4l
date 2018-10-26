#! /bin/bash

##
##  runDrawDists4l2017.sh
##
##  Runs DrawDists4l.cc on all 2017 (data and MC) samples
##

for suffix in "phase_space" "fiducial" "zz_4l"
do
    root.exe -q -b "../macros/DrawDists4l.cc(\"$suffix\")"
done

root.exe -q -b "../macros/DivideDists4l.cc"
