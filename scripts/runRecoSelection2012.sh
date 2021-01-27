#! /bin/bash

##
##  runRecoSelection2012.sh
##
##  Runs RecoSelection.cc and BkgSelection.cc on all 2012 (data and MC) samples
##

for suffix in   "electron_2012" "muon_2012" "zjets_m-50"    "ttbar"     "ttz_2l2nu" "ww_2l2nu" \
                "wz_2l2q"       "wz_3lnu"   "zz_2l2q"       "zz_2l2nu"  "zz_4l"
do
    case $suffix in
        "electron_2012") files=3;;
        "muon_2012") files=4;;
        "zjets_m-50") files=3;;
        "ttbar") files=1;;
        "ttz_2l2nu") files=1;;
        "ww_2l2nu") files=1;;
        "wz_2l2q") files=1;;
        "wz_3lnu") files=1;;
        "zz_2l2q") files=1;;
        "zz_2l2nu") files=1;;
        "zz_4l") files=1;;
    esac

    for i in $(seq 0 $(($files-1)))
    do
        root.exe -q -b "../macros/RecoSelection.cc(\"$suffix\", \"$i\")"
        root.exe -q -b "../macros/BkgSelection.cc(\"$suffix\", \"$i\")"
    done
done
