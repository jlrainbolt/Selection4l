#! /bin/bash

##
##  runRecoSelection2017.sh
##
##  Runs RecoSelection.cc and BkgSelection.cc on all 2017 (data and MC) samples
##

for suffix in   "electron_2017" "muon_2017"     "zjets_m-50"    "ttbar"     "tt_2l2nu"  "ttz_2l2nu" \
                "ww_2l2nu"      "wz_2l2q"       "wz_3lnu"       "zz_2l2nu"  "zz_2l2q"   "zz_4l"     \
                "ggH_zz_4l"     "vbfH_zz_4l"    "wwz_4l2nu"     "wzz_4l2nu" "zzz_4l2nu" "zzg_4l2nu" 
do
    case $suffix in
        "electron_2017") files=6;;
        "muon_2017") files=10;;
        "zjets_m-50") files=18;;
        "ttbar") files=1;;
        "tt_2l2nu") files=1;;
        "ttz_2l2nu") files=1;;
        "ww_2l2nu") files=1;;
        "wz_2l2q") files=2;;
        "wz_3lnu") files=1;;
        "zz_2l2q") files=2;;
        "zz_2l2nu") files=1;;
        "zz_4l") files=1;;
        "ggH_zz_4l") files=1;;
        "vbfH_zz_4l") files=1;;
        "wwz_4l2nu") files=1;;
        "wzz_4l2nu") files=1;;
        "zzz_4l2nu") files=1;;
        "zzg_4l2nu") files=1;;
    esac

    for i in $(seq 0 $(($files-1)))
    do
        root.exe -q -b "../macros/RecoSelection.cc(\"$suffix\", \"$i\")"
        root.exe -q -b "../macros/BkgSelection.cc(\"$suffix\", \"$i\")"
    done
done
