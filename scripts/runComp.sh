#! /bin/sh

suffix="$1"

#root.exe -q -b "genSelection.cc(\"$suffix\")"

root.exe -q -b "drawGenHists.cc(\"$suffix\")"

./haddChannels.sh "PhaseSpace_$suffix"
./haddChannels.sh "Fiducial_$suffix"
./haddChannels.sh "GenSelected_$suffix"
./haddChannels.sh "RecoSelected_$suffix"

root.exe -q -b "overlayHists.cc(\"$suffix\")"
