#! /bin/sh

suffix="$1"

root.exe -q -b ".x genSelection.cc(\"$suffix\")"

root.exe -q -b ".x compHists.cc(\"$suffix\")"

./haddChannels.sh "comp_$suffix"
