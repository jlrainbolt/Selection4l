#! /bin/sh

suffix="$1"

./runHadd.sh trees "$suffix"

root.exe -q -b "drawHists.cc(\"$suffix\")"

./haddChannels.sh "$suffix"
