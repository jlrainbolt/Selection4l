#! /bin/sh

selection="$1"
year="2017"

outFile="$selection""_""$year"".root"

root.exe -q -b "prepareFile.cc(\"\", \"$selection\")"
root.exe -q -b "stackResults.cc(\"$outFile\", \"$selection\")"
echo "Created combined histograms in $outFile".
