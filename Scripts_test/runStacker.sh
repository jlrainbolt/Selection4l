#! /bin/sh

selection="$1"
year="2016"

outFile="$selection""_""$year"".root"

case $selection in
    "mumu")
        lepton="muon"
        ;;
    "ee")
        lepton="electron"
        ;;
esac
dataset="$lepton""_""$year"


# Create empty file to put everything in
root.exe -q -b "prepareFile.cc(\"$selection\")"


# Move histograms from sample files
for suffix in `rootls "$outFile":"Histograms"`
do
    outPath="$outFile":"Histograms/""$suffix"
    inFile="$selection""_""$suffix"".root"
    if rootcp "$inFile":* "$outPath"
    then
        echo "Moved contents from $inFile to $outPath".
#       if rm $inFile
#       then
#           echo "Removed $inFile"
#       fi
    fi
done


# Create canvases and stacks
root.exe -q -b "stackResults.cc(\"$outFile\", \"$selection\")"
echo "Created combined histograms in $outFile".


# Remove crap from file
#rootrm "$outFile":"Histograms/*/color"
#rootrm "$outFile":"Histograms/*/isData"
#rootrm "$outFile":"Histograms/*/lumi"
#rootrm "$outFile":"Histograms/*/xsec"


echo "All done!"
