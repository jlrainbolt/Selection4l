#! /bin/sh

fileName="$1_$2"

# hadd all histograms
if hadd "$fileName".root "$fileName"_*.root
then
    echo ""
    for file in "$fileName"_*.root
    do
        if rm $file
        then
            echo "Removed $file"
        fi
    done
fi
