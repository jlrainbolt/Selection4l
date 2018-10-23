#! /bin/sh

prefix="$1"
suffix="$2"

dir="output"

if hadd ${dir}/${prefix}_${suffix}.root ${dir}/${prefix}_${suffix}_*.root
then
    echo ""
    for file in ${dir}/${prefix}_${suffix}_*.root
    do
        if rm $file
        then
            echo "Removed $file"
        fi
    done
fi
