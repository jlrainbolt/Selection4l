#! /bin/sh

prefix="$1"
suffix="$2"

dir="output"


# Make sure ROOT is loaded on the login node
if [ -z "${ROOTSYS}" ]
then
    echo "Loading ROOT..."
    source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
    source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt/bin/thisroot.sh
    echo "ROOTSYS = " $ROOTSYS
fi

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
