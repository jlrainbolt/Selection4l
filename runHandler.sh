#! /bin/sh

suffix="$1"
id="$2"


# Prep stuff for batch mode
if [ $TEMP ]
then
    echo "date:     " `date`
    echo "host:     " $HOSTNAME
    echo "location: " $TEMP
    echo ""
    echo ""

    source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
    source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt/bin/thisroot.sh
fi


# Run macro
root.exe -q -b "handleSelection.cc(\"$suffix\", \"$id\", \"\")"


# Copy output to scratch dir if needed
if [ $TEMP ]
then
#   cp "$selection"*.root $_CONDOR_SCRATCH_DIR/.


    echo ""
    echo ""
    echo "date:     " `date`
fi
