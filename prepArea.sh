#! /bin/sh

selection="$1"
myROOTSYS="/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt"

echo "date:     " `date`
echo "host:     " $HOSTNAME


# Load ROOT 6 if needed
if [ $ROOTSYS != $myROOTSYS ]
then
    source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
    source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt/bin/thisroot.sh
    echo "Loaded ROOT 6.10.02"
fi

./cleanDirectory.sh

echo "Submitting $selection jobs..."
