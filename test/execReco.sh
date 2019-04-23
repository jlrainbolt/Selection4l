#! /bin/sh

suffix="$1"
id="$2"


echo "Date:     " `date`
echo "Host:     " $HOSTNAME
echo "Location: " $TEMP
echo ""

echo "Loading ROOT..."
source /cvmfs/sft.cern.ch/lcg/views/LCG_93python3/x86_64-slc6-gcc7-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.12.06-05bac/x86_64-slc6-gcc7-opt/bin/thisroot.sh
echo "ROOTSYS = " $ROOTSYS

tar -xzvf source.tar.gz
cd "test"

root.exe -q -b "../macros/RecoSelection.cc(\"$suffix\", \"$id\")"

cp *_${suffix}_*.root ${_CONDOR_SCRATCH_DIR}

echo ""
echo ""
echo "Date:     " `date`
