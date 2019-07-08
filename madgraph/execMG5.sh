#! /bin/sh

model="$1"
finalstate="$2"
nevents="$3"
MUb="$4"
gUbe="$5"
gUbmu="$6"


echo "Date:     " `date`
echo "Host:     " $HOSTNAME
echo "Location: " $TEMP
echo ""

echo "Loading MG5..."
source /cvmfs/sft.cern.ch/lcg/views/LCG_93/x86_64-slc6-gcc62-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/MCGenerators/madgraph5amc/2.6.0-d53b7/x86_64-slc6-gcc62-opt/madgraph5amcenv-genser.sh

tar -xzf source.tar.gz

cd madgraph/${model}_pp_z_${finalstate}

echo "generate_events" >> input.txt
echo "set nevents $nevents" >> input.txt
echo "set iseed 0" >> input.txt

if [[ ${model} == *"U" ]]
then
    echo "set MUb $MUb" >> input.txt
    echo "set WUb Auto" >> input.txt
    echo "set gUbe $gUbe" >> input.txt
    echo "set gUbmu $gUbmu" >> input.txt
fi

./bin/madevent input.txt

echo "Loading ROOT..."
source /cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc62-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.16.00-f8770/x86_64-slc6-gcc62-opt/bin/thisroot.sh
echo "ROOTSYS = " $ROOTSYS

cd Events/run_01
gunzip unweighted_events.lhe.gz

python ${TEMP}/python/LHEtoTTree.py ${nevents}


cp *.root ${_CONDOR_SCRATCH_DIR}

echo ""
echo ""
echo "Date:     " `date`
