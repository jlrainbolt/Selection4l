#! /bin/sh

selection="$1"
nametag="$2"
path="/uscms/home/jrainbol/work/Selection/Scripts/"

if mkdir "$selection""_""$nametag"
then
    cd "$selection""_""$nametag"
    cp 
    mkdir RunReports
    mkdir DagReports
    condor_submit_dag "$path""submitSelection_""$selection"".dag" -usedagdir -outfile_dir DagReports
    cd ..
fi
