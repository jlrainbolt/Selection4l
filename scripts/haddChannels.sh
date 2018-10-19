#! /bin/sh

suffix="$1"
selection=("4m" "2m2e" "2e2m" "4e")

for sel in "${selection[@]}"
do
    rootcp hists_${suffix}.root:${sel}/* "${sel}_${suffix}.root"
done

hadd 4l_${suffix}.root 4m_${suffix}.root 2m2e_${suffix}.root 2e2m_${suffix}.root 4e_${suffix}.root

rootmkdir hists_${suffix}.root:4l
rootcp 4l_${suffix}.root:* hists_${suffix}.root:4l
rm 4l_${suffix}.root


for sel in "${selection[@]}"
do
    rm ${sel}_${suffix}.root
done
