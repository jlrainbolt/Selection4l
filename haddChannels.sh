#! /bin/sh

suffix="$1"
selection=("mumu" "ee" "4m" "2m2e" "2e2m" "4e")

for sel in "${selection[@]}"
do
    rootcp hists_${suffix}.root:${sel}/* "${sel}_${suffix}.root"
done

hadd ll_${suffix}.root mumu_${suffix}.root ee_${suffix}.root

hadd 4l_${suffix}.root 4m_${suffix}.root 2m2e_${suffix}.root 2e2m_${suffix}.root 4m_${suffix}.root

channel=("ll" "4l")

for chan in "${channel[@]}"
do
    rootmkdir hists_${suffix}.root:${chan}
    rootcp ${chan}_${suffix}.root:* hists_${suffix}.root:${chan}
    rm ${chan}_${suffix}.root
done

for sel in "${selection[@]}"
do
    rm ${sel}_${suffix}.root
done
