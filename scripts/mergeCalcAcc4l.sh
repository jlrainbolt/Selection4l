#! /bin/sh

systematics="$1"
hist="$2"
suffix="$3"

prefix="${systematics}${hist}"

./../scripts/mergeOutput.sh ${prefix} ${suffix}

root.exe -q -b "../macros/CalculateAccEff4l.cc(\"${suffix}\", \"${systematics}\", \"${hist}\")"

rm "output/${prefix}_${suffix}.root"
echo "Removed output/${prefix}_${suffix}.root"
