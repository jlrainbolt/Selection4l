#! /bin/sh

systematics="$1"
suffix="$2"

cat output/${systematics}_${suffix}_*.txt > "output/${systematics}_${suffix}.txt"

rm output/${systematics}_${suffix}_*.txt
