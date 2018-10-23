#! /bin/sh

currentDir=${PWD##*/}

if [ "$currentDir" = test ]
then
    tar -czf source.tar.gz --exclude='*.root' --exclude='.git*' ../.
    echo "Created source.tar.gz"

    if [ ! -d "output" ]
    then
        mkdir output
        "Created directory /output"
    fi

    if [ ! -d "output/reports" ]
    then
        mkdir output/reports
        "Created directory /output/reports"
    fi

    mv source.tar.gz output/.
else
    echo "Must submit from /test directory!  source.tar.gz not created."
fi
