#! /bin/sh

currentDir=${PWD##*/}

if [ "$currentDir" = madgraph ]
then
    tar -czf source.tar.gz --exclude='.git*' --exclude='output/*' --exclude='*.tar.gz' ../madgraph ../python
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
    echo "Must submit from /madgraph directory!  source.tar.gz not created."
fi
