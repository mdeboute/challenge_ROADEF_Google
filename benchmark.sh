#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: ./benchmark.sh <dataDir> <timeLimit (sec)> <verbose (int)>"
    exit 1
fi

echo "Experimental Campaign:"
echo "Data directory: $1"
echo "Time limit: $2 seconds"
echo "Verbose: $3"

echo `date` > results.txt

nbFiles=`ls $1 | wc -l`
nbFiles=`expr $nbFiles / 2`
echo "Number of instances: $nbFiles"

modelFile="model"
assignmentFile="assignment"
suffix=`ls $1 | head -1 | cut -d'_' -f2`

for i in `seq 1 $nbFiles`; do
    echo "Solving instance $i"
    echo "Solving instance $i" >> results.txt
    echo "-----------------------------------------------------------------------------" >> results.txt
    python3 src/run.py ./$1/$modelFile"_"$suffix"_"$i.txt ./$1/$assignmentFile"_"$suffix"_"$i.txt $2 $3 >> results.txt
    echo "-----------------------------------------------------------------------------" >> results.txt
done

echo "Experimental Campaign finished!"
exit 0
