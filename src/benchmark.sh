
# Usage: ./benchmark.sh <dataDir> <timeLimit> <verbose>


if [ $# -ne 3 ]; then
    echo "Usage: ./benchmark.sh <dataDir> <timeLimit> <verbose>"
    exit 1
fi

echo "Experimental Campaign:"
echo "Data directory: $1"
echo "Time limit: $2"
echo "Verbose: $3"


cd ..
echo `date` > results.txt





# for i in `seq 1 $nbFiles`
#     for j in `seq 1 4`
#         do
#             echo "Running instance $i$j"
#             echo "Running instance $i$j" >> results.txt
#             if [ $3 -eq 1 ]; then
#                 ./bin/tp2 $1/model_a$i$j.txt $1/assignement_a$i$j.txt $2 >> results.txt
#             else
#                 ./bin/tp2 $1/model_a$i$j.txt $1/assignement_a$i$j.txt $2 >> results.txt 2>&1
#             fi
#         done
#     done

echo "Experimental Campaign finished!"
exit 0