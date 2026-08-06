#!/bin/bash
# this code compile and evaluate uniqueness of the output file (massloading.txt)
# altered on July 23, 2026 because unified version of funcD01 was implemented on July 5.
set -e

nvcc -x cu -DCUDA main.c -o twiceg

start=$(date +%s.%N)

#./twice conf.conf wind.txt topo.txt
./twiceg conf.conf wind.txt ./_topodata/ks_topo_utm1000.txt
#./twiceg conf.conf wind.txt ./_topodata/ks_topo_utm200.txt

end=$(date +%s.%N)

elapsed=$(awk "BEGIN {print $end - $start}")


#if diff massloading.txt massloading_0723.txt > /dev/null; then
if diff massloading.txt massloading_0806.txt > /dev/null; then
    echo "OK: output matches"
else
    echo "ERROR: output differs"
    #diff massloading.txt massloading_0723.txt
    diff massloading.txt massloading_0806.txt
    exit 1
fi


printf "Elapsed time: %.3f sec\n" "$elapsed"
