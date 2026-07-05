#!/bin/bash
set -e

nvcc main.cu -o twice

start=$(date +%s.%N)

./twice conf.conf wind.txt topo.txt

end=$(date +%s.%N)

elapsed=$(awk "BEGIN {print $end - $start}")

if diff decimal_falldriftX_-1.txt decimal_referenece.txt > /dev/null; then
    echo "OK: output matches"
    printf "Elapsed time: %.3f sec\n" "$elapsed"
else
    echo "ERROR: output differs"
    #diff decimal_falldriftX_-1.txt decimal_referenece.txt
    printf "Elapsed time: %.3f sec\n" "$elapsed"
    exit 1
fi