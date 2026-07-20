#!/bin/bash
set -e

nvcc main.cu -o twice

start=$(date +%s.%N)

#./twice conf.conf wind.txt topo.txt
./twice conf.conf wind.txt ./_topodata/ks_topo_utm1000.txt

end=$(date +%s.%N)

elapsed=$(awk "BEGIN {print $end - $start}")

<< EOF
if diff massloading.txt massloading_in_loc.txt > /dev/null; then
    echo "OK: output matches"
else
    echo "ERROR: output differs"
    diff massloading.txt massloading_in_loc.txt
    exit 1
fi
EOF

printf "Elapsed time: %.3f sec\n" "$elapsed"