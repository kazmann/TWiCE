#!/bin/bash

set -e

OUT_SUMMARY="point_number_test.txt"

echo "resolution input_points elapsed_seconds output_points" > "$OUT_SUMMARY"

for res in 1000 500 300 200 100
do
    topo="./_topodata/ks_topo_utm${res}.txt"

    if [ ! -f "$topo" ]; then
        echo "ERROR: $topo が見つかりません"
        exit 1
    fi

    input_points=$(grep -cv '^[[:space:]]*$' "$topo")

    echo "Running resolution ${res}, input points ${input_points}"

    start=$(date +%s.%N)

    ./twice conf.conf wind.txt "$topo"

    end=$(date +%s.%N)

    elapsed=$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.6f", e - s }')

    if [ ! -f "massloading.txt" ]; then
        echo "ERROR: massloading.txt が生成されませんでした"
        exit 1
    fi

    output_points=$(grep -cv '^[[:space:]]*$' massloading.txt)

    new_massloading="massloading.${res}.${input_points}.txt"
    mv massloading.txt "$new_massloading"

    echo "${res} ${input_points} ${elapsed} ${output_points}" >> "$OUT_SUMMARY"

done

echo "Done. Results written to ${OUT_SUMMARY}"