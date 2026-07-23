#!/bin/bash
#
# Benchmark script for CPU and GPU implementations of TWiCE.
#

set -euo pipefail

CC=cc
NVCC=nvcc

SOURCE=main.c
CPU_EXE=twicec
GPU_EXE=twiceg

CONF=conf.conf
WIND=wind.txt
TOPO=./_topodata/ks_topo_utm100.txt

OUTPUT=massloading.txt
LOGFILE=benchmark.txt

# Save all output to both screen and log file
exec > >(tee "$LOGFILE")
exec 2>&1

echo "========================================"
echo "TWiCE CPU/GPU Benchmark"
echo "Started: $(date)"
echo "========================================"
echo

echo "=== Compile CPU version ==="
$CC "$SOURCE" -o "$CPU_EXE" -lm

echo "=== Compile GPU version ==="
$NVCC -x cu -DCUDA "$SOURCE" -o "$GPU_EXE"

echo
echo "=============================="
echo "GPU benchmark (5 runs)"
echo "=============================="

for i in {1..5}
do
    rm -f "$OUTPUT"

    start=$(date +%s.%N)

    ./"$GPU_EXE" "$CONF" "$WIND" "$TOPO" > /dev/null

    end=$(date +%s.%N)

    elapsed=$(awk "BEGIN {print $end - $start}")

    printf "Run %d : %.3f sec\n" "$i" "$elapsed"
done

<< EOF
echo
echo "=============================="
echo "CPU benchmark (5 runs)"
echo "=============================="

for i in {1..5}
do
    rm -f "$OUTPUT"

    start=$(date +%s.%N)

    ./"$CPU_EXE" "$CONF" "$WIND" "$TOPO" > /dev/null

    end=$(date +%s.%N)

    elapsed=$(awk "BEGIN {print $end - $start}")

    printf "Run %d : %.3f sec\n" "$i" "$elapsed"
done

EOF

echo
echo "========================================"
echo "Finished: $(date)"
echo "Results saved to $LOGFILE"
echo "========================================"