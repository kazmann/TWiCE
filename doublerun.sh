#!/bin/bash
#
# Validation script for CPU and GPU implementations of TWiCE.
#
# Before the first use, make this script executable:
#   chmod +x doublerun.sh
#
# Then run:
#   ./doublerun.sh
#

set -euo pipefail

CC=cc
NVCC=nvcc

SOURCE=main.c
CPU_EXE=twicec
GPU_EXE=twiceg

CONF=conf.conf
WIND=wind.txt
TOPO=./_topodata/ks_topo_utm200.txt

OUTPUT=massloading.txt
CPU_OUTPUT=massloading_cpu.txt
GPU_OUTPUT=massloading_gpu.txt


echo "=== Compile CPU version ==="

$CC "$SOURCE" -o "$CPU_EXE" -lm


echo "=== Compile GPU version ==="

$NVCC -x cu -DCUDA "$SOURCE" -o "$GPU_EXE"


rm -f "$OUTPUT" "$CPU_OUTPUT" "$GPU_OUTPUT"


echo "=== Run CPU version ==="

start=$(date +%s.%N)

./"$CPU_EXE" "$CONF" "$WIND" "$TOPO"

end=$(date +%s.%N)

cpu_elapsed=$(awk "BEGIN {print $end - $start}")

if [ ! -f "$OUTPUT" ]; then
    echo "ERROR: CPU version did not produce $OUTPUT"
    exit 1
fi

mv "$OUTPUT" "$CPU_OUTPUT"


echo "=== Run GPU version ==="

start=$(date +%s.%N)

./"$GPU_EXE" "$CONF" "$WIND" "$TOPO"

end=$(date +%s.%N)

gpu_elapsed=$(awk "BEGIN {print $end - $start}")

if [ ! -f "$OUTPUT" ]; then
    echo "ERROR: GPU version did not produce $OUTPUT"
    exit 1
fi

mv "$OUTPUT" "$GPU_OUTPUT"


echo "=== Compare outputs ==="

<< EOF
if diff -q "$CPU_OUTPUT" "$GPU_OUTPUT" > /dev/null; then
    echo "OK: CPU and GPU outputs match exactly"
else
    echo "ERROR: CPU and GPU outputs differ"
    diff "$CPU_OUTPUT" "$GPU_OUTPUT"
    exit 1
fi
EOF

printf "CPU elapsed time: %.3f sec\n" "$cpu_elapsed"
printf "GPU elapsed time: %.3f sec\n" "$gpu_elapsed"