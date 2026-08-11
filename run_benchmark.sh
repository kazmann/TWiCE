#!/usr/bin/env bash
#
# Benchmark TWiCE CUDA performance for combinations of:
#   - chunk_locdim
#   - CUDA blocksize
#
# For each combination:
#   1. Rewrite the two parameters in a temporary source file
#   2. Compile with nvcc
#   3. Run the executable REPEATS times
#   4. Write individual times, mean, and sample standard deviation to CSV
#
# Usage:
#   chmod +x benchmark_cuda.sh
#   ./benchmark_cuda.sh
#
# Edit the SETTINGS section below before running.
#

set -euo pipefail

# ============================================================
# SETTINGS
# ============================================================

SOURCE="main1.c"
NVCC="nvcc"

CONF="conf.conf"
WIND="wind.txt"
LOC="./_topodata/ks_topo_utm200.txt"

REPEATS=25

# Values to test
#CHUNK_SIZES=(256 512 1024 2048 4096 8192 16384 32768 65536)
BLOCK_SIZES=(16 32 64 128 256 512)
CHUNK_SIZES=(16384 32768 65536)
#BLOCK_SIZES=(32 64)

# Output files/directories
RESULT_CSV="cuda_benchmark.csv"
BUILD_DIR=".benchmark_build"

# Optional nvcc flags
NVCC_FLAGS=(-x cu -DCUDA)

# Set to 1 to preserve stdout/stderr from each run as log files.
# Set to 0 to discard normal program output.
SAVE_LOGS=0

# ============================================================
# CHECKS
# ============================================================

command -v "$NVCC" >/dev/null 2>&1 || {
    echo "Error: nvcc was not found." >&2
    exit 1
}

for file in "$SOURCE" "$CONF" "$WIND" "$LOC"; do
    if [[ ! -f "$file" ]]; then
        echo "Error: file not found: $file" >&2
        exit 1
    fi
done

mkdir -p "$BUILD_DIR"

TEMP_SOURCE="$BUILD_DIR/main_benchmark.c"

# The original source is never modified.
cp "$SOURCE" "$TEMP_SOURCE"

# Verify the expected assignments before starting.
chunk_count=$(grep -Ec '^[[:space:]]*int[[:space:]]+chunk_locdim[[:space:]]*=[[:space:]]*[0-9]+[[:space:]]*;' "$TEMP_SOURCE" || true)
block_count=$(grep -Ec '^[[:space:]]*int[[:space:]]+blocksize[[:space:]]*=[[:space:]]*[0-9]+[[:space:]]*;' "$TEMP_SOURCE" || true)

if [[ "$chunk_count" -ne 1 ]]; then
    echo "Error: expected exactly one chunk_locdim assignment, found $chunk_count." >&2
    exit 1
fi

if [[ "$block_count" -lt 1 ]]; then
    echo "Error: no blocksize assignment was found." >&2
    exit 1
fi

# CSV header
{
    printf 'chunk_locdim,blocksize'
    for ((run = 1; run <= REPEATS; run++)); do
        printf ',run_%d_s' "$run"
    done
    printf ',mean_s,stddev_s\n'
} > "$RESULT_CSV"

echo "Benchmark started"
echo "Source      : $SOURCE"
echo "Repeats     : $REPEATS"
echo "CSV output  : $RESULT_CSV"
echo

# ============================================================
# BENCHMARK LOOP
# ============================================================

for chunk in "${CHUNK_SIZES[@]}"; do
    for block in "${BLOCK_SIZES[@]}"; do

        # CUDA permits at most 1024 threads per block on ordinary launches.
        if (( block <= 0 || block > 1024 )); then
            echo "Skip: invalid blocksize=$block" >&2
            continue
        fi

        echo "=== chunk_locdim=$chunk, blocksize=$block ==="

        # Start each combination from the original source.
        cp "$SOURCE" "$TEMP_SOURCE"

        # Replace the single chunk_locdim assignment.
        perl -0pi -e \
            "s/^([ \t]*int[ \t]+chunk_locdim[ \t]*=[ \t]*)[0-9]+([ \t]*;)/\${1}${chunk}\${2}/m" \
            "$TEMP_SOURCE"

        # Replace every blocksize assignment, since the source contains
        # more than one blocksize declaration for different output paths.
        perl -0pi -e \
            "s/^([ \t]*int[ \t]+blocksize[ \t]*=[ \t]*)[0-9]+([ \t]*;)/\${1}${block}\${2}/mg" \
            "$TEMP_SOURCE"

        # Confirm that replacement succeeded.
        if ! grep -Eq "^[[:space:]]*int[[:space:]]+chunk_locdim[[:space:]]*=[[:space:]]*${chunk}[[:space:]]*;" "$TEMP_SOURCE"; then
            echo "Error: failed to set chunk_locdim=$chunk" >&2
            exit 1
        fi

        replaced_blocks=$(grep -Ec "^[[:space:]]*int[[:space:]]+blocksize[[:space:]]*=[[:space:]]*${block}[[:space:]]*;" "$TEMP_SOURCE" || true)
        if [[ "$replaced_blocks" -ne "$block_count" ]]; then
            echo "Error: expected to replace $block_count blocksize assignments, replaced $replaced_blocks." >&2
            exit 1
        fi

        EXE="$BUILD_DIR/twiceg_chunk${chunk}_block${block}"

        echo "Compiling..."
        "$NVCC" "${NVCC_FLAGS[@]}" "$TEMP_SOURCE" -o "$EXE"

        times=()

        for ((run = 1; run <= REPEATS; run++)); do
            log="$BUILD_DIR/run_chunk${chunk}_block${block}_${run}.log"

            # Use nanosecond timestamps to measure complete wall-clock time.
            start_ns=$(date +%s%N)

            if [[ "$SAVE_LOGS" -eq 1 ]]; then
                "$EXE" "$CONF" "$WIND" "$LOC" >"$log" 2>&1
            else
                "$EXE" "$CONF" "$WIND" "$LOC" >/dev/null 2>&1
            fi

            end_ns=$(date +%s%N)

            elapsed=$(awk -v start="$start_ns" -v end="$end_ns" \
                'BEGIN { printf "%.6f", (end-start)/1000000000 }')

            times+=("$elapsed")
            printf '  run %d: %s s\n' "$run" "$elapsed"
        done

        # Calculate arithmetic mean and sample standard deviation (n-1).
        stats=$(printf '%s\n' "${times[@]}" | awk '
            {
                x[NR] = $1
                sum += $1
            }
            END {
                n = NR
                mean = sum / n

                if (n > 1) {
                    for (i = 1; i <= n; i++) {
                        d = x[i] - mean
                        ss += d * d
                    }
                    sd = sqrt(ss / (n - 1))
                } else {
                    sd = 0
                }

                printf "%.6f,%.6f", mean, sd
            }
        ')

        mean=${stats%,*}
        stddev=${stats#*,}

        {
            printf '%d,%d' "$chunk" "$block"
            for t in "${times[@]}"; do
                printf ',%s' "$t"
            done
            printf ',%s,%s\n' "$mean" "$stddev"
        } >> "$RESULT_CSV"

        printf '  mean: %s s\n' "$mean"
        printf '  sd  : %s s\n\n' "$stddev"
    done
done

echo "Benchmark completed."
echo "Results: $RESULT_CSV"
