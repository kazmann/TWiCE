#!/bin/bash
set -e

OUT_SUMMARY="result_comparison.txt"

# ttlmassloading(kg/sq-m) の列番号
# 例: ファイル上で5列目なら 5
MASS_COL=5

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

    new_massloading="massloading.${res}.txt"
    mv massloading.txt "$new_massloading"

    echo "${res} ${input_points} ${elapsed} ${output_points}" >> "$OUT_SUMMARY"
done

echo "Comparing common points against 100m result..."

python3 << PY
import glob

mass_col = int("${MASS_COL}") - 1

files = {}
for path in glob.glob("massloading.*.txt"):
    parts = path.split(".")
    res = int(parts[1])
    files[res] = path

if 100 not in files:
    raise SystemExit("ERROR: 100m massloading file not found")

def read_file(path):
    data = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.split()

            # ヘッダー行を飛ばす
            try:
                val = float(cols[mass_col])
            except ValueError:
                continue

            x = cols[0]
            y = cols[1]

            data[(x, y)] = val

    return data

base = read_file(files[100])

for res in sorted(files):
    if res == 100:
        continue

    outname = f"massloading_compare_{res}m_vs_100m.txt"
    target = read_file(files[res])

    with open(outname, "w") as out:
        out.write("resolution x y mass_100 mass_res abs_diff rel_diff\\n")

        n_common = 0
        max_abs = 0.0
        max_rel = 0.0

        for key, val_res in target.items():
            if key not in base:
                continue

            val_100 = base[key]
            abs_diff = abs(val_100 - val_res)
            rel_diff = abs_diff / max(abs(val_100), 1.0e-30)

            n_common += 1
            max_abs = max(max_abs, abs_diff)
            max_rel = max(max_rel, rel_diff)

            out.write(
                f"{res} {key[0]} {key[1]} "
                f"{val_100:.8e} {val_res:.8e} "
                f"{abs_diff:.8e} {rel_diff:.8e}\\n"
            )

    print(
        f"{res}m: common={n_common}, "
        f"max_abs_diff={max_abs:.8e}, max_rel_diff={max_rel:.8e}"
    )
    print(f"written: {outname}")

print("Comparison files written: massloading_compare_*m_vs_100m.txt")
PY

echo "Done."
echo "Summary written to ${OUT_SUMMARY}"
echo "Comparison files written: massloading_compare_*m_vs_100m.txt"