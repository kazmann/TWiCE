#!/usr/bin/env python3
import math

file1 = "massloading.txt"
file2 = "massloading_in_loc.txt"
outfile = "compare_massloading.txt"

coord_ok = True
max_abs_rel = 0.0
max_info = None

with open(file1) as f1, open(file2) as f2, open(outfile, "w") as out:
    header1 = f1.readline().split()
    header2 = f2.readline().split()

    if header1 != header2:
        coord_ok = False

    out.write("\t".join(header1) + "\n")

    line_no = 1

    for line1, line2 in zip(f1, f2):
        line_no += 1

        a = line1.split()
        b = line2.split()

        if len(a) != len(b):
            coord_ok = False
            continue

        row_out = []

        # columns 1–3: check equality and output original values
        for col in range(3):
            if a[col] != b[col]:
                coord_ok = False
            row_out.append(a[col])

        # columns 6 onward: relative difference
        for col in range(3, len(a)):
            v1 = float(a[col])
            v2 = float(b[col])

            if v2 == 0.0:
                if v1 == 0.0:
                    rel = 0.0
                else:
                    rel = math.inf
            else:
                rel = (v1 - v2) / v2

            row_out.append(f"{rel:.10e}")

            if math.isfinite(rel) and abs(rel) > max_abs_rel:
                max_abs_rel = abs(rel)
                max_info = (line_no, col + 1, header1[col], rel, v1, v2)

        out.write("\t".join(row_out) + "\n")

    # check extra lines
    if f1.readline() or f2.readline():
        coord_ok = False

    out.write("\n")
    out.write("# Summary\n")
    out.write(f"# Coordinate_columns_match: {coord_ok}\n")

    if max_info is None:
        out.write("# Max_relative_difference: 0\n")
    else:
        line_no, col_no, col_name, rel, v1, v2 = max_info
        out.write(f"# Max_abs_relative_difference: {max_abs_rel:.10e}\n")
        out.write(f"# Signed_relative_difference: {rel:.10e}\n")
        out.write(f"# Line: {line_no}\n")
        out.write(f"# Column: {col_no} ({col_name})\n")
        out.write(f"# massloading.txt: {v1:.10e}\n")
        out.write(f"# massloading_in_loc.txt: {v2:.10e}\n")

print(f"Coordinate columns match: {coord_ok}")

if max_info is None:
    print("Max relative difference: 0")
else:
    line_no, col_no, col_name, rel, v1, v2 = max_info
    print(f"Max abs relative difference: {max_abs_rel:.10e}")
    print(f"Signed relative difference: {rel:.10e}")
    print(f"Line: {line_no}")
    print(f"Column: {col_no} ({col_name})")

print(f"Written: {outfile}")