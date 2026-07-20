#!/usr/bin/env python3
"""
Compare all columns in CPU and GPU output files from TWiCE.

The script:

1. Reads the headers of massloading_cpu.txt and massloading_gpu.txt.
2. Stops if the headers are different.
3. Stops if the numbers of data rows are different.
4. Compares every column from left to right.
5. Writes one comparison file for each column.
6. Reports the comparison summary to standard output.

For example, the column

    x(m)

is written to

    compare_x(m).txt

Each output file contains:

    CPU    GPU    log10(CPU/GPU)

First-time setup:
    chmod +x compare_CPUGPU.py

Run:
    ./compare_CPUGPU.py

or:
    python3 compare_CPUGPU.py
"""

import math
import sys
from pathlib import Path


CPU_FILE = Path("massloading_cpu.txt")
GPU_FILE = Path("massloading_gpu.txt")

OUTPUT_PREFIX = "compare_"


def read_table(file_path: Path) -> tuple[list[str], list[list[float]]]:
    """
    Read a whitespace-delimited text file.

    Blank lines and lines beginning with '#' are ignored.
    The first remaining line is treated as the header.
    """

    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")

    with file_path.open("r", encoding="utf-8") as file:
        lines = [
            line.strip()
            for line in file
            if line.strip() and not line.lstrip().startswith("#")
        ]

    if not lines:
        raise ValueError(f"Input file is empty: {file_path}")

    header = lines[0].split()

    if not header:
        raise ValueError(f"No header was found in: {file_path}")

    rows: list[list[float]] = []

    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split()

        if len(fields) != len(header):
            raise ValueError(
                f"Incorrect number of columns in {file_path}, "
                f"line {line_number}.\n"
                f"Expected: {len(header)}\n"
                f"Found:    {len(fields)}\n"
                f"Line:     {line}"
            )

        row: list[float] = []

        for column_index, field in enumerate(fields):
            try:
                value = float(field)
            except ValueError as error:
                column_name = header[column_index]
                raise ValueError(
                    f"Non-numeric value in {file_path}, "
                    f"line {line_number}, column '{column_name}':\n"
                    f"{field}"
                ) from error

            row.append(value)

        rows.append(row)

    return header, rows


def make_output_filename(column_name: str) -> Path:
    """
    Create an output filename from a column name.

    Characters that cannot safely be used in filenames are replaced.
    Parentheses and hyphens are retained.

    Example:
        x(m)                     -> compare_x(m).txt
        ttlmassloading(kg/sq-m)  ->
            compare_ttlmassloading(kg_per_sq-m).txt
    """

    safe_name = column_name

    replacements = {
        "/": "_per_",
        "\\": "_",
        ":": "_",
        "*": "_",
        "?": "_",
        '"': "_",
        "<": "_",
        ">": "_",
        "|": "_",
    }

    for old_character, new_text in replacements.items():
        safe_name = safe_name.replace(old_character, new_text)

    return Path(f"{OUTPUT_PREFIX}{safe_name}.txt")


def calculate_log_ratio(cpu: float, gpu: float) -> float | None:
    """
    Calculate log10(CPU/GPU).

    None is returned if either value is non-finite or not positive.
    """

    if not math.isfinite(cpu) or not math.isfinite(gpu):
        return None

    if cpu <= 0.0 or gpu <= 0.0:
        return None

    return math.log10(cpu / gpu)


def values_match(cpu: float, gpu: float) -> bool:
    """
    Determine whether the CPU and GPU values are exactly equal.

    Two NaN values are not considered equal.
    """

    if math.isnan(cpu) or math.isnan(gpu):
        return False

    return cpu == gpu


def format_value(value: float) -> str:
    """Format a floating-point value in scientific notation."""

    return f"{value:.10e}"


def compare_column(
    column_name: str,
    column_index: int,
    cpu_rows: list[list[float]],
    gpu_rows: list[list[float]],
) -> tuple[
    int,
    int,
    float | None,
    float | None,
    float | None,
    Path,
]:
    """
    Compare one column and write its comparison file.

    Returns:
        number of matching values,
        number of non-matching values,
        mean CPU value,
        mean GPU value,
        maximum absolute log10(CPU/GPU),
        output file path.
    """

    output_path = make_output_filename(column_name)

    match_count = 0
    mismatch_count = 0

    cpu_sum = 0.0
    gpu_sum = 0.0
    cpu_count = 0
    gpu_count = 0

    maximum_absolute_log_ratio: float | None = None

    with output_path.open("w", encoding="utf-8") as output_file:
        output_file.write("CPU\tGPU\tlog10(CPU/GPU)\n")

        for cpu_row, gpu_row in zip(cpu_rows, gpu_rows):
            cpu = cpu_row[column_index]
            gpu = gpu_row[column_index]

            if values_match(cpu, gpu):
                match_count += 1
            else:
                mismatch_count += 1

            if math.isfinite(cpu):
                cpu_sum += cpu
                cpu_count += 1

            if math.isfinite(gpu):
                gpu_sum += gpu
                gpu_count += 1

            log_ratio = calculate_log_ratio(cpu, gpu)

            if log_ratio is None:
                log_ratio_text = "NA"
            else:
                log_ratio_text = format_value(log_ratio)

                absolute_log_ratio = abs(log_ratio)

                if (
                    maximum_absolute_log_ratio is None
                    or absolute_log_ratio > maximum_absolute_log_ratio
                ):
                    maximum_absolute_log_ratio = absolute_log_ratio

            output_file.write(
                f"{format_value(cpu)}\t"
                f"{format_value(gpu)}\t"
                f"{log_ratio_text}\n"
            )

    mean_cpu = cpu_sum / cpu_count if cpu_count > 0 else None
    mean_gpu = gpu_sum / gpu_count if gpu_count > 0 else None

    return (
        match_count,
        mismatch_count,
        mean_cpu,
        mean_gpu,
        maximum_absolute_log_ratio,
        output_path,
    )


def main() -> None:
    """Compare all columns in the CPU and GPU output files."""

    cpu_header, cpu_rows = read_table(CPU_FILE)
    gpu_header, gpu_rows = read_table(GPU_FILE)

    if cpu_header != gpu_header:
        print("ERROR: CPU and GPU headers are different.", file=sys.stderr)
        print(file=sys.stderr)

        print(f"CPU header ({CPU_FILE}):", file=sys.stderr)
        print("  " + " | ".join(cpu_header), file=sys.stderr)
        print(file=sys.stderr)

        print(f"GPU header ({GPU_FILE}):", file=sys.stderr)
        print("  " + " | ".join(gpu_header), file=sys.stderr)

        sys.exit(1)

    if len(cpu_rows) != len(gpu_rows):
        raise ValueError(
            "The CPU and GPU files contain different numbers of data rows:\n"
            f"  CPU: {len(cpu_rows)}\n"
            f"  GPU: {len(gpu_rows)}"
        )

    print("CPU/GPU comparison")
    print(f"CPU file: {CPU_FILE}")
    print(f"GPU file: {GPU_FILE}")
    print(f"Rows:     {len(cpu_rows)}")
    print(f"Columns:  {len(cpu_header)}")
    print()

    print(
        f"{'Column':<32}"
        f"{'Match':>12}"
        f"{'Mismatch':>12}"
        f"{'Mean(CPU)':>18}"
        f"{'Mean(GPU)':>18}"
        f"{'Max |log10(CPU/GPU)|':>24}"
    )
    print("-" * 116)

    for column_index, column_name in enumerate(cpu_header):
        (
            match_count,
            mismatch_count,
            mean_cpu,
            mean_gpu,
            maximum_absolute_log_ratio,
            output_path,
        ) = compare_column(
            column_name,
            column_index,
            cpu_rows,
            gpu_rows,
        )

        mean_cpu_text = (
            "NA" if mean_cpu is None else f"{mean_cpu:.10e}"
        )

        mean_gpu_text = (
            "NA" if mean_gpu is None else f"{mean_gpu:.10e}"
        )

        maximum_text = (
            "NA"
            if maximum_absolute_log_ratio is None
            else f"{maximum_absolute_log_ratio:.10e}"
        )

        print(
            f"{column_name:<32}"
            f"{match_count:>12}"
            f"{mismatch_count:>12}"
            f"{mean_cpu_text:>18}"
            f"{mean_gpu_text:>18}"
            f"{maximum_text:>24}"
        )

        print(f"  -> {output_path}")

    print()
    print("Comparison completed.")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)