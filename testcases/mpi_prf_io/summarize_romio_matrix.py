#!/usr/bin/env python3

"""Aggregate MPI-PRF ROMIO benchmark samples into TSV and Markdown tables."""

from __future__ import annotations

import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def row_value(row: dict[str, str], name: str, legacy_name: str | None = None) -> str:
    """Read the current schema while accepting the original fixture schema."""

    value = row.get(name, "")
    if value:
        return value
    if legacy_name is not None:
        value = row.get(legacy_name, "")
        if value:
            return value
    return "unset"


def numeric_values(rows: list[dict[str, str]], name: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row.get(name, ""))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    return values


def format_stats(values: list[float], digits: int) -> tuple[str, str, str, str]:
    if not values:
        unavailable = "not_available"
        return unavailable, unavailable, unavailable, unavailable
    mean = f"{statistics.mean(values):.{digits}f}"
    stddev = (
        f"{statistics.stdev(values):.{digits}f}"
        if len(values) > 1
        else f"{0:.{digits}f}"
    )
    minimum = f"{min(values):.{digits}f}"
    maximum = f"{max(values):.{digits}f}"
    return mean, stddev, minimum, maximum


def main() -> int:
    if len(sys.argv) not in (2, 3):
        print(
            f"Usage: {Path(sys.argv[0]).name} SUMMARY.tsv [OUTPUT_DIRECTORY]",
            file=sys.stderr,
        )
        return 2

    source = Path(sys.argv[1])
    output_directory = Path(sys.argv[2]) if len(sys.argv) == 3 else source.parent
    output_directory.mkdir(parents=True, exist_ok=True)
    groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    with source.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (
                row.get("ranks", "unset"),
                row.get("profile", "unset"),
                row_value(row, "component"),
                row_value(row, "romio_hints"),
                row_value(row, "hints_source"),
                row_value(row, "cb_buffer_size", "cb_buffer"),
                row_value(row, "romio_ds_write", "ds_write"),
            )
            groups[key].append(row)

    aggregate_path = output_directory / "profile_summary.tsv"
    markdown_path = output_directory / "profile_summary.md"
    fields = (
        "ranks",
        "profile",
        "component",
        "romio_hints",
        "hints_source",
        "cb_buffer_size",
        "romio_ds_write",
        "samples",
        "timed_samples",
        "checkpoint_mean_s",
        "checkpoint_stddev_s",
        "checkpoint_min_s",
        "checkpoint_max_s",
        "throughput_mean_mib_s",
        "write_success",
        "validation",
    )

    aggregate_rows: list[dict[str, str]] = []
    for key, rows in groups.items():
        elapsed = numeric_values(rows, "checkpoint_s")
        throughput = numeric_values(rows, "throughput_mib_s")
        checkpoint_mean, checkpoint_stddev, checkpoint_min, checkpoint_max = format_stats(
            elapsed, 6
        )
        throughput_mean = (
            f"{statistics.mean(throughput):.2f}" if throughput else "not_available"
        )
        write_success = sum(row.get("write_status") == "ok" for row in rows)
        validation = ",".join(
            sorted({row.get("roundtrip_status") or "not_run" for row in rows})
        )
        aggregate_rows.append(
            {
                "ranks": key[0],
                "profile": key[1],
                "component": key[2],
                "romio_hints": key[3],
                "hints_source": key[4],
                "cb_buffer_size": key[5],
                "romio_ds_write": key[6],
                "samples": str(len(rows)),
                "timed_samples": str(len(elapsed)),
                "checkpoint_mean_s": checkpoint_mean,
                "checkpoint_stddev_s": checkpoint_stddev,
                "checkpoint_min_s": checkpoint_min,
                "checkpoint_max_s": checkpoint_max,
                "throughput_mean_mib_s": throughput_mean,
                "write_success": f"{write_success}/{len(rows)}",
                "validation": validation,
            }
        )
    aggregate_rows.sort(key=lambda row: (int(row["ranks"]), row["profile"]))

    with aggregate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(aggregate_rows)

    with markdown_path.open("w", encoding="utf-8") as handle:
        handle.write("# MPI-PRF ROMIO benchmark summary\n\n")
        handle.write(
            "| Ranks | Profile | ROMIO settings | Timed/total | Checkpoint mean ± stddev (s) "
            "| Min–max (s) | Mean MiB/s | Write success | Validation |\n"
        )
        handle.write("| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |\n")
        for row in aggregate_rows:
            settings = (
                f"io={row['component']}; hints={row['romio_hints']}; "
                f"source={row['hints_source']}; cb_buffer_size={row['cb_buffer_size']}; "
                f"romio_ds_write={row['romio_ds_write']}"
            )
            handle.write(
                f"| {row['ranks']} | {row['profile']} | {settings} "
                f"| {row['timed_samples']}/{row['samples']} "
                f"| {row['checkpoint_mean_s']} ± {row['checkpoint_stddev_s']} "
                f"| {row['checkpoint_min_s']}–{row['checkpoint_max_s']} "
                f"| {row['throughput_mean_mib_s']} | {row['write_success']} "
                f"| {row['validation']} |\n"
            )

    print(f"Wrote {aggregate_path}")
    print(f"Wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
