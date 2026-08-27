#!/usr/bin/env python3

"""Aggregate MPI-PRF ROMIO benchmark samples into TSV and Markdown tables."""

from __future__ import annotations

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


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
                row["ranks"],
                row["profile"],
                row["component"],
                row["cb_buffer"],
                row["ds_write"],
            )
            groups[key].append(row)

    aggregate_path = output_directory / "profile_summary.tsv"
    markdown_path = output_directory / "profile_summary.md"
    fields = (
        "ranks",
        "profile",
        "component",
        "cb_buffer",
        "ds_write",
        "samples",
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
        elapsed = [float(row["checkpoint_s"]) for row in rows]
        throughput = [float(row["throughput_mib_s"]) for row in rows]
        write_success = sum(row["write_status"] == "ok" for row in rows)
        validation = ",".join(sorted({row["roundtrip_status"] for row in rows}))
        aggregate_rows.append(
            {
                "ranks": key[0],
                "profile": key[1],
                "component": key[2],
                "cb_buffer": key[3],
                "ds_write": key[4],
                "samples": str(len(rows)),
                "checkpoint_mean_s": f"{statistics.mean(elapsed):.6f}",
                "checkpoint_stddev_s": (
                    f"{statistics.stdev(elapsed):.6f}" if len(elapsed) > 1 else "0.000000"
                ),
                "checkpoint_min_s": f"{min(elapsed):.6f}",
                "checkpoint_max_s": f"{max(elapsed):.6f}",
                "throughput_mean_mib_s": f"{statistics.mean(throughput):.2f}",
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
            "| Ranks | Profile | ROMIO settings | Samples | Checkpoint mean ± stddev (s) "
            "| Min–max (s) | Mean MiB/s | Write success | Validation |\n"
        )
        handle.write("| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |\n")
        for row in aggregate_rows:
            settings = (
                f"io={row['component']}; cb={row['cb_buffer']}; "
                f"ds_write={row['ds_write']}"
            )
            handle.write(
                f"| {row['ranks']} | {row['profile']} | {settings} | {row['samples']} "
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
