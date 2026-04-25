#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_manifest(path: Path):
    data = {}
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            data[key.strip()] = value.strip()
    return data


def validate_dump(dump_dir: Path):
    manifest = read_manifest(dump_dir / "manifest.txt")
    owner_count = int(manifest["q2_owner_count"])
    coarse_elements = int(manifest["coarse_elements"])
    q2_slots = int(manifest["q2_slots_per_coarse"])

    ownership_rows = []
    with (dump_dir / "q2_ownership.csv").open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            ownership_rows.append(row)

    expected_rows = coarse_elements * q2_slots
    if len(ownership_rows) != expected_rows:
        raise SystemExit(f"ownership row count mismatch: got {len(ownership_rows)}, expected {expected_rows}")

    seen_slots = set()
    owner_rows = {}
    dup_hist = {}
    for row in ownership_rows:
        key = (int(row["coarse_elem_global"]), int(row["slot"]))
        if key in seen_slots:
            raise SystemExit(f"duplicate ownership row for coarse/slot {key}")
        seen_slots.add(key)

        owner_index = int(row["owner_index"])
        is_owner = int(row["is_owner"])
        dup_count = int(row["duplicate_count"])
        dup_hist[dup_count] = dup_hist.get(dup_count, 0) + 1
        if is_owner:
            owner_rows.setdefault(owner_index, 0)
            owner_rows[owner_index] += 1

    if len(owner_rows) != owner_count:
        raise SystemExit(f"owner count mismatch: metadata has {len(owner_rows)}, manifest says {owner_count}")

    bad_owner_rows = [idx for idx, count in owner_rows.items() if count != 1]
    if bad_owner_rows:
        raise SystemExit(f"owner rows not unique for owner_index values: {bad_owner_rows[:10]}")

    audit_rows = []
    with (dump_dir / "q2_ownership_audit.csv").open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            audit_rows.append(row)

    if len(audit_rows) != owner_count:
        raise SystemExit(f"audit owner count mismatch: got {len(audit_rows)}, expected {owner_count}")

    q2_fields = sorted(
        p for p in dump_dir.glob("*.csv")
        if p.name not in {"q2_ownership.csv", "q2_ownership_audit.csv", "pressure.csv"}
    )

    payload_counts = {}
    for field in q2_fields:
        with field.open() as f:
            reader = csv.reader(f)
            header = next(reader, None)
            rows = list(reader)
        payload_counts[field.name] = len(rows)
        if len(rows) != owner_count:
            raise SystemExit(f"{field.name}: payload row count mismatch: got {len(rows)}, expected {owner_count}")

    print(f"validated {dump_dir}")
    print(f"  coarse elements      = {coarse_elements}")
    print(f"  q2 slots / coarse    = {q2_slots}")
    print(f"  q2 owner count       = {owner_count}")
    print(f"  duplicate histogram  = {dup_hist}")
    for name, count in payload_counts.items():
        print(f"  {name:<20} = {count}")


def compare_dirs(dir_a: Path, dir_b: Path):
    files_a = sorted(p.name for p in dir_a.iterdir() if p.is_file())
    files_b = sorted(p.name for p in dir_b.iterdir() if p.is_file())
    if files_a != files_b:
        raise SystemExit(f"file set mismatch:\n  {files_a}\n  {files_b}")

    for name in files_a:
        a = (dir_a / name).read_bytes()
        b = (dir_b / name).read_bytes()
        if a != b:
            raise SystemExit(f"file differs: {name}")

    print(f"directories match exactly:\n  {dir_a}\n  {dir_b}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dump_dir", type=Path)
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    validate_dump(args.dump_dir)
    if args.compare:
        compare_dirs(args.dump_dir, args.compare)


if __name__ == "__main__":
    main()
