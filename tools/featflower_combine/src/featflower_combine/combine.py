"""Core logic for combining partitioned FeatFloWer .dmp field files."""

import os
import re
from shutil import copyfile


def _custom_compare(item):
    return item[0]


def read_partitioned_field(elem_entries, file_name):
    """Read a single partitioned .dmp field file and append entries to elem_entries.

    Returns (header_line, header_info).
    """
    header_info = {}
    header_line = ""

    with open(file_name, "r") as f:
        header = f.readline()
        header_line = header
        split_header = header.strip().split(",")
        header_info = {pair[0]: pair[1] for pair in (l.split(":") for l in split_header)}

        while True:
            line = f.readline()
            if not line:
                break

            element_idx = int(line)
            element_components = []
            entry = (element_idx, element_components)

            for i in range(int(header_info["Components"])):
                line = f.readline()
                if not line:
                    break
                values = line.strip().split(" ")
                entry[1].append(values)

            elem_entries.append(entry)

    return header_line, header_info


def write_combined_field(elem_entries, header, components, file_name):
    """Write combined element entries to a single .dmp file."""
    print("Writing combined files: " + str(file_name))

    with open(file_name, "w") as f:
        f.write(header)
        for e in elem_entries:
            f.write(str(e[0]) + "\n")
            for i in range(components):
                f.write(" ".join(e[1][i]) + "\n")


def combine_field(nprocs, field_name, path, out_idx):
    """Read per-processor .dmp files for one field and write a combined output file."""
    element_entries = []

    for i in range(1, nprocs + 1):
        file_to_read = os.path.join(path, "processor_" + str(i), str(out_idx), field_name)
        header_line, header_info = read_partitioned_field(element_entries, file_to_read)

    element_entries.sort(key=_custom_compare)

    del header_info["DofsTotal"]
    header_info["NEL"] = len(element_entries)

    desired_order = ["DofsInElement", "NEL", "Name", "Format", "OutputLevel", "FeSpace", "Version", "Components"]
    header_parts = []
    for key in desired_order:
        if key in header_info:
            header_parts.append(str(key) + ":" + str(header_info[key]))

    header_mod = ",".join(header_parts) + "\n"

    out_dir = os.path.join(path, str(out_idx))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_name = os.path.join(out_dir, field_name)
    write_combined_field(element_entries, header_mod, int(header_info["Components"]), out_name)


def get_dump_fields(nprocs, path, out_idx):
    """Return list of .dmp field file names found in processor_1/<out_idx>/ (excluding time.dmp)."""
    dump_path = os.path.join(path, "processor_1", str(out_idx))
    print("Dump load path: " + dump_path)

    files = []
    for item in os.listdir(dump_path):
        if re.search(r"\w+\.dmp", item) and item != "time.dmp":
            files.append(item)

    return files


def process_dump_dir(nprocs, path, out_idx):
    """Combine all partitioned .dmp fields in path for the given dump index."""
    print("Found " + str(nprocs) + " processor_ directories")
    print("Dump directory: " + str(os.listdir(path)))
    print("Number of processors: " + str(nprocs))

    dump_path = os.path.join(path, "processor_1", str(out_idx))
    print("Backup files: " + str(os.listdir(dump_path)))

    fields = get_dump_fields(nprocs, path, out_idx)

    for f in fields:
        combine_field(nprocs, f, path, out_idx)

    time_file = os.path.join(dump_path, "time.dmp")
    if os.path.exists(time_file):
        dest = os.path.join(path, str(out_idx), "time.dmp")
        print("Writing time file: " + dest)
        copyfile(time_file, dest)


def count_processors(path):
    """Count processor_* subdirectories under path."""
    count = 0
    for entry in os.listdir(path):
        if re.search(r"processor_\d+", entry):
            count += 1
    return count
