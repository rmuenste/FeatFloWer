"""CLI entry point for featflower-combine."""

import argparse
import sys

from featflower_combine.combine import count_processors, process_dump_dir


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="featflower-combine",
        description="Combine partitioned FeatFloWer .dmp field files into single output files.",
        epilog="Example:\n  featflower-combine --dump-path=_dump --idx=2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-d", "--dump-path",
        default="./_dump",
        help="Path to the dump folder (default: ./_dump)",
    )
    parser.add_argument(
        "-i", "--idx",
        default="1",
        help="Index of the dump folder (default: 1)",
    )

    args = parser.parse_args(argv)

    nprocs = count_processors(args.dump_path)
    if nprocs == 0:
        print(f"No processor_* directories found in '{args.dump_path}'", file=sys.stderr)
        return 1

    process_dump_dir(nprocs, args.dump_path, args.idx)
    return 0
