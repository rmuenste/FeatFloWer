import argparse
import sys
from .part_main import mkdir, checkParameters, MainProcess


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="featflower-partition",
        description="Partition a FeatFloWer mesh for parallel execution (METIS 5).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  featflower-partition 4 1 1 NEWFAC testmeshes/2D_FAC/2Dbench.prj"
        ),
    )
    parser.add_argument("NPart",       type=int,  help="Number of partitions (>=1)")
    parser.add_argument("Method",      type=int,  help="Partitioning method (1=recursive, 2/3=kway, negative=axis-based)")
    parser.add_argument("NSubPart",    type=int,  help="Number of subgrids (>=1)")
    parser.add_argument("MeshName",               help="Output mesh name (used as directory key under _mesh/)")
    parser.add_argument("ProjectFile",            help="Path to .prj project file")
    parser.add_argument(
        "--partition-format",
        choices=("legacy", "json"),
        default="legacy",
        help="Partition output format (default: legacy)",
    )

    args = parser.parse_args(argv)

    # Re-use the existing checkParameters convention; build a params list that
    # matches the expected format: [prog, NPart, Method, NSubPart, MeshName, ProjectFile]
    params = [
        parser.prog,
        str(args.NPart),
        str(args.Method),
        str(args.NSubPart),
        args.MeshName,
        args.ProjectFile,
        "--partition-format=%s" % args.partition_format,
    ]

    NPart, PartMethod, NSubPart, MeshName, ProjektFile, partition_format = checkParameters(params)

    mkdir("_mesh")
    MainProcess(NPart, PartMethod, NSubPart, MeshName, ProjektFile,
                partition_format=partition_format)
