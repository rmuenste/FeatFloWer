#!/usr/bin/env python3
"""
Generate sphere positions along a centerline and compute the volume fraction.

Mirrors the C++ functions readVectorsFromFile and generatePointsAlongCenterline
from libs/pe/pe/interface/geometry_utils.h, and the volume fraction computation
from setupATCSerial in libs/pe/pe/interface/sim_setup_serial.h.

Subcommands
-----------
  rings   Full-centerline seeding with margin-based exclusion zones (mirrors C++).
  region  Seeding in a user-defined arc-length fraction window; the number of
          concentric rings is derived automatically from the vessel-radius constraint.
"""

import argparse
import math
import sys
import numpy as np


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_vectors_from_file(filename):
    """Parse a text file of whitespace-separated x y z triplets."""
    vertices = []
    with open(filename) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                print(f"Warning: line {lineno} has fewer than 3 values, skipping.",
                      file=sys.stderr)
                continue
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                vertices.append(np.array([x, y, z], dtype=float))
            except ValueError:
                print(f"Warning: line {lineno} could not be parsed, skipping.",
                      file=sys.stderr)
    return vertices


# ---------------------------------------------------------------------------
# Shared geometry helpers
# ---------------------------------------------------------------------------

def _build_cumulative(vertices):
    """Return (cumulative arc-lengths list, total curve length)."""
    cumulative = [0.0]
    for i in range(len(vertices) - 1):
        cumulative.append(cumulative[-1] + np.linalg.norm(vertices[i + 1] - vertices[i]))
    return cumulative, cumulative[-1]


def _point_on_curve(vertices, cumulative, s):
    """
    Interpolate a position and its tangent direction at arc-length s.
    Returns (point, edge_dir) or (None, None) if s is beyond the curve end.
    """
    num_edges = len(vertices) - 1
    edge_index = 0
    while edge_index < num_edges and s > cumulative[edge_index + 1]:
        edge_index += 1
    if edge_index >= num_edges:
        return None, None

    t = (s - cumulative[edge_index]) / (cumulative[edge_index + 1] - cumulative[edge_index])
    v1, v2 = vertices[edge_index], vertices[edge_index + 1]
    return v1 + (v2 - v1) * t, v2 - v1


def _place_rings(point, edge_dir, sphere_radius, dt, num_rings):
    """
    Return sphere centre positions for all concentric rings around a centerline point.
    Mirrors the inner loop of generatePointsAlongCenterline.
    """
    edge_dir_n = edge_dir / np.linalg.norm(edge_dir)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(edge_dir_n, ref)) > 0.999:
        ref = np.array([0.0, 1.0, 0.0])

    u = np.cross(edge_dir, ref)
    u /= np.linalg.norm(u)
    v = np.cross(edge_dir, u)
    v /= np.linalg.norm(v)

    positions = []
    for j in range(num_rings):
        circle_radius = sphere_radius + dt + j * (2.0 * sphere_radius + dt)
        max_spheres = max(1, int(2.0 * math.pi * circle_radius / (2.0 * sphere_radius)) - 1)
        theta_step = 2.0 * math.pi / max_spheres
        for i in range(max_spheres):
            theta = i * theta_step
            offset = (math.cos(theta) * u + math.sin(theta) * v) * circle_radius
            positions.append(point + offset)
    return positions


# ---------------------------------------------------------------------------
# Generation functions
# ---------------------------------------------------------------------------

def generate_points_along_centerline(vertices, sphere_radius, dt=-1.0,
                                     num_rings=4, num_steps=20, margin_steps=4):
    """
    Generate sphere positions in concentric rings along the full centerline.
    Mirrors the C++ generatePointsAlongCenterline signature exactly.
    """
    if dt < 0.0:
        dt = sphere_radius

    if len(vertices) < 2:
        raise ValueError("Need at least 2 vertices to define a centerline.")

    cumulative, curve_length = _build_cumulative(vertices)
    ds = curve_length / num_steps

    positions = []
    s = margin_steps * ds + 0.2 * ds
    while s <= curve_length - margin_steps * ds:
        point, edge_dir = _point_on_curve(vertices, cumulative, s)
        if point is None:
            break
        positions.extend(_place_rings(point, edge_dir, sphere_radius, dt, num_rings))
        s += ds

    return positions


def generate_points_in_region(vertices, sphere_radius, dt=-1.0,
                               start_frac=0.5, end_frac=0.95,
                               n_cross_sections=10, max_vessel_radius=0.25):
    """
    Seed spheres in a sub-region of the centerline defined by arc-length fractions.

    The number of concentric rings is derived automatically so that the outermost
    ring's centre radius satisfies:

        ring_radius_last  <=  max_vessel_radius - 2 * sphere_radius

    n_cross_sections planes are placed uniformly from start_frac to end_frac of
    the total curve length (both endpoints included).

    Returns
    -------
    positions       : list of np.ndarray(3,) — sphere centre positions
    num_rings       : int   — number of rings actually used
    last_ring_rad   : float — outermost ring centre radius
    curve_length    : float — total centerline length (for informational output)
    """
    if dt < 0.0:
        dt = sphere_radius

    if len(vertices) < 2:
        raise ValueError("Need at least 2 vertices to define a centerline.")
    if not (0.0 <= start_frac < end_frac <= 1.0):
        raise ValueError("Require 0 <= start_frac < end_frac <= 1.")

    cumulative, curve_length = _build_cumulative(vertices)

    s_start = start_frac * curve_length
    s_end   = end_frac   * curve_length

    # Derive num_rings from vessel-radius constraint:
    #   ring_radius_j = sphere_radius + dt + j*(2*sphere_radius + dt) <= max_vessel_radius - 2*r
    max_ring_radius = max_vessel_radius - 2.0 * sphere_radius
    if max_ring_radius <= 0.0:
        raise ValueError(
            f"max_vessel_radius ({max_vessel_radius}) must be > 2 * sphere_radius "
            f"({2 * sphere_radius:.6g}).")

    ring_step = 2.0 * sphere_radius + dt          # radial increment per ring
    first_ring_radius = sphere_radius + dt         # ring 0 centre radius
    if first_ring_radius > max_ring_radius:
        raise ValueError(
            f"Even ring 0 (radius {first_ring_radius:.6g}) exceeds the vessel constraint "
            f"(max_ring_radius = {max_ring_radius:.6g}). Increase max_vessel_radius or "
            f"decrease sphere_radius / dt.")

    num_rings = int((max_ring_radius - first_ring_radius) / ring_step) + 1
    num_rings = max(1, num_rings)

    # Ensure the last ring strictly satisfies the constraint (guard against fp rounding)
    last_ring_rad = first_ring_radius + (num_rings - 1) * ring_step
    if last_ring_rad > max_ring_radius:
        num_rings -= 1
        last_ring_rad = first_ring_radius + (num_rings - 1) * ring_step

    # Evenly spaced cross-sections from s_start to s_end (inclusive)
    ds = (s_end - s_start) / max(n_cross_sections - 1, 1)

    positions = []
    for i in range(n_cross_sections):
        s = s_start + i * ds
        point, edge_dir = _point_on_curve(vertices, cumulative, s)
        if point is None:
            break
        positions.extend(_place_rings(point, edge_dir, sphere_radius, dt, num_rings))

    return positions, num_rings, last_ring_rad, curve_length


# ---------------------------------------------------------------------------
# VTK output
# ---------------------------------------------------------------------------

def write_vtk(filename, positions, sphere_radius, centerline_vertices):
    """
    Write sphere positions and the centerline as a legacy VTK PolyData file.

    Sphere centres are stored as VERTICES (one cell each) with a 'radius' point
    scalar so ParaView's Glyph filter renders each sphere at its true size.
    The centerline is stored as a single LINES polyline for spatial context.
    """
    n_spheres = len(positions)
    n_cl      = len(centerline_vertices)
    n_total   = n_spheres + n_cl

    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Centerline sphere positions\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        f.write(f"POINTS {n_total} float\n")
        for p in positions:
            f.write(f"{p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n")
        for v in centerline_vertices:
            f.write(f"{v[0]:.8f} {v[1]:.8f} {v[2]:.8f}\n")

        f.write(f"\nVERTICES {n_spheres} {2 * n_spheres}\n")
        for i in range(n_spheres):
            f.write(f"1 {i}\n")

        if n_cl >= 2:
            f.write(f"\nLINES 1 {n_cl + 1}\n")
            f.write(str(n_cl))
            for i in range(n_cl):
                f.write(f" {n_spheres + i}")
            f.write("\n")

        f.write(f"\nPOINT_DATA {n_total}\n")
        f.write("SCALARS radius float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for _ in positions:
            f.write(f"{sphere_radius:.8f}\n")
        for _ in centerline_vertices:
            f.write("0.0\n")


# ---------------------------------------------------------------------------
# Shared CLI helpers
# ---------------------------------------------------------------------------

def _add_common_args(p):
    p.add_argument("file",
                   help="Centerline vertices file (x y z per line)")
    p.add_argument("-r", "--sphere-radius", type=float, required=True, metavar="R",
                   help="Sphere radius")
    p.add_argument("-V", "--domain-vol", type=float, required=True, metavar="VOL",
                   help="Domain volume for volume fraction computation")
    p.add_argument("--dt", type=float, default=-1.0,
                   help="Gap from sphere surface to ring centre "
                        "(default: -1.0 → resolved to sphere_radius)")
    p.add_argument("-o", "--output", default=None, metavar="FILE",
                   help="Write sphere positions as x y z text")
    p.add_argument("--vtk", default=None, metavar="FILE.vtk",
                   help="Write VTK PolyData file for ParaView visualization")


def _print_vf(n_particles, sphere_radius, domain_vol):
    part_vol = (4.0 / 3.0) * math.pi * sphere_radius ** 3
    vf = (n_particles * part_vol) / domain_vol * 100.0
    print(f"Particles generated        : {n_particles}")
    print(f"Particle volume            : {part_vol:.6g}")
    print(f"Domain volume              : {domain_vol}")
    print(f"Volume fraction            : {vf:.4f} %")
    return part_vol, vf


def _write_outputs(args, positions, vertices):
    if args.output:
        with open(args.output, "w") as f:
            for p in positions:
                f.write(f"{p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n")
        print(f"Positions written to       : {args.output}")
    if args.vtk:
        write_vtk(args.vtk, positions, args.sphere_radius, vertices)
        print(f"VTK file written to        : {args.vtk}")
        print("  -> In ParaView: Glyph filter, Glyph Type = Sphere,")
        print("                  Scale Array = radius, Scale Factor = 1.0")


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def cmd_rings(args):
    vertices = read_vectors_from_file(args.file)
    if len(vertices) < 2:
        print(f"Error: need at least 2 vertices, got {len(vertices)}.", file=sys.stderr)
        sys.exit(1)

    positions = generate_points_along_centerline(
        vertices,
        sphere_radius=args.sphere_radius,
        dt=args.dt,
        num_rings=args.num_rings,
        num_steps=args.num_steps,
        margin_steps=args.margin_steps,
    )

    dt_eff = args.dt if args.dt >= 0 else args.sphere_radius
    print(f"Centerline vertices loaded : {len(vertices)}")
    print(f"Sphere radius              : {args.sphere_radius}")
    print(f"dt (ring gap)              : {dt_eff}")
    print(f"num_rings                  : {args.num_rings}")
    print(f"num_steps                  : {args.num_steps}")
    print(f"margin_steps               : {args.margin_steps}")
    _print_vf(len(positions), args.sphere_radius, args.domain_vol)
    _write_outputs(args, positions, vertices)


def cmd_region(args):
    vertices = read_vectors_from_file(args.file)
    if len(vertices) < 2:
        print(f"Error: need at least 2 vertices, got {len(vertices)}.", file=sys.stderr)
        sys.exit(1)

    positions, num_rings, last_ring_rad, curve_length = generate_points_in_region(
        vertices,
        sphere_radius=args.sphere_radius,
        dt=args.dt,
        start_frac=args.start_frac,
        end_frac=args.end_frac,
        n_cross_sections=args.n_cross_sections,
        max_vessel_radius=args.max_vessel_radius,
    )

    dt_eff       = args.dt if args.dt >= 0 else args.sphere_radius
    s_start_abs  = args.start_frac * curve_length
    s_end_abs    = args.end_frac   * curve_length
    max_ring_rad = args.max_vessel_radius - 2.0 * args.sphere_radius

    print(f"Centerline vertices loaded : {len(vertices)}")
    print(f"Total centerline length    : {curve_length:.6g}")
    print(f"Sphere radius              : {args.sphere_radius}")
    print(f"dt (ring gap)              : {dt_eff}")
    print(f"start_frac / end_frac      : {args.start_frac} / {args.end_frac}  "
          f"(s = {s_start_abs:.4f} .. {s_end_abs:.4f})")
    print(f"n_cross_sections           : {args.n_cross_sections}")
    print(f"max_vessel_radius          : {args.max_vessel_radius}")
    print(f"Vessel constraint (max r)  : {max_ring_rad:.6g}  "
          f"(= max_vessel_radius - 2*sphere_radius)")
    print(f"Rings derived              : {num_rings}  "
          f"(outermost ring radius = {last_ring_rad:.6g})")
    _print_vf(len(positions), args.sphere_radius, args.domain_vol)
    _write_outputs(args, positions, vertices)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate sphere positions along a centerline and compute volume fraction.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # -- rings subcommand (mirrors C++ generatePointsAlongCenterline) ----------
    p_rings = sub.add_parser(
        "rings",
        help="Full-centerline seeding with margin exclusion zones (mirrors C++).",
        description="Seed spheres along the entire centerline using the same logic as "
                    "generatePointsAlongCenterline() in geometry_utils.h.",
    )
    _add_common_args(p_rings)
    p_rings.add_argument("--num-rings", type=int, default=4,
                         help="Concentric rings per cross-section (default: 4)")
    p_rings.add_argument("--num-steps", type=int, default=20,
                         help="Arc-length divisions along the full curve (default: 20)")
    p_rings.add_argument("--margin-steps", type=int, default=4,
                         help="Stations skipped at each end of the curve (default: 4)")
    p_rings.set_defaults(func=cmd_rings)

    # -- region subcommand ----------------------------------------------------
    p_region = sub.add_parser(
        "region",
        help="Seed spheres in a sub-region of the centerline with vessel-radius constraint.",
        description="Place sphere rings between start_frac and end_frac of the total "
                    "centerline arc-length. The number of rings is derived automatically "
                    "so the outermost ring fits within max_vessel_radius - 2*sphere_radius.",
    )
    _add_common_args(p_region)
    p_region.add_argument("--start-frac", type=float, default=0.5, metavar="FRAC",
                          help="Start of seeding region as fraction of total length (default: 0.5)")
    p_region.add_argument("--end-frac", type=float, default=0.95, metavar="FRAC",
                          help="End of seeding region as fraction of total length (default: 0.95)")
    p_region.add_argument("--n-cross-sections", type=int, default=10, metavar="N",
                          help="Number of cross-section planes in the region (default: 10)")
    p_region.add_argument("--max-vessel-radius", type=float, default=0.25, metavar="R",
                          help="Local vessel radius; outermost ring centre <= R - 2*sphere_radius "
                               "(default: 0.25)")
    p_region.set_defaults(func=cmd_region)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
