#!/usr/bin/env python3
"""compute_I_r.py - Extended with visualization capabilities

Compute the radial index I_r for particle ensembles inside an Archimedes
tube crystallizer (ATC) slug. Supports two modes:

1. Analytic helix centre‑line
   python compute_I_r.py helix points.txt --R 2.835888 --phi0 5.474463 \
                                          --a -0.112848 --z0 -0.678574 --plot

2. Discretised centre‑line (polyline)
   python compute_I_r.py polyline points.txt centre_line_sampled.txt --plot

python3 ./compute_I_r_viz.py polyline point_set.xyz vertices.txt --binning particle_radius --particle_radius 0.0182 --plot --output "particle based" --title "12 RPM 1%wt dp=364"   
"""

import argparse
import math
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.patches as patches


# ----------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------
def read_xyz(filename: str) -> np.ndarray:
    """Read an ASCII file containing x y z columns."""
    return np.loadtxt(filename, dtype=float)


# ----------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------
def frenet_helix(t: float, R: float, phi0: float, a: float, z0: float):
    """Return centre‑line position c(t) and its first/second derivatives."""
    sphi = math.sin(t + phi0)
    cphi = math.cos(t + phi0)
    c = np.array([R * cphi, R * sphi, a * t + z0])
    c1 = np.array([-R * sphi, R * cphi, a])
    c2 = np.array([-R * cphi, -R * sphi, 0.0])
    return c, c1, c2


def newton_project_helix(
    p: np.ndarray,
    R: float,
    phi0: float,
    a: float,
    z0: float,
    tol: float = 1.0e-10,
    max_iter: int = 20,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Project a point p onto the helix using Newton iterations.

    Returns
    -------
    t : float
        Parameter where projection occurs.
    c : np.ndarray
        Centre‑line point at t.
    c1 : np.ndarray
        First derivative (tangent, non‑normalised) at t.
    """
    # Initial guess: angle in xy‑plane minus phase shift
    t = math.atan2(p[1], p[0]) - phi0

    for _ in range(max_iter):
        c, c1, c2 = frenet_helix(t, R, phi0, a, z0)
        g = np.dot(c1, p - c)                       # d/dt |p‑c|^2 / 2
        g_prime = -np.dot(c2, p - c) - np.dot(c1, c1)
        dt = -g / g_prime
        t += dt
        if abs(dt) < tol:
            break

    return t, c, c1


def radial_distance(
    p: np.ndarray, c: np.ndarray, tangent: np.ndarray
) -> float:
    """Distance of p to the centre‑line measured in the cross‑section."""
    T = tangent / np.linalg.norm(tangent)
    r_vec = p - c - np.dot(p - c, T) * T
    return np.linalg.norm(r_vec)


def closest_point_on_segment(
    p: np.ndarray, a: np.ndarray, b: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return closest point on segment AB to P and the segment tangent."""
    ab = b - a
    t = np.dot(p - a, ab) / np.dot(ab, ab)
    t = np.clip(t, 0.0, 1.0)
    c = a + t * ab
    tangent = ab / np.linalg.norm(ab)
    return c, tangent


# ----------------------------------------------------------------------
# Metric computation
# ----------------------------------------------------------------------
def radial_index_from_radii(
    radii: np.ndarray, Nshell: int = 15, R_inner: float | None = None, 
    binning_method: str = "equal_area", particle_radius: float = None
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Compute I_r following the definition in atc_metrics.md.
    
    Parameters
    ----------
    radii : np.ndarray
        Radial distances of particles from centerline
    Nshell : int
        Number of radial shells
    R_inner : float
        Tube inner radius (automatically determined if None)
    binning_method : str
        "equal_area" (default) or "particle_radius"
    particle_radius : float
        Particle radius (required for particle_radius binning)
    """
    if R_inner is None:
        R_inner = radii.max() * 1.001  # small safety margin

    if binning_method == "equal_area":
        edges = np.sqrt(np.linspace(0.0, 1.0, Nshell + 1)) * R_inner
    elif binning_method == "particle_radius":
        if particle_radius is None:
            raise ValueError("particle_radius must be specified for particle_radius binning")
        # Create bins with width = particle_radius
        # Start from particle_radius (particles can't be closer to center than this)
        # End at R_inner - particle_radius (particles can't be closer to wall than this)
        available_radius = R_inner - 2 * particle_radius
        if available_radius <= 0:
            raise ValueError("Particle radius too large for given tube radius")
        
        # Adjust number of shells to fit nicely
        n_bins_fit = int(available_radius / particle_radius)
        actual_bin_width = available_radius / n_bins_fit
        
        edges = np.linspace(particle_radius, R_inner - particle_radius, n_bins_fit + 1)
        print(f"Particle-radius binning: {len(edges)-1} shells, width = {actual_bin_width:.4f} m")
    else:
        raise ValueError("binning_method must be 'equal_area' or 'particle_radius'")

    counts, _ = np.histogram(radii, bins=edges)
    counts = counts.astype(float)
    phi_bulk = counts.mean()
    I_r = 1.0 - 0.5 / phi_bulk * np.sum(np.abs(counts - phi_bulk))
    return I_r, counts, edges


def compute_I_r_helix(
    points_file: str,
    R: float,
    phi0: float,
    a: float,
    z0: float,
    Nshell: int = 15,
    binning_method: str = "equal_area",
    particle_radius: float = None,
) -> Tuple[float, np.ndarray, np.ndarray]:
    P = read_xyz(points_file)
    radii = np.empty(len(P))

    for i, p in enumerate(P):
        _, c, c1 = newton_project_helix(p, R, phi0, a, z0)
        radii[i] = radial_distance(p, c, c1)

    return radial_index_from_radii(radii, Nshell, binning_method=binning_method, 
                                  particle_radius=particle_radius)


def compute_I_r_polyline(
    points_file: str, centerline_file: str, Nshell: int = 15,
    binning_method: str = "equal_area", particle_radius: float = None
) -> Tuple[float, np.ndarray, np.ndarray]:
    P = read_xyz(points_file)
    C = read_xyz(centerline_file)
    radii = np.empty(len(P))

    for i, p in enumerate(P):
        best_c = None
        best_tan = None
        best_d2 = np.inf

        for j in range(len(C) - 1):
            c_try, tan_try = closest_point_on_segment(p, C[j], C[j + 1])
            d2 = np.dot(p - c_try, p - c_try)
            if d2 < best_d2:
                best_d2 = d2
                best_c = c_try
                best_tan = tan_try

        radii[i] = radial_distance(p, best_c, best_tan)

    return radial_index_from_radii(radii, Nshell, binning_method=binning_method, 
                                  particle_radius=particle_radius)


# ----------------------------------------------------------------------
# Visualization functions
# ----------------------------------------------------------------------
def plot_radial_heatmap(counts: np.ndarray, edges: np.ndarray, I_r: float, 
                        output_file: str = None, title: str = None, 
                        binning_method: str = "equal_area", particle_radius: float = None):
    """Create a radial heat map showing particle distribution."""
    
    # Convert counts to percentages
    percentages = 100 * counts / counts.sum()
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 6))
    
    # Subplot 1: Circular heat map
    ax1 = plt.subplot(1, 3, 1)
    ax1.set_aspect('equal')
    
    # Create colormap
    cmap = plt.cm.RdYlBu_r  # Red-Yellow-Blue reversed (blue=low, red=high)
    
    # Normalize colors based on percentages
    vmax = percentages.max()
    vmin = percentages.min()
    
    # Draw concentric annular regions
    for i in range(len(counts)):
        inner_r = edges[i]
        outer_r = edges[i + 1]
        
        # Create annular wedge (full circle)
        wedge = Wedge((0, 0), outer_r, 0, 360, width=outer_r-inner_r, 
                     facecolor=cmap(percentages[i] / vmax), 
                     edgecolor='black', linewidth=0.5, alpha=0.8)
        ax1.add_patch(wedge)
    
    # Add center point
    ax1.plot(0, 0, 'ko', markersize=3)
    
    # Add particle exclusion zones if using particle-radius binning
    if binning_method == "particle_radius" and particle_radius is not None:
        # Inner exclusion zone
        inner_zone = plt.Circle((0, 0), particle_radius, fill=False, 
                               color='red', linestyle='-', linewidth=2, alpha=0.7)
        ax1.add_patch(inner_zone)
#        ax1.text(0, particle_radius + 0.01, 'Particle\nexclusion', 
#                ha='center', va='bottom', fontsize=8, color='red', weight='bold')
        
        # Outer exclusion zone (assuming tube wall)
        max_r = edges[-1] + particle_radius
        outer_zone = plt.Circle((0, 0), max_r, fill=False, 
                               color='red', linestyle='-', linewidth=2, alpha=0.7)
        ax1.add_patch(outer_zone)
#        ax1.text(max_r * 0.7, max_r * 0.7, 'Wall exclusion\nzone', 
#                ha='center', va='center', fontsize=8, color='red', weight='bold')
        plot_max_r = max_r * 1.1
    else:
        plot_max_r = edges[-1] * 1.1
    
    # Add radius labels
    max_r = edges[-1]
    for r in np.linspace(0, max_r, 6)[1:]:  # Skip r=0
        circle = plt.Circle((0, 0), r, fill=False, color='gray', 
                           linestyle='--', alpha=0.5, linewidth=0.5)
        ax1.add_patch(circle)
#        ax1.text(r, 0, f'{r:.2f}m', fontsize=8, ha='left', va='bottom')
    
    ax1.set_xlim(-plot_max_r, plot_max_r)
    ax1.set_ylim(-plot_max_r, plot_max_r)
    ax1.set_title(f'Radial Heat Map\n({binning_method.replace("_", " ").title()})')
    ax1.grid(True, alpha=0.3)
    
    # Subplot 2: Bar chart (histogram)
    ax2 = plt.subplot(1, 3, 2)
    
    # Calculate shell centers for x-axis
    shell_centers = (edges[:-1] + edges[1:]) / 2
    shell_widths = edges[1:] - edges[:-1]
    
    bars = ax2.bar(shell_centers, percentages, width=shell_widths*0.8, 
                   color=[cmap(p / vmax) for p in percentages], 
                   edgecolor='black', linewidth=0.5, alpha=0.8)
    
    # Add average line
    avg_percentage = percentages.mean()
    ax2.axhline(y=avg_percentage, color='red', linestyle='--', linewidth=2, 
                label=f'Average: {avg_percentage:.1f}%')
    
    # Add particle radius indicators if using particle-radius binning
    if binning_method == "particle_radius" and particle_radius is not None:
        ax2.axvline(x=particle_radius, color='red', linestyle='-', linewidth=1, 
                   alpha=0.7, label=f'Particle radius: {particle_radius:.3f}m')
        if len(edges) > 0:
            wall_pos = edges[-1] + particle_radius
            ax2.axvline(x=wall_pos, color='red', linestyle='-', linewidth=1, 
                       alpha=0.7, label=f'Wall position: {wall_pos:.3f}m')
    
    ax2.set_xlabel('Radial Distance (m)')
    ax2.set_ylabel('Particle Percentage (%)')
    ax2.set_title('Radial Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Summary statistics
    ax3 = plt.subplot(1, 3, 3)
    ax3.axis('off')
    
    # Create text summary
    binning_info = f"Binning: {binning_method.replace('_', ' ').title()}"
    if binning_method == "particle_radius" and particle_radius is not None:
        binning_info += f"\nParticle radius: {particle_radius:.4f} m"
        binning_info += f"\nShell width: {shell_widths[0]:.4f} m"
        if len(shell_widths) > 1 and not np.allclose(shell_widths, shell_widths[0]):
            binning_info += f" (avg)"
    
    summary_text = f"""
    Radial Index Analysis
    
    I_r = {I_r:.4f}
    
    {binning_info}
    
    Interpretation:
    • I_r = 1.0: Perfect uniformity
    • I_r = 0.0: Moderate segregation
    • I_r < 0.0: Severe segregation
    
    Distribution Stats:
    • Total particles: {int(counts.sum())}
    • Number of shells: {len(counts)}
    • Max percentage: {percentages.max():.1f}%
    • Min percentage: {percentages.min():.1f}%
    • Average per shell: {avg_percentage:.1f}%
    • Std deviation: {percentages.std():.1f}%
    
    Shell with most particles:
    • Shell {np.argmax(counts) + 1} ({percentages.max():.1f}%)
    • Radius: {edges[np.argmax(counts)]:.3f} - {edges[np.argmax(counts) + 1]:.3f} m
    """
    
    ax3.text(0.05, 0.95, summary_text, transform=ax3.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    
    # Add interpretation color coding
    if I_r >= 0.8:
        status_color = 'green'
        status_text = 'Good mixing'
    elif I_r >= 0.5:
        status_color = 'orange'
        status_text = 'Moderate mixing'
    else:
        status_color = 'red'
        status_text = 'Poor mixing'
    
    ax3.text(0.05, 0.00, f'Status: {status_text}', transform=ax3.transAxes, 
             fontsize=12, fontweight='bold', color=status_color,
             bbox=dict(boxstyle="round,pad=0.3", facecolor=status_color, alpha=0.2))
    
    # Overall title
    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    else:
        fig.suptitle(f'Radial Distribution Analysis (I_r = {I_r:.4f})', 
                     fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=[ax1, ax2], orientation='horizontal', 
                       fraction=0.05, pad=0.1, aspect=30)
    cbar.set_label('Particle Percentage (%)')
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    
    plt.show()


def save_results_csv(counts: np.ndarray, edges: np.ndarray, I_r: float, 
                     output_file: str):
    """Save results to CSV file for further analysis."""
    percentages = 100 * counts / counts.sum()
    shell_centers = (edges[:-1] + edges[1:]) / 2
    
    with open(output_file, 'w') as f:
        f.write("# Radial distribution analysis results\n")
        f.write(f"# I_r = {I_r:.6f}\n")
        f.write(f"# Total particles = {int(counts.sum())}\n")
        f.write("shell_number,inner_radius_m,outer_radius_m,center_radius_m,particle_count,particle_percentage\n")
        
        for i in range(len(counts)):
            f.write(f"{i+1},{edges[i]:.6f},{edges[i+1]:.6f},{shell_centers[i]:.6f},"
                   f"{int(counts[i])},{percentages[i]:.3f}\n")


# ----------------------------------------------------------------------
# Command line interface
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Compute the radial index I_r from particle positions."
    )
    sub = parser.add_subparsers(dest="mode", required=True, metavar="MODE")

    # Helix mode
    p_helix = sub.add_parser("helix", help="Analytic helix centre‑line")
    p_helix.add_argument("points", help="ASCII file with particle positions x y z")
    p_helix.add_argument("--R", type=float, required=True, help="Helix radius [m]")
    p_helix.add_argument("--phi0", type=float, required=True, help="Phase shift [rad]")
    p_helix.add_argument("--a", type=float, required=True, help="Pitch coefficient dz/dt")
    p_helix.add_argument("--z0", type=float, required=True, help="z‑offset [m]")
    p_helix.add_argument("--Nshell", type=int, default=15, help="# equal‑area shells")
    p_helix.add_argument("--plot", action="store_true", help="Generate visualization plots")
    p_helix.add_argument("--output", type=str, help="Output file prefix for plots and data")
    p_helix.add_argument("--title", type=str, help="Custom title for plots")
    p_helix.add_argument("--binning", choices=["equal_area", "particle_radius"], 
                        default="equal_area", help="Binning method")
    p_helix.add_argument("--particle_radius", type=float, 
                        help="Particle radius [m] (required for particle_radius binning)")

    # Polyline mode
    p_poly = sub.add_parser("polyline", help="Discretised centre‑line (polyline)")
    p_poly.add_argument("points", help="ASCII file with particle positions x y z")
    p_poly.add_argument(
        "centerline", help="ASCII file with sampled centre‑line x y z"
    )
    p_poly.add_argument("--Nshell", type=int, default=15, help="# equal‑area shells")
    p_poly.add_argument("--plot", action="store_true", help="Generate visualization plots")
    p_poly.add_argument("--output", type=str, help="Output file prefix for plots and data")
    p_poly.add_argument("--title", type=str, help="Custom title for plots")
    p_poly.add_argument("--binning", choices=["equal_area", "particle_radius"], 
                        default="equal_area", help="Binning method")
    p_poly.add_argument("--particle_radius", type=float, 
                        help="Particle radius [m] (required for particle_radius binning)")

    args = parser.parse_args()

    # Validate particle_radius argument
    if args.binning == "particle_radius" and args.particle_radius is None:
        parser.error("--particle_radius is required when using particle_radius binning")

    if args.mode == "helix":
        I_r, counts, edges = compute_I_r_helix(
            args.points, args.R, args.phi0, args.a, args.z0, args.Nshell,
            binning_method=args.binning, particle_radius=args.particle_radius
        )
    else:
        I_r, counts, edges = compute_I_r_polyline(
            args.points, args.centerline, args.Nshell,
            binning_method=args.binning, particle_radius=args.particle_radius
        )

    # Print basic results
    print(f"I_r = {I_r:.4f}")
    print("Shell edges (m):", edges)
    print("Counts per shell:", counts.astype(int))
    
    # Convert to percentages for display
    percentages = 100 * counts / counts.sum()
    print("Percentages per shell:", [f"{p:.1f}%" for p in percentages])
    
    if args.binning == "particle_radius":
        print(f"Binning method: Particle-radius based (width ≈ {args.particle_radius:.4f} m)")
        print(f"Physical exclusion zones: center < {args.particle_radius:.4f} m, wall > {edges[-1]+args.particle_radius:.4f} m")
    else:
        print("Binning method: Equal-area shells")

    # Generate plots if requested
    if args.plot:
        plot_file = None
        if args.output:
            plot_file = f"{args.output}_{args.binning}_radial_analysis.png"
        
        plot_radial_heatmap(counts, edges, I_r, plot_file, args.title, 
                           args.binning, args.particle_radius)
    
    # Save CSV data if output specified
    if args.output:
        csv_file = f"{args.output}_{args.binning}_radial_data.csv"
        save_results_csv(counts, edges, I_r, csv_file)
        print(f"Data saved to: {csv_file}")


if __name__ == "__main__":
    main()