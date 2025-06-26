#!/usr/bin/env python3
"""compute_metrics.py - Enhanced with radial-vertical visualization

Compute radial index I_r **and** vertical asymmetry A_y for particle sets
in an Archimedes Tube Crystallizer (ATC) slug, with comprehensive visualization
showing both radial distribution and upper/lower asymmetry per shell.

Two modes:

  1. Analytic helix centre‑line
     python compute_metrics.py helix points.txt --R 2.835888 --phi0 5.474463 \
                                                 --a -0.112848 --z0 -0.678574 --plot

  2. Discretised centre‑line (polyline)
     python compute_metrics.py polyline points.txt center_line_sampled.txt --plot

Output:
  I_r,  A_y,  shell edges, counts per shell, radial-vertical analysis
"""

import argparse
import math
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.patches as patches


def read_xyz(fname: str) -> np.ndarray:
    return np.loadtxt(fname, dtype=float)

def frenet_helix(t: float, R: float, phi0: float, a: float, z0: float):
    s, c = math.sin(t + phi0), math.cos(t + phi0)
    c_pos = np.array([R * c, R * s, a * t + z0])
    c1 = np.array([-R * s, R * c, a])
    c2 = np.array([-R * c, -R * s, 0.0])
    return c_pos, c1, c2

def newton_project_helix(p: np.ndarray, R: float, phi0: float, a: float, z0: float, tol: float = 1e-10, max_iter: int = 25):
    t = math.atan2(p[1], p[0]) - phi0
    for _ in range(max_iter):
        c, c1, c2 = frenet_helix(t, R, phi0, a, z0)
        g = np.dot(c1, p - c)
        g_p = -np.dot(c2, p - c) - np.dot(c1, c1)
        dt = -g / g_p
        t += dt
        if abs(dt) < tol:
            break
    T = c1 / np.linalg.norm(c1)
    N_temp = c2 - np.dot(c2, T) * T
    N = N_temp / np.linalg.norm(N_temp)
    B = np.cross(T, N)
    return c, T, B

def closest_point_on_segment(p: np.ndarray, a: np.ndarray, b: np.ndarray):
    ab = b - a
    t = np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0.0, 1.0)
    c = a + t * ab
    T = ab / np.linalg.norm(ab)
    return c, T

def binormal_polyline(C: np.ndarray, idx: int):
    if 0 < idx < len(C) - 1:
        T_prev = C[idx] - C[idx - 1]
        T_next = C[idx + 1] - C[idx]
        T_prev /= np.linalg.norm(T_prev)
        T_next /= np.linalg.norm(T_next)
        dT = T_next - T_prev
        if np.linalg.norm(dT) > 1e-8:
            N = dT / np.linalg.norm(dT)
            T_avg = (T_prev + T_next) / 2.0
            T_avg /= np.linalg.norm(T_avg)
            B = np.cross(T_avg, N)
            if np.linalg.norm(B) > 1e-8:
                return B / np.linalg.norm(B)
    return np.array([0.0, 0.0, 1.0])

def equal_area_edges(R_inner: float, N: int):
    return np.sqrt(np.linspace(0.0, 1.0, N + 1)) * R_inner

def radial_index_with_vertical(radii: np.ndarray, vertical_positions: np.ndarray, R_inner: float, Nshell: int):
    """Compute I_r and detailed radial-vertical distribution."""
    edges = equal_area_edges(R_inner, Nshell)
    
    # Overall radial distribution
    counts, _ = np.histogram(radii, bins=edges)
    counts = counts.astype(float)
    phi_bulk = counts.mean()
    I_r = 1.0 - 0.5 / phi_bulk * np.sum(np.abs(counts - phi_bulk))
    
    # Radial-vertical distribution
    upper_counts = np.zeros(Nshell)
    lower_counts = np.zeros(Nshell)
    
    for i in range(len(radii)):
        # Find which shell this particle belongs to
        shell_idx = np.searchsorted(edges[1:], radii[i])
        if shell_idx >= Nshell:
            shell_idx = Nshell - 1
            
        # Classify as upper or lower
        if vertical_positions[i] > 0:
            upper_counts[shell_idx] += 1
        else:
            lower_counts[shell_idx] += 1
    
    # Calculate asymmetry per shell
    shell_asymmetries = np.zeros(Nshell)
    for i in range(Nshell):
        if upper_counts[i] + lower_counts[i] > 0:
            shell_asymmetries[i] = (upper_counts[i] - lower_counts[i]) / (upper_counts[i] + lower_counts[i])
    
    return I_r, counts, edges, upper_counts, lower_counts, shell_asymmetries

def metrics_helix(points: str, R: float, phi0: float, a: float, z0: float, Nshell: int):
    P = read_xyz(points)
    radii = np.empty(len(P))
    vertical_positions = np.empty(len(P))
    up, down = 0, 0
    
    for i, p in enumerate(P):
        c, T, B = newton_project_helix(p, R, phi0, a, z0)
        r_vec = p - c - np.dot(p - c, T) * T
        radii[i] = np.linalg.norm(r_vec)
        vertical_positions[i] = np.dot(r_vec, B)
        
        if vertical_positions[i] > 0:
            up += 1
        else:
            down += 1
    
    I_r, counts, edges, upper_counts, lower_counts, shell_asymmetries = radial_index_with_vertical(
        radii, vertical_positions, radii.max() * 1.001, Nshell)
    A_y = (up - down) / (up + down)
    
    return I_r, A_y, counts, edges, upper_counts, lower_counts, shell_asymmetries

def metrics_polyline(points: str, center_file: str, Nshell: int):
    P = read_xyz(points)
    C = read_xyz(center_file)
    radii = np.empty(len(P))
    vertical_positions = np.empty(len(P))
    up, down = 0, 0
    
    seg_T = C[1:] - C[:-1]
    seg_T /= np.linalg.norm(seg_T, axis=1)[:, None]
    
    for i, p in enumerate(P):
        best_idx, best_c, best_d2 = None, None, np.inf
        for j in range(len(C) - 1):
            a, b = C[j], C[j + 1]
            c_try, T_try = closest_point_on_segment(p, a, b)
            d2 = np.dot(p - c_try, p - c_try)
            if d2 < best_d2:
                best_d2, best_c, best_idx = d2, c_try, j
        
        T = seg_T[best_idx]
        B = binormal_polyline(C, best_idx + 1)
        r_vec = p - best_c - np.dot(p - best_c, T) * T
        radii[i] = np.linalg.norm(r_vec)
        vertical_positions[i] = np.dot(r_vec, B)
        
        if vertical_positions[i] > 0:
            up += 1
        else:
            down += 1
    
    I_r, counts, edges, upper_counts, lower_counts, shell_asymmetries = radial_index_with_vertical(
        radii, vertical_positions, radii.max() * 1.001, Nshell)
    A_y = (up - down) / (up + down)
    
    return I_r, A_y, counts, edges, upper_counts, lower_counts, shell_asymmetries


def plot_radial_vertical_analysis(counts, edges, upper_counts, lower_counts, shell_asymmetries, 
                                 I_r, A_y, output_file=None, title=None):
    """Create comprehensive radial-vertical analysis visualization."""
    
    # Convert to percentages
    total_particles = counts.sum()
    percentages = 100 * counts / total_particles
    upper_percentages = 100 * upper_counts / total_particles
    lower_percentages = 100 * lower_counts / total_particles
    
    # Create figure with subplots
    fig = plt.figure(figsize=(18, 10))
    
    # Create custom layout
    gs = fig.add_gridspec(3, 4, height_ratios=[2, 2, 1], width_ratios=[1, 1, 1, 1])
    
    # 1. Circular heat map - Total distribution
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_aspect('equal')
    
    cmap = plt.cm.RdYlBu_r
    vmax = percentages.max()
    vmin = percentages.min()
    
    # Draw total distribution
    for i in range(len(counts)):
        inner_r = edges[i]
        outer_r = edges[i + 1]
        wedge = Wedge((0, 0), outer_r, 0, 360, width=outer_r-inner_r, 
                     facecolor=cmap(percentages[i] / vmax), 
                     edgecolor='black', linewidth=0.5, alpha=0.8)
        ax1.add_patch(wedge)
    
    ax1.plot(0, 0, 'ko', markersize=3)
    max_r = edges[-1]
    ax1.set_xlim(-max_r*1.1, max_r*1.1)
    ax1.set_ylim(-max_r*1.1, max_r*1.1)
    ax1.set_title('Total Radial Distribution')
    ax1.grid(True, alpha=0.3)
    
    # 2. Upper half distribution
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_aspect('equal')
    
    upper_vmax = upper_percentages.max() if upper_percentages.max() > 0 else 1
    for i in range(len(counts)):
        inner_r = edges[i]
        outer_r = edges[i + 1]
        # Only show upper semicircle
        wedge = Wedge((0, 0), outer_r, 0, 180, width=outer_r-inner_r, 
                     facecolor=cmap(upper_percentages[i] / upper_vmax), 
                     edgecolor='black', linewidth=0.5, alpha=0.8)
        ax2.add_patch(wedge)
        # Gray out lower half
        wedge_lower = Wedge((0, 0), outer_r, 180, 360, width=outer_r-inner_r, 
                           facecolor='lightgray', edgecolor='black', linewidth=0.5, alpha=0.3)
        ax2.add_patch(wedge_lower)
    
    ax2.plot(0, 0, 'ko', markersize=3)
    ax2.set_xlim(-max_r*1.1, max_r*1.1)
    ax2.set_ylim(-max_r*1.1, max_r*1.1)
    ax2.set_title('Upper Half Distribution')
    ax2.grid(True, alpha=0.3)
    
    # 3. Lower half distribution
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_aspect('equal')
    
    lower_vmax = lower_percentages.max() if lower_percentages.max() > 0 else 1
    for i in range(len(counts)):
        inner_r = edges[i]
        outer_r = edges[i + 1]
        # Gray out upper half
        wedge_upper = Wedge((0, 0), outer_r, 0, 180, width=outer_r-inner_r, 
                           facecolor='lightgray', edgecolor='black', linewidth=0.5, alpha=0.3)
        ax3.add_patch(wedge_upper)
        # Only show lower semicircle
        wedge = Wedge((0, 0), outer_r, 180, 360, width=outer_r-inner_r, 
                     facecolor=cmap(lower_percentages[i] / lower_vmax), 
                     edgecolor='black', linewidth=0.5, alpha=0.8)
        ax3.add_patch(wedge)
    
    ax3.plot(0, 0, 'ko', markersize=3)
    ax3.set_xlim(-max_r*1.1, max_r*1.1)
    ax3.set_ylim(-max_r*1.1, max_r*1.1)
    ax3.set_title('Lower Half Distribution')
    ax3.grid(True, alpha=0.3)
    
    # 4. Shell asymmetry visualization
    ax4 = fig.add_subplot(gs[0, 3])
    ax4.set_aspect('equal')
    
    # Use diverging colormap for asymmetry
    asym_cmap = plt.cm.RdBu_r  # Red = upper dominant, Blue = lower dominant
    asym_max = max(abs(shell_asymmetries.max()), abs(shell_asymmetries.min()), 0.1)
    
    for i in range(len(counts)):
        inner_r = edges[i]
        outer_r = edges[i + 1]
        
        # Normalize asymmetry to [-1, 1] range
        norm_asym = shell_asymmetries[i] / asym_max if asym_max > 0 else 0
        color_val = (norm_asym + 1) / 2  # Map [-1,1] to [0,1]
        
        wedge = Wedge((0, 0), outer_r, 0, 360, width=outer_r-inner_r, 
                     facecolor=asym_cmap(color_val), 
                     edgecolor='black', linewidth=0.5, alpha=0.8)
        ax4.add_patch(wedge)
    
    ax4.plot(0, 0, 'ko', markersize=3)
    ax4.set_xlim(-max_r*1.1, max_r*1.1)
    ax4.set_ylim(-max_r*1.1, max_r*1.1)
    ax4.set_title('Shell Asymmetry\n(Red=Upper, Blue=Lower)')
    ax4.grid(True, alpha=0.3)
    
    # 5. Radial profile with upper/lower breakdown
    ax5 = fig.add_subplot(gs[1, :2])
    
    shell_centers = (edges[:-1] + edges[1:]) / 2
    shell_widths = edges[1:] - edges[:-1]
    
    # Stacked bar chart
    bars_lower = ax5.bar(shell_centers, lower_percentages, width=shell_widths*0.8, 
                        color='lightblue', edgecolor='black', linewidth=0.5, 
                        label='Lower half', alpha=0.8)
    bars_upper = ax5.bar(shell_centers, upper_percentages, width=shell_widths*0.8, 
                        bottom=lower_percentages, color='lightcoral', 
                        edgecolor='black', linewidth=0.5, label='Upper half', alpha=0.8)
    
    # Add total average line
    avg_percentage = percentages.mean()
    ax5.axhline(y=avg_percentage, color='red', linestyle='--', linewidth=2, 
                label=f'Average total: {avg_percentage:.1f}%')
    
    ax5.set_xlabel('Radial Distance (m)')
    ax5.set_ylabel('Particle Percentage (%)')
    ax5.set_title('Radial Distribution: Upper vs Lower Half')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. Shell asymmetry profile
    ax6 = fig.add_subplot(gs[1, 2:])
    
    # Bar chart of asymmetries
    colors = [asym_cmap((asym + asym_max) / (2 * asym_max)) for asym in shell_asymmetries]
    bars_asym = ax6.bar(range(1, len(shell_asymmetries) + 1), shell_asymmetries, 
                       color=colors, edgecolor='black', linewidth=0.5, alpha=0.8)
    
    # Add zero line and overall asymmetry
    ax6.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    ax6.axhline(y=A_y, color='red', linestyle='--', linewidth=2, 
                label=f'Overall A_y: {A_y:.3f}')
    
    ax6.set_xlabel('Shell Number')
    ax6.set_ylabel('Asymmetry Factor A_y')
    ax6.set_title('Vertical Asymmetry per Shell')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    ax6.set_ylim(-1.1, 1.1)
    
    # 7. Summary statistics
    ax7 = fig.add_subplot(gs[2, :])
    ax7.axis('off')
    
    # Calculate additional statistics
    dominant_upper_shells = np.sum(shell_asymmetries > 0.1)
    dominant_lower_shells = np.sum(shell_asymmetries < -0.1)
    symmetric_shells = len(shell_asymmetries) - dominant_upper_shells - dominant_lower_shells
    
    summary_text = f"""
    RADIAL-VERTICAL ANALYSIS SUMMARY
    
    Overall Metrics:  I_r = {I_r:.4f} (radial segregation)  |  A_y = {A_y:.4f} (vertical asymmetry)
    
    Total particles: {int(total_particles)}  |  Upper half: {int(upper_counts.sum())} ({100*upper_counts.sum()/total_particles:.1f}%)  |  Lower half: {int(lower_counts.sum())} ({100*lower_counts.sum()/total_particles:.1f}%)
    
    Shell Analysis:  Upper-dominant shells: {dominant_upper_shells}  |  Lower-dominant shells: {dominant_lower_shells}  |  Symmetric shells: {symmetric_shells}
    
    Physical Interpretation:
    • I_r < 0: Severe radial segregation (annular concentration)
    • A_y > 0: Dean vortices lift particles upward  |  A_y < 0: Gravity dominates, particles settle downward
    • Shell asymmetries reveal radial dependence of vertical forces
    """
    
    ax7.text(0.05, 0.95, summary_text, transform=ax7.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    
    # Overall title
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold')
    else:
        fig.suptitle(f'Radial-Vertical Distribution Analysis (I_r = {I_r:.4f}, A_y = {A_y:.4f})', 
                     fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    
    plt.show()


def save_detailed_results_csv(counts, edges, upper_counts, lower_counts, shell_asymmetries, 
                             I_r, A_y, output_file):
    """Save detailed radial-vertical results to CSV."""
    total_particles = counts.sum()
    percentages = 100 * counts / total_particles
    upper_percentages = 100 * upper_counts / total_particles
    lower_percentages = 100 * lower_counts / total_particles
    shell_centers = (edges[:-1] + edges[1:]) / 2
    
    with open(output_file, 'w') as f:
        f.write("# Radial-vertical distribution analysis results\n")
        f.write(f"# I_r = {I_r:.6f}\n")
        f.write(f"# A_y = {A_y:.6f}\n")
        f.write(f"# Total particles = {int(total_particles)}\n")
        f.write("shell_number,inner_radius_m,outer_radius_m,center_radius_m,")
        f.write("total_count,total_percentage,upper_count,upper_percentage,")
        f.write("lower_count,lower_percentage,shell_asymmetry\n")
        
        for i in range(len(counts)):
            f.write(f"{i+1},{edges[i]:.6f},{edges[i+1]:.6f},{shell_centers[i]:.6f},")
            f.write(f"{int(counts[i])},{percentages[i]:.3f},")
            f.write(f"{int(upper_counts[i])},{upper_percentages[i]:.3f},")
            f.write(f"{int(lower_counts[i])},{lower_percentages[i]:.3f},")
            f.write(f"{shell_asymmetries[i]:.6f}\n")


def main():
    parser = argparse.ArgumentParser(description="Compute I_r and A_y metrics with radial-vertical visualization.")
    sub = parser.add_subparsers(dest="mode", required=True)

    ph = sub.add_parser("helix", help="Analytic helix centre‑line")
    ph.add_argument("points", help="particle positions file (x y z)")
    ph.add_argument("--R", type=float, required=True)
    ph.add_argument("--phi0", type=float, required=True)
    ph.add_argument("--a", type=float, required=True)
    ph.add_argument("--z0", type=float, required=True)
    ph.add_argument("--Nshell", type=int, default=15)
    ph.add_argument("--plot", action="store_true", help="Generate visualization plots")
    ph.add_argument("--output", type=str, help="Output file prefix for plots and data")
    ph.add_argument("--title", type=str, help="Custom title for plots")

    pp = sub.add_parser("polyline", help="Discretised centre‑line")
    pp.add_argument("points", help="particle positions file (x y z)")
    pp.add_argument("centerline", help="centre‑line polyline file (x y z)")
    pp.add_argument("--Nshell", type=int, default=15)
    pp.add_argument("--plot", action="store_true", help="Generate visualization plots")
    pp.add_argument("--output", type=str, help="Output file prefix for plots and data")
    pp.add_argument("--title", type=str, help="Custom title for plots")

    args = parser.parse_args()

    if args.mode == "helix":
        I_r, A_y, counts, edges, upper_counts, lower_counts, shell_asymmetries = metrics_helix(
            args.points, args.R, args.phi0, args.a, args.z0, args.Nshell
        )
    else:
        I_r, A_y, counts, edges, upper_counts, lower_counts, shell_asymmetries = metrics_polyline(
            args.points, args.centerline, args.Nshell
        )

    print(f"I_r = {I_r:.4f}")
    print(f"A_y = {A_y:.4f}")
    print("Shell edges (m):", edges)
    print("Counts per shell:", counts.astype(int))
    print("Upper counts per shell:", upper_counts.astype(int))
    print("Lower counts per shell:", lower_counts.astype(int))
    print("Shell asymmetries:", [f"{a:.3f}" for a in shell_asymmetries])

    # Generate plots if requested
    if args.plot:
        plot_file = None
        if args.output:
            plot_file = f"{args.output}_radial_vertical_analysis.png"
        
        plot_radial_vertical_analysis(counts, edges, upper_counts, lower_counts, 
                                    shell_asymmetries, I_r, A_y, plot_file, args.title)
    
    # Save CSV data if output specified
    if args.output:
        csv_file = f"{args.output}_radial_vertical_data.csv"
        save_detailed_results_csv(counts, edges, upper_counts, lower_counts, 
                                shell_asymmetries, I_r, A_y, csv_file)
        print(f"Data saved to: {csv_file}")

if __name__ == "__main__":
    main()