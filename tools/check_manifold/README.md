# Triangle Mesh Validator

A CGAL-based tool for validating triangle meshes in OFF format. This tool checks for manifold properties and various mesh quality issues.

## Overview

The `validate_triangle_mesh` tool performs comprehensive validation of triangle meshes, checking for:

- **Manifold properties**: Whether the mesh represents a valid 2-manifold surface
- **Degenerate triangles**: Triangles with zero area or collinear vertices
- **Small area triangles**: Triangles below a specified area threshold
- **Extreme angle triangles**: Triangles with very acute or obtuse angles
- **High aspect ratio triangles**: Triangles that are excessively elongated
- **High circumradius/inradius ratio**: Triangles with poor geometric properties

## Building

The tool requires CGAL and is built as part of the main FeatFloWer build system:

```bash
# Configure with CGAL support
cmake -DUSE_CGAL=ON -DCMAKE_BUILD_TYPE=Release ..

# Build the project
make
```

The executable will be built in the `tools/check_manifold` directory within your build folder.

## Usage

### Basic Usage (Production Mode)

```bash
./validate_triangle_mesh <mesh_file.off>
```

This provides a clean, summary-only output suitable for production validation workflows.

### Command Line Options

- `--enable-repair-check`: Enable additional manifold validation that includes auto-repair checking
- `--verbose`: Show detailed face-by-face information for debugging

### Usage Patterns

#### Quick Production Check (Default)
```bash
# Fast validation with summary output only
./validate_triangle_mesh mesh.off
```

#### Detailed Debugging Mode
```bash
# Show detailed face-by-face triangle information
./validate_triangle_mesh mesh.off --verbose
```

#### Enable Auto-Repair Validation
```bash
# Include manifold check with potential auto-repair
./validate_triangle_mesh mesh.off --enable-repair-check
```

#### Full Debug Mode
```bash
# Both verbose output and auto-repair checking
./validate_triangle_mesh mesh.off --verbose --enable-repair-check
```

### Custom Thresholds

Thresholds can be specified after the options:

```bash
./validate_triangle_mesh mesh.off [OPTIONS] [area_threshold] [min_angle] [max_angle] [aspect_ratio] [circumradius_ratio]
```

### Threshold Parameters

- `area_threshold`: Minimum triangle area (default: 1e-12)
- `min_angle`: Minimum angle in degrees (default: 1.0)
- `max_angle`: Maximum angle in degrees (default: 179.0)
- `aspect_ratio`: Maximum aspect ratio (default: 30.0)
- `circumradius_ratio`: Maximum circumradius/inradius ratio (default: 10.0)

### Complete Examples

```bash
# Production validation with custom thresholds
./validate_triangle_mesh my_mesh.off 1e-6 2.0 178.0 25.0 8.0

# Debug mode with strict quality criteria
./validate_triangle_mesh my_mesh.off --verbose --enable-repair-check 1e-6 1.0 179.0 30.0 10.0
```

## Output

The tool provides different levels of output depending on the options used:

### Default Output (Production Mode)

```
CGAL Triangle Mesh Validator
CGAL version: 5.6
Loaded mesh with 8 vertices and 12 faces

=== Manifold Check ===
Raw manifold check (no auto-repair): MANIFOLD
Auto-repair check disabled (use --enable-repair-check to enable)
Final verdict: Mesh is MANIFOLD

=== Triangle Quality Analysis ===

=== Summary ===
Total triangles: 12
Degenerate triangles: 0
Small area triangles: 0
Extreme angle triangles: 0
High aspect ratio triangles: 0
High circumradius/inradius triangles: 0
```

### Verbose Output (Debug Mode)

With `--verbose` flag, detailed face-by-face information is shown:

```
Configuration:
  Repair check enabled: NO
  Verbose output: YES
  Thresholds: area=1e-12, minAng=1°, maxAng=179°, edgeRatio=30, RrRatio=10

=== Triangle Quality Analysis ===
Area threshold: 1e-12
Min angle threshold: 1°
Max angle threshold: 179°
Aspect ratio threshold: 30
Circumradius/inradius threshold: 10

[Individual triangle issues would be listed here if found]
Extreme angle triangle (face 5): min = 0.5°, max = 179.2°
High aspect ratio triangle (face 8): ratio = 45.2

=== Summary ===
...
```

### With Auto-Repair Check

With `--enable-repair-check` flag:

```
=== Manifold Check ===
Raw manifold check (no auto-repair): NON-MANIFOLD
Loaded mesh manifold check: MANIFOLD
Final verdict: Mesh is NON-MANIFOLD
```

## Workflow Separation

The tool is designed to serve different user workflows:

- **Production Users**: Quick validation with clean, summary-only output
- **Debug/Repair Users**: Detailed analysis with face-by-face information and auto-repair checking

## Test Suites

Two test suites are provided to validate the tool functionality:

### 1. Bash Test Suite (`run_tests.sh`)

```bash
# Run all tests
./run_tests.sh

# Test with verbose output
./run_tests.sh --verbose

# Test with auto-repair checking
./run_tests.sh --enable-repair-check
```

Features:
- Colored output for easy result interpretation
- Tests various mesh types (clean, degenerate, non-manifold)
- Supports testing both production and debug modes
- Automatic build if needed
- Comprehensive test result summary

### 2. Python Test Suite (`test_validator.py`)

```bash
# Run Python test suite
python3 test_validator.py

# Test specific modes
python3 test_validator.py --test-verbose
python3 test_validator.py --test-repair-check
```

Features:
- Structured test cases with expected results
- Tests both production and debug output modes
- Detailed output parsing and validation
- Automatic build if needed
- Comprehensive error reporting

## Test Files

The `test_files/` directory contains various test meshes:

- `clean_tetrahedron.off`: Perfect tetrahedron (should pass all tests)
- `cube.off`: Simple cube mesh
- `zero_area_triangle.off`: Contains degenerate triangles
- `sliver_triangle.off`: Contains very thin triangles
- `nonmanifold_edge.off`: Non-manifold mesh example
- `open_cylinder.off`: Mesh with boundary edges

## Integration with FeatFloWer

The tool is automatically built when:
- Building on non-Windows systems
- CGAL support is enabled (`-DUSE_CGAL=ON`)
- The CGAL external project dependency is satisfied

## Troubleshooting

### Build Issues

- Ensure CGAL is properly installed or the external project is configured
- Check that C++17 support is available
- Verify CMAKE_BUILD_TYPE is set appropriately

### Runtime Issues

- Verify the input file is in valid OFF format
- Check file permissions and path accessibility
- Ensure the mesh file contains valid geometric data

## Dependencies

- **CGAL**: Computational Geometry Algorithms Library
- **C++17**: Modern C++ standard support
- **CMake 3.12+**: Build system requirement