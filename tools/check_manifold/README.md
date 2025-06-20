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

### Basic Usage

```bash
./validate_triangle_mesh <mesh_file.off>
```

### With Custom Thresholds

```bash
./validate_triangle_mesh <mesh_file.off> [area_threshold] [min_angle] [max_angle] [aspect_ratio] [circumradius_ratio]
```

### Parameters

- `area_threshold`: Minimum triangle area (default: 1e-12)
- `min_angle`: Minimum angle in degrees (default: 0.1)
- `max_angle`: Maximum angle in degrees (default: 179.9)
- `aspect_ratio`: Maximum aspect ratio (default: 100.0)
- `circumradius_ratio`: Maximum circumradius/inradius ratio (default: 100.0)

### Example

```bash
# Use strict quality criteria
./validate_triangle_mesh my_mesh.off 1e-6 1.0 179.0 30.0 10.0
```

## Output

The tool provides detailed reports including:

```
=== Triangle Mesh Validation Report ===
Mesh file: my_mesh.off
Vertices: 8
Triangles: 12

=== Manifold Check ===
Mesh is MANIFOLD

=== Quality Analysis ===
Degenerate triangles: 0
Small area triangles: 0
Extreme angle triangles: 0
High aspect ratio triangles: 0
High circumradius/inradius triangles: 0

Final verdict: Mesh is MANIFOLD and passes all quality checks
```

## Test Suites

Two test suites are provided to validate the tool functionality:

### 1. Bash Test Suite (`run_tests.sh`)

```bash
# Run all tests
./run_tests.sh
```

Features:
- Colored output for easy result interpretation
- Tests various mesh types (clean, degenerate, non-manifold)
- Automatic build if needed
- Comprehensive test result summary

### 2. Python Test Suite (`test_validator.py`)

```bash
# Run Python test suite
python3 test_validator.py
```

Features:
- Structured test cases with expected results
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