# HashGrid Acceleration Test Application

## Purpose

This test application verifies that the HashGrid-accelerated point-inside queries produce identical results to the baseline linear search implementation.

## Test Setup

- **Particles**: 1000 spheres arranged in a 10×10×10 grid
- **Domain**: [0,1]³ unit cube
- **Sphere radius**: 0.01
- **Sphere spacing**: 0.1 (center-to-center)
- **Material**: Iron (PE standard material)
- **PE Mode**: Serial mode (no MPI within PE library)

## Test Procedure

The test performs a single CFD timestep with only alpha field computation (no NS solve). This exercises the point-inside query system for all mesh nodes against all 1000 spheres.

### Success Criteria

Both methods (accelerated ON/OFF) should produce **identical solid/fluid node counts**.

## Build Instructions

### 1. Build Baseline Version (no acceleration)

```bash
cd /path/to/ff-accel
mkdir build_test_baseline
cd build_test_baseline

cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_PE=ON \
      -DUSE_PE_SERIAL_MODE=ON \
      -DUSE_ACCELERATED_POINT_QUERY=OFF \
      -DBUILD_APPLICATIONS=ON \
      ..

make -j8 q2p1_hashgrid_test
```

### 2. Build Accelerated Version

```bash
cd /path/to/ff-accel
mkdir build_test_accelerated
cd build_test_accelerated

cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_PE=ON \
      -DUSE_PE_SERIAL_MODE=ON \
      -DUSE_ACCELERATED_POINT_QUERY=ON \
      -DBUILD_APPLICATIONS=ON \
      ..

make -j8 q2p1_hashgrid_test
```

## Running the Test

### Prerequisites

You need a mesh file and corresponding partition. The test expects:
- `_data/q2p1_param.dat` - Parameter file (already provided)
- A mesh file (`.tri` format) - copy from another application or create
- Partition files for the mesh

Example setup (assuming 4 MPI processes):

```bash
cd applications/q2p1_hashgrid_test/_data

# Copy a simple mesh from another application
cp ../../q2p1_fc2/_data/cube.tri .

# Update the mesh filename in your simulation parameters if needed
```

### Run Baseline

```bash
cd build_test_baseline/applications/q2p1_hashgrid_test
mpirun -np 4 ./q2p1_hashgrid_test > ../../../output_baseline.txt 2>&1
```

### Run Accelerated

```bash
cd build_test_accelerated/applications/q2p1_hashgrid_test
mpirun -np 4 ./q2p1_hashgrid_test > ../../../output_accelerated.txt 2>&1
```

### Compare Results

```bash
cd /path/to/ff-accel

# Extract and compare solid node counts
grep "Total dofs inside" output_baseline.txt
grep "Total dofs inside" output_accelerated.txt

# Full difference
diff output_baseline.txt output_accelerated.txt
```

## Expected Output

The test should print:

```
--HASHGRID TEST SETUP
Domain: [0, 1]^3
Spheres: 1000 (10x10x10 grid)
Radius: 0.01
Spacing: 0.1
Material: Iron (density=7.874)

Number of particles: 1000

HASHGRID TEST: Computing alpha field
> Total dofs inside: XXXXX
> Dofs per Particle: YYY

Test complete.
```

The key metric is **"Total dofs inside"** - this count must be identical between baseline and accelerated runs.

## Verification Checklist

- [ ] PE setup executed (see "HASHGRID TEST SETUP" in output)
- [ ] 1000 particles created (see "Number of particles: 1000")
- [ ] Alpha field computed (see "Total dofs inside" message)
- [ ] No errors or warnings
- [ ] Counts match between baseline and accelerated builds

## Troubleshooting

### HashGrid Not Activating

The HashGrid activates automatically when the particle count exceeds a threshold. With 1000 spheres, it should activate. To verify:

Add diagnostic output in `libs/pe/src/interface/object_queries.cpp` to print which method is being used.

### Different Counts

If the counts differ:

1. Check that both builds are using PE_SERIAL_MODE
2. Verify the CMake configuration (check for `USE_ACCELERATED_POINT_QUERY`)
3. Enable debug output in the point query code to trace the discrepancy
4. Try with fewer spheres (e.g., 5×5×5 = 125) to test below the activation threshold

### Build Errors

- Ensure `USE_PE=ON` is set before enabling `USE_PE_SERIAL_MODE`
- Use the two-step CMake configuration as shown above
- Check that all PE library headers are up to date

## Files

- `q2p1_hashgrid_test.f90` - Main program
- `app_init.f90` - Initialization (mesh loading, PE setup)
- `CMakeLists.txt` - Build configuration
- `_data/q2p1_param.dat` - Parameter file
- `README.md` - This file
