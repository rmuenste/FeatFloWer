# FeatFloWer Runtime Process Control

**Version**: 1.0
**Last Updated**: 2025-01-28
**Source**: `source/ProcCtrl.f90`

## Introduction

The **ProcessControl** system provides interactive runtime control for FeatFloWer simulations. This allows you to modify simulation parameters, create checkpoints, generate output, and reload configurations **without stopping or restarting** the solver. This is particularly valuable for long-running HPC simulations where you want to adjust settings based on intermediate results or create backups during execution.

## How It Works

### Basic Mechanism

1. **During simulation**, FeatFloWer periodically checks for a control file: `_data/ProcCtrl.txt`
2. **If the file exists**, the solver reads and executes commands from it
3. **After execution**, the file is automatically cleared to prevent re-execution
4. **All MPI processes** synchronize to ensure consistent state across parallel domains

### Control Flow

```
Simulation running → Check _data/ProcCtrl.txt
                            ↓
                     File exists?
                            ↓ Yes
                     Read commands
                            ↓
                     Execute commands
                            ↓
                     MPI barrier sync
                            ↓
                     Clear file
                            ↓
                     Continue simulation
```

### File Format

Commands in `_data/ProcCtrl.txt` use a simple `Command = value` format:

```
Command_Name = value
```

**Example:**
```bash
Dump_Out = 05
```

Multiple commands can be issued at once (one per line):
```bash
Dump_Out = 10
GMV_Out = 1234
Reload_Velo = _data/new_solver_params.dat
```

---

## Available Commands

### 1. Output Commands

#### `Dump_Out` - Save Solution State

**Syntax:**
```
Dump_Out = NN
```

**Parameters:**
- `NN`: File number (00-99)

**Description:**
Creates a checkpoint by saving the current solution state (velocity, pressure fields) and Fictitious Boundary Method (FBM) data to numbered files. Use this for backup during long simulations or to restart from a specific state later.

**Output files created:**
- Solution data: `sol.NN` (or similar naming scheme)
- FBM data: corresponding FBM state files

**Example:**
```bash
# Save checkpoint to file number 05
echo "Dump_Out = 05" > _data/ProcCtrl.txt
```

**Use cases:**
- Create periodic backups during overnight runs
- Save state before testing parameter changes
- Checkpoint before potentially unstable phases

---

#### `Dump_In` - Load Solution State

**Syntax:**
```
Dump_In = filename
```

**Parameters:**
- `filename`: Path to solution file to load

**Description:**
Loads a previously saved solution state from file. The simulation continues from this loaded state. Useful for restarting from checkpoints or testing different parameter settings from the same initial condition.

**Example:**
```bash
# Load from checkpoint file
echo "Dump_In = sol.05" > _data/ProcCtrl.txt
```

**Use cases:**
- Restart simulation from checkpoint
- Branch simulation to test different scenarios
- Recover from unstable iterations

**⚠️ Warning:** Ensure the loaded state is compatible with current mesh and decomposition.

---

#### `GMV_Out` - Generate Visualization Output

**Syntax:**
```
GMV_Out = NNNN
```

**Parameters:**
- `NNNN`: File number (0000-9999)

**Description:**
Generates visualization output in GMV format (or configured format) for the current solution state. Creates output files for post-processing with ParaView or other visualization tools.

**Output files:**
- VTU/VTK files (depending on configuration)
- Mesh and solution data at current time step

**Example:**
```bash
# Generate output file 1234
echo "GMV_Out = 1234" > _data/ProcCtrl.txt
```

**Use cases:**
- Generate intermediate visualization during simulation
- Create output at specific time points of interest
- Monitor simulation progress visually

---

### 2. Parameter Reload Commands

These commands allow you to modify simulation parameters **during runtime** without restarting. The new parameters take effect immediately.

#### `Reload_Velo` - Reload Velocity Solver Parameters

**Syntax:**
```
Reload_Velo = parameter_file.dat
```

**Parameters:**
- `parameter_file.dat`: Path to parameter file containing `Velo@` parameters

**Description:**
Reloads velocity solver configuration from the specified file. This includes multigrid settings, solver tolerances, iteration counts, and smoothing parameters. Changes take effect immediately for subsequent time steps.

**Reloadable Velo@ parameters:**
- Solver type (Jacobi, BiCGSTAB)
- Convergence tolerances (`defCrit`, `MinDef`)
- Multigrid levels and cycle types
- Smoother iterations and relaxation parameters
- Coarse solver settings

**Example:**
```bash
# Switch to tighter convergence criteria
cat > _data/tighter_velo.dat << EOF
Velo@defCrit = 1d-8
Velo@MinDef = 1d-10
Velo@MGMaxIterCyc = 50
EOF

echo "Reload_Velo = _data/tighter_velo.dat" > _data/ProcCtrl.txt
```

**Use cases:**
- Tighten solver tolerances if solution diverges
- Adjust multigrid parameters for better performance
- Switch solver types during different flow regimes
- Experiment with solver settings without full restart

---

#### `Reload_Pres` - Reload Pressure Solver Parameters

**Syntax:**
```
Reload_Pres = parameter_file.dat
```

**Parameters:**
- `parameter_file.dat`: Path to parameter file containing `Pres@` parameters

**Description:**
Reloads pressure solver configuration. Similar to `Reload_Velo` but for the pressure Poisson equation solver.

**Reloadable Pres@ parameters:**
- Multigrid levels and iteration counts
- Convergence criteria
- Smoother settings
- Relaxation parameters
- Coarse solver type

**Example:**
```bash
# Increase pressure solver iterations
cat > _data/pres_params.dat << EOF
Pres@MGMaxIterCyc = 100
Pres@MGCriterion1 = 1d-7
EOF

echo "Reload_Pres = _data/pres_params.dat" > _data/ProcCtrl.txt
```

**Use cases:**
- Adjust pressure solver for incompressibility constraint
- Fine-tune performance on specific problems
- Debug solver convergence issues

---

#### `Reload_SimPar` - Reload Simulation Parameters

**Syntax:**
```
Reload_SimPar = parameter_file.dat
```

**Parameters:**
- `parameter_file.dat`: Path to parameter file containing `SimPar@` parameters

**Description:**
Reloads general simulation parameters including time stepping, output settings, and mesh levels. This is the most comprehensive reload option.

**Reloadable SimPar@ parameters:**
- Time step size (`TimeStep`)
- Time discretization scheme (`TimeScheme`)
- Output frequencies (`OutputFreq`, `BackUpFreq`)
- Mesh adaptation settings
- Matrix renewal strategy
- Physical model switches (tracer, viscoelastic, etc.)

**Example:**
```bash
# Reduce time step for stability
cat > _data/simpar.dat << EOF
SimPar@TimeStep = 0.001d0
SimPar@OutputFreq = 0.01d0
EOF

echo "Reload_SimPar = _data/simpar.dat" > _data/ProcCtrl.txt
```

**Use cases:**
- Adapt time step based on CFL condition
- Change output frequency
- Enable/disable physical models dynamically
- Adjust backup frequency

**⚠️ Important:** Some changes (like mesh levels) may require consistent solution state. Test carefully.

---

#### `Reload_Prop` - Reload Physical Properties

**Syntax:**
```
Reload_Prop = parameter_file.dat
```

**Parameters:**
- `parameter_file.dat`: Path to parameter file containing `Prop@` parameters

**Description:**
Reloads physical properties such as fluid density, viscosity, gravity, and material parameters. Useful for simulating changing physical conditions or parameter studies.

**Reloadable Prop@ parameters:**
- Fluid density (`Density`)
- Dynamic viscosity (`Viscosity`)
- Gravity vector (`Gravity`)
- Surface tension (`Sigma`)
- Diffusion coefficients (`DiffCoeff`)
- Non-Newtonian parameters (`PowerLawExp`)
- Viscoelastic parameters (`ViscoLambda`, `ViscoAlphaImp`, etc.)

**Example:**
```bash
# Change to different fluid properties
cat > _data/water_props.dat << EOF
Prop@Density = 1000d0,1000d0
Prop@Viscosity = 0.001d0,0.001d0
Prop@Gravity = 0d0,0d0,-9.81d0
EOF

echo "Reload_Prop = _data/water_props.dat" > _data/ProcCtrl.txt
```

**Use cases:**
- Simulate heating/cooling (changing viscosity)
- Test sensitivity to physical parameters
- Model multi-phase flows with property transitions
- Gradually ramp gravity or other forcing terms

---

## Practical Usage Guide

### Basic Workflow

1. **Start your simulation** as normal:
   ```bash
   mpirun -np 4 ./applications/q2p1_fc2/q2p1_fc2
   ```

2. **Monitor the simulation** (check log files, protocol files)

3. **Issue commands** by creating `_data/ProcCtrl.txt`:
   ```bash
   echo "Dump_Out = 10" > _data/ProcCtrl.txt
   ```

4. **Wait for next ProcessControl check** (typically each time step or output cycle)

5. **Commands execute and file is cleared automatically**

6. **Verify results** in output/log files

### Common Scenarios

#### Scenario 1: Create Regular Checkpoints

During a long simulation, create checkpoints every few hours:

```bash
# Script to run alongside simulation
while true; do
  sleep 3600  # Wait 1 hour
  TIMESTAMP=$(date +%H)
  echo "Dump_Out = ${TIMESTAMP}" > _data/ProcCtrl.txt
  echo "Created checkpoint at hour ${TIMESTAMP}"
done
```

#### Scenario 2: Adjust Time Step on the Fly

If you notice stability issues:

```bash
# Create reduced time step parameter file
cat > _data/smaller_dt.dat << EOF
SimPar@TimeStep = 0.0001d0
EOF

# Apply the change
echo "Reload_SimPar = _data/smaller_dt.dat" > _data/ProcCtrl.txt
```

#### Scenario 3: Generate Intermediate Visualization

Produce visualization output without waiting for scheduled output:

```bash
echo "GMV_Out = 9999" > _data/ProcCtrl.txt
```

#### Scenario 4: Parameter Sweep During Single Run

Test different viscosities sequentially:

```bash
#!/bin/bash
# Viscosity parameter sweep script

for visc in 0.001 0.01 0.1; do
  # Create parameter file
  cat > _data/visc_${visc}.dat << EOF
Prop@Viscosity = ${visc}d0,${visc}d0
EOF

  # Reload properties
  echo "Reload_Prop = _data/visc_${visc}.dat" > _data/ProcCtrl.txt

  # Wait for change to take effect and simulation to stabilize
  sleep 300  # 5 minutes

  # Save checkpoint
  echo "Dump_Out = ${visc//./}" > _data/ProcCtrl.txt

  # Generate output
  sleep 60
  echo "GMV_Out = ${visc//./}" > _data/ProcCtrl.txt

  sleep 60
done
```

#### Scenario 5: Emergency Recovery

If simulation becomes unstable:

```bash
# 1. Load last good checkpoint
echo "Dump_In = sol.10" > _data/ProcCtrl.txt

# Wait for load to complete (check logs)
sleep 30

# 2. Apply more conservative settings
cat > _data/conservative.dat << EOF
SimPar@TimeStep = 0.0005d0
Velo@defCrit = 1d-8
Pres@MGMaxIterCyc = 100
EOF

echo "Reload_SimPar = _data/conservative.dat" > _data/ProcCtrl.txt
sleep 10
echo "Reload_Velo = _data/conservative.dat" > _data/ProcCtrl.txt
sleep 10
echo "Reload_Pres = _data/conservative.dat" > _data/ProcCtrl.txt
```

---

## Best Practices

### ✅ Do's

1. **Use absolute or relative paths** consistently for parameter files
   ```bash
   echo "Reload_Velo = _data/params.dat" > _data/ProcCtrl.txt
   ```

2. **Check log files** after issuing commands to verify execution
   - Look for messages like "Reloading Velocity parameters from file..."

3. **Test parameter changes** on small problems first
   - Verify compatibility before applying to production runs

4. **Keep backup parameter files** organized
   ```
   _data/
   ├── q2p1_param.dat          # Original
   ├── params_conservative.dat # Tested alternatives
   ├── params_aggressive.dat
   └── params_test.dat
   ```

5. **Document your changes** in simulation logs
   ```bash
   echo "# Changed to smaller dt at $(date)" >> simulation_log.txt
   echo "Reload_SimPar = _data/smaller_dt.dat" > _data/ProcCtrl.txt
   ```

6. **Use version control** for parameter files
   ```bash
   git add _data/params_*.dat
   git commit -m "Parameter configurations for run XYZ"
   ```

### ❌ Don'ts

1. **Don't modify ProcCtrl.txt while solver is reading it**
   - Wait until file is cleared before issuing next command

2. **Don't reload incompatible states**
   - Mesh topology/decomposition must match for `Dump_In`

3. **Don't change too many parameters at once**
   - Makes debugging difficult if something goes wrong

4. **Don't forget MPI synchronization time**
   - Large parallel jobs may take longer to process commands

5. **Don't rely on timing**
   - ProcessControl checks happen at solver-determined intervals
   - Not suitable for real-time control with precise timing

6. **Don't use for initial setup**
   - Set initial parameters in `q2p1_param.dat`
   - Use ProcessControl only for runtime adjustments

---

## Advanced Usage

### Conditional Command Execution

Use shell scripts to issue commands based on simulation state:

```bash
#!/bin/bash
# Monitor and adapt based on protocol file

while true; do
  # Check if residual is too high (example metric)
  residual=$(tail -n 1 _data/prot.txt | awk '{print $5}')

  if (( $(echo "$residual > 1e-3" | bc -l) )); then
    echo "High residual detected, tightening solver..."
    echo "Reload_Velo = _data/tight_tolerance.dat" > _data/ProcCtrl.txt
  fi

  sleep 60
done
```

### Automated Checkpoint Strategy

Implement smart checkpointing based on simulation progress:

```bash
#!/bin/bash
# Checkpoint at specific time values

TARGET_TIMES=(0.1 0.5 1.0 2.0 5.0 10.0)
CHECKPOINT_NUM=0

while true; do
  # Read current simulation time from protocol file
  current_time=$(tail -n 1 _data/prot.txt | awk '{print $2}')

  for target in "${TARGET_TIMES[@]}"; do
    if (( $(echo "$current_time >= $target" | bc -l) )) && \
       [ ! -f "_data/.checkpoint_${target}_done" ]; then

      echo "Dump_Out = ${CHECKPOINT_NUM}" > _data/ProcCtrl.txt
      echo "GMV_Out = ${CHECKPOINT_NUM}" >> _data/ProcCtrl.txt
      touch "_data/.checkpoint_${target}_done"

      echo "Created checkpoint at t=${target}"
      CHECKPOINT_NUM=$((CHECKPOINT_NUM + 1))
    fi
  done

  sleep 30
done
```

### Parameter Optimization Loop

Automated parameter tuning during simulation:

```bash
#!/bin/bash
# Gradually increase Reynolds number

RE_VALUES=(100 200 500 1000 2000)

for RE in "${RE_VALUES[@]}"; do
  VISC=$(echo "scale=10; 1.0 / $RE" | bc)

  cat > _data/re_${RE}.dat << EOF
Prop@Viscosity = ${VISC}d0,${VISC}d0
EOF

  echo "Reload_Prop = _data/re_${RE}.dat" > _data/ProcCtrl.txt

  # Wait for flow to reach steady state (monitor residuals)
  echo "Testing Re=${RE}..."
  sleep 600  # 10 minutes

  # Save state
  echo "Dump_Out = ${RE}" > _data/ProcCtrl.txt
  sleep 60
  echo "GMV_Out = ${RE}" > _data/ProcCtrl.txt
  sleep 60
done
```

---

## Troubleshooting

### Commands Not Executing

**Symptom:** Commands in `ProcCtrl.txt` are ignored

**Possible causes:**
1. File permissions issue - ensure readable by MPI processes
2. File located in wrong directory - must be `_data/ProcCtrl.txt`
3. ProcessControl not called in simulation loop
4. Syntax error in command

**Solution:**
```bash
# Check file exists and has correct permissions
ls -la _data/ProcCtrl.txt

# Verify format (no spaces around =)
cat _data/ProcCtrl.txt

# Check protocol file for error messages
tail -f _data/prot.txt
```

### Parameter Reload Has No Effect

**Symptom:** Parameters reloaded but behavior unchanged

**Possible causes:**
1. Wrong parameter file or section
2. Parameter not actually reloadable (hardcoded)
3. Parameter cached/computed once at startup
4. Typo in parameter name

**Solution:**
```bash
# Verify file contains correct section
grep "Velo@" _data/params.dat

# Check protocol file for "Reloading..." message
grep "Reloading" _data/prot.txt

# Try with known working parameter
echo "Reload_Velo = _data/q2p1_param.dat" > _data/ProcCtrl.txt
```

### Dump_In Fails

**Symptom:** Loading checkpoint causes crash or errors

**Possible causes:**
1. Mesh incompatibility (different refinement level)
2. MPI decomposition mismatch (different number of processes)
3. File corruption
4. Wrong file format

**Solution:**
```bash
# Verify checkpoint file exists and is readable
ls -lh sol.*

# Check number of processes matches
# Must use same number as when checkpoint was created

# Try loading on same configuration that created it
```

### File Keeps Getting Recreated

**Symptom:** `ProcCtrl.txt` appears again after being cleared

**Possible causes:**
1. Multiple scripts writing to the file
2. Timing issue with file system
3. Automated monitoring script running

**Solution:**
```bash
# Check for other processes writing to file
lsof _data/ProcCtrl.txt

# Ensure only one control script is running
ps aux | grep ProcCtrl
```

---

## Technical Implementation Details

### Source Code
- **File:** `source/ProcCtrl.f90`
- **Main routine:** `ProcessControl(MFILE,MTERM)`
- **Finalization:** `Finalize()` and `Finalize_Particles()`

### Call Frequency
ProcessControl is typically called:
- Once per time step in explicit schemes
- Once per outer iteration in implicit schemes
- After each output cycle
- Frequency depends on application implementation

### MPI Synchronization
- Only `showid` process reads the file
- Commands broadcast to all processes via `ShareValueC_myMPI`
- All processes execute commands in sync
- `Barrier_myMPI()` ensures consistency

### File Handling
- File unit: 989
- Action: `READWRITE`
- Auto-clear: File rewritten with "Nothing to be done..." after execution

---

## Related Documentation

- **Parameter Reference:** See `PARAMETER_REFERENCE.md` for all `Velo@`, `Pres@`, `Prop@`, and `SimPar@` parameters
- **Checkpoint Format:** See application-specific documentation for `SolToFile`/`SolFromFile` format
- **Output Formats:** See `Output_Profiles` documentation for GMV/VTU format details

---

## Summary

The ProcessControl system provides powerful interactive control over running simulations:

| Command | Purpose | Use Case |
|---------|---------|----------|
| `Dump_Out` | Save checkpoint | Backup, recovery points |
| `Dump_In` | Load checkpoint | Restart, branching scenarios |
| `GMV_Out` | Generate visualization | Intermediate output |
| `Reload_Velo` | Update velocity solver | Adjust convergence, performance |
| `Reload_Pres` | Update pressure solver | Adjust pressure solution |
| `Reload_SimPar` | Update simulation settings | Time step, output frequency |
| `Reload_Prop` | Update physical properties | Material properties, gravity |

**Key benefits:**
- ✅ No simulation restart required
- ✅ Interactive parameter tuning
- ✅ Flexible checkpointing strategy
- ✅ On-demand visualization output
- ✅ MPI-safe parallel execution

**Perfect for:**
- Long HPC runs that need mid-flight adjustments
- Parameter studies without multiple runs
- Adaptive simulations responding to flow features
- Creating backups during critical phases
- Debugging solver convergence issues

Start with simple commands (`Dump_Out`, `GMV_Out`) and gradually explore parameter reloading as you become familiar with the system!
