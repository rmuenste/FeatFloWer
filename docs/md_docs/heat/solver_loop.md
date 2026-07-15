# heat Solver Loop

## Entry Points

- `applications/heat/heat.f90`
- `source/src_LinSc/LinSc_transport_extensions.f90`
- `source/src_LinSc/LinSc_def.f90`
- `source/src_LinSc/LinSc_user.f90`

## Main Time Loop

After initialization, `heat.f90` runs:

```text
do itns = 1, nitns
    itnsr = 0
    timnsh = timens
    dt = tstep
    timens = timens + dt
    inonln_u = 2
    inonln_t = 2

    Transport_LinScalar_EWIKON(
        Boundary_LinSc_Val_EWIKON,
        AddSource_EWIKON,
        ufile,
        inonln_t)

    postprocessing_app_heat
    print_time
    handle_statistics
    istep_ns = istep_ns + 1

    exit on final time, DivergedSolution, or ConvergedSolution
end do
```

The velocity fields passed into scalar assembly are inherited from shared Q2/P1
state. The heat app itself advances only the Q1 scalar temperature field.

## EWIKON Scalar Transport

`Transport_LinScalar_EWIKON` performs one time step:

```text
NLMAX = NLMAX + 1
thstep = 0.5 * tstep
Tracer%prm%AFC = .false.

on worker ranks:
    save old scalar solution
    Create_RhoCpConvMat_Ewikon
    clear defect
    Matdef_LinScalar_EWIKON(Tracer, 1, 0)
    AddSource_EWIKON
    AddBoundaryHeatFlux(1)
    Boundary_LinSc_Def
    store constant RHS
    Matdef_LinScalar_EWIKON(Tracer, -1, 1)
    AddBoundaryHeatFlux(2)
    E011Sum
    Boundary_LinSc_Def
    Resdfk_General_LinScalar

compute scalar defect criterion
loop nonlinear scalar iterations:
    Solve_General_LinScalar
    rebuild defect
    Matdef_LinScalar_EWIKON(Tracer, -1, 0)
    E011Sum
    Boundary_LinSc_Def
    Resdfk_General_LinScalar
    Protocol_linScalar
    stop when defect criterion and NLmin are satisfied

sum heat flux/source integrals
IntegrateOutputQuantities
CreateSensorOutputs
copy Tracer solution to Temperature
update time-averaged Temperature_AVG after 10 percent of run
NLMAX = NLMAX - 1
```

## Boundary Values

`Boundary_LinSc_Val_EWIKON` maps `Tracer%knpr` markers to temperature values.
Some markers are hard-coded legacy values, while markers `3` and `4` use the
segment initial temperature:

```text
iSeg = myHeatObjects%Segment(i)
Tracer%val(NLMAX)%x(i) = mySigma%mySegment(iSeg)%InitTemp
```

Changing boundary semantics usually requires checking both marker creation
(`LinSc_Knpr`/geometry classification) and this callback.

## Source Terms

`AddSource_EWIKON` maps each scalar node to its heat segment:

```text
iSeg = myHeatObjects%Segment(i)
iMat = segment material index or air material 0
if ObjectType is WIRE or MELT:
    dSource = 1e10 * UseHeatSource / segment Volume
    Tracer%def(i) += MLmat(i) * dSource * tstep
```

`UseHeatSource` is initialized from `VolumetricHeatSourceMax` and can be changed
by the regulation/convergence machinery.

## Heat Fluxes And Output Quantities

`Transport_LinScalar_EWIKON` adds explicit and implicit boundary heat fluxes via
`AddBoundaryHeatFlux(1)` and `AddBoundaryHeatFlux(2)`. At the end of each step
it reduces:

- conductive area/flux
- convective area/flux
- radiative area/flux
- total heat source

These are printed with the `Time[s]_{Conductive,Convective,Radiative}...`
protocol line.
