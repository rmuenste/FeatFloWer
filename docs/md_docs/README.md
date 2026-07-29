# FeatFloWer Documentation Index

This directory contains implementation references, operational guides, and
focused engineering analyses for FeatFloWer. Use this index to find the relevant
document quickly. Documents marked as analysis, roadmap, or experimental should
not be treated as canonical behavior without checking the current source.

## Build and Operational Guides

- [containerization.md](containerization.md): Build and run FeatFloWer with the bundled Docker setup.
- [low_noise_build.md](low_noise_build.md): Configure low-noise builds while retaining complete logs and diagnostics.
- [guide_01_q2p1_fc_ext_cylinder_benchmark_from_scratch.md](guide_01_q2p1_fc_ext_cylinder_benchmark_from_scratch.md): Full build, partition, run, and validation workflow for the `q2p1_fc_ext` cylinder benchmark.
- [guide_02_q2p1_bench_sedimentation_pe_serial_from_scratch.md](guide_02_q2p1_bench_sedimentation_pe_serial_from_scratch.md): Build and stage the PE-serial sedimentation benchmark.
- [guide_03_q2p1_sse_tse_gendie_from_scratch.md](guide_03_q2p1_sse_tse_gendie_from_scratch.md): Configure, stage, and submit `q2p1_sse` SSE/TSE workflows.
- [guide_04_q2p1_atc_pe_serial_fbm_kvel_from_scratch.md](guide_04_q2p1_atc_pe_serial_fbm_kvel_from_scratch.md): Build `q2p1_ATC` with PE serial mode, CGAL, and FBM acceleration.
- [guide_05_agent_new_application_workflow.md](guide_05_agent_new_application_workflow.md): Checklist for creating and integrating a new Q2/P1 application, including PE bridges.
- [guide_06_q2p1_bench_fluidization_pe_serial_fbm_kvel_gcc14_from_scratch.md](guide_06_q2p1_bench_fluidization_pe_serial_fbm_kvel_gcc14_from_scratch.md): GCC 14 workflow for the PE-serial fluidization benchmark.

## Tools and Testing

- [automated_test_system.md](automated_test_system.md): Current architecture, data model, and remaining work for the `featflower-test` runner.
- [featflower_test_usage_guide.md](featflower_test_usage_guide.md): Installation and practical CLI usage for YAML-driven local and SLURM tests.
- [featflower_partitioner_usage_guide.md](featflower_partitioner_usage_guide.md): Installation and usage of the maintained Python partitioner package.
- [documentation_audit_status.md](documentation_audit_status.md): Record of the June 2026 documentation cleanup and removed historical notes.

## Configuration and Runtime Control

- [parameter_reference.md](parameter_reference.md): Reference for `q2p1_param.dat` categories and commonly used parameters.
- [process_control.md](process_control.md): Runtime control-file commands for output, checkpointing, and parameter updates.

## Mesh and Parallel Decomposition

- [mesh_structure.md](mesh_structure.md): Mesh types, multilevel layout, connectivity arrays, and access conventions.
- [domain_decomposition.md](domain_decomposition.md): Partitioning strategy, master-worker layout, and parallel mesh workflow.
- [communication_structures.md](communication_structures.md): MPI communication structures produced by `PARENTCOMM` and related routines.
- [dof_mapping.md](dof_mapping.md): Analysis of `NDFGL` and hierarchical Q2 degree-of-freedom numbering.

## QuadSc Matrices and Solvers

- [quadsc_current_implementation.md](quadsc_current_implementation.md): Current module organization, CSR conventions, assembly ownership, and coarse solvers.
- [quadsc_refactoring_status.md](quadsc_refactoring_status.md): Completed structural refactoring and remaining modernization work.
- [matrix_structures_guide.md](matrix_structures_guide.md): Detailed matrix families, allocation patterns, and level-dependent storage.
- [matrix_assembly_guide.md](matrix_assembly_guide.md): Q2/P1 matrix assembly workflow, generic helpers, and parallel handling.
- [hypre_csr_analysis.md](hypre_csr_analysis.md): Pressure Schur-complement CSR construction and HYPRE conversion.
- [hypre_gpu_stage1_tardis.md](hypre_gpu_stage1_tardis.md): External Hypre CUDA build, managed-memory interface, and Tardis cylinder validation.
- [solver_baseline_phase0.md](solver_baseline_phase0.md): FF_TIMING instrumentation, CPU baseline job harness, and summarizer for the solver-library evaluation.
- [solver_phase1_direct_solvers.md](solver_phase1_direct_solvers.md): External SuiteSparse UMFPACK for coarse types 2/4 and MKL PARDISO as new type 9.
- [solver_phase2_velocity_coarse.md](solver_phase2_velocity_coarse.md): External GCC-compatible MUMPS (drops the Intel requirement) and the velocity coarse-solver sweep.
- [solver_phase3_fine_level_pressure.md](solver_phase3_fine_level_pressure.md): Fine-level pressure solve via Hypre GMRES+BoomerAMG (Pres@MGCrsSolverType = 10).
- [fetchcontent_dependencies.md](fetchcontent_dependencies.md): Provider pattern for hypre/SuiteSparse/MKL, pinned FetchContent versions and hashes, and offline/cluster pre-staging.
- [solver_type_validation.md](solver_type_validation.md): Startup validation of Pres/Velo@MGCrsSolverType against compiled-in solver libraries (MPI_Abort with a clear message instead of late STOPs or silent no-ops).
- [umfpack_hypre_analysis.md](umfpack_hypre_analysis.md): Comparison of UMFPACK and HYPRE coarse-solver interfaces.
- [ff_theta_solution_scheme.md](ff_theta_solution_scheme.md): Analytical walkthrough of the fractional-step theta transport routine; verify details against current source.

## Finite Element Machinery

- [fem_machinery.md](fem_machinery.md): Q2 element evaluation, derivative switches, mapping, and recurring FEM implementation conventions.
- [velocity_midpoint_evaluation.md](velocity_midpoint_evaluation.md): Q2 velocity interpolation at element centers with reference implementation notes.
- [strain_rate_dissipation_calculation.md](strain_rate_dissipation_calculation.md): Strain-rate dissipation integral derivation and example implementation.

## Fictitious Boundary Method

- [fbm_implementation_report.md](fbm_implementation_report.md): Map from FBM concepts to the geometry, constraint, solve, and force paths in the code.
- [hydrodynamic_force_computation.md](hydrodynamic_force_computation.md): Derivation and implementation analysis of alpha-gradient hydrodynamic force integration.
- [fbm_acceleration_usage.md](fbm_acceleration_usage.md): Compile-time and runtime control of HashGrid and KVEL acceleration.
- [hashgrid_verification.md](hashgrid_verification.md): Runtime comparison of accelerated point queries with the baseline implementation.
- [kvel_force_acceleration.md](kvel_force_acceleration.md): KVEL/KEEL/KAAL candidate-element acceleration for particle force integration.
- [particle_reynolds_evaluation.md](particle_reynolds_evaluation.md): Center, interface, and extended-interface particle Reynolds-number methods.

## PE Integration

- [pe_usage_overview.md](pe_usage_overview.md): When to enable PE, build modes, and the high-level coupling flow.
- [pe_initialization.md](pe_initialization.md): Detailed Fortran/C++ initialization flow, communicator handling, and setup hooks.
- [mapping.md](mapping.md): Mapping between Fortran particle-query interfaces and C++ implementations.

## Euler-Lagrange and Restart Development

- [euler_lagrange_development.md](euler_lagrange_development.md): Development analysis and staged roadmap for unresolved Euler-Lagrange coupling; requires manual review before use as current guidance.
- [el_frozen_field_application_parameters.md](el_frozen_field_application_parameters.md): Runtime parameter reference for `q2p1_el_frozen_trace`.
- [repartition_dump_workflow_notes.md](repartition_dump_workflow_notes.md): Investigation of legacy dump merging and ambiguous Q2 reconstruction across partitions.
- [provenance_dump_format.md](provenance_dump_format.md): Experimental ownership-aware dump format, reader/writer behavior, and validation limits.

## Domain-Specific Analysis

- [atc_metrics.md](atc_metrics.md): Interpretation of qualitative suspension-flow-map regimes and possible ATC metrics.
