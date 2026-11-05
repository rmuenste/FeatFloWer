# MANIFEST — every file, its origin and a one-line description

Origin paths are absolute paths in the FeatFloWer checkout (`/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer`, branch `feature/dns-validation`, 2026-10-02) or in the mesh repository / website checkout named. "written"/"generated" = created for this hand-off. Sizes in bytes.

| file | size | origin | description |
|---|---|---|---|
| `CASES.md` | 11239 | `written` | Per-case parameters, derived numbers, experiment, FeatFloWer results per level/dt, run-length recommendations |
| `LITERATURE.md` | 5758 | `written` | Full citations and what each source is used for |
| `PITFALLS.md` | 18116 | `written` | Every campaign difficulty as symptom/cause/check, library-agnostic |
| `README.md` | 19890 | `written for this hand-off` | Hand-off entry point: benchmark, cases, FeatFloWer configuration, certified results, reference discipline, source conflicts |
| `campaign_docs/dns_practitioners_guide.md` | 35881 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/docs/md_docs/dns_practitioners_guide.md (v2.5, 2026-10-02)` | Full practitioner guide (sections 1-7, 10 are the ones cited) |
| `campaign_docs/plan_extracts.md` | 9074 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/dns-validation-campaign-plan.md (line-range extracts)` | Plan sections D0.4, D1.3/D1.4, section 11 risks, section 13 literature |
| `literature/brenner_1961.pdf` | 809551 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/literature/brenner_1961.pdf` | Brenner, Chem. Eng. Sci. 16 (1961) 242-251 (copyrighted, internal use) |
| `literature/ten_cate_nieuwstad_derksen_vandenakker_2002_phys_fluids_14_4012.pdf` | 1176319 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_piv.pdf` | ten Cate et al., Phys. Fluids 14 (2002) 4012-4025 (copyrighted, internal use) |
| `literature/uhlmann_dusek_2014.pdf` | 7262967 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/literature/uhlmann_dusek_2014.pdf` | Uhlmann & Dusek, IJMF 59 (2014) 221-243 (copyrighted, internal use) |
| `mesh/MESH.md` | 7333 | `written` | Both meshes, refinement rule, element counts, sphere placement, partitioning |
| `mesh/fullbox_ten_cate_mesh_v1/box.tri` | 27428 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/box.tri` | Full-box coarse mesh 6x6x9, 324 hex / 490 vertices (not used for certified rows) |
| `mesh/fullbox_ten_cate_mesh_v1/box.vtk` | 29037 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/box.vtk` | VTK view of the full-box coarse mesh |
| `mesh/fullbox_ten_cate_mesh_v1/file.prj` | 62 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/file.prj` | Project file |
| `mesh/fullbox_ten_cate_mesh_v1/generator/hex_ex.py` | 10535 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/quadextrude/hex_ex.py` | Hex-extrusion mesh generator (hex_ex.py and its mesh package) that exported box.tri |
| `mesh/fullbox_ten_cate_mesh_v1/generator/mesh/__init__.py` | 812 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/quadextrude/mesh/__init__.py` | Hex-extrusion mesh generator (hex_ex.py and its mesh package) that exported box.tri |
| `mesh/fullbox_ten_cate_mesh_v1/generator/mesh/mesh.py` | 49370 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/quadextrude/mesh/mesh.py` | Hex-extrusion mesh generator (hex_ex.py and its mesh package) that exported box.tri |
| `mesh/fullbox_ten_cate_mesh_v1/generator/mesh/mesh_io.py` | 24706 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/quadextrude/mesh/mesh_io.py` | Hex-extrusion mesh generator (hex_ex.py and its mesh package) that exported box.tri |
| `mesh/fullbox_ten_cate_mesh_v1/generator/simple_ex.py` | 6921 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/quadextrude/simple_ex.py` | Hex-extrusion mesh generator (hex_ex.py and its mesh package) that exported box.tri |
| `mesh/fullbox_ten_cate_mesh_v1/generator/tri2vtk_converter.py` | 4603 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/quadextrude/tri2vtk_converter.py` | Hex-extrusion mesh generator (hex_ex.py and its mesh package) that exported box.tri |
| `mesh/fullbox_ten_cate_mesh_v1/xmax.par` | 296 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/xmax.par` | Boundary patch of the full box (all Wall) |
| `mesh/fullbox_ten_cate_mesh_v1/xmin.par` | 291 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/xmin.par` | Boundary patch of the full box (all Wall) |
| `mesh/fullbox_ten_cate_mesh_v1/ymax.par` | 297 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/ymax.par` | Boundary patch of the full box (all Wall) |
| `mesh/fullbox_ten_cate_mesh_v1/ymin.par` | 286 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/ymin.par` | Boundary patch of the full box (all Wall) |
| `mesh/fullbox_ten_cate_mesh_v1/zmax.par` | 228 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/zmax.par` | Boundary patch of the full box (all Wall) |
| `mesh/fullbox_ten_cate_mesh_v1/zmin.par` | 166 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/ten_cate_mesh_v1/zmin.par` | Boundary patch of the full box (all Wall) |
| `mesh/quarterbox_benchSym/bench.prj` | 57 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/bench.prj` | Project file listing mesh + par files |
| `mesh/quarterbox_benchSym/bot.par` | 235 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/bot.par` | Boundary patch of the quarter box |
| `mesh/quarterbox_benchSym/grid.vtu` | 110619 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/grid.vtu` | VTK view of the quarter-box coarse mesh |
| `mesh/quarterbox_benchSym/mesh.tri` | 67799 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/mesh.tri (== repo benchSym/mesh.tri == rundir _mesh/NEWFAC/GRID.tri)` | Quarter-box coarse mesh of the certified runs, 876 hex / 1239 vertices |
| `mesh/quarterbox_benchSym/top.par` | 279 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/top.par` | Boundary patch of the quarter box |
| `mesh/quarterbox_benchSym/x.par` | 1286 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/x.par` | Boundary patch of the quarter box |
| `mesh/quarterbox_benchSym/xwall.par` | 342 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/xwall.par` | Boundary patch of the quarter box |
| `mesh/quarterbox_benchSym/y.par` | 1091 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/y.par` | Boundary patch of the quarter box |
| `mesh/quarterbox_benchSym/ywall.par` | 342 | `/data/warehouse17/rmuenste/work/MESH/mesh_repo/benchSym/ywall.par` | Boundary patch of the quarter box |
| `reference/README.md` | 4887 | `written` | Provenance of every data file under reference/ |
| `reference/datasheet_rows_tencate.csv` | 29656 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/docs/md_docs/dns_validation_datasheet.csv (grep by row id)` | Verbatim datasheet rows cited in the hand-off (header + 49 rows) |
| `reference/experiment/case_E1_h.csv` | 1826 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/, == ff-redesign site source-data)` | Digitised experimental trajectory E1: t [s], h/d |
| `reference/experiment/case_E2_h.csv` | 1144 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/, == ff-redesign site source-data)` | Digitised experimental trajectory E2: t [s], h/d |
| `reference/experiment/case_E3_h.csv` | 1112 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/, == ff-redesign site source-data)` | Digitised experimental trajectory E3: t [s], h/d |
| `reference/experiment/case_E4_h.csv` | 1171 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/, == ff-redesign site source-data)` | Digitised experimental trajectory E4: t [s], h/d |
| `reference/experiment/ref_E1.dat` | 1049 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/)` | Digitised experimental settling velocity E1: t [s], v_z [m/s] |
| `reference/experiment/ref_E2.dat` | 1082 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/)` | Digitised experimental settling velocity E2: t [s], v_z [m/s] |
| `reference/experiment/ref_E3.dat` | 889 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/)` | Digitised experimental settling velocity E3: t [s], v_z [m/s] |
| `reference/experiment/ref_E4.dat` | 644 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip -> tenCateData/ (== tc-ref/)` | Digitised experimental settling velocity E4: t [s], v_z [m/s] |
| `reference/experiment/tc-ref_README.md` | 2405 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tc-ref/README.md` | Repo README for the digitised data incl. the 2026-08-01 digitisation audit |
| `reference/experiment/tenCateData.zip` | 4715 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tenCateData.zip` | Digitised PIV curves, original archive (2026-07-31) |
| `reference/featflower/decks/e1_l3/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/cube.json` | Rigid bodies (sphere + bottom plane) of run e1_l3 |
| `reference/featflower/decks/e1_l3/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/example.json` | PE/particle configuration of run e1_l3 |
| `reference/featflower/decks/e1_l3/job.sbatch` | 611 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/job.sbatch` | Slurm job script of run e1_l3 (binary, ranks, memory) |
| `reference/featflower/decks/e1_l3/q2p1_param.dat` | 2982 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/_data/q2p1_param.dat` | CFD deck of run e1_l3 |
| `reference/featflower/decks/e1_l3/start/data.TXT` | 1732 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/start/data.TXT` | PE boilerplate init file read unconditionally by the application (values superseded by example.json) |
| `reference/featflower/decks/e1_l3/start/sampleRigidBody.xml` | 1247 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/start/sampleRigidBody.xml` | PE boilerplate init file read unconditionally by the application (values superseded by example.json) |
| `reference/featflower/decks/e1_l3_dt0p5_sync/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_dt0p5_sync/cube.json` | Rigid bodies (sphere + bottom plane) of run e1_l3_dt0p5_sync |
| `reference/featflower/decks/e1_l3_dt0p5_sync/example.json` | 1294 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_dt0p5_sync/example.json` | PE/particle configuration of run e1_l3_dt0p5_sync |
| `reference/featflower/decks/e1_l3_dt0p5_sync/job.sbatch` | 834 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_dt0p5_sync/job.sbatch` | Slurm job script of run e1_l3_dt0p5_sync (binary, ranks, memory) |
| `reference/featflower/decks/e1_l3_dt0p5_sync/q2p1_param.dat` | 2906 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_dt0p5_sync/_data/q2p1_param.dat` | CFD deck of run e1_l3_dt0p5_sync |
| `reference/featflower/decks/e1_l3_g3def/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_g3def/cube.json` | Rigid bodies (sphere + bottom plane) of run e1_l3_g3def |
| `reference/featflower/decks/e1_l3_g3def/example.json` | 1331 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_g3def/example.json` | PE/particle configuration of run e1_l3_g3def |
| `reference/featflower/decks/e1_l3_g3def/job.sbatch` | 630 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_g3def/job.sbatch` | Slurm job script of run e1_l3_g3def (binary, ranks, memory) |
| `reference/featflower/decks/e1_l3_g3def/q2p1_param.dat` | 2982 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_g3def/_data/q2p1_param.dat` | CFD deck of run e1_l3_g3def |
| `reference/featflower/decks/e1_l4/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l4/cube.json` | Rigid bodies (sphere + bottom plane) of run e1_l4 |
| `reference/featflower/decks/e1_l4/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l4/example.json` | PE/particle configuration of run e1_l4 |
| `reference/featflower/decks/e1_l4/job.sbatch` | 610 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l4/job.sbatch` | Slurm job script of run e1_l4 (binary, ranks, memory) |
| `reference/featflower/decks/e1_l4/q2p1_param.dat` | 2982 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l4/_data/q2p1_param.dat` | CFD deck of run e1_l4 |
| `reference/featflower/decks/e2_l2/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l2/cube.json` | Rigid bodies (sphere + bottom plane) of run e2_l2 |
| `reference/featflower/decks/e2_l2/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l2/example.json` | PE/particle configuration of run e2_l2 |
| `reference/featflower/decks/e2_l2/job.sbatch` | 611 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l2/job.sbatch` | Slurm job script of run e2_l2 (binary, ranks, memory) |
| `reference/featflower/decks/e2_l2/q2p1_param.dat` | 2981 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l2/_data/q2p1_param.dat` | CFD deck of run e2_l2 |
| `reference/featflower/decks/e2_l3/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3/cube.json` | Rigid bodies (sphere + bottom plane) of run e2_l3 |
| `reference/featflower/decks/e2_l3/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3/example.json` | PE/particle configuration of run e2_l3 |
| `reference/featflower/decks/e2_l3/job.sbatch` | 611 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3/job.sbatch` | Slurm job script of run e2_l3 (binary, ranks, memory) |
| `reference/featflower/decks/e2_l3/q2p1_param.dat` | 2982 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3/_data/q2p1_param.dat` | CFD deck of run e2_l3 |
| `reference/featflower/decks/e2_l3_g3def/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3_g3def/cube.json` | Rigid bodies (sphere + bottom plane) of run e2_l3_g3def |
| `reference/featflower/decks/e2_l3_g3def/example.json` | 1331 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3_g3def/example.json` | PE/particle configuration of run e2_l3_g3def |
| `reference/featflower/decks/e2_l3_g3def/job.sbatch` | 630 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3_g3def/job.sbatch` | Slurm job script of run e2_l3_g3def (binary, ranks, memory) |
| `reference/featflower/decks/e2_l3_g3def/q2p1_param.dat` | 2982 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3_g3def/_data/q2p1_param.dat` | CFD deck of run e2_l3_g3def |
| `reference/featflower/decks/e2_l4/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l4/cube.json` | Rigid bodies (sphere + bottom plane) of run e2_l4 |
| `reference/featflower/decks/e2_l4/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l4/example.json` | PE/particle configuration of run e2_l4 |
| `reference/featflower/decks/e2_l4/job.sbatch` | 610 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l4/job.sbatch` | Slurm job script of run e2_l4 (binary, ranks, memory) |
| `reference/featflower/decks/e2_l4/q2p1_param.dat` | 2981 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l4/_data/q2p1_param.dat` | CFD deck of run e2_l4 |
| `reference/featflower/decks/e3_l2/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l2/cube.json` | Rigid bodies (sphere + bottom plane) of run e3_l2 |
| `reference/featflower/decks/e3_l2/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l2/example.json` | PE/particle configuration of run e3_l2 |
| `reference/featflower/decks/e3_l2/job.sbatch` | 611 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l2/job.sbatch` | Slurm job script of run e3_l2 (binary, ranks, memory) |
| `reference/featflower/decks/e3_l2/q2p1_param.dat` | 2981 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l2/_data/q2p1_param.dat` | CFD deck of run e3_l2 |
| `reference/featflower/decks/e3_l3/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l3/cube.json` | Rigid bodies (sphere + bottom plane) of run e3_l3 |
| `reference/featflower/decks/e3_l3/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l3/example.json` | PE/particle configuration of run e3_l3 |
| `reference/featflower/decks/e3_l3/job.sbatch` | 611 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l3/job.sbatch` | Slurm job script of run e3_l3 (binary, ranks, memory) |
| `reference/featflower/decks/e3_l3/q2p1_param.dat` | 2982 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l3/_data/q2p1_param.dat` | CFD deck of run e3_l3 |
| `reference/featflower/decks/e3_l4/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l4/cube.json` | Rigid bodies (sphere + bottom plane) of run e3_l4 |
| `reference/featflower/decks/e3_l4/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l4/example.json` | PE/particle configuration of run e3_l4 |
| `reference/featflower/decks/e3_l4/job.sbatch` | 610 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l4/job.sbatch` | Slurm job script of run e3_l4 (binary, ranks, memory) |
| `reference/featflower/decks/e3_l4/q2p1_param.dat` | 2981 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l4/_data/q2p1_param.dat` | CFD deck of run e3_l4 |
| `reference/featflower/decks/e4_l2/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l2/cube.json` | Rigid bodies (sphere + bottom plane) of run e4_l2 |
| `reference/featflower/decks/e4_l2/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l2/example.json` | PE/particle configuration of run e4_l2 |
| `reference/featflower/decks/e4_l2/job.sbatch` | 611 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l2/job.sbatch` | Slurm job script of run e4_l2 (binary, ranks, memory) |
| `reference/featflower/decks/e4_l2/q2p1_param.dat` | 2981 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l2/_data/q2p1_param.dat` | CFD deck of run e4_l2 |
| `reference/featflower/decks/e4_l3/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3/cube.json` | Rigid bodies (sphere + bottom plane) of run e4_l3 |
| `reference/featflower/decks/e4_l3/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3/example.json` | PE/particle configuration of run e4_l3 |
| `reference/featflower/decks/e4_l3/job.sbatch` | 611 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3/job.sbatch` | Slurm job script of run e4_l3 (binary, ranks, memory) |
| `reference/featflower/decks/e4_l3/q2p1_param.dat` | 2981 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3/_data/q2p1_param.dat` | CFD deck of run e4_l3 |
| `reference/featflower/decks/e4_l3_dt0p25_sync/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p25_sync/cube.json` | Rigid bodies (sphere + bottom plane) of run e4_l3_dt0p25_sync |
| `reference/featflower/decks/e4_l3_dt0p25_sync/example.json` | 1295 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p25_sync/example.json` | PE/particle configuration of run e4_l3_dt0p25_sync |
| `reference/featflower/decks/e4_l3_dt0p25_sync/job.sbatch` | 838 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p25_sync/job.sbatch` | Slurm job script of run e4_l3_dt0p25_sync (binary, ranks, memory) |
| `reference/featflower/decks/e4_l3_dt0p25_sync/q2p1_param.dat` | 2906 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p25_sync/_data/q2p1_param.dat` | CFD deck of run e4_l3_dt0p25_sync |
| `reference/featflower/decks/e4_l3_dt0p5_sync/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p5_sync/cube.json` | Rigid bodies (sphere + bottom plane) of run e4_l3_dt0p5_sync |
| `reference/featflower/decks/e4_l3_dt0p5_sync/example.json` | 1294 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p5_sync/example.json` | PE/particle configuration of run e4_l3_dt0p5_sync |
| `reference/featflower/decks/e4_l3_dt0p5_sync/job.sbatch` | 835 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p5_sync/job.sbatch` | Slurm job script of run e4_l3_dt0p5_sync (binary, ranks, memory) |
| `reference/featflower/decks/e4_l3_dt0p5_sync/q2p1_param.dat` | 2905 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p5_sync/_data/q2p1_param.dat` | CFD deck of run e4_l3_dt0p5_sync |
| `reference/featflower/decks/e4_l4/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4/cube.json` | Rigid bodies (sphere + bottom plane) of run e4_l4 |
| `reference/featflower/decks/e4_l4/example.json` | 1293 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4/example.json` | PE/particle configuration of run e4_l4 |
| `reference/featflower/decks/e4_l4/job.sbatch` | 610 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4/job.sbatch` | Slurm job script of run e4_l4 (binary, ranks, memory) |
| `reference/featflower/decks/e4_l4/q2p1_param.dat` | 2981 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4/_data/q2p1_param.dat` | CFD deck of run e4_l4 |
| `reference/featflower/decks/e4_l4_dt0p5_sync/cube.json` | 183 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4_dt0p5_sync/cube.json` | Rigid bodies (sphere + bottom plane) of run e4_l4_dt0p5_sync |
| `reference/featflower/decks/e4_l4_dt0p5_sync/example.json` | 1294 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4_dt0p5_sync/example.json` | PE/particle configuration of run e4_l4_dt0p5_sync |
| `reference/featflower/decks/e4_l4_dt0p5_sync/job.sbatch` | 834 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4_dt0p5_sync/job.sbatch` | Slurm job script of run e4_l4_dt0p5_sync (binary, ranks, memory) |
| `reference/featflower/decks/e4_l4_dt0p5_sync/q2p1_param.dat` | 2905 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4_dt0p5_sync/_data/q2p1_param.dat` | CFD deck of run e4_l4_dt0p5_sync |
| `reference/featflower/ff_E1_L3_dt0p5.csv` | 679770 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E1 level L3 dt 0.5 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E1_L3_dt1.csv` | 340058 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E1 level L3 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E1_L3_dt1_lubON.csv` | 340098 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E1 level L3 dt 1 ms, lubrication add-on ON (row d22_g3_tencate): t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E1_L4_dt1.csv` | 340058 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E1 level L4 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E2_L2_dt1.csv` | 213452 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E2 level L2 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E2_L3_dt1.csv` | 213450 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E2 level L3 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E2_L3_dt1_lubON.csv` | 213540 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E2 level L3 dt 1 ms, lubrication add-on ON (row d22_g3_tencate): t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E2_L4_dt1.csv` | 213658 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E2 level L4 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E3_L2_dt1.csv` | 142390 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E3 level L2 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E3_L3_dt1.csv` | 142395 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E3 level L3 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E3_L4_dt1.csv` | 142561 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E3 level L4 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E4_L2_dt1.csv` | 102968 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E4 level L2 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E4_L3_dt0p25.csv` | 410985 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E4 level L3 dt 0.25 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E4_L3_dt0p5.csv` | 205635 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E4 level L3 dt 0.5 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E4_L3_dt1.csv` | 102974 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E4 level L3 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E4_L4_dt0p5.csv` | 205595 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E4 level L4 dt 0.5 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_E4_L4_dt1.csv` | 103059 | `generated from raw_logs/particle_force_*.log` | FeatFloWer settling curve E4 level L4 dt 1 ms: t, z_center, gap, h/d, u_z, F_z |
| `reference/featflower/ff_peaks_summary.csv` | 2421 | `generated from raw_logs/ (extract_ff_curves.py)` | Peak, timing, touchdown, rest, final gap/force per certified run, with datasheet cross-check |
| `reference/featflower/raw_logs/particle_force_e1_l3.log` | 915946 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/particle_force.log` | Raw per-step particle log of certified run e1_l3 |
| `reference/featflower/raw_logs/particle_force_e1_l3_dt0p5_sync.log` | 1831846 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_dt0p5_sync/particle_force.log` | Raw per-step particle log of certified run e1_l3_dt0p5_sync |
| `reference/featflower/raw_logs/particle_force_e1_l3_g3def.log` | 915946 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_g3def/particle_force.log` | Raw per-step particle log of certified run e1_l3_g3def |
| `reference/featflower/raw_logs/particle_force_e1_l4.log` | 915946 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l4/particle_force.log` | Raw per-step particle log of certified run e1_l4 |
| `reference/featflower/raw_logs/particle_force_e2_l2.log` | 575146 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l2/particle_force.log` | Raw per-step particle log of certified run e2_l2 |
| `reference/featflower/raw_logs/particle_force_e2_l3.log` | 575146 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3/particle_force.log` | Raw per-step particle log of certified run e2_l3 |
| `reference/featflower/raw_logs/particle_force_e2_l3_g3def.log` | 575146 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l3_g3def/particle_force.log` | Raw per-step particle log of certified run e2_l3_g3def |
| `reference/featflower/raw_logs/particle_force_e2_l4.log` | 575146 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e2_l4/particle_force.log` | Raw per-step particle log of certified run e2_l4 |
| `reference/featflower/raw_logs/particle_force_e3_l2.log` | 383446 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l2/particle_force.log` | Raw per-step particle log of certified run e3_l2 |
| `reference/featflower/raw_logs/particle_force_e3_l3.log` | 383446 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l3/particle_force.log` | Raw per-step particle log of certified run e3_l3 |
| `reference/featflower/raw_logs/particle_force_e3_l4.log` | 383446 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l4/particle_force.log` | Raw per-step particle log of certified run e3_l4 |
| `reference/featflower/raw_logs/particle_force_e4_l2.log` | 276946 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l2/particle_force.log` | Raw per-step particle log of certified run e4_l2 |
| `reference/featflower/raw_logs/particle_force_e4_l3.log` | 276946 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3/particle_force.log` | Raw per-step particle log of certified run e4_l3 |
| `reference/featflower/raw_logs/particle_force_e4_l3_dt0p25_sync.log` | 1107646 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p25_sync/particle_force.log` | Raw per-step particle log of certified run e4_l3_dt0p25_sync |
| `reference/featflower/raw_logs/particle_force_e4_l3_dt0p5_sync.log` | 553846 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l3_dt0p5_sync/particle_force.log` | Raw per-step particle log of certified run e4_l3_dt0p5_sync |
| `reference/featflower/raw_logs/particle_force_e4_l4.log` | 276946 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4/particle_force.log` | Raw per-step particle log of certified run e4_l4 |
| `reference/featflower/raw_logs/particle_force_e4_l4_dt0p5_sync.log` | 553846 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e4_l4_dt0p5_sync/particle_force.log` | Raw per-step particle log of certified run e4_l4_dt0p5_sync |
| `reference/figures/d21_brenner_crossover.png` | 155206 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/docs/md_docs/dns_figures/d21_brenner_crossover.png` | Resolved FBM wall-approach force vs Brenner: the gap <~ 2h crossover |
| `reference/figures/d22_g3_tencate.png` | 115328 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/docs/md_docs/dns_figures/d22_g3_tencate.png` | E1/E2 bottom approach, FBM with/without lubrication add-on vs PIV (Fig. 13 analogue) |
| `reference/figures/e1_l3_vs_piv.png` | 62641 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3/e1_l3_vs_piv.png` | E1 L3 velocity and trajectory vs digitised PIV (compare_tencate.py layout) |
| `reference/reference_values.csv` | 2852 | `written from paper Table I/II + derived` | Reference numbers with per-row source/formula |
| `reference/site_lubrication_export/approach_E1_base.csv` | 49449 | `/data/warehouse17/rmuenste/code/FF-EL/ff-redesign/scripts/source-data/sedimentation/lubrication/approach_E1_base.csv` | Benchmark-website source data (read-only copy): approach-window curves / bands / case table |
| `reference/site_lubrication_export/approach_E1_lub.csv` | 68978 | `/data/warehouse17/rmuenste/code/FF-EL/ff-redesign/scripts/source-data/sedimentation/lubrication/approach_E1_lub.csv` | Benchmark-website source data (read-only copy): approach-window curves / bands / case table |
| `reference/site_lubrication_export/approach_E2_base.csv` | 30241 | `/data/warehouse17/rmuenste/code/FF-EL/ff-redesign/scripts/source-data/sedimentation/lubrication/approach_E2_base.csv` | Benchmark-website source data (read-only copy): approach-window curves / bands / case table |
| `reference/site_lubrication_export/approach_E2_lub.csv` | 42320 | `/data/warehouse17/rmuenste/code/FF-EL/ff-redesign/scripts/source-data/sedimentation/lubrication/approach_E2_lub.csv` | Benchmark-website source data (read-only copy): approach-window curves / bands / case table |
| `reference/site_lubrication_export/brenner_bands.csv` | 225 | `/data/warehouse17/rmuenste/code/FF-EL/ff-redesign/scripts/source-data/sedimentation/lubrication/brenner_bands.csv` | Benchmark-website source data (read-only copy): approach-window curves / bands / case table |
| `reference/site_lubrication_export/cases.csv` | 127 | `/data/warehouse17/rmuenste/code/FF-EL/ff-redesign/scripts/source-data/sedimentation/lubrication/cases.csv` | Benchmark-website source data (read-only copy): approach-window curves / bands / case table |
| `tools/compare_tencate.py` | 4549 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/compare_tencate.py` | Comparison protocol: peak, timing, RMS vs PIV, h/d, touchdown |
| `tools/d22_g3_tencate_analysis.py` | 3769 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/d22_g3_tencate_analysis.py` | Lubrication-ON vs OFF approach analysis (Fig. 13 analogue) |
| `tools/stage_tencate_case.py` | 8000 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/stage_tencate_case.py` | How the campaign staged E1-E4 rundirs (fluid properties on both sides, dt sync, step counts) |
| `tools/tencate_error_decomposition.py` | 4274 | `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/tencate_error_decomposition.py` | Additive S(h)+T(dt) error fit over the ladder points |

Total: 165 files, 25.6 MB.
