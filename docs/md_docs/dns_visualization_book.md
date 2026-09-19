# DNS validation campaign — Visualization book, v1.0

Status: **v1.0**, 2026-09-19. Inventory of the six rendered cases (ten Cate
sedimentation, DKT, hindered settling, numerical viscometer, Oberbeck
spheroid, Jeffery orbit) as they stand after the gallery went live on
`ff-redesign` master (PR #4, 2026-09-16; PR #5 and commit `8cfbe0a`,
2026-09-17). Scripts collected into `tools/viz/` the same day.

## 0. Purpose and how to use

This book is the ledger for visualizations as `dns_validation_datasheet.md`
is for numbers: every asset of the campaign — ParaView still, hand-off
texture, wall-note plot, Blender render, website ladder — has an entry that
says what it shows, which run and frame it was made from, how it was made,
at what quality level it stands and what the next step is. Entries are
dated and are never silently rewritten: a superseding entry names the one
it replaces (as the datasheet's erratum rows do), so the history of a
picture stays readable. The datasheet carries no `viz`/`gallery` rows
(checked 2026-09-19); the per-case Log sections here play that role. Add a
case by copying §4's section skeleton; move a case up the ladder (§2) by
adding a dated Log line and updating "Current level" — the evidence named
there must exist on disk or in a commit.

## 1. Pipeline and environment

**Headless rendering.** Only the Kitware binary works on the warehouse login
node: `/sfw/paraview/6.0.1/bin/pvpython --mesa --force-offscreen-rendering
<script> <args>` (bundled llvmpipe; ~12 s start-up). The `paraview/5.13.0`
and `5.11.2` osmesa modules do not work (they link a `libOSMesa.so.8` that
is absent). Plot scripts are plain `python3` + matplotlib (`Agg`).

**Review loop.** Render a PNG into the scratchpad → view it (Read) → deliver
without file copies: `SendUserFile` for a quick look, or a data-URI figure in
a private review artifact (16 MB page cap; the six-case board is artifact
`f910a9d8` "Sedimentation Stills", 21 figures, template kept in the
scratchpad as `viz/sedimentation_stills.template.html`). The owner works
over ssh and reviews stills in the browser, not on the machine.

**Blender hand-off.** One folder per case, `blender_viz/<case>/`, at the repo
root; listed in `.git/info/exclude`, never committed (also excluded:
`d62_v1b_frames/`, `d61_frames/`, the copied Fritz frames). The owner copies
the folder (`blender_viz/website_assets.zip`, 24 MB, 2026-09-16 22:31 is the
transport bundle) and Blender agents build the scene there. Contract:

- `README.md` — scene brief: what is shown, coordinate frame table, file
  table, texture mapping, colour scale, composition suggestions, provenance
  (run, frame, fluid, mesh, render date). `scene_spec.json` — the same,
  machine-readable, plus a dated `decision` block for owner choices.
- Units: SI (metres) where the run is SI (ten Cate); otherwise the run's
  box units with a **recommended scene scale** stated once (1 d = 1 cm for
  DKT, hindered settling, viscometer, Jeffery; 1 cell = 10 cm for
  Oberbeck). z is up in every case; the camera side is given so that +x is
  to the right.
- Plane textures: `*_color.png` (colour field only), `*_overlay.png`
  (isolines + arrows on transparent), `*_combined.png`; the pixel frame
  equals the data plane exactly, no margins, so the UV rule is a plain
  rectangle (u = 0 at the left coordinate, v = 0 at the bottom).
- Alpha rule: alpha = 0 outside the domain and inside every body footprint,
  cut **8 % inside** the true section so the Blender body always covers the
  hole (viscometer: holes from the logged centres; ten Cate: the sphere
  footprint; spheroids: the analytic ellipse). The body must be present in
  the scene or the hole shows.
- `colorbar*.png`, 500 × 1600 px, transparent — for compositing, not the 3-D
  scene. Optional geometry: `isolines*.obj` (polylines in scene units) and
  `arrows*.ply` (triangle meshes) for relief and shadows.
- Wall-note plots (the "note on the lab wall" a scene can show on a
  monitor): dark style, 300 dpi, opaque, see §3.
- Per-body data as CSV with a comment header stating run, time and units.

**What is tracked.** Committed: analysis figures that substantiate datasheet
rows and their generator tools (`docs/md_docs/dns_figures/*.png`,
`tools/*.py`), the Q0 setup-scene tool `tools/render_setup_scene.py`, and now
the hand-off producers in `tools/viz/`. Untracked: every render, every
hand-off texture, `.pvsm` state files (the three setup-scene `.pvsm` files in
`dns_figures/` are the exception the owner accepted with
`setup_scenes_units.md`), the owner's Blender scenes and PNG masters.

**Website rule.** `ff-redesign` ships WebP width ladders only under
`public/benchmark-assets/<id>/media/gallery/`; PNG masters stay on the
owner's machine (`blender_viz/website_assets/<case>/`). The registry
`src/data/gallery.ts` holds id, title, caption, provenance line, family,
aspect, focal point, ladder and a `kind: "image" | "video"` field ready for
animations. Captions are plain language with no "how this image was made"
text; the site carries clean results only.

## 2. Quality ladder

| level | name | "done" means | evidence recorded in the Log |
|---|---|---|---|
| Q0 | ParaView still for review | one headless still of the chosen frame, viewed and accepted by the owner as the right frame/plane/quantity | script + args, PNG name and size, review artifact or SendUserFile date |
| Q1 | hand-off folder complete | `blender_viz/<case>/` holds textures (color/overlay/combined) or per-body CSVs, colourbar, geometry where useful, wall-note plot(s), README and `scene_spec.json`; pixel frames verified against the README | folder date, file count, texture pixel sizes, README provenance block |
| Q2 | Blender still rendered | the owner's PNG master exists at full resolution and was reviewed | master file name, resolution, date |
| Q3 | hero render on the website | WebP ladder derived from the master, registry entry with caption and provenance, cross-link figure on the benchmark page, gallery order fixed | gallery id, ladder rungs, `ff-redesign` commit/PR |
| Q4 | refined hero | at least one recorded revision of the published hero (composition, materials, colour treatment, labels, resolution) with the old render named | new master + ladder, commit, what changed and why |
| Q5 | animation | a `kind: "video"` gallery item with poster; needs a dense-output re-run (frame cadence the production runs never wrote) | re-run rundir, frame cadence, encoder settings, poster |

Lower levels stay valid when a higher one is reached; a case can sit at Q3
for one asset and Q1 for a companion (e.g. the viscometer φ = 0.05 cup).

## 3. Conventions

**Colour maps.** Data planes: **Viridis (matplotlib), filled bands** — the
owner's decision for ten Cate and Jeffery (2026-09-16/17) and the default
of the later exports (viscometer, Oberbeck). Band counts: ten Cate 25 bands
of width 0.04 over 0–1 (the paper's own filled-contour step); Jeffery 20
bands over 0–0.4 (edges every 0.02); viscometer 20 bands over 0–1 (isolines
every 0.1); Oberbeck 20 bands over 0–1 (isolines 0.2, 0.4, 0.6, 0.8, 0.9,
0.95). Isolines sit on band edges unless listed. The first treatment, "TU
Petrol Sequential" (`tu-colormaps.json`, stops 0.0 #c7dae2, 0.2 #a0bfcd,
0.4 #78a4b7, 0.6 #5a91a7, 0.8 #3d7e98, 1.0 #326a82, continuous), stays in
the ten Cate and Jeffery folders as reference; Turbo and Inferno banded
variants exist for ten Cate only. Per-body colour (hindered settling):
Viridis over the settling ratio 0–1.2. Fixed body colours (DKT): leader cyan
`#4fd8ff`, trailer amber `#ffb84a`.

**Normalisations.** ten Cate |u|/u_∞, u_∞ = 0.03845 m/s (Re ν/d_p); Jeffery
|u′|/(γ̇ a) with γ̇ a = 0.1 and u′ = u − γ̇ z x̂; viscometer |u|/U_wall,
U_wall = 0.5; Oberbeck |u|/u_max of each run (u_max 2.5952e-3 parallel,
2.1109e-3 perpendicular); hindered settling −v_z/u_t, u_t = 0.4061.

**Line colour on Viridis.** Black isolines/arrows vanish on the dark-violet
far field; the Jeffery `_lightlines` set uses light grey (0.93). Direction
arrows are unit-length (direction only), masked in the body and where the
normalised speed < 0.01.

**Coordinate frames and scales.** z up everywhere; camera on the −y side
sees +x to the right (viscometer mid plane: seen from +z, +y up). Scene
scales as in §1. Periodic cells (Oberbeck) may tile the texture along x.

**Wall-note plot style** (all six cases): ground `#2a2d31`, text `#f2f2f2`,
grid `#5a5f66`; series cyan `#4fd8ff`, amber `#ffb84a`, lime `#8dff6a`, pink
`#ff7ab6`; reference curves white/dashed; 300 dpi, opaque; font size 12–13,
titles 13–14. The dark versions are the chosen ones (owner, 2026-09-16); the
earlier white/transparent versions stay as reference only.

**Captions.** On the website: one sentence of what the reader sees, one
provenance line (regime, moment, run family), no method text, no defect
history. In this book: full provenance including script and arguments.

## 4. Cases

### 4.1 ten Cate sedimentation (E1) — `blender_viz/ten_cate_sedimentation/`

**Data.** Run `q2p1_dns_rundir_e1_l3_g3def` (Re 1.5, level 3, D/h ≈ 12 on
the output mesh, lubrication variant inactive at this gap; datasheet rows
`d22_g3_tencate`, `e1_l3`), frame `_vtk/main.03399.pvtu`, t = 3.400 s,
sphere centre z = 0.015398 m from `particle_force.log` (h/d_p = 0.53, the
paper's Fig. 9 moment). E2 wall note from `q2p1_dns_rundir_e2_l3_g3def`'s
log (no E2 frame used).

**What is shown.** Plane y = 0 through the sphere over the full column
(x ∈ [−0.05, 0.05], z ∈ [0, 0.16] m; quarter domain mirrored across x = 0
with u_x flipped): |u|/u_∞ colour (peak 0.69 in frame), isolines, direction
arrows on a 34 × 54 grid. Wall notes: settling velocity and gap height vs
time, DNS vs PIV, for E1 and the E2 sibling.

**Assets.**

| asset | file(s) | size / format | date |
|---|---|---|---|
| Q0 still | `reference_still_fig9_style.png` | 2000 × 1750 RGB | 2026-09-16 |
| plane textures, 4 treatments | `slice_y0_{color,overlay,combined}[_viridis_banded|_turbo_banded|_inferno_banded].png` | 2000 × 3200 RGBA, 50 µm/px | 2026-09-16 |
| colourbars | `colorbar[_variant].png` | 500 × 1600 RGBA | 2026-09-16 |
| geometry | `isolines[_variant].obj` (24 lines banded, 12 petrol), `arrows.ply` | OBJ / PLY, metres | 2026-09-16 |
| wall notes (chosen) | `e1_velocity_trajectory_plot_dark_300dpi.png`, `e2_…dark_300dpi.png` | 3300 × 1200 RGBA | 2026-09-16 |
| wall notes (reference) | `e{1,2}_…{white_150dpi,transparent_300dpi}.png` | 1540 × 560 / 3300 × 1200 | 2026-09-16 |
| trajectory | `sphere_trajectory.csv` (4300 rows) | CSV, SI | 2026-09-16 |
| Blender master (owner) | `website_assets/ten_cate_sedimentation/gallery/ten_cate_gallery_master.png` | 2400 × 1800 RGBA, 5.1 MB | 2026-09-16 22:23 |
| site ladder | `sedimentation/media/gallery/ten_cate_gallery_1200.webp` (1200 w), `ten_cate_gallery.webp` (2400 w) | WebP | 2026-09-16 |

Website: gallery id `sedimentation`, family "single body", aspect 4:3, focal
50 % 60 %. Caption: "A single sphere half a diameter above the tank floor:
the squeeze flow spreads along the bottom while the wake still follows the
sphere down." Provenance line: "Re = 1.5, t = 3.40 s, the moment of ten
Cate's flow-field comparison".

**Pipeline.** `tools/viz/tencate_fig9.py RUNDIR 03399 0.03845 0.015398
OUT.png "label"` (Q0); `tools/viz/export_blender_slice.py RUNDIR 03399
0.03845 0.015398 OUTDIR ["Viridis (matplotlib)" 25 _viridis_banded]` (also
run with Turbo/Inferno and with the petrol default); wall notes
`tools/viz/tencate_plot_dark.py RUN_LOG E1 TCREF OUT.png` around
`tools/compare_tencate.py`. Frame and log present on disk 2026-09-19.

**Current level.** **Q3** since 2026-09-16 (commit `bb259f4`/`527b518`, PR
#4 merged `9edeb00`). Evidence: master 2400 × 1800, ladder 1200/2400.

**Known limitations.** Output mesh 12 cells across the diameter — the
isolines carry that texture. The top of every colour map is unused (peak
0.69), deliberately. The E2 frame `main.03399` does not exist on disk
(`e2_l3_g3def/_vtk`), so an E2 plane would need the right E2 frame index.
`scene_spec.json` still lists the petrol map as `color_scale.colormap` with
the Viridis choice only under `recommended_variant`.

**Next quality step.** E2 (Re 4.1) companion plane at its own h/d_p ≈ 0.5
frame so the second "screen" has a matching column; then a Q4 pass on the
hero (the plane's far-field violet against the lab lighting was the open
point in the README). Animation of the whole settling (Q5) needs the E1
frames at a regular cadence — the frame cadence the g3def run actually
wrote to `_vtk` was not checked; check it before planning.

**Log.**
- 2026-09-16 — Q0 Fig. 9-style still and petrol textures; E1/E2 plots in
  white and transparent.
- 2026-09-16 — Turbo/Inferno/Viridis banded variants added; owner chooses
  `_viridis_banded` and the dark plot style.
- 2026-09-16 — owner's Blender master 2400 × 1800; registered as
  `sedimentation` in PR #4.

### 4.2 Drafting, kissing, tumbling — `blender_viz/dkt_tumbling/`

**Data.** Run `q2p1_dns_rundir_dkt_nofric_long` (job 137965; box 6 × 6 × 24
in d, level 3, D/h = 8, dt = 0.005, t = 0–40, frictionless contact),
`particle_force.log` only; datasheet rows `dkt_nofric_long`, `dkt_nofric`,
`dkt_tkiss_correction`. The README states no field output beyond t = 0
(`_vtk` holds 117 entries; not re-verified which times).

**What is shown.** Geometry-only stroboscopic still: eight labelled stages
(Release t = 0.005, Drafting 10, Kissing 18.04, Rolling contact 25, Tumbling
28, Separation 29.64, Side by side 32.92, Role exchange 40) plus ghost pairs
every 2 t.u., leader cyan / trailer amber; no fluid plane (owner decision).
Wall note: heights, centre distance and tilt vs time with the stages marked.

**Assets.**

| asset | file(s) | size / format | date |
|---|---|---|---|
| stages | `stages.csv`, `stages.json` (with label text) | 8 rows | 2026-09-16 |
| ghosts | `ghosts_every_2tu.csv` (21 pairs) | CSV | 2026-09-16 |
| trajectory | `trajectory_full.csv` (8000 rows, 0.005 t.u.) | CSV | 2026-09-16 |
| wall note | `dkt_phases_plot_dark_300dpi.png` | 3300 × 2550 RGBA | 2026-09-16 |
| Q0 preview | `strobe_preview_dark.png` | 1300 × 2600 RGBA | 2026-09-16 |
| Blender master (owner) | `website_assets/dkt_tumbling/renders/dkt_strobe_1500x2000.png` | 1500 × 2000 RGBA, 3.3 MB | 2026-09-16 22:23 |
| site ladder | `dkt/media/gallery/dkt_strobe_{600,1000,1500}w.webp` | WebP | 2026-09-16 |

Website: id `dkt`, family "pairs", aspect 3:4, focal 50 % 45 %. Caption:
"Two spheres released one above the other, shown at eight moments: the
trailer drafts and catches up, the pair touches, tumbles over and separates
with the roles exchanged." Provenance: "Frictionless contact, D/h = 8,
t = 0 to 40".

**Pipeline.** `tools/viz/dkt_export.py RUNDIR OUTDIR` (python3; produces
everything above). Campaign tools `tools/dkt_plot_trajectory.py`,
`tools/dkt_export_series.py` are the website-figure producers (copies in the
folder are identical to the tracked ones).

**Current level.** **Q3** since 2026-09-16 (PR #4). Evidence: master
1500 × 2000, ladder 600/1000/1500.

**Known limitations.** No fluid data at any stage; the tilt arcs and labels
are Blender-side design. README table quotes the Release stage at t = 0
while `stages.csv` carries the first log record t = 0.005.

**Next quality step.** Add a fluid vorticity/velocity plane at the kissing
and tumbling stages once a dense-output re-run exists (the geometry here
would not change). Q5 animation candidate (prerequisite: the same re-run at
a frame cadence of ~0.25 t.u.).

**Log.**
- 2026-09-16 — export, dark plot, preview; owner decision: geometry only.
- 2026-09-16 — owner's Blender master; registered as `dkt` in PR #4.

### 4.3 Hindered settling — `blender_viz/hindered_settling/`

**Data.** Cloud run `q2p1_dns_rundir_d32_n120_s1` (N = 120, seed 1, walled
6 d column, D/h = 8, ρ_p/ρ_f = 1.14, frictionless), reference
`q2p1_dns_rundir_d32_ut` (u_t = 0.4061); snapshot t = 20 in the settled
window 15–25; datasheet rows `d32_phi_ladder`, `d32_ut_ref`. Logs only.

**What is shown.** 120 spheres coloured by −v_z/u_t (Viridis 0–1.2; mean
0.597), 3 t.u. trails (13 instants, 0.25 apart), and the lone reference
sphere at (0, 0, 13.97) in a second identical column. Wall note: U/u_t vs
cloud φ for N = 20…120 × 3 seeds, fit n = 4.58, unbounded RZ band 2.7–3.0.

**Assets.**

| asset | file(s) | size / format | date |
|---|---|---|---|
| cloud | `cloud_t20.csv` (120 rows: id, x, y, z, v, ratio) | CSV | 2026-09-16 |
| trails | `trails_t20.csv` (1560 rows) | CSV | 2026-09-16 |
| snapshot | `snapshot.json` | JSON | 2026-09-16 |
| wall note | `hindered_rz_plot_dark_300dpi.png` | 2700 × 1950 RGBA | 2026-09-16 |
| Q0 preview | `cloud_preview_dark.png` | 1400 × 2400 RGBA | 2026-09-16 |
| Blender master (owner) | `website_assets/hindered_settling/renders/hindered_gallery_1600x2000.png` | 1600 × 2000 RGBA, 3.7 MB | 2026-09-16 22:23 |
| site ladder | `hindered-settling/media/gallery/hindered_gallery_{600,1000,1500}w.webp` | WebP | 2026-09-16 |

Website: id `hindered-settling`, family "collective", aspect 4:5, focal
50 % 50 %. Caption: "A cloud of 120 spheres settles at 60 % of the speed of
the lone sphere in the neighbouring column; colour is each sphere's own
settling speed." Provenance: "N = 120, walled 6 d column, t = 20".

**Pipeline.** `tools/viz/hindered_export.py CLOUD_RUNDIR UT_RUNDIR OUTDIR
20` (python3). Ladder numbers from `tools/d32_ladder_analysis.py --window
15 25 --ut 0.4061` (copy in folder identical to tracked).

**Current level.** **Q3** since 2026-09-16 (PR #4). Evidence: master
1600 × 2000, ladder 600/1000/1500.

**Known limitations.** No fluid field (only t = 0 frames exist). The
colour legend for the settling ratio is not supplied as a colourbar file
(the README states the range in words).

**Next quality step.** Add a settling-ratio colourbar PNG to the hand-off
(same 500 × 1600 format) so the render can carry a legend; consider a second
snapshot (t = 15 or 25) to show the cloud's spreading; the wide-column
family `d32w_*` (row `d32_wide_attribution`) has no picture yet and would
make the confinement argument visible side by side.

**Log.**
- 2026-09-16 — export, dark RZ plot, preview; owner decision: no fluid plane.
- 2026-09-16 — owner's Blender master; registered as `hindered-settling` in
  PR #4.

### 4.4 Numerical viscometer — `blender_viz/numerical_viscometer/`

**Data.** Hero `q2p1_dns_rundir_d52_v23` (φ = 0.20, 900 spheres, level 3,
D/h = 9.15, plateau t ≥ 235, η_r = 1.7143; row `d52_v23_phi20`), frame
`_vtk/main.14000.pvtu` (t = 250, 1.2 M cells); companion
`q2p1_dns_rundir_d52_v21` (φ = 0.05, 225 spheres, η_r = 1.1062; row
`d52_v21_einstein`), frame `main.13000.pvtu`. Sphere centres from the last
`particle_force.log` record. Plot rows also `d52_v22_phi10`, `d52_v20_baseline`.

**What is shown.** Annular Couette cell (r_i = 5 rotating at Ω = 0.1,
r_a = 10 fixed, height 10): mid plane z = 5 and vertical plane y = 0 with
|u|/U_wall in 20 Viridis bands, isolines every 0.1, exact sphere holes; all
spheres placed from the log. Wall note: η_r vs φ, DNS vs Einstein,
Batchelor, Krieger–Dougherty composites (+0.60 / −0.78 / −0.94 %).

**Assets.**

| asset | file(s) | size / format | date |
|---|---|---|---|
| mid-plane textures | `phi020_mid_{color,overlay,combined}.png`, `phi005_mid_*.png` | 4000 × 4000 RGBA, 200 px/unit | 2026-09-16 |
| vertical textures | `phi020_vert_*.png`, `phi005_vert_*.png` | 4000 × 2000 RGBA | 2026-09-16 |
| colourbar | `colorbar_visco.png` | 500 × 1600 RGBA | 2026-09-16 |
| spheres | `phi020_spheres_t250.csv` (900), `phi005_spheres_t250.csv` (225) | CSV, box units | 2026-09-16 |
| side files | `_centres_phi020.npy`, `_centres_phi005.npy` | npy (script intermediates) | 2026-09-16 |
| wall note | `viscometer_eta_plot_dark_300dpi.png` | 2700 × 1950 RGBA | 2026-09-16 |
| Blender masters (owner) | `renders/viscometer_lab_2000x1500.png` (v1, all spheres opaque), `renders/viscometer_lab_ghost_2000x1500.png` (v2, upper half faint) | 2000 × 1500 RGBA, 13 MB each | 2026-09-16 23:02 / 2026-09-17 10:38 |
| owner's test renders | `renders/test_v{1,2,3}.png`, `viscometer_lab_test_v3.png` | 600 × 450 / 1000 × 750 | 2026-09-16 23:01 |
| owner's ladders | `viscometer_lab_{480…2000}w.webp` (8 rungs, v1); `viscometer_lab_ghost_{600,1000,1500,2000}w.webp` (v2) | WebP | 2026-09-16 / 17 |
| site ladder | `numerical-viscometer/media/gallery/viscometer_lab_{600,1000,1500,2000}w.webp` — **these are the ghost renders** (byte sizes match `_ghost_*`), shipped under the v1 file names | WebP | 2026-09-17 |

Website: id `numerical-viscometer`, family "collective", aspect 4:3, focal
50 % 55 %. Caption: "A suspension of 900 neutrally buoyant spheres between a
rotating bob and a fixed cup; the data plane shows how the particles bend
the Couette profile." Provenance: "phi = 0.20 at the torque plateau,
t = 250".

**Pipeline.** Write the sphere CSV first (the export reads
`OUTDIR/<TAG>_spheres_t250.csv`; the CSV producer was an inline step, not a
tracked script — see §6), then
`tools/viz/export_blender_viscometer.py PVTU OUTDIR phi020` (and `phi005`);
holes via a ProgrammableFilter on min_i |x − x_i|²/R². Wall note:
`tools/viz/wallnote_plots.py OUT_VISCO.png OUT_OBERBECK.png`. Probe:
`tools/viz/probes/probe_visco.py PVTU`.

**Current level.** **Q4** since 2026-09-17: hero revised once after
publication — the first render (all spheres opaque, PR #4, 2026-09-16) was
replaced by the ghost render (upper half faint so the mid plane stays
visible; PR #5, commit `8e62adb`, merged `5916579`). Evidence: both masters
on the owner's machine, site ladder dated 2026-09-17 13:06.

**Known limitations.** The render shows the lower half of the cup with the
spheres above the plane as ghosts (owner decision) — the vertical plane is
not used on the site. The φ = 0.05 companion sits at Q1 (textures exist, no
render). The export script writes its `.npy` side file to the session
scratch path (hard-coded).

**Next quality step.** Exact hole masks are already in place; a top-half
fade of the plane or a second cup with the φ = 0.05 companion are the
composition options left open in the README. Q5 candidate (the rotating
suspension) needs a dense-output re-run of v23 over one bob revolution
(2π/Ω ≈ 63 t.u.). The D/h = 16 rung (`d52_v25f`, closed 2026-09-19, row
`d52_v25f_l4_verdict`, φ = 0.05 only) could supply a higher-resolution
companion frame; the φ = 0.20 hero has no level-4 counterpart.

**Log.**
- 2026-09-16 — textures, CSVs, colourbar, dark η_r plot; owner decision:
  hero φ = 0.20, companion φ = 0.05.
- 2026-09-16 — owner's first Blender master (`viscometer_lab_2000x1500.png`)
  registered in PR #4.
- 2026-09-17 — ghost-sphere revision (`viscometer_lab_ghost_2000x1500.png`)
  replaces it on the site (PR #5 `8e62adb`); file names on the site unchanged.

### 4.5 Oberbeck spheroid drag — `blender_viz/oberbeck_spheroid/`

**Data.** Runs `q2p1_dns_rundir_d61_v4a` (parallel) and `_v4b`
(perpendicular), level-4 rung (2b/h = 19), Fritz jobs 4171401/02; frames
`d61_frames/v4{a,b}/_vtk/main.00350.pvtu` (t = 3.5, forces on the plateau
slope, written one level below the solve, 46 656 cells; copied back from
Fritz). Plot rows `d61_v123_l3`, `d61_v4_resolution`, `d61_v5_halfsize`,
`d61_review_corrections`. Level-3 frames of V1, V2 and the 45° V3 are also
under `d61_frames/` (5832 cells).

**What is shown.** Periodic unit cell, plane y = 0.5 through the body for
both orientations: |u|/u_max in 20 Viridis bands, isolines, white
streamlines seeded across the bottom edge, direction arrows. Wall note:
anisotropy ratio over the resolution/size ladder vs Oberbeck 1.14532 with
±2 % band, and absolute drags within ±3 %.

**Assets.**

| asset | file(s) | size / format | date |
|---|---|---|---|
| plane textures | `{parallel,perpendicular}_plane_{color,overlay,combined,streamlines}.png` | 4000 × 4000 RGBA | 2026-09-16 |
| colourbar | `colorbar_oberbeck.png` | 500 × 1600 RGBA | 2026-09-16 |
| wall note | `oberbeck_drag_plot_dark_300dpi.png` | 3600 × 1680 RGBA | 2026-09-16 |
| Blender master (owner) | `website_assets/oberbeck_spheroid/renders/oberbeck_gallery_2000x1200.png` | 2000 × 1200 RGBA, 2.7 MB | 2026-09-16 22:23 |
| site ladder | `oberbeck-spheroid-drag/media/gallery/oberbeck_gallery_{600,1000,1500,2000}w.webp` | WebP | 2026-09-16 |

Website: id `oberbeck-spheroid-drag`, family "non-spherical", aspect 5:3,
focal 45 % 50 %. Caption: "The same spheroid held along and across a slow
flow: the drag across the axis is 1.145 times the drag along it."
Provenance: "Stokes flow, aspect ratio 2, level-4 runs".

**Pipeline.** `tools/viz/export_blender_oberbeck.py PVTU 0 0 1 OUTDIR
parallel` and `… 1 0 0 OUTDIR perpendicular`; wall note from
`tools/viz/wallnote_plots.py` (numbers hard-coded from the rows above);
analysis `tools/d61_oberbeck_analysis.py`. Probe `tools/viz/probes/probe_d61.py`.

**Current level.** **Q3** since 2026-09-16 (PR #4). Evidence: master
2000 × 1200, ladder 600/1000/1500/2000.

**Known limitations.** Each orientation is normalised by its own u_max, so
the two planes are not on a common colour axis (stated in the README). The
periodic array is translated into "channel with a sting" in the scene — an
honest but not literal rendering of the setup.

**Next quality step.** Streamline density and seeding (one seed line at the
bottom edge; the wake side is sparse) — re-export with a second seed line
or a stream-tracer on a grid; add the V3 45° orientation as a third body
from the existing level-3 frame (`d61_frames/v3`) to show the lateral
force (row `d61_v3b_transverse`).

**Log.**
- 2026-09-16 — frames copied from Fritz, textures with streamlines, dark
  ladder plot.
- 2026-09-16 — owner's Blender master; registered as
  `oberbeck-spheroid-drag` in PR #4.

### 4.6 Jeffery orbit (D6.2 V1b) — `blender_viz/jeffery_orbit/`

**Data.** Fritz rundir `q2p1_dns_rundir_d62_v1b` (H = 8 box, level 4, 2b/h
= 10.4 per `d62_resolution_pinned`, ρ_r = 10, Re_a = 0.05, three chained
segments); frames `d62_v1b_frames/_vtk/main.09000.pvtu` (t = 90, chosen) and
`main.10000.pvtu` (t = 100, reference), written one level below the solve
(430 080 hexahedra, 20 cells across the body length; no indicator field —
the body is masked analytically from the logged axis). Axis trace from the
`DNS_PART_AXIS` records of the three segment logs (seams 40.01, 80.01;
`d62_v1b_frames/v1b_axis.log`). Datasheet row `d62_v1b_orbit`; committed
analysis figure `docs/md_docs/dns_figures/d62_jeffery_v1b.png`.

**What is shown.** Planar Couette cell (x ∈ [−4, 4], y ∈ [−3, 3], z ∈ [−4,
4], plates ±0.8 x̂, γ̇ = 0.2), spheroid a = 0.5, b = c = 0.25 at the origin,
axis (0.8537, 0, −0.5208) at t = 90 (Blender Euler Y = 31.39°). Plane y = 0
with the disturbance field |u′|/(γ̇ a) (peak 0.47), full and zoom windows;
20 Viridis bands over 0–0.4 with light-grey lines in the chosen set. Wall
note: axis components, rotation rate vs time and vs orientation, DNS vs
Jeffery (T γ̇ = 15.708 predicted, 15.755 measured, +0.30 %).

**Assets.**

| asset | file(s) | size / format | date |
|---|---|---|---|
| petrol textures | `t{090,100}_{full,zoom}_dist_{color,overlay,combined}.png`, `t{090,100}_full_total_color.png` | 4000 × 4000 RGBA (full 500 px/unit, zoom 1000 px/unit) | 2026-09-16 |
| Viridis banded, black lines | `t090_{full,zoom}_dist_*_viridis_banded.png`, `t090_full_total_color_viridis_banded.png` | 4000 × 4000 RGBA | 2026-09-17 11:05 |
| Viridis banded, light lines | `t090_*_viridis_banded_lightlines.png` | 4000 × 4000 RGBA | 2026-09-17 11:06 |
| colourbars | `colorbar_dist.png`, `colorbar_dist_viridis_banded[_lightlines].png` | 500 × 1600 RGBA | 2026-09-16 / 17 |
| geometry | `t{090,100}_zoom_isolines.obj`, `t{090,100}_zoom_arrows.ply` and the `_viridis_banded[_lightlines]` twins (arrows identical) | OBJ / PLY, box units | 2026-09-16 / 17 |
| wall notes | `jeffery_plot_dark_300dpi.png` (chosen style); `jeffery_plot_transparent_300dpi.png`, `jeffery_plot_white_150dpi.png` (reference) | 2550 × 3150 / 1275 × 1575 RGBA | 2026-09-16 |
| traces | `axis_trace.csv` (12 001 rows, 0.01 t.u.), `axis_strobe_t080_t100.csv` (11 orientations) | CSV | 2026-09-16 |
| Blender masters (owner) | `website_assets/jeffery_orbit/gallery/jeffery_gallery_master.png` (petrol plane), `jeffery_gallery_viridis_master.png` (Viridis plane) | 2400 × 1800 RGBA, 4.2 / 4.3 MB | 2026-09-17 11:20 |
| owner's ladders | `jeffery_gallery[_1200].webp` (petrol), `jeffery_gallery_viridis[_1200].webp` | WebP 1200 / 2400 w | 2026-09-17 |
| site ladder | `jeffery-orbit/media/gallery/jeffery_gallery_1200.webp`, `jeffery_gallery.webp` — **the Viridis renders** (byte sizes match `_viridis`), shipped under the petrol file names | WebP | 2026-09-17 13:06 |

Website: id `jeffery-orbit`, family "non-spherical", aspect 4:3, focal
40 % 50 %. Caption: "A spheroid tumbling in simple shear between two moving
plates; the plane shows the four-lobed flow the body itself induces."
Provenance: "t = 90, aspect ratio 2, wall clearance 8 semi-axes".

**Pipeline.** `tools/viz/concat_axis.py > v1b_axis.log` (Fritz paths);
`tools/viz/export_blender_jeffery.py PVTU 0.8537 0 -0.5208 OUTDIR t090`
(petrol), `… t090 "Viridis (matplotlib)" 20 _viridis_banded` and
`… _viridis_banded_lightlines 0.93,0.93,0.93`; t100 with axis (−0.2153, 0,
−0.9766). Wall note `tools/viz/jeffery_plot_dark.py AXIS_LOG OUT.png` around
`tools/d62_jeffery_analysis.py`.

**Current level.** **Q4** since 2026-09-17: the published hero (petrol
plane, PR #4) was replaced by the Viridis-plane render on master
(`8cfbe0a`, direct to master). Evidence: both masters on the owner's
machine, site ladder dated 2026-09-17 13:06. Whether the Viridis master
used the black-line or the `_lightlines` set is not recorded in the site
commit or the folder (not verified).

**Known limitations.** Two frames only; the strobe (scene idea 1) and the
"instrument" plates are Blender-side. The body has no indicator field in
the frames, so the hole is analytic (8 % inside the ellipse).
`scene_spec.json` still names the petrol map as `color_scale.colormap`,
the Viridis choice under `color_variants`, and `plot_assets.opaque` points
to the white plot although the dark one is the chosen style. The gallery
caption's "wall clearance 8 semi-axes" is the H = 8 box (plate at 4 = 8 a
from the centre) — correct, but phrased differently from the README's
"H = 8".

**Next quality step.** Stroboscopic ghosts from `axis_strobe_t080_t100.csv`
in the hero (the uneven spacing is the physics); the H = 4 clearance rung
(`d62_v2_clearance`, figure `d62_jeffery_v2_h4.png`) as a companion still
if the site ever tells the wall story. Q5: one full orbit (T ≈ 78.5 t.u.)
needs a dense-output re-run at the L3 resolution (the "Jeffery L3" item of
the artifact backlog), frame cadence ≈ 0.5 t.u.

**Log.**
- 2026-09-16 — petrol textures for t = 90 and t = 100, geometry, axis traces,
  plots (white, transparent, dark); owner chooses t = 90 and the dark plot.
- 2026-09-16 — owner's Blender master (petrol plane) registered in PR #4.
- 2026-09-17 — Viridis banded set (20 bands) and `_lightlines` twin added on
  owner request; owner renders `jeffery_gallery_viridis_master.png`; site
  swapped to it (`8cfbe0a`).

## 5. Backlog

**Animations (Q5) — none started.** Each needs a dense-output re-run at a
frame cadence the production runs never wrote; the gallery registry's
`kind: "video"` and `poster` fields are ready.

| candidate | prerequisite run | cadence idea | poster |
|---|---|---|---|
| DKT tumble | `dkt_nofric_long` re-run with field output every ~0.25 t.u. over t = 15–35 | 80 frames | the strobe still |
| Jeffery orbit | V1b at L3 (cheaper) over one period, frames every 0.5 t.u. | ~160 frames | t = 90 still |
| viscometer | v23 continued over one bob revolution (63 t.u.) at the plateau, frames every 0.5 t.u. | ~125 frames | ghost still |
| ten Cate settling | E1 g3def frames at ≥ 20 fps of physical time | check existing `_vtk` cadence first | Fig. 9 still |

**Cases without a visualization entry** (datasheet families checked
2026-09-19; the only pictures that exist are the tracked analysis figures
and the Q0 setup scenes):

| family / case | what exists | one-line idea |
|---|---|---|
| D1.1 Hasimoto cell (`d11_l4`) | Q0 setup scene `dns_figures/d11_hasimoto_setup.png` + `.pvsm` (`tools/render_setup_scene.py`) | periodic cell with the sphere and its images tiled, mid-plane |u|/U_sup; the a_eff = a − 0.14h story as an inset of the interface band |
| D1.2 noise floor, D1 dt/level ladders | none | not a picture case; a ladder plot exists in the guide |
| D1.3 ten Cate E2–E4 | E2 wall-note plot only | E4 (Re 32) plane at h/d_p ≈ 0.5 — the wake is the visual payoff missing from E1 |
| D2.1 Brenner crossover | analysis figure `d21_brenner_crossover.png`, `d21_g2_lubrication.png` | sphere-wall approach with the film band highlighted on a zoomed plane (prescribed-motion run, frames if any) |
| D2.2 lubrication ladder G0–G3 | analysis figure `d22_g3_tencate.png` | E1 lubrication ON vs OFF final-approach planes side by side |
| D2.3 DKT frictional vs frictionless (`d23_result`, `d23_omegafix_rerun`) | none | two-column strobe: frictional (repaired binary) vs frictionless — same eight stages |
| D3.1 Beetstra arrays | Q0 setup scenes `d31_beetstra_p020_setup.png`, `d31_beetstra_p010_re30_setup.png` (+ `.pvsm`) | random array cell with images, mid-plane speed; Re 27 vs Stokes at the same φ |
| D3.2 wide column (`d32w_*`) | none | 12 d column beside the 6 d one — the confinement exponent made visible |
| D5.1 (`d51_*`, superseded by D5.2) | none | none needed |
| D5.2 lubrication pairs (`d52_v22L/v23L`), D/h = 16 rung (`d52_v25f`, closed 2026-09-19) | none | v25f final-plateau frame (t = 320, Fritz `_vtk`) as a φ = 0.05 companion at twice the resolution; lubrication ON vs OFF mid planes of v23L/v23 side by side |
| D6.1 V3 45° (`d61_v3b_transverse`) | level-3 frame under `d61_frames/v3` | third body in the Oberbeck scene with its lateral force arrow |
| D6.2 V2 clearance (`d62_v2_clearance`) | analysis figure `d62_jeffery_v2_h4.png` | H = 4 cell beside the H = 8 one at the same phase |
| legacy site pages rb2/rb3/fac3 | none in this campaign | out of scope until the owner asks |

## 6. Change log of the book

- **v1.0, 2026-09-19** — first edition. Six cases inventoried from the
  hand-off folders, the owner's `website_assets/`, `ff-redesign`
  `src/data/gallery.ts` and the site's `media/gallery/` ladders; quality
  ladder Q0–Q5 defined; conventions pinned from the READMEs and
  `scene_spec.json` files; scripts collected into `tools/viz/` (see its
  README). Open items found while inventorying, for the lead: (a) the site
  ships the viscometer ghost and the Jeffery Viridis renders under the
  earlier file names, so the file name no longer says which render it is —
  the registry entry or this book is the only record; (b) the ten Cate and
  Jeffery `scene_spec.json` files still carry the petrol map as the primary
  `color_scale` with the Viridis decision in a sub-key; (c) the viscometer
  sphere-CSV producer (the step before `export_blender_viscometer.py`) was
  not a saved script; (d) ten Cate `README` says `arrows.ply` has 30 k
  vertices, `scene_spec.json` says 20 700.
