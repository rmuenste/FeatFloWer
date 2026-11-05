# tools/viz — DNS campaign visualization scripts

Producers of the Blender hand-off folders (`blender_viz/<case>/`, untracked) and
of the review stills, collected here 2026-09-19 from the session scratchpad.
Ledger and conventions: `docs/md_docs/dns_visualization_book.md`.

Every copy is byte-identical to the one inside the hand-off folder it served
(verified with `diff` on 2026-09-19); no versions diverged. Scripts are kept
verbatim, including their hard-coded paths (listed below) — fix those before
running from another checkout.

pvpython scripts run only with the headless command
`/sfw/paraview/6.0.1/bin/pvpython --mesa --force-offscreen-rendering <script> <args>`;
plot scripts run with plain `python3` (matplotlib, numpy).

| script | runs under | purpose | arguments | case |
|---|---|---|---|---|
| `export_blender_slice.py` | pvpython | plane textures (colour / overlay / combined), colourbar, `isolines.obj`, `arrows.ply` for the y = 0 plane of a ten Cate frame; mirrors the quarter domain | `RUNDIR FRAME U_INF PZ OUTDIR [CMAP BANDS SUFFIX]` — `CMAP` preset name (default "TU Petrol Sequential"), `BANDS` filled bands (0 = continuous), `SUFFIX` file-name suffix | ten Cate sedimentation |
| `tencate_fig9.py` | pvpython | ParaView still in the paper's Fig. 9 half-plane layout (Q0 review image `reference_still_fig9_style.png`) | `RUNDIR FRAME U_INF PZ OUT.png "label"` | ten Cate sedimentation |
| `tencate_plot_dark.py` | python3 | dark-style wall note; wraps `tools/compare_tencate.py` (absolute path hard-coded) and restyles it | `RUN_LOG CASE TCREF OUT.png` | ten Cate E1 and E2 |
| `export_blender_jeffery.py` | pvpython | disturbance-field textures u' = u − γ̇ z x̂ on the full and zoom y = 0 planes, total-velocity colour, colourbar, zoom isolines/arrows geometry | `PVTU AX AY AZ OUTDIR TAG [CMAP BANDS SUFFIX LINECOL]` — `AX AY AZ` body axis at the frame, `LINECOL` "r,g,b" for light isolines/arrows | Jeffery orbit |
| `jeffery_plot_dark.py` | python3 | dark-style orbit plot; wraps `tools/d62_jeffery_analysis.py` (absolute path hard-coded) with `--gammadot 0.2 --tmin 0.5 --seams 40.01,80.01` | `AXIS_LOG OUT.png` | Jeffery orbit |
| `concat_axis.py` | python3 | concatenates the `DNS_PART_AXIS` records of the three V1b segment logs into one axis log (input of the plot); Fritz path hard-coded | none (edit `D` and `segs`) | Jeffery orbit |
| `dkt_export.py` | python3 | stages CSV/JSON, ghost sequence, full trajectory, dark phases plot, 2-D strobe preview from `particle_force.log` | `RUNDIR OUTDIR` | DKT |
| `export_blender_viscometer.py` | pvpython | mid (z = 5) and vertical (y = 0) plane textures with exact sphere holes from the logged centres, colourbar; reads `OUTDIR/<TAG>_spheres_t250.csv` (must exist first) and writes `_centres_<TAG>.npy` to a hard-coded scratch path | `PVTU OUTDIR TAG` (`TAG` = phi020 or phi005) | numerical viscometer |
| `export_blender_oberbeck.py` | pvpython | y = 0.5 plane textures (colour / overlay / combined / streamlines) and colourbar for one orientation | `PVTU AX AY AZ OUTDIR TAG` (`TAG` = parallel or perpendicular) | Oberbeck spheroid |
| `wallnote_plots.py` | python3 | the two dark wall notes: viscometer η_r(φ) and Oberbeck ratio/absolute ladders; numbers hard-coded from the datasheet rows named in its docstring | `OUT_VISCO.png OUT_OBERBECK.png` | numerical viscometer, Oberbeck spheroid |
| `hindered_export.py` | python3 | cloud snapshot CSV, trails CSV, `snapshot.json`, dark RZ plot, 2-D preview; u_t = 0.4061 hard-coded | `CLOUD_RUNDIR UT_RUNDIR OUTDIR T_SNAP` | hindered settling |
| `probes/smoke.py` | pvpython | renders one sphere — the headless-pipeline smoke test (the scratch `smokeA.py`/`smokeB.py` differed only in the output file name and were not copied) | none | environment |
| `probes/probe.py` | pvpython | lists Reflect/PLYWriter properties and exporters of the installed ParaView | none | environment |
| `probes/probe2.py` | pvpython | checks colour-map preset names (Viridis/Inferno/Turbo) and `Discretize`/`NumberOfTableValues` | none | environment |
| `probes/probe_d61.py` | pvpython | prints cell count, bounds, velocity range of a pvtu | `PVTU` | Oberbeck spheroid |
| `probes/probe_visco.py` | pvpython | prints cell count, bounds, Mixer/velocity ranges and |u| at r ≈ 5, 7.5, 10 of a viscometer pvtu | `PVTU` | numerical viscometer |

Hard-coded absolute paths (all in this checkout): `tu-colormaps.json` at the repo
root (`export_blender_slice.py`, `export_blender_jeffery.py`,
`export_blender_viscometer.py`, `export_blender_oberbeck.py`, `tencate_fig9.py`);
`tools/compare_tencate.py` (`tencate_plot_dark.py`); `tools/d62_jeffery_analysis.py`
(`jeffery_plot_dark.py`); the session scratch dir for the `.npy` side file
(`export_blender_viscometer.py`); the Fritz workspace (`concat_axis.py`).

Analysis tools that the hand-off folders also carry as reference copies are NOT
duplicated here — their tracked originals are the source: `tools/compare_tencate.py`,
`tools/d62_jeffery_analysis.py`, `tools/d32_ladder_analysis.py`,
`tools/d61_oberbeck_analysis.py`, `tools/dkt_plot_trajectory.py`,
`tools/dkt_export_series.py` (hand-off copies verified identical 2026-09-19).
The Q0 setup-scene renderer for the D1.1 / D3.1 cells is `tools/render_setup_scene.py`
(pvbatch, argparse; see `docs/md_docs/dns_figures/setup_scenes_units.md`).

Not copied: `v25f_torque.py` (scratch, D5.2 v25f torque-log extraction, not a
visualization script) and the `.png`/`.npy`/`.html` artefacts of the scratch dir.
