# Pipeline: from the source data to the figures

```
 public data: Landsat C2 L2, Sentinel-2 L2A, Copernicus DEM, WorldCover, AgERA5 ─┐
 Beijing-3 imagery (commercial, restricted) ─────────────────────────────────────┤  A. upstream steps
 station polygons with construction years (restricted) ───────────────────────────┘     (research repository; not in this repository)
                                   │
                                   ▼
 frozen inputs of the v1.2 analysis (restricted; project data tree)
                                   │  B. analysis/ via scripts/10_units.py … 13_secondary.py
                                   ▼
 v1.2 results (<project>/data/eco/releases/fixed_ew_v2_orientation/analysis/v12; restricted: unit identifiers and years)
                                   │  C. tools/export_derived.py (de-identification)
                                   ▼
 data/derived (this repository)
                                   │  D. scripts/30_tables.py, 31_numbers.py, 40_figures.py, 90_verify.py
                                   ▼
 Tables 1–2 and S1, numbers of the text, Figs. 1–8 and S1–S12, verification against the published values
```

Steps B to D are in this repository and were run for the verification (`VERIFICATION.md`). Step A produced the
frozen inputs of the analysis; it is documented below with the producing scripts of the research repository (branch
`paper/en-pv-ecology`, commit `d065677`) so that every input can be traced, but its code is not included: it needs
the raw imagery and manual steps of the authors (drawing and dating the station polygons, reviewing mask failures),
and its outputs enter the analysis as frozen files (`<output>/analysis/v12/inputs/manifest.json` lists the source
and SHA-256 of the inputs that the analysis snapshots).

Paths in this document are relative to the project root of the research repository unless they start with a folder
of this repository. `REL` = `data/eco/releases/fixed_ew_v2_orientation`. Restricted = contains Beijing-3 pixels,
station polygons or construction years, or unit coordinates.

## A. Upstream steps (research repository)

### A1. Radiometric calibration of the Beijing-3 imagery (Section 2.2.1, Fig. S1)

| Product | Producer | Inputs |
| --- | --- | --- |
| `data/eco/bj3cal/{coef_per_scene.csv, coefficients.json, cache/}` | `eco/bj3_calibration.py`, stages `s2 → s2check → prep → cal` (the stages `trace`, `overlap`, `prof`, `yoy`, `figs` are diagnostics) | Beijing-3 L4 scenes; 2 m panel fraction of the row-band maps; Sentinel-2 C1 L2A B04/B08/SCL from Earth Search (public); Copernicus DEM (trace stage) |
| `data/eco/bj3cal/s2_bg/coefficients_bg.json` (blue, green; BJ3A1) | `eco/bj3ndvi/s2_bg.py fetch → prep → cal` | Sentinel-2 B02/B03 (public), BJ3A1 bands 1–2 |
| `data/eco/bj3ndvi/<tag>/` (calibrated 2 m red/NIR, QA, classes, 30 m sums) | `eco/bj3ndvi/run.py [config] → flags → product → support → qc → regcheck`; `bgprod.py` for the blue/green sums | the scenes, the panel masks, the coefficients |
| `data/eco/chrono/celltab_<tag>_{raw,fe_red,fe_nir}.parquet` (2 m sums per 150 m unit, scene offsets) | `eco/chrono/a0_tables.py 06 08 2026_08` | 30 m sums of `data/eco/bj3ndvi/<tag>/`, scene table, seam offsets |

Fig. S1 reads the calibration diagnostics (`pv2mveg.data`, tier "restricted").

### A2. Row-band mapping and footprints (Section 3.1)

| Product | Producer | Inputs |
| --- | --- | --- |
| 0.5 m panel masks of each acquisition (`outputs/06`, `outputs/08`, `data/eco/mask2026/2026_08`) | U-Net segmentation and row constraints: `src/talatan`, `scripts/`, `eco/mask2026/run_pipeline.ps1` | Beijing-3 panchromatic scenes |
| common grids `data/eco/grid`, `data/eco/grid2m/<tag>` (panel fraction, array interior, distance to arrays, station id, orientation, pitch, elevation, slope) | `eco/build_frac.py` → `build_zones.py` → `build_config.py` → `grid2m.py` → `build_dem.py` | masks, labelled station outlines, Copernicus DEM |
| fixed-array release `data/eco/releases/fixed_ew_v1` | `experiments/bright_geometry/fixed_release.py` → `compose_fixed_release.py` → `audit_fixed_release.py` → `fixed_core.py --core-policy {legacy,exact,clearance1p5,clearance1p0}` | masks, grids, 2 m products |
| repaired release `REL/<tag>/{rowband,certified_ew}_0p5m.tif`, 2 m clearance layers `REL/ecology_1p5_v3/<tag>/`, native footprint sums `REL/native/<tag>/native_sums.npz` | `orientation_repair.py --tags 06,08,2026_08` → `verify_orientation_repair.py <tag>` → `orientation_release.py --compose --measure` (all in `experiments/bright_geometry/`) | the v1 release, panchromatic scenes, quality exclusions |
| subset measurements of the original cohort `REL/analysis/measurements/<tag>/` | `orientation_diagnostics.py --stage cache`, `--stage measure` | release layers, 2 m products, blue/green bands |

### A3. Landsat time series, installation timing and event study (Sections 2.2.2, 3.2, 3.4; Figs. 3, S3)

| Product | Producer | Inputs |
| --- | --- | --- |
| Landsat scenes and annual midsummer composites `data/eco/ts/ls_annual/` | `eco/fetch_landsat.py`, `eco/ts/fetch_2026.py`, `eco/ts/s1_harmonize.py` → `s2_composite.py` | Landsat C2 L2 from Planetary Computer (public) |
| detected installation year `data/eco/ts/detect_30m.tif`, unit series `data/eco/ts/cell_y.npz`, unit pool `data/eco/ts/cells_pool.parquet` (east–west: `ts/ew/cells_pool.parquet`) | `eco/ts/s0_worldcover.py`, `s3_detect.py`, `s4_units.py`, `s5_pool.py`, `ew0_cells.py --thr 0.9` (`TS_SCOPE=ew`) | composites, WorldCover (public), grids, exclusion lists |
| matched event study `data/eco/ts/ew/es_raw_sum_<index>_det_treated/{curves.csv, matched.npz}`, placebo `ts/ew/placebo_sum_ndvi_det/` | `eco/ts/s6_event.py --timing det --outcome sum_<index> --B 400` and `s7_placebo.py --timing det --R 150 --B 200` (`TS_SCOPE=ew`, `run_ew.ps1`) | unit pool, unit series, AgERA5 table |
| estimator comparison `ts/ew/review2/X1_same_target.csv`, rainfall regression `ts/ew/review2/X5_precip_null.csv` | `eco/ts/ew_extra.py` (`TS_SCOPE=ew`) | the same, and the precipitation anomaly |
| relative-magnitudes sensitivity `ts/ew_honest/` | `eco/ts/h1_honest_boot.py` → `h2_honest_run.py` → `h3_honest_validate.py` | event-study bootstrap draws |
| AgERA5 growing-season table `data/eco/climate/agera5_annual.csv` | `eco/climate_table.py` | AgERA5 daily data (public) |

The event-study curves, the estimator table, the unit counts and the band effects enter this repository as
`data/derived/ts/ew/` and `data/derived/legacy/fig2_*.csv` (written by the earlier figure script
`paper_eco/figcode/fig2.py`); the rainfall regression as `data/derived/frozen/rainfall_frozen.csv`.

### A4. Station polygons, construction years and analysis units (Section 3.2)

| Product | Producer | Inputs |
| --- | --- | --- |
| hand-drawn station polygons with construction years `data/eco/anomaly/v4/hand_stations_4543.gpkg`; unit exclusion lists `data/eco/ts/ew/exclude_cells*.txt` | `eco/anomaly/h1_entries.py` … `h5_apply_decisions.py` (drawn and dated by the first author) | 0.5 m imagery, annual Landsat images |
| units of the coverage expansion `data/eco/coverage_expansion/units/cells_expansion.parquet` | `experiments/coverage_expansion/x01_units.py` | unit pool, polygons, detected years, orientation |
| matched controls on the construction-year window `coverage_expansion/match/matched_{parent_hand,ew_hy}.npz` | `x03_match.py` (the matching of `eco/ts/es.py`, Mahalanobis with distance penalty) | unit pool, unit series, AgERA5 table |
| native footprint sums and gap profiles of the phased units, per-unit estimates of the expansion | `x02_native.py --set ew_hy`, `m01_gap_profile.py --set ew`, `x05_estimate.py --set {ew_hy,parent_hand}` | release layers, 2 m products |

`scripts/10_units.py` copies these into `<output>/analysis/v12/inputs/` and checks that the per-unit contrasts of
the merged units equal the per-unit estimates of the expansion.

### A5. Imaging kernels, panchromatic profiles and red edge (Section 3.6)

| Product | Producer | Inputs |
| --- | --- | --- |
| band-response kernels, held-out tests, shadow parameters `data/eco/adjacency/{k1_*,k2_*}`, `ew_only/{k3_*,a4_shadow_params.csv}` | `eco/adjacency/k1_kernel.py` → `k2_summary.py` → `k3_holdout.py`, `a4_ew_forward.py` (`run_ew_only.ps1`) | panchromatic and multispectral scenes, masks |
| per-unit forward models and optical scenario composition `REL/analysis/optical_scenario_composition.csv` | `orientation_diagnostics.py --stage a2`, `--stage parameters`, `--stage a4`, `--stage optics-summary` | kernels, scenes, masks, 2 m products |
| 1.2 m six-band profiles `REL/analysis/1p5m/adjacency/spectral/e1_<date>_{prof,ctrl}.parquet` | `orientation_diagnostics.py --stage auxiliary`, `--stage aux-summary`; registration `eco/adjacency/b1_bj3n.py` | Beijing-3 1.2 m imagery, release row bands |

### A6. Unmixing tests (Fig. S2)

The panel end-member and unmixing tests of Fig. S2 come from an earlier experiment (`data/eco/purify/`); the rows
shown are in `data/derived/purify/`.

## B. The v1.2 analysis (this repository)

`analysis/` holds the code (`analysis/README.md`); `scripts/10_units.py` … `13_secondary.py` run it.

| Step | Function | Main inputs | Outputs (`<output>/analysis/v12/`) |
| --- | --- | --- | --- |
| `10_units.py` | units, timing, controls; snapshot of A4 | A4, `ts/cells_pool.parquet`, `ts/cell_y.npz` | `inputs/`, `units_validation.json` |
| `11_measure.py` | lanes, distance bands, gap sides, 0.5 m boundary perturbations, visible bands of the phased units | A2, `inputs/` | `measurements/<tag>/` |
| `12_estimate.py` `estimate` | matched cross-scene contrasts of inter-row ground and mixed pixel at 1.5 (primary), 1 and 2 m; 600 joint resamples of stations and control blocks plus calibration variance | `chrono/celltab_<tag>_*.parquet` (A1), native sums (A2, A4), `inputs/` | `cross_scene/{effects.csv, cells_<tag>.csv}`, `1m/`, `2m/` |
| `buffer` | unmatched buffer baseline and its variants, paired resamples | A2, `ts/cells_pool.parquet` | `buffer_primary/` |
| `lanes`, `spectral`, `distance` | service lanes, band differences (blue/green on BJ3A1), distance bands | `measurements/`, A1 | `cross_scene/{lanes,spectral,distance}.csv` |
| `aux` | red edge of the 1.2 m imagery at three clearances | A5 | `cross_scene/aux/`, `secondary/aux/` |
| `optics` | share reproducible by imaging scenarios | A5 | `cross_scene/optics.csv` |
| `annotate`, `terrain`, `origin` | origin and construction year of the units, terrain of the buffer cells | `inputs/` | unit tables, `buffer_primary/terrain.json`, `cross_scene/origin.csv` |
| `13_secondary.py` | balance and held-out checks, leave-one-control-out, robustness, leave-one-station-out, regressions, mixed-pixel bias, common footprints, boundary perturbations, forward/backward checks | the outputs above, `REL/analysis/crossperiod/` caches, `ts/ew/anchored/` | `cross_scene/secondary/` |

## C. De-identified derived data

`tools/export_derived.py` (run with the project layout) copies the tables that steps D read from the v1.2 results
and the upstream products into `data/derived/`, dropping coordinates, construction and detection years and station
labels and replacing unit, station and block identifiers by order-preserving pseudonyms (`data/README.md`).
`data/derived/MANIFEST.json` lists the source of each file and the SHA-256 of source and copy.

## D. Tables, numbers and figures

| Script | Output | Reads |
| --- | --- | --- |
| `scripts/30_tables.py` | Tables 1–2 (`tables/table1.md`, `table2.md`, CSV), Table S1, the held-out checks read by the numbers | `v12/cross_scene` (effects, unit tables, optics), `v12/buffer_primary/cells_08.csv`, `v12/cross_scene/secondary/1p5m/statistics/holdout.csv` |
| `scripts/31_numbers.py` | `numbers/numbers.json`: every number of the text | all of `data/derived` (`pv2mveg.numbers`) |
| `scripts/40_figures.py` | `figures/<name>.png/.pdf`, plotted values in `figdata/` | per figure: `pv2mveg.catalog` and README.md |
| `scripts/90_verify.py` | `verification.json` | `reference/` (published values) |

With the project layout, the same scripts read the full project tree and also draw the six restricted figures
(Figs. 1, 5, 6, S1, S3, S8); Figs. 6 and S8 then write the two summaries that `31_numbers.py` reads (in the
repository layout these are shipped in `data/derived/summaries/`).
