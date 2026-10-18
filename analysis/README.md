# Analysis code of the v1.2 analysis (mirrored)

The results reported in the manuscript (version 2.0) come from the v1.2 analysis of the authors' research
repository; versions 1.3 to 2.0 changed the text and the figures, not the analysis. This folder holds the code of
that analysis, mirrored from the research repository (branch `paper/en-pv-ecology`, commit `d065677`, the v2.0
manuscript). It needs the restricted inputs (Beijing-3 products, station polygons and construction years, unit
coordinates; see `../data/README.md`) and runs through `../scripts/10_units.py` to `../scripts/13_secondary.py`.

## What was mirrored, and how

The four analysis scripts and the 44 modules that they execute were mirrored, in their original folder layout:

| Script | Stage of the analysis | Wrapper |
| --- | --- | --- |
| `experiments/bright_geometry/v12_units.py` | units, timing, matched controls; snapshot of the frozen inputs | `scripts/10_units.py` |
| `experiments/bright_geometry/v12_measure.py` | native 2 m subset measurements of the phased units | `scripts/11_measure.py` |
| `experiments/bright_geometry/v12_run.py` | main estimates, buffer baseline, lanes, spectra, distance bands, red edge, imaging scenarios | `scripts/12_estimate.py` |
| `experiments/bright_geometry/v12_secondary.py` | supplementary diagnostics | `scripts/13_secondary.py` |

The modules were selected by running the four scripts on the project data with an audit hook that recorded every
imported module, every file read or written and every executed function. The mirror is mechanical:

* **Pruned.** Top-level functions and classes that were not executed and are not referenced by executed code were
  removed; the mirrored modules have about 5,300 lines, half of the 10,515 lines of the original modules. Module and
  function names, the folder layout and the remaining code are unchanged, because the code imports its modules by
  flat names after `sys.path` inserts, loads `eco/common.py` by path, and patches functions at run time from their
  source text (`orientation_diagnostics.rewrite`, `v12_run.optics`), which needs the original function bodies.
* **Roots configured.** After its docstring, each module has a five-line header (a provenance comment and four lines
  of code) that imports `pv2mveg.legacy` (`../src/pv2mveg/legacy.py`), which provides the roots of the
  configuration (`../config/paths.example.toml`):

  | Name | Meaning | Replaces |
  | --- | --- | --- |
  | `PROJECT` | project data root (`project_root`; read only) | the absolute project path of the research machine |
  | `CODE` | this folder | the code root used for module search paths |
  | `V12_OUT` | `<output>/analysis/v12`, where everything is written | `data/eco/releases/fixed_ew_v2_orientation/analysis/v12` |
  | `FONT_FILE` | figure font of the diagnostic plots (configured font, else STIX) | the absolute font path |
  | `SELECTED_STATIONS` | labelled station outlines of the earlier installation-year check | the absolute path beside the project |

  The lines that used these paths now use the names above (marked `pv2mveg configuration`). The analysis writes
  only below `V12_OUT`: the verification runs blocked every write into the project data.
* **Imports trimmed.** Imports of modules that the v1.2 analysis does not execute were removed, notably the
  release-building chain (`build_candidate` → `array_constraints` → `reconstruct`, `round2`), which needs PyTorch and
  OpenCV and produced the row-band release that the analysis only reads; unused third-party imports were removed too.
* **Comments in English.** Chinese comments and docstrings were translated; the syntax tree of every module, with
  docstrings blanked, is identical before and after the translation, so no code changed. Module docstrings describe
  the original modules, including command-line uses and outputs that the v1.2 analysis does not use.
* **Kept verbatim.** 24 lines keep Chinese string literals: run-time log messages, labels of two diagnostic plots
  (`eco/adjacency/e2_spec_analysis.py`, not used in the manuscript), the value `'是'` ("yes") of an exclusion list,
  scene-geometry labels of `eco/rowband/rb_common.py`, and the patterns that `orientation_diagnostics.rewrite`
  must find in the source of the functions it patches. The most frequent message,

  ```
  [cc.scope_filter] SCOPE=ew tag=06：排除单元清单 主清单 3514；合计 3514 个单元，本表命中 3514，其中原本在口径内被剔除 282
  ```

  reads "exclusion lists of units: main list 3514; 3514 units in total, 3514 present in this table, 282 of which
  were inside the east–west scope and are removed".

## Running it

```bash
python scripts/10_units.py       # about 20 s
python scripts/11_measure.py     # about 8 min
python scripts/12_estimate.py    # about 30 min
python scripts/13_secondary.py   # about 10 min, one process per stage
python tools/compare_analysis.py # file-by-file comparison with the archived results
```

Times are those of the authors' workstation (32 logical processors). The wrappers need `layout = "project"` and
`project_root` (`../config/paths.example.toml`), start every script, and every stage of the supplementary
diagnostics, in a fresh process with UTF-8 mode on and single-threaded BLAS, and stop at the first failure. One process per stage matters:
several diagnostics resample with a module-level random generator (`eco/chrono/j_review.py`), so their confidence
intervals depend on the sequence of calls within a process; the archived results were produced stage by stage.

The outputs have the layout of the archived v1.2 folder: `inputs/` (frozen inputs with a SHA-256 manifest),
`units_validation.json`, `measurements/<tag>/`, `cross_scene/` (primary design at 1.5 m; `1m/`, `2m/`; `aux/`
red edge; `secondary/` diagnostics) and `buffer_primary/`. In the legacy-format tables the prefix `same_` and the flag
`valid_same` denote the primary design, and `station` holds the hand-drawn station number.

## Map of the modules

| Function (manuscript section) | Modules |
| --- | --- |
| Units, timing and controls (3.2, 3.3) | `experiments/bright_geometry/v12_units.py`; matching and event-time tools `eco/ts/es.py`; Landsat unit series `eco/ts/common.py`; scope, exclusion lists and control means `eco/chrono/cc.py` |
| Footprints and subset measurements (3.1) | `v12_measure.py`, `revision_measure.py`, `orientation_diagnostics.py` (adapters of the earlier code to the repaired row-band release), `fixed_core.py`, `core_geometry.py`, `dataio.py`, `ecology_revision.py`; row-band and module-rule layers `eco/rowband/rb_common.py`, `eco/rowband/b_component.py`; 2 m product configuration `eco/bj3ndvi/cfg.py`, `eco/bj3_calibration.py` (scene list) |
| Matched estimates and resampling (3.3, 3.5) | `cross_scene.py` (primary cross-scene design), `orientation_statistics.py`, `orientation_release.py`, `eco/chrono/a1_chrono.py`; station-clustered bootstrap `eco/chrono/j_review.py` |
| Unmatched buffer baseline (3.3) | `buffer_primary.py`, `buffer_baseline.py`, `buffer_variants.py` |
| Spectral and red-edge diagnostics (3.6) | `cross_scene_mechanisms.py` (`spectral`, `distance`), blue and green indices `eco/chrono/k_bg_greenness.py`, red edge of the 1.2 m imagery `eco/adjacency/e2_spec_analysis.py`, `revision_aux_summary.py` |
| Imaging scenarios (3.6) | `cross_scene_mechanisms.py` (`optics`), `revision_optics.py`, `revision_pan.py`, `eco/adjacency/a2_cells.py`, `eco/adjacency/adjc.py`, `eco/adjacency/a5_model_set.py`, panchromatic profiles `eco/panprof/p1_extract.py`, `eco/panprof/pp.py`, texture `eco/adjacency/e3_pan_texture.py`, `eco/adjacency/e6_texture_analysis.py` |
| Supplementary diagnostics (supplementary tables and figures) | `v12_secondary.py` with `revision_statistics.py`, `revision_balance.py`, `revision_placebo.py`, `revision_subsets.py`, `revision_crossperiod.py`, `revision_temporal.py`, `revision_forward.py` |
| Shared grid definitions | `eco/common.py`, `fixed_ecology.py`, `fixed_release.py` |

Modules listed without a folder are in `experiments/bright_geometry/`. `../docs/PIPELINE.md` describes the inputs of
each stage and the upstream steps that produced them.
