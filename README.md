# PV2mVeg: vegetation between solar-panel rows from sub-metre satellite imagery

Code and derived data behind the manuscript

> **Quantifying vegetation greenness between photovoltaic panel rows using sub-metre satellite imagery: A case study
> of the Gonghe Talatan solar park, China**

This repository is prepared for the review of the manuscript and is not yet public. With the data shipped here it
reproduces Tables 1–2, the computed supplementary tables, every number quoted in the text and 14 of the 20 figures,
identical to the published versions; the remaining six figures and the analysis itself need restricted inputs that
the authors provide on request (see [Data availability](#data-availability)).

## Version and source

| | |
| --- | --- |
| Manuscript | version 2.0 |
| Source | the authors' research repository, branch `paper/en-pv-ecology`, commit `d065677` (2026-10-16), from which the v2.0 manuscript, supplement and figures were built |
| Later commits | `bf9bf41`, `2d59183` and `36aa0a6` (2026-10-18) check the reference list and change the AgERA5 citation to version 2.0; they change no code, data or number used here. The results here use the AgERA5 series of the archived analysis (see Public data) |
| Analysis | the v1.2 analysis of the research repository, unchanged up to v2.0 (`analysis/README.md`) |

## Quick start

Python 3.12 is required (the results were produced with CPython 3.12 on Windows 11).

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt      # Linux/macOS: .venv/bin/python
.venv/Scripts/python scripts/30_tables.py                    # Tables 1-2, Table S1
.venv/Scripts/python scripts/31_numbers.py                   # every number of the text -> output/numbers/numbers.json
.venv/Scripts/python scripts/40_figures.py                   # the 14 figures of tier "repo" -> output/figures/
.venv/Scripts/python scripts/90_verify.py                    # comparison with the published values (reference/)
```

`90_verify.py` ends with `VERIFY PASSED` when the numbers and tables are identical to the published ones, the
plotted-value tables of the figures are identical to those of the release, and every recomputed cell of the
supplementary tables agrees. The figures are pixel-identical to the published PNG files when the figure font of the
publication is configured (see [Figure font](#figure-font)); without it they are drawn with a fallback font and
only checked for being rendered. `python -m pytest` runs the same checks in a temporary folder (a few minutes;
`-m "not slow"` leaves out the figures).

## Repository structure

```
README.md            this file
LICENSE              placeholder: the licence is to be decided by the authors
requirements.txt     pinned Python packages (the environment of the published results)
pyproject.toml       package metadata of src/pv2mveg
config/              paths.example.toml: data layout, data and output folders, figure font
src/pv2mveg/         configuration, data readers, figure style, tables, numbers, verification, analysis runner
figures/             one script per figure (fig01_*.py ... fig08_*.py, figS01_*.py ... figS12_*.py)
scripts/             numbered entry points: 10-13 analysis (restricted inputs), 30 tables, 31 numbers,
                     40 figures, 90 verification
analysis/            the v1.2 analysis code, mirrored from the research repository (analysis/README.md)
tools/               export_derived.py (writes data/derived), compare_analysis.py, make_reference.py
data/derived/        de-identified derived tables (about 28 MB) with MANIFEST.json (data/README.md)
data/restricted/     placeholder for the inputs provided on request (not tracked)
reference/           published values: numbers.json, Tables 1-2, supplementary tables, figure pixel hashes,
                     plotted-value tables of the release
docs/                PIPELINE.md (from the source data to the figures), VERIFICATION.md (results of the checks)
tests/               pytest: numbers, tables, supplement cells, figures, data integrity, no local paths
```

## What reproduces what

Tier **repo**: reproduced from `data/derived` in this repository. Tier **restricted**: the code is here and was
verified on the authors' project data; it needs the restricted inputs (layout `project`).

| Item | Script | Inputs | Tier |
| --- | --- | --- | --- |
| Table 1 (units, contrasts by age group and baseline) | `scripts/30_tables.py` | `v12/cross_scene/effects.csv`, `cells_<tag>.csv` | repo |
| Table 2 (share reproducible by imaging) | `scripts/30_tables.py` | `v12/cross_scene/optics.csv` | repo |
| Table S1 (units by satellite) | `scripts/30_tables.py` | `v12/buffer_primary/cells_08.csv` | repo |
| Tables S5, S8–S13 (every printed estimate, bound and count) | `scripts/90_verify.py` (`pv2mveg.supplement`) | `v12/cross_scene`, `v12/buffer_primary`, `v12/cross_scene/secondary` | repo |
| Numbers of the text, abstract, highlights and captions (3,243 values) | `scripts/31_numbers.py` (`pv2mveg.numbers`) | all of `data/derived` | repo |
| Fig. 1 study area, unit and cross-section | `figures/fig01_study_area.py` | Beijing-3 true colour and panchromatic scene, station polygons, unit coordinates, row-band and clearance layers | restricted |
| Fig. 2 workflow | `figures/fig02_workflow.py` | `numbers.json` (two counts) | repo |
| Fig. 3 mixed-pixel event study | `figures/fig03_event_study.py` | `ts/ew/` event-study curves, `legacy/fig2_*.csv` | repo |
| Fig. 4 time axis, two records, three surfaces | `figures/fig04_time_axis.py` | `v12/cross_scene` (unit tables, effects, lanes), `secondary/1p5m/chrono/` (block labels), `ts/ew/` curves | repo |
| Fig. 5 matched and unmatched baselines | `figures/fig05_baselines.py` | Beijing-3 true colour, unit coordinates, matched controls and Landsat series of the units | restricted |
| Fig. 6 spatial pattern and station means | `figures/fig06_spatial_pattern.py` | Beijing-3 true colour, unit coordinates, station polygons | restricted |
| Fig. 7 spectral, red-edge, distance and clearance evidence | `figures/fig07_spectral_evidence.py` | `v12/cross_scene/{spectral,distance,effects}.csv`, `1m/`, `2m/`, `red_edge/e2_indices.csv` | repo |
| Fig. 8 share reproducible by imaging alone | `figures/fig08_imaging_share.py` | `v12/cross_scene/optics.csv` | repo |
| Fig. S1 radiometric calibration | `figures/figS01_calibration.py` | calibration samples and metrics of the Beijing-3 scenes | restricted |
| Fig. S2 unmixing tests | `figures/figS02_unmixing_tests.py` | `purify/` | repo |
| Fig. S3 installation-year detection | `figures/figS03_installation_year.py` | Landsat annual composites, detected years, station grid of the park | restricted |
| Fig. S4 unit-level distributions | `figures/figS04_unit_distributions.py` | `v12/cross_scene/cells_<tag>.csv` | repo |
| Fig. S5 common-footprint changes | `figures/figS05_common_footprints.py` | `secondary/crossperiod/` | repo |
| Fig. S6 mixed minus inter-row contrast | `figures/figS06_mixed_minus_interrow.py` | `secondary/1p5m/statistics/mixed_bias.csv` | repo |
| Fig. S7 single-year cohorts | `figures/figS07_single_year_cohorts.py` | as Fig. 4 | repo |
| Fig. S8 coverage decomposition | `figures/figS08_coverage.py` | Beijing-3 true colour, unit coordinates, station polygons | restricted |
| Fig. S9 visible and NIR spectral contrasts | `figures/figS09_spectral_contrasts.py` | `secondary/measurements/<tag>/subset_summary.csv` | repo |
| Fig. S10 distance bands | `figures/figS10_distance_bands.py` | `v12/cross_scene/distance.csv` | repo |
| Fig. S11 imaging and boundary sensitivity | `figures/figS11_imaging_boundary.py` | `v12/cross_scene/optics.csv`, `secondary/measurements/<tag>/boundary_paired.csv` | repo |
| Fig. S12 mixed-pixel rainfall diagnostic | `figures/figS12_rainfall.py` | `frozen/rainfall_frozen.csv` | repo |
| v1.2 analysis results (`v12/`) | `scripts/10_units.py` … `13_secondary.py` (`analysis/`) | frozen inputs of the analysis (`docs/PIPELINE.md`) | restricted |

Paths in the Inputs column are relative to `data/derived/` (`secondary/` = `v12/cross_scene/secondary/`). Each
figure script also writes the values it plots to `output/figdata/`; `pv2mveg.catalog` lists them, and
`90_verify.py` compares them with those of the release (`reference/figdata/`). Not reproduced here: Tables S2–S3
(blind interpretation of earlier mask designs, counts from the annotation sheets), Table S4 (rules), and Tables
S6–S7 (identification tests and robust inference of the mixed-pixel event study, outputs of the event-study
release; their producers are listed in `docs/PIPELINE.md`).

## Run order

**Repository data (layout `repo`, the default):** `30_tables.py` → `31_numbers.py` → `40_figures.py` →
`90_verify.py`, as in [Quick start](#quick-start). `40_figures.py` draws each figure in its own process
(`-j N` in parallel; names such as `Fig4 FigS7` select figures).

**Project data (layout `project`, restricted inputs):** set `layout = "project"` and `project_root` in
`config/paths.toml` (copy of `config/paths.example.toml`), then

```bash
python scripts/30_tables.py
python scripts/40_figures.py Fig6 FigS8     # write the two summaries that the numbers read
python scripts/31_numbers.py
python scripts/40_figures.py                # all 20 figures
python scripts/90_verify.py
```

**Analysis (layout `project`):** `10_units.py` → `11_measure.py` → `12_estimate.py` → `13_secondary.py`
(about 50 minutes; `analysis/README.md`), then `python tools/compare_analysis.py` compares the results, file by
file, with the archived v1.2 results of the project data. The figure, table and number scripts read the archived
results; `docs/PIPELINE.md` shows how the analysis results become `data/derived`.

## Configuration

`src/pv2mveg/config.py` reads, in order, the defaults, `config/paths.toml`, a file named by `PV2MVEG_CONFIG` and the
environment variables `PV2MVEG_LAYOUT`, `PV2MVEG_DATA_ROOT`, `PV2MVEG_PROJECT_ROOT`, `PV2MVEG_OUTPUT` and
`PV2MVEG_FONT`. Outputs go to `output/` (`figures/`, `figdata/`, `numbers/`, `tables/`, `checks/`, and
`analysis/v12/` for the analysis). No path of the authors' machines is hard-coded.

## Environment

`requirements.txt` pins every package imported by the verification runs, at the versions of the authors'
environment (Python 3.12, numpy 2.5, pandas 3.0, matplotlib 3.11, rasterio 1.5, geopandas 1.1). Pixel identity of
the figures was verified with exactly these versions; other matplotlib or FreeType versions can move anti-aliased
pixels without changing the plotted values. The analysis scripts use the same environment; they run with single-
threaded BLAS (`OMP_NUM_THREADS=1`) and UTF-8 mode, both set by `pv2mveg.analysis`.

## Data availability

**In this repository** (`data/derived`, about 28 MB): unit-level and aggregated result tables of the analysis,
the mixed-pixel event-study curves, the rows of the unmixing tests shown in Fig. S2, the frozen rainfall regression
and two summaries of the restricted map figures. They are de-identified: no coordinates, construction or detection
years of units, or station labels; unit, station and block identifiers are order-preserving pseudonyms; array age is
kept (`data/README.md`). `data/derived/MANIFEST.json` gives the source file and SHA-256 of each table.

**Restricted, on reasonable request to the authors:** the Beijing-3 imagery (commercial; 0.5 m panchromatic, 2 m
and 1.2 m multispectral) and everything derived at pixel level from it (calibrated products, row-band maps,
footprint layers, image chips), the hand-drawn station polygons with their construction years, the unit coordinates
and the intermediate tables that hold them. Use of the imagery is subject to the provider's licence. These inputs
are delivered in the folder layout of the project data and used with the layout `project`.

**Public data** used upstream are not redistributed; the research repository downloaded them as follows
(`docs/PIPELINE.md` names the scripts):

| Data | Source |
| --- | --- |
| Landsat 5/7/8/9 Collection 2 Level-2, day of year 121–288 of 2000–2026, cloud cover < 80% (674 scenes) | Microsoft Planetary Computer STAC, collection `landsat-c2-l2` |
| Sentinel-2 Level-2A, 2019–2026 (354 scenes); red, NIR, SCL, and blue and green for the visible calibration | Earth Search v1 STAC (Element 84), collections `sentinel-2-c1-l2a` and `sentinel-2-l2a` |
| Copernicus GLO-30 DEM, tiles N35_E100 and N36_E100 | `copernicus-dem-30m` open bucket (AWS) |
| ESA WorldCover 10 m 2021 v200, tiles N33E099 and N36E099 | `esa-worldcover` open bucket (AWS) |
| AgERA5 daily agrometeorological indicators, point series at 100.30°E, 36.01°N, 2000–2025 (version 2.0) | Copernicus Climate Data Store |

## De-identification and its limits

The analysis units are cells of a fixed 150 m grid; their identifiers, coordinates and construction years would
locate and date the stations. The repository tables therefore carry neither coordinates nor years, and identifiers
are replaced by pseudonyms. Three limits remain, which the authors should weigh before making the repository
public:

* **Array age is kept**, because every result is stratified by it; with the acquisition year it gives the
  construction year of a unit, though not its place.
* **Pseudonyms preserve order** (dense ranks of the original identifiers), so that every sort, grouping and
  resampling of the code gives the published results. The original cell numbers run row by row over the grid, so
  the order of the pseudonyms still reflects the north–south order of the cells.
* **Three original cell numbers remain**: the scripts of Fig. 1 and Fig. 5 select their example units by their
  original numbers, which the plotted-value tables of the release also record (`reference/figdata/fig1_counts.json`,
  `fig4_summary.json`). With the grid defined in `analysis/eco/common.py` and `analysis/eco/rowband/rb_common.py`
  they give the place of these three units, which Figs. 1 and 5 show on the image anyway. The other files of
  `reference/` and `data/derived` hold no coordinates, years or original identifiers (checked by key and by value
  range).

## Verification

`docs/VERIFICATION.md` reports the checks run for this repository: the published numbers, tables, supplementary
cells and plotted values reproduced exactly from `data/derived` and from the project data; the 20 figures
pixel-identical in the layouts that can draw them; and the mirrored analysis code reproducing the archived v1.2
results file by file.

## Licence and citation

The licence is to be decided by the authors (`LICENSE`). Please cite the manuscript; a citation entry will be added
on publication.
