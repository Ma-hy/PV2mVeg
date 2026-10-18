# Verification

This record lists the checks run on 2026-10-18 for this repository and their results. All runs used Windows 11
(x86-64, 32 logical processors) and CPython 3.12.13. The published values were taken from the research repository
at commit `d065677` (`tools/make_reference.py`, stored in `reference/`): the numbers of the text
(`numbers.json`), Tables 1–2 (Markdown and CSV), the pipe tables of the supplement, the size and SHA-256 of the
RGBA pixels of the 20 published figures, and the plotted-value tables written by the figure scripts of the release.

## 1. Fresh clone and fresh environment

The repository at commit `7d4e492` was cloned into an empty folder, and a new environment was created from
`requirements.txt` with CPython 3.12.13 (`uv venv`, `uv pip install -r requirements.txt`; 39 packages, exactly
the pinned set). All checks below that mention the clone ran in this clone and environment.

| Check | Result |
| --- | --- |
| `python -m pytest` (publication font configured) | 10 passed in 82 s: derived data match `MANIFEST.json`, no identifying columns, no local paths, mirrored analysis code parses and has English comments, the analysis refuses the repository layout, numbers, tables and supplementary cells identical, the 14 repository figures pixel-identical with identical plotted values |

## 2. Repository data (layout `repo`)

| Check | Publication font | Fallback font (STIX) |
| --- | --- | --- |
| `numbers.json` (3,243 values) | byte-identical | byte-identical |
| Tables 1–2 (Markdown, CSV) | byte-identical | byte-identical |
| Table S1 | identical cells | identical cells |
| Supplementary tables S5, S8–S13 (276 cells or groups) | all agree | all agree |
| Plotted-value tables of the 14 figures (17 tables) | byte-identical | byte-identical |
| Figures of tier "repo" (14) | pixel-identical | rendered; pixels differ (expected) |

Publication font: the clone, `30_tables.py` → `31_numbers.py` → `40_figures.py` → `90_verify.py` (`VERIFY PASSED`).
Fallback font: the same scripts in the working copy (the code of figures, tables and numbers of `7d4e492`) without
`PV2MVEG_FONT` (`VERIFY PASSED`; the figures are drawn with STIXGeneral and were inspected visually).

## 3. Project data (layout `project`, restricted inputs)

The clone, with `layout = "project"` and the authors' project data as `project_root`, ran `30_tables.py` →
`40_figures.py Fig6 FigS8` → `31_numbers.py` → `40_figures.py` → `90_verify.py` with the publication font
(`VERIFY PASSED`):

| Check | Result |
| --- | --- |
| `numbers.json` (3,243 values) | byte-identical |
| Tables 1–2 (Markdown, CSV), Table S1 | byte-identical; identical cells |
| Supplementary tables S5, S8–S13 (276 cells or groups) | all agree |
| Figures (20: Figs. 1–8, S1–S12) | all pixel-identical |
| Plotted-value tables (23) | byte-identical |

An audit hook refused every write into the project data during this run; no write was attempted.

## 4. The analysis

The archived v1.2 results (`<project>/data/eco/releases/fixed_ew_v2_orientation/analysis/v12`, 17,157 files)
were regenerated three times and compared file by file with `tools/compare_analysis.py`: byte identity, or, for the
stage records (`complete.json`), equality apart from run times and dates.

| Run | Code | Environment | Result |
| --- | --- | --- | --- |
| a | original code of the research repository at `d065677`: the four scripts, the supplementary diagnostics one process per stage | authors' environment, audit hook blocking writes to the project data | 17,157 of 17,157 identical |
| b | mirrored code (`analysis/`) before the translation of the comments, same stages | authors' environment, audit hook blocking writes to the project data | 17,157 of 17,157 identical |
| c | committed repository (`scripts/10_units.py` … `13_secondary.py`) | the fresh clone and environment of section 1, audit hook blocking writes to the project data | 17,157 of 17,157 identical; 40 min (units 8 s, measurements 7 min, estimates 26 min, diagnostics 7 min) |

The re-runs also write three side files that the archive does not hold (`1p5m/adjacency/log_e2.txt` and two
diagnostic plots of the red-edge script, `1p5m/figs/e2_*.png`).

Two observations from these runs:

* **One process per stage.** Running the supplementary diagnostics in one process
  (`v12_secondary.py tables stats … forward`, as the notes of the research repository suggest) changes the
  confidence bounds of some diagnostics, because several of them resample with one module-level random generator
  (`analysis/eco/chrono/j_review.py`) whose state then depends on the preceding stages. The archived results were
  produced stage by stage; `scripts/13_secondary.py` does the same.
* **Writes.** Under the audit hook, the analysis wrote only below its output folder (`<output>/analysis/v12`); in
  run c (53 processes) no write into the project data was attempted.

## 5. Export of the de-identified data

`tools/export_derived.py`, run in the clone on the project data (after the figure run above, which writes the two
summaries it copies), wrote 72 files (27.5 MB) that are byte-identical to `data/derived`, including
`MANIFEST.json`: the export is deterministic, and the repository data are exactly what the tool derives from the
project data. The pseudonym maps covered 13,622 cells, 150 stations and 149 control blocks.

## 6. Not verified here

* The upstream steps of `docs/PIPELINE.md` (section A: calibration, row-band mapping, Landsat time series and event
  study, station polygons, imaging kernels) were not re-run; the analysis reads their outputs as frozen files of the
  project data (it snapshots the inputs of the coverage expansion with their SHA-256,
  `<output>/analysis/v12/inputs/manifest.json`).
* Pixel identity was checked on Windows with the pinned packages only. Other operating systems or library versions
  can move anti-aliased pixels; the plotted-value tables are the platform-independent check.
* Tables S2–S4, S6 and S7 are not recomputed (README.md, "What reproduces what").
