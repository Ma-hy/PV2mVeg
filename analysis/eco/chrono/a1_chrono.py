"""A: age series (space for time): 2 m inter-row core ground − the same set of matched controls (pure ground >300 m outside arrays), minus the pre-construction Landsat baseline difference.

Usage: python eco/chrono/a1_chrono.py --tag 06 --timing det --thr 50 --B 400
  --tag     acquisition (06, 08, 2026_08)
  --timing  det = detection year (year of panel appearance, detect_year of ts), lab = label year; matched controls from es_raw_sum_ndvi_<timing>_treated/matched.npz of ts
  --thr     minimum number of 2 m inter-row core pixels of a treated unit
Variants:
  same    only the matched controls with the modal scene of the treated unit (>=2), no scene correction (main result)
  all     all matched controls, no correction (reference: includes scene offsets)
  fe_red  all matched controls, scene offsets of the seams subtracted from Σ first (offset assumed to be all in red)
  fe_nir  the same (offset assumed to be all in NIR); both assumptions give almost the same NDVI, NIRv/MSAVI take the two as an interval
Per treated unit: Δ = I(inter-row core) − mean_j I(control j); D_pre = mean over e<=−2 of the Landsat midsummer composite (treated − the same control subset); Δadj = Δ − D_pre;
raw2025 = Landsat 2025 midsummer (treated − the same control subset) − D_pre (the estimator of the ts difference method), bias = raw2025 − Δadj.
Uncertainty: joint bootstrap of 2 km blocks (one set of block multiplicities for treated units and controls, as in ts), B draws.
Output: data/eco/chrono/a1_<tag>_<timing>_thr<thr>_{cells.csv, curves.csv}
"""
# Mirrored from eco/chrono/a1_chrono.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import numpy as np
import cc

C, es = cc.C, cc.es
VARIANTS = ['same', 'all', 'fe_red', 'fe_nir']
CTRL_MIN_N = 2000      # minimum number of 2 m pure-ground pixels >300 m from arrays in a control unit (of 5625 per unit)
SAME_MIN = 2           # minimum number of same-scene controls
PF_BINS = [(0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 1.01)]
CFG_BINS = [('EW 10–11.5 m', 1, 10, 11.5), ('EW 9–11 m', 1, 9, 11), ('EW ≥11 m', 1, 11, 99), ('EW <9 m', 1, 0, 9), ('NS ≥11 m', 2, 11, 99), ('NS <11 m', 2, 0, 11)]
CFG_AGES = [('0–1', 0, 1), ('2–3', 2, 3), ('4–5', 4, 5), ('6–9', 6, 9), ('≥10', 10, 99)]


QTY = ['dadj', 'd', 'dpre', 'raw', 'bias', 'ic', 'it']
KS = cc.IDX + ['red', 'nir']   # indices + single bands (band decomposition, section H)


def per_cell(D, v, k, wpool=None):
    t = D['tabs'][D['tab'][v]]
    It = t[f'ir_{k}'].values[D['tidx']]
    Icells = np.where(t.ctrl_n.values >= CTRL_MIN_N, t[f'ctrl_{k}'].values, np.nan)
    Icm, _ = cc.ctrl_mean(Icells, D['M'], D['sel'][v], wpool)
    Dpre, raw = cc.landsat_base(D['YL'][k], D['tidx'], D['T'], D['M'], D['sel'][v], wpool, D['year'])
    d = It - Icm
    dadj = d - Dpre
    out = dict(dadj=dadj, d=d, dpre=Dpre, raw=raw, bias=raw - dadj, ic=Icm, it=It)
    ok = D['valid'][v] & np.isfinite(dadj)
    return out, ok
