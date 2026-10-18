"""J: additions for the review: bootstrap clustered by station × detection year, total uncertainty (with the scene-level systematic error), full tables of same-age units of different batches, composition of the >=3 m subset,
label year (BJ3A1 rule), configuration relations, main tests and multiple comparisons, full table of the one-year change of old arrays. Uses only existing per-unit results (a1_*_cells.csv, e_yoy_*_cells.csv).

Clustered bootstrap: treated units are resampled with replacement by "cluster", cluster = station × detection year (units without a station label by 2 km block, the grouping of s6_event.py);
the control mean and the pre-construction baseline of each treated unit are drawn with it (controls "linked"; control blocks are not resampled separately).
The 2 km joint block bootstrap intervals of the original a1/e_yoy (one set of block multiplicities for treated units and controls) are given alongside. Groups with <6 clusters report point estimates only.
Total uncertainty: bootstrap SE and the scene-level calibration systematic error (1σ of the difference of two group means in one scene: NDVI 0.007; converted from the red error: NIRv 0.002, MSAVI 0.0043, red 0.0029, NIR 0.003)
combined in quadrature, ±1.96σ; relative quantities add a systematic error of the same order for the denominator (control mean).
Output: data/eco/chrono/j_*.csv
"""
# Mirrored from eco/chrono/j_review.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import warnings
import numpy as np
import cc

warnings.simplefilter('ignore')
O = cc.out_dir()
es = cc.es
B = 2000
SYS = {'ndvi': 0.007, 'nirv': 0.002, 'msavi': 0.0043, 'red': 0.0029, 'nir': 0.003}
PER = [('2025-06', '06', ''), ('2025-08 (A1)', '08', '_BJ3A1'), ('2026-08', '2026_08', '')]
KS = ['ndvi', 'nirv', 'msavi', 'red', 'nir']
SUMG = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('≥10', 10, 99)]
AGRP = [('0', 0, 0), ('1', 1, 1), ('2', 2, 2), ('3', 3, 3), ('4–5', 4, 5), ('6–7', 6, 7), ('8–9', 8, 9), ('10–13', 10, 13)]
MIN_CL = 6
rng = np.random.default_rng(20261006)


def cboot(x, clus, masks, den=None):
    """Clustered bootstrap. x (n,) unit values; masks {name: bool}; den (n,) optional (denominator of relative quantities: ratio of the equally weighted means of one group).
    Returns {name: dict(cell, cell_lo, cell_hi, cell_se, stn, stn_lo, stn_hi, stn_se, n, ncl, rel...)}."""
    f = np.isfinite(x) & (np.isfinite(den) if den is not None else True)
    uc, inv = np.unique(clus[f], return_inverse=True)
    xv = x[f]
    dv = den[f] if den is not None else None
    names = list(masks)
    M = np.stack([masks[k][f] for k in names], 1).astype(float)          # (n, m)
    nc = len(uc)
    S = np.zeros((nc, len(names)))
    N = np.zeros((nc, len(names)))
    np.add.at(S, inv, M * xv[:, None])
    np.add.at(N, inv, M)
    if dv is not None:
        SD = np.zeros((nc, len(names)))
        np.add.at(SD, inv, M * dv[:, None])
    has = N > 0
    Cm = np.where(has, S / np.where(has, N, 1), 0)
    K = np.stack([np.bincount(rng.integers(0, nc, nc), minlength=nc) for _ in range(B)]).astype(float)   # (B, nc)
    with np.errstate(invalid='ignore', divide='ignore'):
        cell_b = (K @ S) / (K @ N)
        stn_b = (K @ Cm) / (K @ has)
        cell_p = S.sum(0) / N.sum(0)
        stn_p = Cm.sum(0) / has.sum(0)
        if dv is not None:
            rel_b = (K @ S) / (K @ SD)
            rel_p = S.sum(0) / SD.sum(0)
            dmean = SD.sum(0) / N.sum(0)
    out = {}
    for j, k in enumerate(names):
        ncl = int(has[:, j].sum())
        r = dict(n=int(N[:, j].sum()), ncl=ncl, cell=cell_p[j], stn=stn_p[j])
        for nm, bb in (('cell', cell_b[:, j]), ('stn', stn_b[:, j])):
            g = bb[np.isfinite(bb)]
            ok = ncl >= MIN_CL and g.size > 50
            r[f'{nm}_se'] = g.std() if ok else np.nan
            r[f'{nm}_lo'], r[f'{nm}_hi'] = (np.percentile(g, [2.5, 97.5]) if ok else (np.nan, np.nan))
        if dv is not None:
            g = rel_b[:, j][np.isfinite(rel_b[:, j])]
            r['rel'] = rel_p[j]
            r['rel_se'] = g.std() if (ncl >= MIN_CL and g.size > 50) else np.nan
            r['den'] = dmean[j]
        out[k] = r
    return out
