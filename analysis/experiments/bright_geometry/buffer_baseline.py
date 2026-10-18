"""Unmatched buffer baseline for the primary inter-row contrasts (v1.1).

The primary analysis compares each treated 150 m unit with its original
same-scene matched controls. This script replaces the matched controls by an
unselected local buffer: every 150 m cell within R km of the treated unit whose
ground at 300-1,000 m from the mapped arrays of the same acquisition contains at
least 2,000 native 2 m pixels from the unit's frozen reference scene, and which
passes the land-cover quality rules of the control pool (R1, R3-R6 of
eco/ts/s5_pool.py: valid Landsat record, no station polygon, no Landsat
darkening step, WorldCover cropland/built-up/water/wetland/tree limits, no
seasonal water). No covariate matching or common-support rule (R7, R8) is used.

For unit i and index k
    Ic_buf(i)   = equal-cell mean of k over the buffer cells (reflectances of the
                  buffer-band pixels averaged within each cell before the index)
    Dpre_buf(i) = mean over Landsat midsummer years with e <= -2 (>= 3 years) of
                  [Y_k(i) - equal-cell mean of Y_k over the same buffer cells]
    dadj_buf(i) = I_ground(i) - Ic_buf(i) - Dpre_buf(i)
so the estimand differs from the primary one only in the comparison ground.

Uncertainty: the 600 joint resamples of the primary analysis with the same seed
and draw sequence (whole-station treated clusters; independent 2 km blocks of
comparison cells), so buffer and matched replicates are paired; the frozen
open-ground calibration SE is added in quadrature. Sensitivity variants (radius
3 and 8 km; ground beyond 1 km) use 2,000 whole-station resamples with the
comparison ground fixed, the convention of the paper's secondary diagnostics.

Outputs: data/eco/releases/fixed_ew_v2_orientation/analysis/buffer/
  effects_buffer.csv    group estimates, joint and joint+calibration intervals,
                        matched estimates on the same units and the paired
                        buffer - matched difference
  sensitivity.csv       variant estimates with conditional station intervals
  cells_<tag>.csv       unit-level matched and buffer contrasts
  complete.json         audit
Run: .venv/Scripts/python.exe experiments/bright_geometry/buffer_baseline.py [--B 600]
"""
# Mirrored from experiments/bright_geometry/buffer_baseline.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for _k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(_k, '4')
import sys
import warnings
from pathlib import Path
import numpy as np

warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_statistics as N

E = N.E
cc = E.cc
es = cc.es
J = E.J
TAGS = ['06', '08', '2026_08']
METRICS = ['ndvi', 'nirv', 'msavi', 'red', 'nir']
GROUPS = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('≥10', 10, 99), ('all', 0, 99)]
OUT = N.OUT / 'buffer'
SEED = 20261010           # the primary joint-resampling seed (orientation_statistics.analyze)
MIN_PIX = 2000            # the matched-control pixel minimum
MIN_CELLS = 5
VARIANTS = {               # name: (band, radius m)
    'ring300_1000_R5km': ('o300_1000', 5000),
    'ring300_1000_R3km': ('o300_1000', 3000),
    'ring300_1000_R8km': ('o300_1000', 8000),
    'beyond1000_R8km': ('ogt1000', 8000),
}
MAIN = 'ring300_1000_R5km'


def quality(cells):
    """Land-cover quality rules R1, R3-R6 of the control pool, without R2, R7 and R8."""
    Y = np.load(os.path.join(cc.TS, 'cell_y.npz'))
    nvalid = np.isfinite(Y['sum_ndvi']).sum(1)
    c = cells
    return (c.in_aoi.values & (nvalid >= 20) & (c.station_any.values == 0) & (c.strict_n.values == 0)
            & (c.wc40.values < 0.05) & (c.wc50.values < 0.05) & (c.wc80.values < 0.01)
            & (c.wc90.values < 0.05) & (c.wc10.values < 0.1) & (c.water_frac_max.values < 0.05))


class Buffer:
    """Vectorised buffer means of the 2 m index and the Landsat pre-construction difference."""

    def __init__(self, d, W, band):
        self.d, self.W = d, W
        raw = d['tabs']['raw']
        I = cc.index_from_sums(raw[f'{band}_sr'].values, raw[f'{band}_sn'].values, raw[f'{band}_n'].values)
        self.I = {k: I[k] for k in METRICS}
        self.Y = {k: d['YL'][k] for k in METRICS}
        self.e = es.YEARS[None, :] - d['T'][:, None]
        self.pre = self.e <= -2
        self.yi = d['year'] - es.YEARS[0]
        self.n = np.asarray(W.sum(1)).ravel()
        cols = np.unique(W.indices)
        self.cols = cols
        self.Wc = W[:, cols]
        # stacked Landsat matrix of the buffer columns: (ncols, 27 * metrics)
        self.Ystack = np.hstack([np.nan_to_num(self.Y[k][cols]) for k in METRICS])
        self.Fstack = np.hstack([np.isfinite(self.Y[k][cols]).astype(float) for k in METRICS])
        self.Istack = np.stack([np.nan_to_num(self.I[k][cols]) for k in METRICS], 1)
        self.FI = np.stack([np.isfinite(self.I[k][cols]).astype(float) for k in METRICS], 1)

    def compute(self, wp=None):
        w = np.ones(len(self.cols)) if wp is None else wp[self.cols]
        num = self.Wc @ (w[:, None] * self.Istack * self.FI)
        den = self.Wc @ (w[:, None] * self.FI)
        with np.errstate(invalid='ignore', divide='ignore'):
            Ic = np.where(den > 0, num / den, np.nan)
        ynum = self.Wc @ (w[:, None] * self.Ystack * self.Fstack)
        yden = self.Wc @ (w[:, None] * self.Fstack)
        with np.errstate(invalid='ignore', divide='ignore'):
            Yc = np.where(yden > 0, ynum / yden, np.nan)
        out = {}
        t = self.d['tidx']
        for j, k in enumerate(METRICS):
            Yt = self.Y[k][t]
            D = Yt - Yc[:, 27 * j:27 * (j + 1)]
            bm = self.pre & np.isfinite(D)
            nb = bm.sum(1)
            b = np.where(bm, D, 0).sum(1) / np.maximum(nb, 1)
            b[nb < 3] = np.nan
            rawy = D[:, self.yi] - b if 0 <= self.yi < D.shape[1] else np.full(len(t), np.nan)
            out[k] = dict(ic=Ic[:, j], dpre=b, raw=rawy)
        return out


def wmean(x, w):
    ok = np.isfinite(x) & (w > 0)
    return float(np.sum(x[ok] * w[ok]) / w[ok].sum()) if ok.any() else np.nan


def station_boot(x, clusters, valid, masks, B=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    u, inv = np.unique(clusters, return_inverse=True)
    res = {g: [] for g in masks}
    for _ in range(B):
        wt = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u)).astype(float)[inv]
        for g, m in masks.items():
            res[g].append(wmean(x, wt * valid * m))
    return {g: np.asarray(v) for g, v in res.items()}
