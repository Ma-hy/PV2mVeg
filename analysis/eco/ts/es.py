"""Core functions of the event study: matching, estimators (difference method / regression counterfactual / Callaway–Sant'Anna staggered DID), aggregation, block bootstrap.

Notation: Y (n_cell, 27) unit outcomes 2000–2026; T event year of the unit (e = year − T; e=−1, 0 construction period; e<=−2 pre-construction).
"""
# Mirrored from eco/ts/es.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import numpy as np
import common as C

YEARS = np.arange(2000, 2027)
NYR = len(YEARS)
EV = np.arange(-12, 13)
BLOCK = 2000.0


# ---------------------------------------------------------------- tools
def nan_slope(x, y):
    """Row-wise NaN-aware OLS slope. x (m,) or (n,m); y (n,m)."""
    x = np.broadcast_to(x, y.shape)
    v = np.isfinite(y) & np.isfinite(x)
    n = v.sum(1)
    xm = np.where(v, x, 0).sum(1) / np.maximum(n, 1)
    ym = np.where(v, y, 0).sum(1) / np.maximum(n, 1)
    dx = np.where(v, x - xm[:, None], 0)
    dy = np.where(v, y - ym[:, None], 0)
    sxx = (dx * dx).sum(1)
    with np.errstate(invalid='ignore', divide='ignore'):
        b = (dx * dy).sum(1) / sxx
    b[(n < 3) | (sxx <= 0)] = np.nan
    return b, ym - b * xm, n


def block_id(x, y):
    bx = np.floor((np.asarray(x) - C.X0) / BLOCK).astype(np.int64)
    by = np.floor((C.Y1 - np.asarray(y)) / BLOCK).astype(np.int64)
    return by * 1000 + bx


def to_event(D, T, ev=EV):
    """D (n,27) by year -> (n, len(ev)) by event time."""
    n = D.shape[0]
    out = np.full((n, len(ev)), np.nan, np.float64)
    for j, e in enumerate(ev):
        yi = (T + e) - YEARS[0]
        ok = (yi >= 0) & (yi < NYR)
        out[ok, j] = D[np.flatnonzero(ok), yi[ok]]
    return out


def year_mask(T, cond):
    """(n,27) boolean: cond(e) at the event time of each cell."""
    e = YEARS[None, :] - np.asarray(T)[:, None]
    return cond(e)


# ---------------------------------------------------------------- matching
FEAT_NAMES = ['ndvi_mean', 'ndvi_trend', 'ndvi_precip_sens', 'red_mean', 'elev', 'slope']


def features(Yo, red, cells, idx, T, pz, offset=6):
    """Features over the matching window 2000..T−offset. Returns (len(idx), 6) and the number of valid years in the window."""
    wy = YEARS <= T - offset
    y = Yo[idx][:, wy]
    r = red[idx][:, wy]
    nv = np.isfinite(y).sum(1)
    f = np.full((len(idx), 6), np.nan)
    f[:, 0] = np.nanmean(y, 1)
    f[:, 1] = nan_slope(YEARS[wy].astype(float), y)[0]
    f[:, 2] = nan_slope(pz[wy], y)[0]
    f[:, 3] = np.nanmean(r, 1)
    f[:, 4] = cells.elev.values[idx]
    f[:, 5] = cells.slope.values[idx]
    return f, nv


def est_diff(Yt, Yc, T, base=lambda e: e <= -2):
    D = Yt - Yc
    bm = year_mask(T, base) & np.isfinite(D)
    b = np.where(bm, D, 0).sum(1) / np.maximum(bm.sum(1), 1)
    b[bm.sum(1) < 3] = np.nan
    return D - b[:, None]


def est_reg(Yt, Yc, T, fit=lambda e: e <= -2):
    """Regression counterfactual: Yt = a + b·Yc fitted over the fit years, effect = Yt − (a + b·Yc)."""
    fm = year_mask(T, fit)
    yy = np.where(fm, Yt, np.nan)
    xx = np.where(fm, Yc, np.nan)
    b, a, n = nan_slope(xx, yy)
    b[n < 5] = np.nan
    return Yt - (a[:, None] + b[:, None] * Yc), b


def holdout(Yt, Yc, T, offset=6):
    """Hold-out differences on T−5..T−2, (n,4) ×2, with only the years <=T−offset as baseline (difference method) and for the fit (regression method)."""
    D = est_diff(Yt, Yc, T, base=lambda e: e <= -offset)
    R, _ = est_reg(Yt, Yc, T, fit=lambda e: e <= -offset)
    ev = np.arange(-5, -1)
    return to_event(D, T, ev), to_event(R, T, ev)


def multiplicities(cluster_ids, rng):
    """Resample clusters with replacement; returns the number of draws of each cluster (vectorised mapping in dict form)."""
    u, inv = np.unique(cluster_ids, return_inverse=True)
    draw = rng.integers(0, len(u), len(u))
    m = np.bincount(draw, minlength=len(u)).astype(float)
    return m[inv]
