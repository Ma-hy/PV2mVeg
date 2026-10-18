"""Shared functions of eco/chrono: acquisition configuration, 30 m sums of the 2 m products -> aggregation to the 150 m units of ts, scene offsets, matched controls, pre-construction Landsat baseline, block bootstrap.

Every script starts with `import cc` (this module adds src/ and eco/ts to the path and imports talatan; it reuses eco/ts/common.py and es.py).
Acquisitions are registered in PERIODS; a new acquisition (such as 2026_08) only needs its product path and acquisition year here; readiness is checked by ready().
"""
# Mirrored from eco/chrono/cc.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
import json

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
import pv2mveg  # noqa: F401,E402  (PROJ/GDAL data of the rasterio wheel, no WMI)

sys.path.insert(0, os.path.join(CODE, 'eco', 'ts'))
import common as C  # noqa: E402
import es  # noqa: E402,F401

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BJ = os.path.join(ROOT, 'data', 'eco', 'bj3ndvi')
TS = os.path.join(ROOT, 'data', 'eco', 'ts')
OUT = os.path.join(ROOT, 'data', 'eco', 'chrono')
FIG = os.path.join(OUT, 'figs')
for _d in (OUT, FIG):
    os.makedirs(_d, exist_ok=True)

# acquisitions: product folder tag, acquisition year (array age = year − year of panel appearance), readiness flag file (None = ready when the product exists)
PERIODS = {
    '06': dict(year=2025, label='2025-06', done=None),
    '08': dict(year=2025, label='2025-08', done=None),
    # 2026-08: can be re-run directly once the mask and the 2 m product are ready (data/eco/mask2026/DONE.json, and eco/bj3ndvi run.py product 2026_08 done)
    '2026_08': dict(year=2026, label='2026-08', done=os.path.join(ROOT, 'data', 'eco', 'mask2026', 'DONE.json')),
}
IDX = ['ndvi', 'nirv', 'msavi']
CTRL_BANDS = ['out_300_1000', 'out_gt1000']          # control units: pure ground >300 m outside arrays
FE_BANDS = ['out_100_300', 'out_300_1000', 'out_gt1000']  # pure ground outside arrays used to estimate scene offsets


def prod_path(tag, kind='30m'):
    return os.path.join(BJ, tag, f'bj3ndvi_{tag}_{kind}.tif')


def bands_json(tag):
    with open(os.path.join(BJ, tag, f'bj3ndvi_{tag}_bands.json'), encoding='utf-8') as f:
        return json.load(f)


def read_bands(path, names=None):
    """Read a multi-band 30 m raster; returns {band name: (H,W) float64}."""
    import rasterio
    with rasterio.open(path) as r:
        desc = list(r.descriptions)
        assert r.height == C.H and r.width == C.W and abs(r.transform.c - C.X0) < 1e-6 and abs(r.transform.f - C.Y1) < 1e-6
        want = desc if names is None else names
        out = {}
        for nm in want:
            out[nm] = r.read(desc.index(nm) + 1).astype(np.float64)
    return out


def cell_sum(a):
    """(H,W) -> (ncell,) 5×5 sums (the partition of the 150 m units of ts, from the upper left)."""
    v = C.cell_view(a)
    return v.sum(axis=(-3, -1)).reshape(-1)


def index_from_sums(sr, sn, n):
    """Indices from mean reflectances (averaged first). Returns dict(ndvi, nirv, msavi, red, nir)."""
    with np.errstate(invalid='ignore', divide='ignore'):
        red = np.where(n > 0, sr / np.maximum(n, 1), np.nan)
        nir = np.where(n > 0, sn / np.maximum(n, 1), np.nan)
    return dict(red=red, nir=nir, ndvi=C.ndvi(red, nir), nirv=C.nirv(red, nir), msavi=C.msavi(red, nir))


def add_index(T, k):
    I = index_from_sums(T[f'{k}_sr'].values, T[f'{k}_sn'].values, T[f'{k}_n'].values)
    for nm, v in I.items():
        T[f'{k}_{nm}'] = v
    return T


# ------------------------------------------------------------------ matched controls and Landsat baseline
def matched(timing='det', kind='treated'):
    d = os.path.join(TS, f'es_raw_sum_ndvi_{timing}_{kind}')
    m = np.load(os.path.join(d, 'matched.npz'))
    return m['tidx'], m['T'], m['M']


def landsat_Y():
    Y = np.load(os.path.join(TS, 'cell_y.npz'))
    return {k: Y[f'sum_{k}'].astype(np.float64) for k in IDX + ['red', 'nir']}


def landsat_base(Yk, tidx, T, M, sel, wpool=None, year=2025):
    """Pre-construction baseline difference D_pre of the Landsat midsummer composite (mean over e<=−2, >=3 years) and the raw effect of the given year.
    sel (n,k) boolean: controls used (e.g. the same-scene subset); wpool: bootstrap multiplicities. Returns (D_pre, raw_year_effect)."""
    YM = Yk[np.clip(M, 0, None)]  # (n,k,27)
    w = sel.astype(float) * (M >= 0)
    if wpool is not None:
        w = w * wpool[np.clip(M, 0, None)]
    w3 = w[:, :, None] * np.isfinite(YM)
    sw = w3.sum(1)
    with np.errstate(invalid='ignore', divide='ignore'):
        Yc = np.where(sw > 0, (np.nan_to_num(YM) * w3).sum(1) / np.maximum(sw, 1e-12), np.nan)
    Yt = Yk[tidx]
    D = Yt - Yc
    e = es.YEARS[None, :] - T[:, None]
    bm = (e <= -2) & np.isfinite(D)
    nb = bm.sum(1)
    b = np.where(bm, D, 0).sum(1) / np.maximum(nb, 1)
    b[nb < 3] = np.nan
    yi = year - es.YEARS[0]
    raw = D[:, yi] - b if 0 <= yi < len(es.YEARS) else np.full(len(T), np.nan)
    return b, raw


def ctrl_mean(Ic, M, sel, wpool=None):
    """(Weighted) mean of the control-unit index Ic (ncell,) over the matched-control subset of each treated unit."""
    V = Ic[np.clip(M, 0, None)]
    w = sel.astype(float) * (M >= 0) * np.isfinite(V)
    if wpool is not None:
        w = w * wpool[np.clip(M, 0, None)]
    sw = w.sum(1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(sw > 0, (np.nan_to_num(V) * w).sum(1) / np.maximum(sw, 1e-12), np.nan), sw


AGE_GROUPS = [('0', 0, 0), ('1', 1, 1), ('2', 2, 2), ('3', 3, 3), ('4–5', 4, 5), ('6–7', 6, 7), ('8–9', 8, 9), ('10–13', 10, 13)]
AGE_SUMMARY = [('0–1', 0, 1), ('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('≥10', 10, 99)]


# ------------------------------------------------------------------ analysis scope (section K: east-west fixed arrays only; optional station exclusion list)
# with CHRONO_SCOPE=ew the rules are identical to eco/ts/ew0_cells.py (TS_SCOPE=ew, threshold 0.9):
#   fixed-array 30 m pixel = orient 1 (east-west rows) and orient_purity >= 0.9; east-west share of a unit = fixed-array pixels / pixels with rows (orient>0);
#   treated units are kept with a share >= 0.9; ring units are kept when >= 0.8 of their "nearest reliable array pixels" are fixed-array pixels.
#   stations of the exclusion lists are also removed: data/eco/anomaly/exclude_stations.csv (column station; optional column exclude)
#   and data/eco/ts/ew/exclude_stations.txt (list of the ts side, one station_id per line); the union of both.
# outputs go to data/eco/chrono/ew/. The default 'all' is the original rule.
SCOPE = os.environ.get('CHRONO_SCOPE', 'all')
EXCLUDE_CSV = os.path.join(ROOT, 'data', 'eco', 'anomaly', 'exclude_stations.csv')
EXCLUDE_TXT = os.path.join(ROOT, 'data', 'eco', 'ts', 'ew', 'exclude_stations.txt')
#   and the unit exclusion list data/eco/ts/ew/exclude_cells.txt (one 150 m unit number per line, lines starting with # are comments; written by eco/anomaly/h5_apply_decisions.py).
#   removed = station list ∪ unit list; without a unit list, or with an empty one, the behaviour is unchanged.
EXCLUDE_CELLS = os.path.join(ROOT, 'data', 'eco', 'ts', 'ew', 'exclude_cells.txt')
#   per-acquisition removal: data/eco/ts/ew/exclude_cells_<tag>.txt (same format; e.g. units affected by thin cloud / cloud shadow in 2026-08 are removed from 2026_08 only).
#   scope_filter(cells, tag) also removes these units when tag is given and the file exists and is not empty; tag may be one acquisition or several (e.g. the union of the two acquisitions of a year-to-year difference).
#   without the file, with an empty file, or without tag, the behaviour is unchanged.
EXCLUDE_CELLS_PERIOD = {t: os.path.join(ROOT, 'data', 'eco', 'ts', 'ew', f'exclude_cells_{t}.txt') for t in PERIODS}
EW_SHARE, EW_PURITY, RING_SHARE = 0.9, 0.9, 0.8
#   sensitivity analyses (main results unaffected): CHRONO_EXCLUDE_EXTRA = paths of additional unit exclusion lists (separated by ;, removed from every acquisition),
#   CHRONO_OUT_NAME = name of the output subfolder (default SCOPE; e.g. ew_sens_mask -> data/eco/chrono/ew_sens_mask/). Without both, the behaviour is unchanged.
EXCLUDE_CELLS_EXTRA = [p for p in os.environ.get('CHRONO_EXCLUDE_EXTRA', '').split(';') if p.strip()]
OUT_NAME = os.environ.get('CHRONO_OUT_NAME', '').strip()
_CACHE = {}


def out_dir():
    d = os.path.join(OUT, OUT_NAME) if OUT_NAME else (OUT if SCOPE == 'all' else os.path.join(OUT, SCOPE))
    os.makedirs(d, exist_ok=True)
    return d


def excluded_stations():
    ex = set()
    if os.path.exists(EXCLUDE_CSV):
        d = pd.read_csv(EXCLUDE_CSV)
        if 'station' in d.columns:
            if 'exclude' in d.columns:
                d = d[d['exclude'].astype(str).str.strip().str.lower().isin(['1', 'true', 'yes', 'y', '是'])]
            ex |= set(pd.to_numeric(d['station'], errors='coerce').dropna().astype(int))
    if os.path.exists(EXCLUDE_TXT):
        ex |= {int(x) for x in open(EXCLUDE_TXT, encoding='utf-8').read().split() if x.strip().isdigit()}
    return ex


def excluded_cells(path=None):
    """Unit exclusion list (one cell per line; empty lines and lines starting with # are ignored). Empty set when the file does not exist or is empty."""
    path = path or EXCLUDE_CELLS
    if not os.path.exists(path):
        return set()
    out = set()
    for ln in open(path, encoding='utf-8'):
        s = ln.strip()
        if s and not s.startswith('#') and s.isdigit():
            out.add(int(s))
    return out


def ew_share():
    if 'ew' not in _CACHE:
        o = C.grid('orient')
        pu = np.nan_to_num(C.grid('orient_purity'))
        ew = ((o == 1) & (pu >= EW_PURITY)).astype(float)
        has = (o > 0).astype(float)
        _CACHE['ew'] = (C.cell_pixels(ew).sum(1) / np.maximum(C.cell_pixels(has).sum(1), 1), C.cell_pixels(has).sum(1))
    return _CACHE['ew']


def scope_filter(cells, tag=None):
    """(ncell,) boolean: treated units kept under the current SCOPE (with SCOPE=all only the exclusion lists apply).
    tag: acquisition (key of PERIODS, e.g. '06', '08', '2026_08') or a list / tuple of them; when given, the units of EXCLUDE_CELLS_PERIOD[tag] are also removed (if the file exists and is not empty)."""
    keep = np.ones(len(cells), bool)
    if SCOPE == 'ew':
        sh, nrow = ew_share()
        keep &= (sh >= EW_SHARE) & (nrow > 0)
    ex = excluded_stations()
    if ex:
        keep &= ~(np.isin(cells.station.values, list(ex)) & (cells.station.values > 0))
    exc = excluded_cells()
    msg = [f'主清单 {len(exc)}']
    if tag is not None:
        for t in ([tag] if isinstance(tag, str) else list(tag)):
            p = EXCLUDE_CELLS_PERIOD.get(t)
            if p is not None:
                ep = excluded_cells(p)
                if ep:
                    msg.append(f'按期 {t} {len(ep)}（其中不在主清单 {len(ep - exc)}）')
                exc = exc | ep
    for p in EXCLUDE_CELLS_EXTRA:
        ex2 = excluded_cells(p)
        msg.append(f'额外（CHRONO_EXCLUDE_EXTRA）{os.path.basename(p)} {len(ex2)}（其中新增 {len(ex2 - exc)}）')
        exc = exc | ex2
    if exc:
        cid = cells.cell.values if 'cell' in cells.columns else np.arange(len(cells))
        hit = np.isin(cid, list(exc))
        print(f'[cc.scope_filter] SCOPE={SCOPE} tag={tag}：排除单元清单 ' + '，'.join(msg) + f'；合计 {len(exc)} 个单元，'
              f'本表命中 {int(hit.sum())}，其中原本在口径内被剔除 {int((hit & keep).sum())}', flush=True)
        keep &= ~hit
    return keep
