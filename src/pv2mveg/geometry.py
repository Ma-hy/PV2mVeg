"""Grids, projections and Beijing-3 scene registry of the project data (restricted-data figures, analysis).

Acquisitions are processed on their native grids: 2025-06 and 2025-08 in EPSG:4543, 2026-08 in EPSG:4542 (CGCS2000
3-degree Gauss-Kruger zones). Analysis units are cells of 150 m (5 x 5 Landsat pixels of 30 m) of the Landsat grid in
EPSG:32647, numbered row-major from the upper-left corner (ts_cell_index). The 2 m computation grid of each
acquisition is given by grid_of; 2 m pixel centres are projected to EPSG:32647 with exact knots every 32 native pixels
and bilinear interpolation in between (utm_grid_coords), with the knots of the product tiles preserved (coords).

All paths are relative to the project root (config: P.project) and need the project layout.
"""
import glob
import json
import os
import re
from functools import lru_cache

import numpy as np

from pv2mveg.config import P

PER = {
    '06': dict(label='2025-06', year=2025, epsg=4543, raw='06', mask='outputs/06/talatan_pv_mask_0p5m_06.tif',
               source='outputs/06/talatan_pv_source_0p5m_06.tif', summary='outputs/06/summary_06.json'),
    '08': dict(label='2025-08', year=2025, epsg=4543, raw='08', mask='outputs/08/talatan_pv_mask_0p5m_08.tif',
               source='outputs/08/talatan_pv_source_0p5m_08.tif', summary='outputs/08/summary_08.json'),
    '2026_08': dict(label='2026-08', year=2026, epsg=4542, raw='2026_08', mask='data/eco/mask2026/2026_08/talatan_pv_mask_0p5m_2026_08.tif',
                    source='data/eco/mask2026/2026_08/talatan_pv_source_0p5m_2026_08.tif',
                    summary='data/eco/mask2026/2026_08/summary_2026_08.json'),
}

G30 = dict(x0=596985.0, y0=4014015.0, res=30.0, w=1935, h=1567)   # Landsat grid (EPSG:32647)

# inter-row footprint layers of the row-band release (data/eco/releases/fixed_ew_v2_orientation)
CLEARANCES = ['1m', '1p5m', '2m']
DIRS = {'1m': 'ecology_1m_v3', '1p5m': 'ecology_1p5_v3', '2m': 'ecology_exact_v2'}
LAYERS = dict(zip(CLEARANCES, ['pure_interrow_2m', 'clearance_1p5m_2m', 'core_2m_2m']))


def path(tag, key):
    P.need_project('Beijing-3 products')
    return str(P.project / PER[tag][key])


def _key(p):
    name = os.path.basename(p).replace('_ORTHO_MS', '').replace('_ORTHO_PAN', '')
    return re.sub(r'_(MS|PAN)(\d?)_', r'_X\2_', name)[:-4]


def scenes(tag):
    """Scene list of an acquisition in the order of the scene index of the 0.5 m source raster, with the paths of the
    panchromatic and multispectral Level-4 files (data/raw/<tag>)."""
    keys = json.load(open(path(tag, 'summary'), encoding='utf8'))['scenes']
    files = glob.glob(os.path.join(str(P.project), 'data', 'raw', PER[tag]['raw'], '*.tif'))
    pans = {_key(f): f for f in files if re.search(r'_PAN\d?_', os.path.basename(f))}
    mss = {_key(f): f for f in files if re.search(r'_MS\d?_', os.path.basename(f))}
    return [dict(idx=i, key=k, pan=pans.get(k), ms=mss.get(k)) for i, k in enumerate(keys)]


_TR = {}


def to_utm(tag, x, y):
    from pyproj import Transformer
    if tag not in _TR:
        _TR[tag] = Transformer.from_crs(PER[tag]['epsg'], 32647, always_xy=True)
    return _TR[tag].transform(x, y)


def from_utm(tag, x, y):
    from pyproj import Transformer
    k = tag + '_inv'
    if k not in _TR:
        _TR[k] = Transformer.from_crs(32647, PER[tag]['epsg'], always_xy=True)
    return _TR[k].transform(x, y)


def utm_grid_coords(tag, x0, y0, h, w, res=2.0, k=32):
    """EPSG:32647 coordinates of the cell centres of a native window (upper left x0, y0; pixel size res): exact
    projection every k pixels, bilinear interpolation in between (as the 2 m products)."""
    from scipy.interpolate import RegularGridInterpolator
    rr = np.r_[np.arange(0, h, k), h - 1] if (h - 1) % k else np.arange(0, h, k)
    cc = np.r_[np.arange(0, w, k), w - 1] if (w - 1) % k else np.arange(0, w, k)
    if len(rr) < 2:
        rr = np.array([0, max(h - 1, 1)])
    if len(cc) < 2:
        cc = np.array([0, max(w - 1, 1)])
    Rg, Cg = np.meshgrid(rr, cc, indexing='ij')
    X, Y = to_utm(tag, x0 + (Cg + 0.5) * res, y0 - (Rg + 0.5) * res)
    pr, pc = np.meshgrid(np.arange(h), np.arange(w), indexing='ij')
    q = np.stack([pr.ravel(), pc.ravel()], 1)
    ux = RegularGridInterpolator((rr, cc), X, bounds_error=False, fill_value=None)(q).reshape(h, w)
    uy = RegularGridInterpolator((rr, cc), Y, bounds_error=False, fill_value=None)(q).reshape(h, w)
    return ux, uy


def ts_cell_index(ux, uy):
    """EPSG:32647 coordinates -> 150 m unit number (5 x 5 Landsat pixels, row-major from the upper left); -1 outside."""
    CS, CH, CW = 5, 1567 // 5, 1935 // 5
    c30 = np.floor((ux - G30['x0']) / 30).astype(np.int64)
    r30 = np.floor((G30['y0'] - uy) / 30).astype(np.int64)
    ok = (c30 >= 0) & (r30 >= 0) & (c30 < CW * CS) & (r30 < CH * CS)
    return np.where(ok, (r30 // CS) * CW + c30 // CS, -1)


@lru_cache(maxsize=3)
def grid_of(tag):
    """2 m computation grid of an acquisition (upper-left corner and size)."""
    if tag in ('06', '08'):
        return dict(x0=329890.0, y0=4013270.0, width=26146, height=20420)
    import rasterio
    with rasterio.open(path(tag, 'mask')) as d:
        t = d.transform
        return dict(x0=t.c, y0=t.f, width=int(np.ceil(d.width / 4)), height=int(np.ceil(d.height / 4)))


def canonical_centres(tag, gr, r0, c0, h, w):
    """The product interpolates projection knots every 32 native pixels, with a shortened last interval at its
    1024-pixel tile edge; smaller tiles keep those knots (ending at pixel 255 instead of knot 256 can shift a centre
    by <1 mm across a UTM aggregation edge)."""
    rh = h + int(r0 + h < gr['height'] and (r0 + h) % 1024 != 0)
    cw = w + int(c0 + w < gr['width'] and (c0 + w) % 1024 != 0)
    ux, uy = utm_grid_coords(tag, gr['x0'] + c0 * 2, gr['y0'] - r0 * 2, rh, cw)
    return ux[:h, :w], uy[:h, :w]


def canonical_tile(tag, r, c):
    g = grid_of(tag); h = min(256, g['height'] - r); w = min(256, g['width'] - c)
    return canonical_centres(tag, g, r, c, h, w)


def coords(tag, r, c, h, w):
    """EPSG:32647 centres of any 2 m window, with the frozen projection knots of the 256-pixel tiles."""
    x = np.zeros((h, w)); y = x.copy()
    for rr in range(r // 256 * 256, (r + h - 1) // 256 * 256 + 1, 256):
        for cc in range(c // 256 * 256, (c + w - 1) // 256 * 256 + 1, 256):
            a, b = canonical_tile(tag, rr, cc)
            rlo = max(r, rr); rhi = min(r + h, rr + a.shape[0]); clo = max(c, cc); chi = min(c + w, cc + a.shape[1])
            dst = np.s_[rlo - r:rhi - r, clo - c:chi - c]; src = np.s_[rlo - rr:rhi - rr, clo - cc:chi - cc]
            x[dst] = a[src]; y[dst] = b[src]
    return x, y


def preview(src, bounds, tf, crs, shape, fill=255, resampling=None):
    """Display-only resampling of a raster window (statistics always use the full-resolution grid)."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.vrt import WarpedVRT
    resampling = Resampling.nearest if resampling is None else resampling
    nodata = 0 if not np.isfinite(fill) else fill
    with rasterio.open(src) as ds:
        with WarpedVRT(ds, crs=crs, transform=tf, width=shape[1], height=shape[0], nodata=nodata, resampling=resampling) as vrt:
            return vrt.read(1, masked=True, out_dtype='float32').filled(fill)


def stretch(a):
    """1-99 % linear stretch to 8 bit (display)."""
    a = np.asarray(a, float)
    vals = a[np.isfinite(a)]
    if not len(vals):
        return np.zeros(a.shape, np.uint8)
    lo, hi = np.percentile(vals, [1, 99])
    return (np.nan_to_num(np.clip((a - lo) / max(hi - lo, 1), 0, 1)) * 255).astype(np.uint8)


def unit_xy(tag, cell):
    """EPSG:32647 centre of a 150 m unit (read from the unit table of the v1.1 analysis, restricted)."""
    import pandas as pd
    t = pd.read_csv(P.analysis / '1p5m' / 'chrono' / f'a1_{tag}_det_thr50_cells.csv', usecols=['cell', 'x', 'y']).set_index('cell')
    return float(t.loc[cell, 'x']), float(t.loc[cell, 'y'])


def native_example(tag, cell):
    """One 150 m unit in native imagery: 0.5 m panchromatic (8-bit stretch), the mapped 0.5 m row band and the 2 m
    footprint selections at the three clearances, on a 192 m window aligned to the 2 m grid. Every footprint count
    of the unit must equal the release sums (the example shows the complete unit)."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import from_origin
    from rasterio.windows import Window
    rel = P.release
    grid = grid_of(tag)
    ux0, uy0 = unit_xy(tag, cell)
    x, y = from_utm(tag, ux0, uy0)
    h = w = 96; r = int((grid['y0'] - y) / 2) - 48; c = int((x - grid['x0']) / 2) - 48
    x0 = grid['x0'] + c * 2; y0 = grid['y0'] - r * 2
    bounds = (x0, y0 - h * 2, x0 + w * 2, y0); tf = from_origin(x0, y0, .5, .5); shape = (h * 4, w * 4)
    epsg = PER[tag]['epsg']
    src = preview(path(tag, 'source'), bounds, tf, epsg, shape, fill=255)
    pan = np.full(shape, np.nan, np.float32)
    sc = scenes(tag)
    for sid in np.unique(src[src != 255]):
        p = preview(sc[int(sid)]['pan'], bounds, tf, epsg, shape, fill=np.nan, resampling=Resampling.bilinear)
        pan[src == sid] = p[src == sid]
    mask = preview(rel / tag / 'rowband_0p5m.tif', bounds, tf, epsg, shape, fill=255)
    ux, uy = coords(tag, r, c, h, w); inside = ts_cell_index(ux, uy) == cell
    selections = {}
    for p in CLEARANCES:
        f = rel / DIRS['1p5m'] / tag / (LAYERS[p] + '.tif')
        with rasterio.open(f) as ds:
            selections[p] = ds.read(1, window=Window(c, r, w, h), boundless=True, fill_value=255) == 1
        sums = np.load(rel / DIRS[p] / tag / 'core_sums.npz')['sums'][2, 2, cell]
        assert int((selections[p] & inside).sum()) == int(sums), (tag, cell, p, 'example must show full cell')
    assert not (selections['2m'] & ~selections['1p5m']).any()
    assert not (selections['1p5m'] & ~selections['1m']).any()
    return stretch(pan), mask, inside, selections
