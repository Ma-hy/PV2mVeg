"""Shared configuration of eco/rowband: image row-band rule vs module-projection rule.

Acquisitions (tag): '06' = 2025-06-16 (BJ3A1), '08' = 2025-08 (the analysis uses BJ3A1 scenes only; A4/N3 scenes only for the park-wide module mask),
'2026_08' = 2026-08-04 (BJ3A1, two orbits). Each acquisition is processed on its native grid: EPSG:4543 for 2025, EPSG:4542 for 2026.
Every script starts with `import rb_common as R` (this module adds src/ to the path and imports talatan). Existing products are only read; new files go to
data/eco/rowband/, data/eco/grid2m/<tag>/component/ (file names with _comp) and data/eco/manual_check/.
"""
from __future__ import annotations
# Mirrored from eco/rowband/rb_common.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401

import glob
import json
import math
import os
import re

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
import pv2mveg  # noqa: F401,E402  (PROJ/GDAL data of the rasterio wheel, no WMI)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

FONT = FONT_FILE
OUT = os.path.join(ROOT, 'data', 'eco', 'rowband')
GRID = os.path.join(ROOT, 'data', 'eco', 'grid')
GRID2M = os.path.join(ROOT, 'data', 'eco', 'grid2m')
BJ = os.path.join(ROOT, 'data', 'eco', 'bj3ndvi')
MANUAL = os.path.join(ROOT, 'data', 'eco', 'manual_check')
os.makedirs(OUT, exist_ok=True)

PER = {
    '06': dict(label='2025-06', year=2025, epsg=4543, raw='06',
               mask='outputs/06/talatan_pv_mask_0p5m_06.tif', source='outputs/06/talatan_pv_source_0p5m_06.tif',
               summary='outputs/06/summary_06.json', arrays='outputs/06/talatan_pv_arrays_06.gpkg',
               config=GRID, pf10='panel_frac_06_10m.tif', filled=None),
    '08': dict(label='2025-08', year=2025, epsg=4543, raw='08',
               mask='outputs/08/talatan_pv_mask_0p5m_08.tif', source='outputs/08/talatan_pv_source_0p5m_08.tif',
               summary='outputs/08/summary_08.json', arrays='outputs/08/talatan_pv_arrays_08.gpkg',
               config=GRID, pf10='panel_frac_08_10m.tif', filled='outputs/08/talatan_pv_filled_08.gpkg'),
    '2026_08': dict(label='2026-08', year=2026, epsg=4542, raw='2026_08',
                    mask='data/eco/mask2026/2026_08/talatan_pv_mask_0p5m_2026_08.tif',
                    source='data/eco/mask2026/2026_08/talatan_pv_source_0p5m_2026_08.tif',
                    summary='data/eco/mask2026/2026_08/summary_2026_08.json',
                    arrays='data/eco/mask2026/2026_08/talatan_pv_arrays_2026_08.gpkg',
                    config=os.path.join(BJ, '2026_08', 'config'), pf10=None, filled=None,
                    cloud='data/eco/mask2026/2026_08/talatan_pv_cloud_2m_2026_08.tif'),
}
for _t, _c in PER.items():
    for _k in ('mask', 'source', 'summary', 'arrays', 'filled', 'cloud'):
        if _c.get(_k):
            _c[_k] = os.path.join(ROOT, _c[_k])
TAGS = ['06', '08', '2026_08']

# construction batches (build_year 10 m layer) and row-pitch classes
BATCHES = [('2013–15', 2013, 2015), ('2016–19', 2016, 2019), ('2020–22', 2020, 2022), ('2023', 2023, 2023),
           ('2024', 2024, 2024), ('2025', 2025, 2025)]
PITCH = [('<7 m', 0, 7), ('7–9.5 m', 7, 9.5), ('9.5–11 m', 9.5, 11), ('≥11 m', 11, 99)]


# ------------------------------------------------------------------ scenes
def _key(path):
    name = os.path.basename(path).replace('_ORTHO_MS', '').replace('_ORTHO_PAN', '')
    return re.sub(r'_(MS|PAN)(\d?)_', r'_X\2_', name)[:-4]


def scenes(tag):
    """Scene table of the acquisition (order = pixel values of the source map), with panchromatic/multispectral paths, sensor, orbit, sun and satellite geometry."""
    keys = json.load(open(PER[tag]['summary'], encoding='utf8'))['scenes']
    files = glob.glob(os.path.join(ROOT, 'data', 'raw', PER[tag]['raw'], '*.tif'))
    pans = {_key(f): f for f in files if re.search(r'_PAN\d?_', os.path.basename(f))}
    mss = {_key(f): f for f in files if re.search(r'_MS\d?_', os.path.basename(f))}
    geo = _geometry(tag)
    out = []
    for i, k in enumerate(keys):
        sensor = k.split('_')[0]
        cam = k.split('_')[1]
        m = re.search(r'_(\d{14})_L4_([0-9A-Z]+)_SC_(\d+)', k)
        g = geo.get(k, {})
        out.append(dict(idx=i, key=k, pan=pans.get(k), ms=mss.get(k), sensor=sensor, camera=cam, time=m.group(1),
                        strip=m.group(2), sc=int(m.group(3)), short=f'{sensor}-{cam}-{m.group(2)}-SC{int(m.group(3))}', **g))
    return out


def sun_position(y, mo, d, h, mi, s, lat=36.055, lon=100.402):
    """Michalsky (1988) approximation (as a06_sun_position.py of the earlier paper sources). Returns (elevation, azimuth) in degrees."""
    if mo <= 2:
        y -= 1; mo += 12
    a = y // 100; b = 2 - a + a // 4
    jd = int(365.25 * (y + 4716)) + int(30.6001 * (mo + 1)) + d + b - 1524.5 + (h + mi / 60 + s / 3600) / 24
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360
    g = math.radians((357.528 + 0.9856003 * n) % 360)
    lam = math.radians(L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g))
    eps = math.radians(23.439 - 0.0000004 * n)
    ra = math.degrees(math.atan2(math.cos(eps) * math.sin(lam), math.cos(lam))) % 360
    dec = math.asin(math.sin(eps) * math.sin(lam))
    ut = h + mi / 60 + s / 3600
    gmst = (6.697375 + 0.0657098242 * n + ut) % 24
    ha = math.radians(((gmst + lon / 15) * 15 - ra + 180) % 360 - 180)
    phi = math.radians(lat)
    el = math.asin(math.sin(phi) * math.sin(dec) + math.cos(phi) * math.cos(dec) * math.cos(ha))
    az = math.atan2(-math.sin(ha), math.tan(dec) * math.cos(phi) - math.sin(phi) * math.cos(ha))
    return math.degrees(el), math.degrees(az) % 360


def _geometry(tag):
    """{scene key: dict(sun_el, sun_az, sat_el, sat_az, geo_src)}. No metadata for June: the sun position is computed from the time in the scene name; the satellite geometry is unknown (NaN)."""
    out = {}
    if tag == '2026_08':
        g = pd.read_csv(os.path.join(ROOT, 'data', 'eco', 'mask2026', 'check', 'scene_geometry_2026.csv'))
        for _, r in g.iterrows():
            out[r.key] = dict(sun_el=r.sun_elevation_center, sun_az=r.sun_azimuth_center, sat_el=r.satellite_elevation_center,
                              sat_az=r.satellite_azimuth_center, geo_src='L1 元数据（景中心）')
        return out
    g = pd.read_csv(os.path.join(ROOT, 'ancillary', 'scene_geometry.csv'))
    for _, r in g[g.month.astype(str).str.zfill(2) == tag].iterrows():
        if np.isfinite(r.sun_elevation):
            out[r.key] = dict(sun_el=r.sun_elevation, sun_az=r.sun_azimuth, sat_el=r.sat_elevation, sat_az=r.sat_azimuth,
                              geo_src='L1 元数据（左上角）')
        else:
            t = re.search(r'_(\d{14})_L4_', r.key).group(1)
            el, az = sun_position(int(t[:4]), int(t[4:6]), int(t[6:8]), int(t[8:10]), int(t[10:12]), int(t[12:14]))
            out[r.key] = dict(sun_el=el, sun_az=az, sat_el=np.nan, sat_az=np.nan, geo_src='太阳按时刻计算；无卫星几何')
    return out


# ------------------------------------------------------------------ coordinates and sampling of the 10 m layers
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


G10 = dict(x0=597000.0, y0=4014000.0, res=10.0, w=5802, h=4698)
G30 = dict(x0=596985.0, y0=4014015.0, res=30.0, w=1935, h=1567)


def read10(name, cfg_dir=GRID):
    import rasterio
    with rasterio.open(os.path.join(cfg_dir, f'{name}_10m.tif')) as d:
        return d.read(1)


def sample10(arr, ux, uy, fill=0):
    c = np.floor((np.asarray(ux) - G10['x0']) / 10).astype(np.int64)
    r = np.floor((G10['y0'] - np.asarray(uy)) / 10).astype(np.int64)
    ok = (c >= 0) & (c < G10['w']) & (r >= 0) & (r < G10['h'])
    out = np.full(np.shape(ux), fill, dtype=arr.dtype if np.issubdtype(arr.dtype, np.floating) or fill == 0 else np.float64)
    out[ok] = arr[r[ok], c[ok]]
    return out


def utm_grid_coords(tag, x0, y0, h, w, res=2.0, k=32):
    """EPSG:32647 coordinates of the cell centres (upper left x0,y0, pixel size res): exact projection every k cells, bilinear interpolation in between (as eco/bj3ndvi/product._utm_coords)."""
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
    """EPSG:32647 coordinates -> number of the 150 m unit of eco/ts (5×5 cells of the Landsat 30 m grid, from the upper left); -1 outside the grid."""
    CS, CH, CW = 5, 1567 // 5, 1935 // 5
    c30 = np.floor((ux - G30['x0']) / 30).astype(np.int64)
    r30 = np.floor((G30['y0'] - uy) / 30).astype(np.int64)
    ok = (c30 >= 0) & (r30 >= 0) & (c30 < CW * CS) & (r30 < CH * CS)
    return np.where(ok, (r30 // CS) * CW + c30 // CS, -1)
