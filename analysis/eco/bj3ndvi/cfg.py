"""Beijing-3 L4 2 m inter-row ground vegetation-index product: acquisition configuration, scene tables, calibration coefficients and shared constants.

Acquisitions (tag): '06' = 2025-06-16, '08' = 2025-08-13/14/17, '2026_08' = 2026-08-04.
Each acquisition needs: the 0.5 m panel mask of the month, the 0.5 m scene-source map (pixel value = index into summary.scenes), the scene order (scenes of summary.json),
the array-area polygons, optional polygons of the "no-image fill" area, the MS scene folder, the 2 m computation grid and the coordinate system (from the mask file).
Calibration coefficients are always read from data/eco/bj3cal/coefficients.json (output of the cal stage of eco/bj3_calibration.py); there is no other rule.
"""
from __future__ import annotations
# Mirrored from eco/bj3ndvi/cfg.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401

import json
import os
import sys

import numpy as np

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
sys.path.insert(0, os.path.join(CODE, 'eco'))

import bj3_calibration as BC  # noqa: E402  reuses its scene_list / _sid_from_key / s2_load / fit_huber / ndvi
from common import G10, G30, GRID_DIR, set_mpl_font  # noqa: E402,F401

OUT = os.path.join(ROOT, 'data', 'eco', 'bj3ndvi')
GRID2M = os.path.join(ROOT, 'data', 'eco', 'grid2m')
BJ3CAL = os.path.join(ROOT, 'data', 'eco', 'bj3cal')
FONT = FONT_FILE

# the 2 m grid shared by the two 2025 acquisitions = the grid of data/eco/grid2m (union extent of the June and August 2 m products, custom EPSG:4543 WKT)
UNION2025 = dict(x0=329890.0, y0=4013270.0, width=26146, height=20420)

PERIODS = {
    '06': dict(label='2025-06', cal_period='2025-06', raw='06',
               mask='outputs/06/talatan_pv_mask_0p5m_06.tif', frac='outputs/06/talatan_pv_fraction_2m_06.tif',
               source='outputs/06/talatan_pv_source_0p5m_06.tif', summary='outputs/06/summary_06.json',
               arrays='outputs/06/talatan_pv_arrays_06.gpkg', filled=None, grid=UNION2025),
    '08': dict(label='2025-08', cal_period='2025-08', raw='08',
               mask='outputs/08/talatan_pv_mask_0p5m_08.tif', frac='outputs/08/talatan_pv_fraction_2m_08.tif',
               source='outputs/08/talatan_pv_source_0p5m_08.tif', summary='outputs/08/summary_08.json',
               arrays='outputs/08/talatan_pv_arrays_08.gpkg', filled='outputs/08/talatan_pv_filled_08.gpkg',
               grid=UNION2025),
    # 2026: the mask was produced by another agent (data/eco/mask2026/DONE.json), on the native EPSG:4542 grid (same grid and scenes as MS, no resampling).
    # clouds/cloud shadows must be removed; the configuration layers (orient/pitch) are rebuilt from the row polygons of this acquisition with the algorithm of eco/build_config.py (own_config, see config_layers.py).
    '2026_08': dict(label='2026-08', cal_period='2026-08', raw='2026_08',
                    mask='data/eco/mask2026/2026_08/talatan_pv_mask_0p5m_2026_08.tif',
                    frac='data/eco/mask2026/2026_08/talatan_pv_fraction_2m_2026_08.tif',
                    source='data/eco/mask2026/2026_08/talatan_pv_source_0p5m_2026_08.tif',
                    summary='data/eco/mask2026/2026_08/summary_2026_08.json',
                    arrays='data/eco/mask2026/2026_08/talatan_pv_arrays_2026_08.gpkg', filled=None,
                    cloud='data/eco/mask2026/2026_08/talatan_pv_cloud_2m_2026_08.tif', own_config=True, cache_mask_borrowed=True, grid=None),
}

# pixel classes (2 m)
CLS_NODATA, CLS_PANEL, CLS_EDGE, CLS_INTERROW, CLS_LANE, CLS_OUT, CLS_OUTNEAR = 0, 1, 2, 3, 4, 5, 6
CLS_NAMES = {0: 'nodata', 1: 'panel', 2: 'edge', 3: 'core_interrow', 4: 'core_lane', 5: 'outside_pure', 6: 'outside_near_panel'}
# quality bits (2 m qa)
QA_ANOM = 1        # this 10 m calibration cell was flagged as anomalous by the "residual >5× robust SD" rule (clouds/cloud shadows/heliostat field/construction etc.)
QA_ANOM_BUF = 2    # 1-cell (10 m) buffer around anomalous cells
QA_UNTESTED = 4    # S2 not clear here or the BJ-3 10 m cell incomplete: no anomaly test possible (kept, not removed)
QA_FILLED = 8      # no-image fill area of August (taken from the June product): removed
QA_SEAM = 16       # the 16 0.5 m sub-pixels of the 2 m cell come from different scenes (scene seam)
QA_NOMS = 32       # multispectral invalid (DN <=1) or no source scene
QA_NEARNOD = 64    # no-image area of the mask within 2 m (the distance is an upper bound): removed
QA_CLOUD = 128     # clouds/cloud shadows (with a 10 m buffer; cloud_2m of the 2026 mask agent), treated as no image: removed
QA_EXCLUDE = QA_ANOM | QA_ANOM_BUF | QA_FILLED | QA_NOMS | QA_NEARNOD | QA_CLOUD

CORE_MIN = 2.0     # core ground: >=2 m from the nearest panel edge
RELAX_MIN = 1.0    # relaxed rule for north-south arrays (comparison only, not in the main product)
DIST_BINS = [(2, 3), (3, 4), (4, 6), (6, 10), (10, 1e9)]
DIST_BIN_NAMES = ['d2_3', 'd3_4', 'd4_6', 'd6_10', 'd10p']
DIRS = {1: 'N', 2: 'E', 3: 'S', 4: 'W'}
OUT_BANDS = [(0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e9)]
OUT_BAND_NAMES = ['out_0_30', 'out_30_100', 'out_100_300', 'out_300_1000', 'out_gt1000']
GAP_TOL = 0.0      # inter-row criterion: clearance across the rows <= local row pitch pitch_m + GAP_TOL


def pcfg(tag):
    c = dict(PERIODS[tag])
    c['tag'] = tag
    for k in ('mask', 'frac', 'source', 'summary', 'arrays', 'filled', 'cloud'):
        if c.get(k):
            c[k] = os.path.join(ROOT, c[k])
    if c['grid'] is None:  # 2 m grid derived from the mask extent (4×4 aggregation of the 0.5 m mask)
        import rasterio
        with rasterio.open(c['mask']) as d:
            t = d.transform
            c['grid'] = dict(x0=t.c, y0=t.f, width=int(np.ceil(d.width / 4)), height=int(np.ceil(d.height / 4)))
    import rasterio
    with rasterio.open(c['mask']) as d:
        c['crs_wkt'] = d.crs.to_wkt()
    c['out'] = os.path.join(OUT, tag)
    c['g2m'] = os.path.join(GRID2M, tag)
    c['work'] = os.path.join(OUT, tag, 'work')
    c['config_dir'] = os.path.join(OUT, tag, 'config') if c.get('own_config') else GRID_DIR
    return c


def scenes(c):
    """Scene table of the acquisition in the order of summary.scenes (= pixel values of the source map), with MS path, sid, calibration coefficients, S2 date and co-registration shift."""
    keys = json.load(open(c['summary'], encoding='utf8'))['scenes']
    lst = {s['sid']: s for s in BC.scene_list() if s['dir'] == c['raw']}
    coefs = json.load(open(os.path.join(BJ3CAL, 'coefficients.json')))
    out = []
    for k, key in enumerate(keys):
        sid = BC._sid_from_key(key)
        if sid not in lst:
            raise KeyError(f'no MS file for {key} -> {sid}')
        cf = coefs[sid]
        if c.get('coef_refit'):   # coefficients refitted with the mask of this acquisition (output of flags.py), see README
            cf = json.load(open(os.path.join(c['out'], 'coefficients_refit.json')))[sid]
        z = np.load(os.path.join(BJ3CAL, 'cache', sid + '.npz'))
        out.append(dict(idx=k, key=key, sid=sid, path=lst[sid]['path'], sensor=lst[sid]['sensor'], date=lst[sid]['date'],
                        s2=cf['s2'], a_red=cf['red'][0], b_red=cf['red'][1], a_nir=cf['nir'][0], b_nir=cf['nir'][1],
                        shift_dy=int(z['shift_dy_px']), shift_dx=int(z['shift_dx_px'])))
    return out
