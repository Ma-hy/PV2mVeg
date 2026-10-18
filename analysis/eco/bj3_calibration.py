"""Feasibility test of the scene-by-scene cross-calibration of Beijing-3 L4 2 m multispectral data against Sentinel-2 L2A.

Stages (python eco/bj3_calibration.py <stage>):
  s2       download the S2 C1 L2A closest to each BJ-3 acquisition date (B04/B08/SCL, park window, UTM 47N 10 m)
  prep     per scene, reproject BJ-3 2 m red/NIR to a 2 m UTM grid aligned with S2, co-register, aggregate to 10 m, save npz
  cal      per-scene regression, spatial block cross-validation, NDVI comparison, residual diagnostics, mixed-pixel test
  overlap  NDVI difference in the overlaps of adjacent scenes after calibration
  prof     2 m distance-to-panel × direction profiles of typical arrays (native GK grid)
  yoy      year-to-year difference of stable pixels, 2025-08 -> 2026-08
All outputs in data/eco/bj3cal/. Reads only outputs/ and data/raw/.
"""
from __future__ import annotations
# Mirrored from eco/bj3_calibration.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[2] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401

import glob
import os

import pv2mveg  # noqa: F401  (PROJ/GDAL data of the rasterio wheel, no WMI)
from rasterio.transform import Affine

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
OUT = os.path.join(ROOT, 'data', 'eco', 'bj3cal')
S2DIR = os.path.join(OUT, 's2')
CACHE = os.path.join(OUT, 'cache')
FONT = FONT_FILE
PROXY = 'http://127.0.0.1:10808'
UTM = 'EPSG:32647'
# park window (UTM 47N)
WX0, WY0, WX1, WY1 = 597000, 3967020, 655020, 4014000
W10 = Affine(10, 0, WX0, 0, -10, WY1)
NC10, NR10 = (WX1 - WX0) // 10, (WY1 - WY0) // 10

# BJ-3 acquisition date -> S2 dates (first choice + alternative), by the tile cloud cover in stac_candidates.json; see README
S2_DATES = ['20250612', '20250615', '20250812', '20250814', '20250817', '20250819', '20260802', '20260804']


# ---------------------------------------------------------------- scene list
MS_DIRS = {'06': ('2025-06', 'outputs/06/talatan_pv_fraction_2m_06.tif'),
           '08': ('2025-08', 'outputs/08/talatan_pv_fraction_2m_08.tif'),
           '2026_08': ('2026-08', 'outputs/08/talatan_pv_fraction_2m_08.tif')}  # no mask for 2026: borrows 2025-08 (see README)
SCL_OK = (4, 5, 7)
SCL_BAD = (0, 1, 2, 3, 8, 9, 10, 11)
PSF_SIGMA_PX = 2.4  # 2 m pixels; Nyquist MTF of the S2 10 m bands ≈0.3 -> total Gaussian σ≈4.9 m, minus about 1 m of BJ-3 itself


def scene_list():
    out = []
    for d, (period, frac) in MS_DIRS.items():
        for f in sorted(glob.glob(os.path.join(ROOT, 'data', 'raw', d, '*_MS*.tif'))):
            name = os.path.basename(f)
            parts = name.split('_')
            ts = parts[2]
            short = f"{parts[0]}_{parts[1]}_{ts[:8]}_{ts[8:14]}_SC{parts[6].replace('.tif', '')}"
            out.append(dict(dir=d, period=period, path=f, sid=short, date=ts[:8], utc=ts[8:14],
                            sensor=parts[0], frac=os.path.join(ROOT, frac)))
    return out


# ---------------------------------------------------------------- L1/L2/L4/S2 traceability
L1L2 = os.path.join(OUT, 'l1l2')
TRACE_CASES = [  # (label, L1/L2 file, L4 file of the same scene, S2 date)
    ('2025-06 L1', 'BJ3A1_MS1_20250616042256_L1_105DC1_SC_002_0120250618010001_001',
     'data/raw/06/BJ3A1_MS1_20250616042256_L4_105DC1_SC_002_0120250618010001_001.tif', '20250615'),
    ('2025-06 L2', 'BJ3A1_MS1_20250616042256_L2_105DC1_SC_002_0120250618010001_001',
     'data/raw/06/BJ3A1_MS1_20250616042256_L4_105DC1_SC_002_0120250618010001_001.tif', '20250615'),
    ('2025-08 L1', 'BJ3A1_MS1_20250817042257_L1_10617A_SC_002_1121250905002001_001',
     'data/raw/08/BJ3A1_MS1_20250817042257_L4_10617A_SC_002_1121250905002001_001_ORTHO_MS.tif', '20250817'),
    ('2026-08 L1', 'BJ3A1_MS1_20260804041228_L1_107835_SC_002_1121260601007003_001',
     'data/raw/2026_08/BJ3A1_MS1_20260804041228_L4_107835_SC_002.tif', '20260804'),
]


# ---------------------------------------------------------------- 2 m profiles (typical arrays)
def _sid_from_key(key):
    p = key.split('_')
    cam = p[1][1:]
    return f"{p[0]}_MS{cam}_{p[2][:8]}_{p[2][8:14]}_SC{p[6]}"
