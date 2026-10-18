"""Shared functions of eco/ts: paths, 30 m grid, Landsat scene reading and QA, sensor harmonisation, indices, aggregation to 150 m units, figure font.

Every script starts with `import common` (this module adds src/ to the path and imports talatan).
"""
# Mirrored from eco/ts/common.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
import pv2mveg  # noqa: F401,E402  (PROJ/GDAL data of the rasterio wheel, no WMI)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

LS_DIR = os.path.join(ROOT, 'data', 'eco', 'landsat')
LS_INDEX = os.path.join(LS_DIR, 'index.csv')
GRID = os.path.join(ROOT, 'data', 'eco', 'grid')
ANN = os.path.join(ROOT, 'data', 'eco', 'ls_annual')
OUT = os.path.join(ROOT, 'data', 'eco', 'ts')
FIG = os.path.join(OUT, 'figs')
CLIMATE = os.path.join(ROOT, 'data', 'eco', 'climate', 'agera5_annual.csv')
for _d in (ANN, OUT, FIG):
    os.makedirs(_d, exist_ok=True)

# analysis scope: TS_SCOPE empty = full sample (results in data/eco/ts/); TS_SCOPE=ew = fixed-array (east-west) scope:
# unit tables and all results go to data/eco/ts/ew/ (the shared inputs cell_y.npz, annual composites, S2 unit composites and detections are still read from data/eco/ts/).
SCOPE = os.environ.get('TS_SCOPE', '').strip()
ROOTR = os.path.join(OUT, SCOPE) if SCOPE else OUT  # root of the results (not to be confused with the pixel size RES)
CELLS = os.path.join(ROOTR, 'cells_pool.parquet')
FIGR = os.path.join(ROOTR, 'figs')
for _d in (ROOTR, FIGR):
    os.makedirs(_d, exist_ok=True)

H, W = 1567, 1935
X0, Y1, RES = 596985.0, 4014015.0, 30.0
CRS = 'EPSG:32647'
YEARS = list(range(2000, 2027))

BANDS = ['blue', 'green', 'red', 'nir08', 'swir16', 'swir22']
SHORT = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
SCALE, OFFSET = 0.0000275, -0.2

# QA_PIXEL bits
B_FILL, B_DILATED, B_CIRRUS, B_CLOUD, B_SHADOW, B_SNOW, B_CLEAR, B_WATER = 0, 1, 2, 3, 4, 5, 6, 7
BAD = (1 << B_FILL) | (1 << B_DILATED) | (1 << B_CLOUD) | (1 << B_SHADOW) | (1 << B_SNOW)

# Roy et al. (2016, RSE 185:57-70) Table 2, ETM+ -> OLI surface reflectance (OLI = a + b * ETM+)
ROY_RMA = {'blue': (-0.0095, 0.9785), 'green': (-0.0016, 0.9542), 'red': (-0.0022, 0.9825),
           'nir': (-0.0021, 1.0073), 'swir1': (-0.0030, 1.0171), 'swir2': (0.0029, 0.9949)}
ROY_OLS = {'blue': (0.0003, 0.8474), 'green': (0.0088, 0.8483), 'red': (0.0061, 0.9047),
           'nir': (0.0412, 0.8462), 'swir1': (0.0254, 0.8937), 'swir2': (0.0172, 0.9071)}
HARM_JSON = os.path.join(OUT, 'harmonization.json')

FONT_PATH = FONT_FILE


# ------------------------------------------------------------------ grid layers
def grid(name):
    import rasterio
    with rasterio.open(os.path.join(GRID, f'{name}_30m.tif')) as r:
        a = r.read(1)
    return a


# ------------------------------------------------------------------ indices
def ndvi(red, nir):
    with np.errstate(divide='ignore', invalid='ignore'):
        return (nir - red) / (nir + red)


def nirv(red, nir):
    return ndvi(red, nir) * nir


def msavi(red, nir):
    with np.errstate(invalid='ignore'):
        return (2 * nir + 1 - np.sqrt((2 * nir + 1) ** 2 - 8 * (nir - red))) / 2


INDEX_FUN = {'ndvi': ndvi, 'nirv': nirv, 'msavi': msavi}

# ------------------------------------------------------------------ 150 m units (5x5 Landsat pixels, aligned with the grid, from the upper left)
CS = 5
CH, CW = H // CS, W // CS   # 313 x 387


def cell_view(a):
    """(..., H, W) -> (..., CH, CS, CW, CS) view (incomplete rows and columns of fewer than 5 dropped)."""
    a = a[..., :CH * CS, :CW * CS]
    return a.reshape(a.shape[:-2] + (CH, CS, CW, CS))


def cell_pixels(a):
    """(H,W) -> (CH*CW, 25) pixel values of each unit."""
    v = cell_view(a)
    return np.moveaxis(v, -3, -2).reshape(CH * CW, CS * CS)


# ------------------------------------------------------------------ climate
def climate():
    c = pd.read_csv(CLIMATE).set_index('year')
    return c
