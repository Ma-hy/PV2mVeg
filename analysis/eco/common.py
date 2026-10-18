"""Shared definitions of the common grid of the ecological analysis (EPSG:32647).

Grid conventions
- G10: 10 m, pixel edges at multiples of 10 (Sentinel-2 convention). x 597000–655020, y 3967020–4014000, 5802 columns × 4698 rows.
- G30: 30 m, pixel centres at multiples of 30, edges at 15+30k (Landsat convention). Smallest grid covering the AOI:
  x 596985–655035, y 3967005–4014015, 1935 columns × 1567 rows (upper-left corner 596985, 4014015).
- G5: 5 m counting grid, edges at multiples of 5, extent as G30; a 10 m cell = 2×2 G5 cells (offset 3 cells), a 30 m cell = 6×6 G5 cells.
- G25: 2.5 m distance-transform grid, pixel centres at multiples of 2.5; the 10 m and 30 m pixel centres fall on G25 pixel centres.
"""
from __future__ import annotations
# Mirrored from eco/common.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[2] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401

import os

import numpy as np
from rasterio.crs import CRS
from rasterio.transform import from_origin

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
OUT = os.path.join(ROOT, 'outputs')
GRID_DIR = os.path.join(ROOT, 'data', 'eco', 'grid')
GRID2M_DIR = os.path.join(ROOT, 'data', 'eco', 'grid2m')
WORK = os.path.join(ROOT, 'data', 'eco', 'work')
STATIONS = str(SELECTED_STATIONS)

CRS_UTM = CRS.from_epsg(32647)
CRS_GK = CRS.from_epsg(4543)

AOI = (597000.0, 3967020.0, 655020.0, 4014000.0)  # l, b, r, t


class Grid:
    def __init__(self, name, x0, y0, res, width, height):
        self.name, self.x0, self.y0, self.res, self.width, self.height = name, x0, y0, res, width, height
        self.transform = from_origin(x0, y0, res, res)

    @property
    def shape(self):
        return self.height, self.width

    @property
    def bounds(self):
        return self.x0, self.y0 - self.height * self.res, self.x0 + self.width * self.res, self.y0

    def centers(self):
        xs = self.x0 + (np.arange(self.width) + 0.5) * self.res
        ys = self.y0 - (np.arange(self.height) + 0.5) * self.res
        return xs, ys


G10 = Grid('10m', 597000.0, 4014000.0, 10.0, 5802, 4698)
G30 = Grid('30m', 596985.0, 4014015.0, 30.0, 1935, 1567)
G5 = Grid('5m', 596985.0, 4014015.0, 5.0, 11610, 9402)
G25 = Grid('2.5m', 596983.75, 4014016.25, 2.5, 23221, 18805)
GRIDS = {'10m': G10, '30m': G30}

# offset of G10 in G5 (row, column)
G10_OFF = (int(round((G5.y0 - G10.y0) / 5)), int(round((G10.x0 - G5.x0) / 5)))
assert G10_OFF == (3, 3)


def set_mpl_font():
    import matplotlib
    from matplotlib import font_manager
    fp = FONT_FILE  # figure font (pv2mveg configuration)
    font_manager.fontManager.addfont(fp)
    name = font_manager.FontProperties(fname=fp).get_name()
    matplotlib.rcParams['font.family'] = name
    matplotlib.rcParams['axes.unicode_minus'] = False
    return name
