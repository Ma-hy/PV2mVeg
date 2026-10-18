"""B/C layers of the module-projection rule and propagation of boundary errors (per acquisition, whole park in tiles).

Module-projection mask (east-west arrays): pixels of the 0.5 m row-band mask with orient_10m = 1 (east-west); per column, the north end of each north-south run is eroded by e metres,
e = the distance "north edge of the mask to the dark-strip/module boundary" measured in section A (eco/rowband/a/erosion_table.csv, per scene × construction batch; the scene from the 0.5 m source map,
the batch from build_year_10m). When e is not a multiple of 0.5 m, the columns are dithered in order (golden-ratio sequence) so that the expected eroded area matches e;
quantisation error of a single column ±0.25 m. North-south and unoriented pixels are not eroded (module rule = row-band rule).
0.5 m result data/eco/rowband/b/<tag>/talatan_pv_mask_0p5m_<tag>_comp.tif: 1 = module, 2 = dark strip (eroded part of the row band), 0 = not panel, 255 = no image.

2 m layers (data/eco/grid2m/<tag>/component/, on the grid of the eco/bj3ndvi/<tag> products, single-acquisition mask, algorithm of eco/bj3ndvi/product.py):
  cover_2m_comp       module cover % (4×4 sub-pixels), 255 = no image
  dark_2m_comp        dark-strip cover % (eroded part of the row band)
  dist_ctr_2m_comp    distance from the cell centre to the nearest module pixel centre (m, capped at 100), NaN = no image; to the module edge = this value − 0.25
  dir_2m_comp         direction of the nearest module pixel: 1 north 2 east 3 south 4 west, 0 = the centre sub-pixel is module, 9 = none within 100 m, 255 = no image
  dist_rb_2m_comp     the same algorithm to the nearest row-band pixel (= the original rule; on the north side, to the outer edge of the dark strip) (m, capped at 100)
  dist_dark_2m_comp   distance to the nearest dark-strip pixel (m, capped at 100)
  dir_rb_2m_comp      direction of the nearest row-band pixel (defined as dir_2m). 3 = row band to the south -> the cell is on the north side of a row (next to the dark strip); 1 = on the south side of a row (next to the modules)
  class_2m_comp       pixel classes redrawn with the module rule (codes as class_2m): 0 no data, 1 module (module cover >0), 2 panel edge (ground inside the array, <2 m from the module edge), 3 inter-row core (>=2 m from the module edge, with modules within one row pitch on both sides across the rows),
                      4 lane core, 5 pure ground outside arrays, 6 near modules outside arrays; dark-strip cells (row-band cover >0, module cover = 0) are classed as ground 2–6,
                      distinguished by bit 16 of core_flags or by dark_2m_comp; QA rules unchanged (removal still follows qa_2m)
  core_flags_2m_comp  bits: 1 = inter-row core of the row-band rule (recomputed, should equal class_2m == 3), 2 = inter-row core of the module rule,
                      4 = core of the module rule with a footprint >=0.5 m from the dark strip/row band (>=1.5 m from the row-band edge), 8 = on the north side of a row (dir_rb = 3),
                      16 = the cell contains dark-strip sub-pixels, 32 = east-west (orient 1), 64 = passed QA (counted)
10 m / 30 m (EPSG:32647): 0.5 m sub-pixel counts (module, row band, with image) summed by the 2 m cell centres falling in each cell,
  data/eco/grid/panel_frac_<tag>_{10m,30m}_comp.tif (module fraction), with the row-band fractions, east-west shares etc. recomputed under data/eco/rowband/b/grid/.
C (propagation): only on tiles with east-west treated units (cells_pool + scope_filter with CHRONO_SCOPE=ew), for 3 pixel selections × 11 boundary cases,
  Σred, Σnir and n summed per 150 m unit (calibrated reflectance from eco/bj3ndvi/<tag>/red_2m, nir_2m; qa, scene and array area as in the product),
  result data/eco/rowband/b/<tag>/c_cells.npz.
Usage: python eco/rowband/b_component.py --tag 06 [--workers 12] [--tiles r0,c0;...]
"""
# Mirrored from eco/rowband/b_component.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rb_common as R  # noqa: E402

import numpy as np  # noqa: E402


CORE = 1024        # 2 m cells
HALO = 200         # 0.5 m pixels = 100 m
CAP = 100.0
CLOUD_BUF = 5
QA_EXCLUDE = 1 | 2 | 8 | 32 | 64 | 128
PHI = 0.6180339887498949
# boundary cases (applied to the row-band mask; the module rule is eroded again on the perturbed row band); dark± changes only the module erosion width
SCEN = ['base', 'wid-0.5', 'wid+0.5', 'shiftN0.5', 'shiftS0.5', 'Nedge-0.5', 'Nedge+0.5', 'Sedge-0.5', 'Sedge+0.5', 'dark-0.5', 'dark+0.5']
VARS = ['rb', 'cp', 'cpx']   # row band >=2 m; module >=2 m; module >=2 m and >=1.5 m from the row-band edge
G = {}


def up(a):      # content moved north by 1 pixel: new[i] = a[i+1]
    o = np.zeros_like(a); o[:-1] = a[1:]; return o


def dn(a):      # content moved south by 1 pixel: new[i] = a[i-1]
    o = np.zeros_like(a); o[1:] = a[:-1]; return o


def perturb(pan, s):
    if s in ('base', 'dark-0.5', 'dark+0.5'):
        return pan
    if s == 'wid-0.5':
        return pan & dn(pan) & up(pan)
    if s == 'wid+0.5':
        return pan | dn(pan) | up(pan)
    if s == 'shiftN0.5':
        return up(pan)
    if s == 'shiftS0.5':
        return dn(pan)
    if s == 'Nedge-0.5':        # north edge moved in (narrower row)
        return pan & dn(pan)
    if s == 'Nedge+0.5':        # north edge moved out
        return pan | up(pan)
    if s == 'Sedge-0.5':
        return pan & up(pan)
    if s == 'Sedge+0.5':
        return pan | dn(pan)
    raise KeyError(s)


def gaps(pan, nod, h, w, axis):
    """As eco/bj3ndvi/product.py: clearance (m) between the nearest panels on both sides across the rows, through the cell centre; both sides must end at a panel."""
    stop = pan | nod
    BIG = 1 << 30
    if axis == 'ns':
        cols = HALO + np.r_[np.arange(w) * 4 + 1, np.arange(w) * 4 + 2]
        S = stop[:, cols]; P = pan[:, cols]
        idx = np.arange(S.shape[0], dtype=np.int32)[:, None]
        last = np.maximum.accumulate(np.where(S, idx, -1), axis=0)
        nxt = np.minimum.accumulate(np.where(S, idx, BIG)[::-1], axis=0)[::-1]
        ra = HALO + np.arange(h) * 4 + 1
        L, N = last[ra], nxt[ra + 1]
        jj = np.broadcast_to(np.arange(S.shape[1])[None, :], L.shape)
        ok = (L >= 0) & P[np.clip(L, 0, None), jj] & (N < BIG) & P[np.clip(N, None, S.shape[0] - 1), jj]
        g = np.where(ok, (N - L - 1) * 0.5, np.nan)
        both = ok[:, :w] & ok[:, w:]
        return np.where(both, np.fmax(g[:, :w], g[:, w:]), np.nan)
    rows = HALO + np.r_[np.arange(h) * 4 + 1, np.arange(h) * 4 + 2]
    S = stop[rows, :]; P = pan[rows, :]
    idx = np.arange(S.shape[1], dtype=np.int32)[None, :]
    last = np.maximum.accumulate(np.where(S, idx, -1), axis=1)
    nxt = np.minimum.accumulate(np.where(S, idx, BIG)[:, ::-1], axis=1)[:, ::-1]
    ca = HALO + np.arange(w) * 4 + 1
    L, N = last[:, ca], nxt[:, ca + 1]
    ii = np.broadcast_to(np.arange(S.shape[0])[:, None], L.shape)
    ok = (L >= 0) & P[ii, np.clip(L, 0, None)] & (N < BIG) & P[ii, np.clip(N, None, S.shape[1] - 1)]
    g = np.where(ok, (N - L - 1) * 0.5, np.nan)
    both = ok[:h] & ok[h:]
    return np.where(both, np.fmax(g[:h], g[h:]), np.nan)


_MX0 = {}


SPEC = [('cover', 'cover_2m_comp.tif', 'uint8', 255), ('dark', 'dark_2m_comp.tif', 'uint8', 255),
        ('dist', 'dist_ctr_2m_comp.tif', 'float32', float('nan')), ('dir', 'dir_2m_comp.tif', 'uint8', 255),
        ('dist_rb', 'dist_rb_2m_comp.tif', 'float32', float('nan')), ('dist_dark', 'dist_dark_2m_comp.tif', 'float32', float('nan')),
        ('dir_rb', 'dir_rb_2m_comp.tif', 'uint8', 255), ('cls', 'class_2m_comp.tif', 'uint8', None),
        ('flags', 'core_flags_2m_comp.tif', 'uint8', None)]


def grid_of(tag):
    """2 m computation grid (as eco/bj3ndvi/cfg.pcfg)."""
    if tag in ('06', '08'):
        return dict(x0=329890.0, y0=4013270.0, width=26146, height=20420)
    import rasterio
    with rasterio.open(R.PER[tag]['mask']) as d:
        t = d.transform
        return dict(x0=t.c, y0=t.f, width=int(np.ceil(d.width / 4)), height=int(np.ceil(d.height / 4)))
