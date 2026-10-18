"""K-BG: is the higher inter-row NDVI green vegetation? Ratio-type greenness indices of the blue and green bands (BJ3A1 scenes only; east-west fixed-array treated units only).

Usage: python eco/chrono/k_bg_greenness.py
Scope: CHRONO_SCOPE=ew (forced by this script), i.e. cc.scope_filter (the rule of eco/ts/ew0_cells.py: fixed east-west pixels >=90% of the pixels with rows in the unit, stations of the exclusion list removed); all north-south (tracking) arrays removed.
Acquisitions: 2025-06 (all A1), 2025-08 (only units whose modal inter-row scene is BJ3A1; the blue/green sums themselves contain only A1 pixels), 2026-08 (all A1).
Data: bj3ndvi_<tag>_30m_bg.tif of eco/bj3ndvi bgprod.py (calibrated blue, green, red and NIR Σ and n of one pixel set).
Treated units, matched controls and scene matching exactly as in a1_chrono.py: the 150 m units of ts and their 10 matched controls (detection-year timing),
only controls with the modal scene of the treated unit (>=2; scenes from celltab_<tag>_raw.parquet); control value = pure ground >300 m outside arrays (>=2000 pixels),
treated value = inter-row core (>=50 pixels); the lane core is computed separately (>=50 pixels).
Indices (from mean reflectances): GRVI=(G−R)/(G+R); VARI=(G−R)/(G+R−B); ExG=2g−r−b (chromatic coordinates g=G/(R+G+B) etc.); G/R; NDVI.
Also the reflectance differences and relative differences (difference/control mean) of blue, green, red and NIR.
No pre-construction baseline (the Landsat baseline has red/NIR only), so Δ is the direct difference "inter-row − same-scene matched controls" (NDVI also gives Δadj for comparison).
Intervals: bootstrap clustered by station × detection year (j_review.cboot, B=2000, linked controls); the net difference (old − new) from the difference distribution of the same cluster draws.
Output: data/eco/chrono/ew/k_bg_{cells,groups,net,panel}.csv
"""
# Mirrored from eco/chrono/k_bg_greenness.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
import warnings

os.environ['CHRONO_SCOPE'] = 'ew'
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np  # noqa: E402
import cc  # noqa: E402
import j_review as J  # noqa: E402

warnings.simplefilter('ignore')
es = cc.es
PER = [('2025-06', '06'), ('2025-08 (A1)', '08'), ('2026-08', '2026_08')]
THR, CTRL_MIN, SAME_MIN = 50, 2000, 2
IDXS = ['grvi', 'vari', 'exg', 'g_r', 'ndvi']
BANDS = ['blue', 'green', 'red', 'nir']
AGES = [('0–1', 0, 1), ('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('≥10', 10, 99)]
NETS = [('≥6 − 0–1', '≥6', '0–1'), ('3–5 − 0–1', '3–5', '0–1'), ('≥10 − 0–1', '≥10', '0–1')]
rng = np.random.default_rng(20261006)


def indices(b, g, r, n):
    with np.errstate(invalid='ignore', divide='ignore'):
        s = b + g + r
        return dict(blue=b, green=g, red=r, nir=n, grvi=(g - r) / (g + r), vari=(g - r) / (g + r - b),
                    exg=(2 * g - r - b) / s, g_r=g / r, ndvi=(n - r) / (n + r))


def cell_vals(B, groups):
    n = sum(cc.cell_sum(B[f'{gp}_n']) for gp in groups)
    m = {}
    for b in BANDS:
        s = sum(cc.cell_sum(B[f'{gp}_sum_{b}']) for gp in groups)
        with np.errstate(invalid='ignore', divide='ignore'):
            m[b] = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return indices(m['blue'], m['green'], m['red'], m['nir']), n
