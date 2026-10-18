"""Additional comparison-ground variants of the unmatched buffer baseline (v1.1 review).

Adds, with 2,000 whole-station resamples and comparison ground fixed:
  ring_cell500_R5km      the main ring, restricted to cells whose 25 Landsat pixels are
                         all >= 500 m from the array union (pool rule R2, cells_pool.dist_min)
  ring300_1000_R5km_pix  the main ring with pixel-count weights instead of equal cells
  beyond1000_R8km_supp   ground beyond 1 km within 8 km, restricted to the pool's slope,
                         elevation and pre-construction support rules (R7, R8)
  pool_R5km              every control-pool cell (all rules R1-R8) within 5 km, same scene,
                         using its >300 m ground: the matching population without matching
Writes analysis/buffer/sensitivity_extra.csv and appends a note to complete.json.
"""
# Mirrored from experiments/bright_geometry/buffer_variants.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
import warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import buffer_baseline as BB

N, E, cc, J = BB.N, BB.E, BB.cc, BB.J
GROUPS = BB.GROUPS


def support_rules(cells):
    Y = np.load(os.path.join(cc.TS, 'cell_y.npz'))
    nd, red = Y['sum_ndvi'], Y['sum_red']
    t = cells[cells.treated_det | cells.treated_lab]
    pre_nd = np.nanmean(nd[:, :11], 1); pre_red = np.nanmean(red[:, :11], 1)
    lo_e, hi_e = np.percentile(t.elev, [1, 99]); s99 = np.percentile(t.slope, 99)
    tn, tr = pre_nd[t.index], pre_red[t.index]
    r7 = (cells.slope.values <= s99) & (cells.elev.values >= lo_e - 50) & (cells.elev.values <= hi_e + 50)
    r8 = ((pre_nd >= np.nanpercentile(tn, .5)) & (pre_nd <= np.nanpercentile(tn, 99.5)) &
          (pre_red >= np.nanpercentile(tr, .5)) & (pre_red <= np.nanpercentile(tr, 99.5)))
    return r7 & r8
