"""A5 artefact range over the model set (east-west fixed arrays): joint sensitivity to kernel × panel end-member alternative × shadow scheme, with each kernel labelled "rejected / not rejected" by the k3 hold-out validation.

End-member (panel) schemes: p10 (main: row centre >=1.0 m), p05, p15 (row-centre width), swap (June and August exchanged: median end-member of the same camera in the other acquisition; 2026-08 exchanged with 2025-06),
  x0.8, x1.2 (both bands scaled), r+n− (red ×1.2, NIR ×0.8: lower panel NDVI), r−n+ (red ×0.8, NIR ×1.2).
  Inversion: p = (p_obs − g(1 − w_int))/w_int (same kernel); for units without one, the median of the same camera and acquisition.
Shadow: none (dark strip merged with the panel, main rule); shadow (dark strip as a third component, reflectance = panchromatic dark-strip/ground ratio s × control ground; panel end-member inverted from row-centre pixels without dark strip);
  shadow_blue (as shadow, with the NIR of the dark strip ×0.85 to mimic blue-shifted skylight illumination: lower dark-strip NDVI). Computed only for main, wide, redabs, myoung, npkr.
Hold-out labels of the kernels (rejected_ref of k3_holdout_summary.csv): red component rejected on both red|pan|out and red|pan|young = rejected, on one only = partly rejected (still in the range);
  NIR component rejected on nir|mask|young = rejected. s08/s12 (width ±20%) are not tested separately.
Intervals: bootstrap clustered by station × detection year (grp of the a1 table), B = 1000.
Output: data/eco/adjacency/ew_only/a5_{joint,shadow,modelset}.csv
"""
# Mirrored from eco/adjacency/a5_model_set.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import adjc as A  # noqa: E402
import a2_cells as B  # noqa: E402

EW = os.path.join(A.OUT, 'ew_only')
AGEG = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99)]
SWAP = {'06': '08', '08': '06', '2026_08': '06'}
RED_MAP = {'box': 'box', 'main': 'main', 'equal': 'main', 'shiftN': 'main', 'shiftS': 'main', 'wide': 'wide', 'sharp': 'sharp',
           'npk': 'npk', 'npkr': 'npk', 'myoung': 'longtail', 'mold': 'longtail_old', 'redabs': 'redabs', 's08': None, 's12': None}
NIR_MAP = {'box': 'box', 'main': 'main', 'equal': 'equal', 'shiftN': 'main', 'shiftS': 'main', 'wide': 'narrow', 'sharp': 'sharp',
           'npk': 'npk', 'npkr': None, 'myoung': 'main', 'mold': 'longtail_old', 'redabs': 'main', 's08': None, 's12': None}
