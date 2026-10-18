"""E6 summary of the panchromatic texture (output of e3): inter-row core − same-scene matched controls, by age group; 95% bootstrap clustered by station × detection year.

Two comparisons:
  raw     difference of the unit median textures (cv16, hp, strip), and the panchromatic brightness ratio (inter-row/control − 1)
  bmatch  brightness matching: per scene, the control pixel sample is split into 20 brightness quantile bins; the median texture of each bin gives the curve of the "expected texture of ground of the same brightness";
          residual of each pixel = texture − expectation (linear interpolation in brightness), unit value = median residual; inter-row unit residual − mean residual of its control units.
          This excludes the bias "dark surfaces have larger relative noise, hence a naturally higher relative texture" (uniformly darkened soil should leave no positive residual after brightness matching).
Output: data/eco/adjacency/ew_only/e6_texture.csv; figure figs/e6_texture.png
"""
# Mirrored from eco/adjacency/e6_texture_analysis.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import adjc as A  # noqa: E402

EW = os.path.join(A.OUT, 'ew_only')
CELLS_CSV = {'06': 'a1_06_det_thr50_cells.csv', '08': 'a1_08_det_thr50_BJ3A1_cells.csv', '2026_08': 'a1_2026_08_det_thr50_cells.csv'}
AGEG = [('0–1', 0, 1), ('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('≥10', 10, 99)]
MET = ['cv16', 'hp', 'strip']
