"""Shared part of eco/panprof: acquisitions, treated units of the main a1 rule, same-scene matched control pairs, 2 km blocks, age groups and batches.

Reuses the paths, scene tables, window reading and grid offsets of eco/adjacency/adjc.py (read-only calls, no changes).
Rules identical to section K of eco/chrono (main a1 rule):
  treated units = units of data/eco/chrono/ew/a1_<tag>_det_thr50[_BJ3A1]_cells.csv with valid_same and finite same_ndvi_dadj
            (east-west fixed-array selection, exclusion list of anomalous arrays and per-acquisition exclusion lists already applied; in 2025-08 only units whose modal scene is BJ3A1);
  controls  = the matched controls of the unit in data/eco/ts/es_raw_sum_ndvi_det_treated/matched.npz
            whose modal scene (ctrl_scene) equals the modal scene of the treated unit (ir_scene) and that have >=2000 2 m pure-ground pixels >300 m from arrays (the same rule of a1).
"""
# Mirrored from eco/panprof/pp.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
sys.path.insert(0, os.path.join(CODE, 'eco', 'adjacency'))
sys.path.insert(0, os.path.join(CODE, 'eco', 'ts'))
import adjc as A  # noqa: E402  (import talatan: PROJ_LIB, WMI workaround)


OUT = os.path.join(ROOT, 'data', 'eco', 'panprof')
os.makedirs(OUT, exist_ok=True)
CHR_EW = os.path.join(ROOT, 'data', 'eco', 'chrono', 'ew')
CHR = os.path.join(ROOT, 'data', 'eco', 'chrono')
TS = os.path.join(ROOT, 'data', 'eco', 'ts')
TAGS = ['06', '08', '2026_08']
LABEL = {'06': '2025-06', '08': '2025-08', '2026_08': '2026-08'}
YEAR = {'06': 2025, '08': 2025, '2026_08': 2026}
CELLS_CSV = {'06': 'a1_06_det_thr50_cells.csv', '08': 'a1_08_det_thr50_BJ3A1_cells.csv', '2026_08': 'a1_2026_08_det_thr50_cells.csv'}
CTRL_MIN_N = 2000
SAME_MIN = 2
DBIN = 0.5            # width of the distance bins (m)
NBIN = 12             # 0–6 m
SIDES = ['N', 'S', 'C']   # N = pixel on the north side of a row (nearest row band to the south, next to the north edge of the modules / shadow); S = pixel on the south side of a row; C = equidistant (row middle)
AGE_GROUPS = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99)]
BATCHES = [('2013–14', 2012, 2014), ('2015–17', 2015, 2017), ('2018–21', 2018, 2021), ('2022–26', 2022, 2026)]


_ES = {}


def log(msg, fn='log_panprof.txt'):
    import time
    s = time.strftime('%H:%M:%S') + ' ' + str(msg)
    print(s, flush=True)
    with open(os.path.join(OUT, fn), 'a', encoding='utf-8') as f:
        f.write(s + '\n')


CLEAR_CLASSES = [('clear ≤5.0 m', 0, 5.0), ('clear 5.5–6.0 m', 5.5, 6.0), ('clear ≥6.5 m', 6.5, 99)]
