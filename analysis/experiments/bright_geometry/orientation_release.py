"""Versioned, evidence-gated direction repair and three-clearance native replay.

The frozen v1 masks, parent cohort, projection knots and radiometry are inputs.
Only previously uncertain direction pixels passing the separate repair evidence
can enter the new qualification. All other exclusions retain their old meaning.
"""
# Mirrored from experiments/bright_geometry/orientation_release.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for _k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(_k, '1')
import dataio as D  # Initialise the project PROJ/GDAL environment before rasterio.
import fixed_core as F

OLD_ROOT = D.ROOT/'data/eco/releases/fixed_ew_v1'
NEW_ROOT = D.ROOT/'data/eco/releases/fixed_ew_v2_orientation'
REPORT = D.ROOT/'paper_eco/revisions/orientation_repair_v1/core'
TAGS = ['06', '08', '2026_08']
CLEARANCES = ['1m', '1p5m', '2m']
POLICIES = dict(zip(CLEARANCES, ['clearance1p0', 'clearance1p5', 'exact']))
DIRS = {'1m':'ecology_1m_v3', '1p5m':'ecology_1p5_v3', '2m':'ecology_exact_v2'}
LAYERS = dict(zip(CLEARANCES, ['pure_interrow_2m', 'clearance_1p5m_2m', 'core_2m_2m']))
G = {}
HANDLES=[]


def old_ecology(clearance, tag):
    return OLD_ROOT/DIRS[clearance]/tag
