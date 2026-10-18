"""Paired ecological reanalysis with frozen matched controls and Landsat baseline."""
# Mirrored from experiments/bright_geometry/fixed_ecology.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
import dataio as D
from fixed_release import ROOT
from fixed_core import ecology_dir,POLICIES
os.environ['CHRONO_SCOPE']='ew'
sys.path.insert(0,str(CODE/'eco/chrono'))
sys.path.insert(0,str(CODE/'eco/ts'))
import cc
import a1_chrono as A1
import j_review as J

METRICS=['ndvi','nirv','msavi','red','nir']
GROUPS=[(str(i),i,i) for i in range(15)]+[('0–2',0,2),('3–5',3,5),('≥6',6,99),('≥10',10,99),('all',0,99)]
