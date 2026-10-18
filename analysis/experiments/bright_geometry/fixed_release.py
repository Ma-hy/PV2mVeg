"""Versioned fixed-array release; original products and NS experiments stay separate."""
# Mirrored from experiments/bright_geometry/fixed_release.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import dataio as D

ROOT=D.ROOT/'data/eco/releases/fixed_ew_v1'
REVIEW=D.OUT/'candidate_2026_v2/human_review/user_review_606619e321fd3196.json'
