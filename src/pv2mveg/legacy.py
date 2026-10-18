"""Roots of the analysis code mirrored from the research repository (analysis/ of this repository).

The mirrored modules import this module before anything else:

  PROJECT            project root of the configuration (layout "project"): every data path of the analysis code is
                     relative to it (data/eco/..., outputs/..., ancillary/...); read only
  CODE               the mirrored code tree (analysis/), for the module search paths that the code sets itself
  V12_OUT            output folder of the v1.2 analysis (<output>/analysis/v12): everything the analysis scripts
                     write, including the snapshot of their frozen inputs (inputs/)
  FONT_FILE          figure font of the diagnostic plots that the analysis writes alongside its tables (the
                     configured font, else STIXGeneral of matplotlib)
  SELECTED_STATIONS  labelled station outlines of the earlier installation-year check (restricted)
"""
import importlib.util
from pathlib import Path

from pv2mveg.config import P, REPO, SETTINGS

P.need_project('The analysis scripts (analysis/, scripts/1x_*.py)')
PROJECT = P.project
CODE = REPO / 'analysis'
V12_OUT = P.out / 'analysis' / 'v12'


def _font_file():
    if SETTINGS.font is not None and SETTINGS.font.exists():
        return str(SETTINGS.font)
    # STIXGeneral of matplotlib, found without importing matplotlib: the mirrored modules import this module before
    # they set the thread counts of the numerical libraries, so numpy must not be loaded here
    return str(Path(importlib.util.find_spec('matplotlib').origin).parent / 'mpl-data' / 'fonts' / 'ttf' / 'STIXGeneral.ttf')


FONT_FILE = _font_file()
SELECTED_STATIONS = P.selected_stations
