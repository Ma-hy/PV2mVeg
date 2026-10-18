"""Paths and settings.

Every script finds its inputs and outputs through the object `P`:

    from pv2mveg.config import P
    effects = P.cs12 / 'effects.csv'

Settings are read in this order (later sources win):
  1. defaults: layout "repo", data root <repository>/data/derived, output folder <repository>/output;
  2. config/paths.toml (copy config/paths.example.toml; the copy is not tracked by git);
  3. a TOML file named by PV2MVEG_CONFIG;
  4. environment variables PV2MVEG_LAYOUT, PV2MVEG_DATA_ROOT, PV2MVEG_PROJECT_ROOT, PV2MVEG_OUTPUT and
     PV2MVEG_FONT.

Layouts
  repo      the de-identified derived tables shipped in data/derived. They reproduce every figure, table and
            number marked "repository data" in README.md.
  project   the data tree of the research project (<project_root>/data and <project_root>/paper_eco), which
            also holds the restricted inputs: Beijing-3 products, hand-drawn station polygons and construction
            years, unit coordinates. Needed for the figures marked "restricted" and for the analysis scripts.

Both layouts expose the same logical locations (attributes of P); a file has the same name and the same
columns in both, except the identifying columns removed from the repository copies (data/README.md).
"""
from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TAGS = ['06', '08', '2026_08']                       # acquisitions 2025-06, 2025-08, 2026-08
LAYOUTS = ('repo', 'project')


def _read_toml(path: Path) -> dict:
    with open(path, 'rb') as f:
        return tomllib.load(f)


def _resolve(p: str | os.PathLike | None, base: Path = REPO) -> Path | None:
    if p is None or str(p).strip() == '':
        return None
    q = Path(os.path.expandvars(os.path.expanduser(str(p))))
    return q if q.is_absolute() else (base / q)


@dataclass
class Settings:
    layout: str = 'repo'
    data_root: Path = REPO / 'data' / 'derived'
    project_root: Path | None = None
    output: Path = REPO / 'output'
    font: Path | None = None
    sources: list = field(default_factory=list)

    @classmethod
    def load(cls) -> 'Settings':
        s = cls()
        files = [REPO / 'config' / 'paths.toml']
        if os.environ.get('PV2MVEG_CONFIG'):
            files.append(Path(os.environ['PV2MVEG_CONFIG']))
        for f in files:
            if f.exists():
                c = _read_toml(f)
                d, o, g = c.get('data', {}), c.get('output', {}), c.get('figures', {})
                s.layout = d.get('layout', s.layout)
                s.data_root = _resolve(d.get('data_root')) or s.data_root
                s.project_root = _resolve(d.get('project_root')) or s.project_root
                s.output = _resolve(o.get('root')) or s.output
                s.font = _resolve(g.get('font')) or s.font
                s.sources.append(str(f))
        env = os.environ
        s.layout = env.get('PV2MVEG_LAYOUT', s.layout)
        s.data_root = _resolve(env.get('PV2MVEG_DATA_ROOT')) or s.data_root
        s.project_root = _resolve(env.get('PV2MVEG_PROJECT_ROOT')) or s.project_root
        s.output = _resolve(env.get('PV2MVEG_OUTPUT')) or s.output
        s.font = _resolve(env.get('PV2MVEG_FONT')) or s.font
        if s.layout not in LAYOUTS:
            raise ValueError(f'layout must be one of {LAYOUTS}, not {s.layout!r}')
        if s.layout == 'project' and s.project_root is None:
            raise ValueError('layout "project" needs project_root (config/paths.toml or PV2MVEG_PROJECT_ROOT)')
        return s


class Paths:
    """Logical locations of inputs and outputs for the active layout."""

    def __init__(self, s: Settings):
        self.settings = s
        self.layout = s.layout
        if s.layout == 'repo':
            d = s.data_root
            self.v12 = d / 'v12'                      # v1.2 analysis results
            self.ts = d / 'ts'                        # Landsat unit series and the mixed-pixel event study
            self.legacy = d / 'legacy'                # tables of earlier figure code used as inputs (Fig. 3)
            self.legacy_v110 = d / 'legacy' / 'v110'
            self.frozen = d / 'frozen'                # frozen result tables of earlier releases (Fig. S12)
            self.purify = d / 'purify'                # panel end-member and unmixing tests (Fig. S2)
            self.summaries = d / 'summaries'          # side outputs of the restricted map figures, used by numbers
            self.units = d / 'units'                  # de-identified unit table (origin, timing offsets)
            self.project = None
        else:
            r = s.project_root
            self.project = r
            self.eco = r / 'data' / 'eco'
            self.release = self.eco / 'releases' / 'fixed_ew_v2_orientation'
            self.analysis = self.release / 'analysis'  # v1.1 analysis root (frozen inputs of v1.2)
            self.v12 = self.analysis / 'v12'
            self.ts = self.eco / 'ts'
            self.legacy = r / 'paper_eco' / 'figdata'
            self.legacy_v110 = self.legacy / 'v110'
            self.frozen = r / 'paper_eco' / 'revisions' / 'full_v4'
            self.purify = self.eco / 'purify'
            self.bj3cal = self.eco / 'bj3cal'
            self.summaries = None                     # produced by the map figures into P.figdata
            self.units = None
            self.stations = self.eco / 'anomaly' / 'v2' / 'hand_stations_32647.gpkg'
            self.rgb_cache = r / 'work' / 'maps_v110'
            # the labelled station outlines of the earlier installation-year check (Fig. S3), kept beside the project
            self.selected_stations = (Path(os.environ['PV2MVEG_SELECTED_STATIONS']) if os.environ.get('PV2MVEG_SELECTED_STATIONS')
                                      else r.parent / 'S2' / 'full_data' / 'SelectedStations.shp')
        self.cs12 = self.v12 / 'cross_scene'
        self.bp12 = self.v12 / 'buffer_primary'
        self.sec12 = self.cs12 / 'secondary'
        # red-edge indices of the auxiliary 1.2 m imagery: folder 'aux' in the project data, which is a reserved
        # device name on Windows and cannot be checked out there; 'red_edge' in the repository
        aux = 'red_edge' if s.layout == 'repo' else 'aux'
        self.aux12 = self.cs12 / aux
        self.sec_aux12 = self.sec12 / aux
        self.es = self.ts / 'ew'
        self.out = s.output
        self.figures = self.out / 'figures'
        self.figdata = self.out / 'figdata'
        self.numbers = self.out / 'numbers'
        self.tables = self.out / 'tables'
        self.checks = self.out / 'checks'
        if self.summaries is None:
            self.summaries = self.figdata

    def need_project(self, what: str):
        """Stop with a clear message when a step needs the restricted project data."""
        if self.layout != 'project':
            raise SystemExit(f'{what} needs restricted inputs: run with layout "project" and project_root set '
                             f'(see README.md, "Data availability").')

    def mkdirs(self):
        for p in (self.figures, self.figdata, self.numbers, self.tables, self.checks):
            p.mkdir(parents=True, exist_ok=True)

    def __repr__(self):
        return f'Paths(layout={self.layout!r}, v12={self.v12}, out={self.out})'


SETTINGS = Settings.load()
P = Paths(SETTINGS)
