"""Remeasure three-period core ground with frozen radiometry and parent cells.

Legacy replay, old-with-quality, and new-with-identical-quality are accumulated
on the same native 2m centres. No formal product or ecological table is changed.
"""
# Mirrored from experiments/bright_geometry/fixed_core.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import dataio as D
from fixed_release import ROOT
import b_component as B
from rasterio.windows import Window,from_bounds
from core_geometry import exact_center_distance

CORE=256;HALO=200;NVAR=3;NSC=64
NAMES=['legacy_replay','old_qualified','new_qualified']
GRID_SCHEME='original_1024_knots_complete_tile_coverage_v3'
POLICIES=['legacy','exact','clearance1p5','clearance1p0','ground']
THRESHOLDS={'legacy':2.,'exact':2.,'clearance1p5':1.5,'clearance1p0':1.,'ground':0.}
G={}


def read_bounds(ds,bounds,fill,shape):
    win=from_bounds(*bounds,ds.transform).round_offsets().round_lengths()
    a=ds.read(1,window=win,boundless=True,fill_value=fill)
    assert a.shape==shape,(ds.name,a.shape,shape)
    return a


def canonical_centres(tag,gr,r0,c0,h,w):
    # Original product interpolates projection knots every 32 native pixels,
    # with a shortened last interval at its 1024-pixel tile edge. Preserve
    # those knots when measuring smaller tiles: ending at pixel 255 instead
    # of knot 256 can shift a centre by <1 mm across a UTM aggregation edge.
    rh=h+int(r0+h<gr['height'] and (r0+h)%1024!=0)
    cw=w+int(c0+w<gr['width'] and (c0+w)%1024!=0)
    ux,uy=D.R.utm_grid_coords(tag,gr['x0']+c0*2,gr['y0']-r0*2,rh,cw)
    return ux[:h,:w],uy[:h,:w]


def ecology_dir(policy):
    assert policy in POLICIES
    return ROOT/{'legacy':'ecology','exact':'ecology_exact_v2','clearance1p5':'ecology_1p5_v3','clearance1p0':'ecology_1m_v3','ground':'ecology_ground_v3'}[policy]
